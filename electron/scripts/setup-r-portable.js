#!/usr/bin/env node

// Build a relocatable, minimal R runtime for the macOS Electron application.
// The official CRAN framework embeds absolute /Library/Frameworks paths, so a
// plain copy is not portable. This script copies only the packages used by the
// API, rewrites every Mach-O dependency to @loader_path and verifies the result
// with an isolated environment before electron-builder is allowed to run.

const { execFileSync, spawnSync } = require('child_process');
const fs = require('fs');
const os = require('os');
const path = require('path');

const electronDir = path.resolve(__dirname, '..');
const root = path.resolve(electronDir, '..');
const portableDir = path.join(electronDir, 'R-portable');
const apiDir = path.join(root, 'api');
const requiredPackages = ['plumber', 'jsonlite'];
const basePackages = [
  'base', 'compiler', 'datasets', 'graphics', 'grDevices', 'grid', 'methods',
  'parallel', 'splines', 'stats', 'stats4', 'tools', 'utils',
];
const frameworkPattern = /^\/Library\/Frameworks\/R\.framework\/(?:Versions\/[^/]+\/)?Resources\/?/;

function section(label) {
  process.stdout.write(`\n${label}\n`);
}

function run(file, args, options = {}) {
  const result = spawnSync(file, args, {
    encoding: 'utf8',
    maxBuffer: 64 * 1024 * 1024,
    ...options,
  });
  if (result.status !== 0) {
    const detail = [result.stdout, result.stderr].filter(Boolean).join('\n').trim();
    throw new Error(`${file} ${args.join(' ')} failed${detail ? `:\n${detail}` : ''}`);
  }
  return (result.stdout || '').trim();
}

function rEval(expression) {
  return run('Rscript', ['-e', expression]);
}

function listFiles(directory) {
  const files = [];
  const visit = (current) => {
    for (const entry of fs.readdirSync(current, { withFileTypes: true })) {
      const full = path.join(current, entry.name);
      if (entry.isDirectory()) visit(full);
      else if (entry.isFile()) files.push(full);
    }
  };
  visit(directory);
  return files;
}

function machoDependencies(file) {
  const result = spawnSync('otool', ['-L', file], { encoding: 'utf8' });
  if (result.status !== 0) return null;
  return result.stdout
    .split(/\r?\n/)
    .slice(1)
    .map((line) => line.trim().match(/^(\S+)\s+\(/)?.[1])
    .filter(Boolean);
}

function portableReference(file, original) {
  const suffix = original.replace(frameworkPattern, '');
  const relativeRoot = path.relative(path.dirname(file), portableDir).split(path.sep).join('/');
  return `@loader_path/${relativeRoot ? `${relativeRoot}/` : ''}${suffix}`;
}

function patchMachOFiles() {
  section('Relocating Mach-O dependencies');
  const candidates = listFiles(portableDir).filter((file) => {
    const extension = path.extname(file).toLowerCase();
    if (extension === '.so' || extension === '.dylib') return true;
    return (fs.statSync(file).mode & 0o111) !== 0;
  });
  let machoCount = 0;
  let changedCount = 0;

  for (const file of candidates) {
    const dependencies = machoDependencies(file);
    if (!dependencies) continue;
    machoCount += 1;
    let changed = false;

    for (const dependency of dependencies) {
      if (!frameworkPattern.test(dependency)) continue;
      const replacement = portableReference(file, dependency);
      run('install_name_tool', ['-change', dependency, replacement, file]);
      changed = true;
    }

    const identity = spawnSync('otool', ['-D', file], { encoding: 'utf8' });
    const currentId = identity.status === 0
      ? identity.stdout.split(/\r?\n/).slice(1).find((line) => frameworkPattern.test(line.trim()))
      : null;
    if (currentId) {
      const suffix = currentId.trim().replace(frameworkPattern, '');
      run('install_name_tool', ['-id', `@rpath/${suffix}`, file]);
      changed = true;
    }

    if (changed) {
      run('codesign', ['--force', '--sign', '-', file]);
      changedCount += 1;
    }
  }

  const unresolved = [];
  const external = [];
  for (const file of candidates) {
    const dependencies = machoDependencies(file);
    if (!dependencies) continue;
    for (const dependency of dependencies) {
      if (frameworkPattern.test(dependency)) unresolved.push(`${file}: ${dependency}`);
      if (/^\/(?:opt|usr\/local)\//.test(dependency)) external.push(`${file}: ${dependency}`);
    }
  }
  if (unresolved.length) {
    throw new Error(`Non-relocatable R.framework references remain:\n${unresolved.join('\n')}`);
  }
  if (external.length) {
    throw new Error(`Unbundled third-party library references remain:\n${external.join('\n')}`);
  }
  console.log(`  Inspected ${machoCount} Mach-O files; relocated ${changedCount}.`);
}

function writeLaunchers() {
  const launcher = `#!/bin/sh
R_HOME="$(CDPATH= cd -- "$(dirname -- "$0")/.." && pwd)"
export R_HOME
export R_SHARE_DIR="$R_HOME/share"
export R_INCLUDE_DIR="$R_HOME/include"
export R_DOC_DIR="$R_HOME/doc"
export R_LIBS="$R_HOME/library"
export R_LIBS_USER="$R_HOME/library"
exec "$R_HOME/bin/exec/R" --no-echo --no-restore "$@"
`;
  const rscript = `#!/bin/sh
R_HOME="$(CDPATH= cd -- "$(dirname -- "$0")/.." && pwd)"
export R_HOME
export R_SHARE_DIR="$R_HOME/share"
export R_INCLUDE_DIR="$R_HOME/include"
export R_DOC_DIR="$R_HOME/doc"
export R_LIBS="$R_HOME/library"
export R_LIBS_USER="$R_HOME/library"
if [ "$#" -gt 0 ] && [ "\${1#-}" = "$1" ]; then
  script="$1"
  shift
  exec "$R_HOME/bin/exec/R" --no-echo --no-restore --file="$script" --args "$@"
fi
exec "$R_HOME/bin/exec/R" --no-echo --no-restore "$@"
`;
  fs.writeFileSync(path.join(portableDir, 'bin', 'R'), launcher, { mode: 0o755 });
  fs.writeFileSync(path.join(portableDir, 'bin', 'Rscript'), rscript, { mode: 0o755 });
}

function isolatedEnvironment() {
  const library = path.join(portableDir, 'library');
  return {
    PATH: '/usr/bin:/bin:/usr/sbin:/sbin',
    HOME: os.tmpdir(),
    LANG: process.env.LANG || 'en_US.UTF-8',
    LC_ALL: process.env.LC_ALL || '',
    R_HOME: portableDir,
    R_LIBS: library,
    R_LIBS_USER: library,
    R_DEFAULT_PACKAGES: 'datasets,utils,grDevices,graphics,stats,methods',
  };
}

function verifyPortableRuntime() {
  section('Verifying the relocated runtime');
  const rscript = path.join(portableDir, 'bin', 'Rscript');
  const expression = [
    'stopifnot(normalizePath(R.home()) == normalizePath(Sys.getenv("R_HOME")))',
    'library(plumber)',
    'library(jsonlite)',
    'cat(R.version.string, "\\n")',
    'cat("R_HOME=", R.home(), "\\n", sep="")',
    'cat("PORTABLE_R_OK\\n")',
  ].join('; ');
  const output = run(rscript, ['-e', expression], { env: isolatedEnvironment() });
  if (!output.includes('PORTABLE_R_OK')) throw new Error('Portable R verification marker missing.');
  console.log(output);

  const apiExpression = [
    `setwd(${JSON.stringify(apiDir)})`,
    'source("simulation.R.dtd-backup")',
    'source("ddm_inspector.R")',
    'source("helpers.R")',
    'source("templates.R")',
    'source("self_tests.R")',
    'stopifnot(length(DDM.self_test_catalog()) == 6)',
    'cat("PORTABLE_API_OK\\n")',
  ].join('; ');
  const apiOutput = run(rscript, ['-e', apiExpression], { env: isolatedEnvironment() });
  if (!apiOutput.includes('PORTABLE_API_OK')) throw new Error('Portable API source verification failed.');
  console.log(apiOutput);
}

function main() {
  if (process.platform !== 'darwin') {
    throw new Error('setup-r-portable.js builds the macOS runtime. Use setup-r-portable-win.js on Windows.');
  }
  console.log('DDM-UI — relocatable R runtime setup');

  const missing = requiredPackages.filter((pkg) => {
    const result = spawnSync('Rscript', ['-e', `quit(status=if(requireNamespace('${pkg}', quietly=TRUE)) 0 else 1)`]);
    return result.status !== 0;
  });
  if (missing.length) {
    section(`Installing missing build packages: ${missing.join(', ')}`);
    rEval(`install.packages(c(${missing.map((pkg) => `'${pkg}'`).join(',')}), repos='https://cloud.r-project.org')`);
  }

  const rHomeReported = run('R', ['RHOME']);
  const rHome = fs.realpathSync(rHomeReported);
  console.log(`System R: ${rHome}`);

  const direct = requiredPackages.map((pkg) => `'${pkg}'`).join(',');
  const packageLines = rEval([
    `direct <- c(${direct})`,
    'deps <- tools::package_dependencies(direct, db=installed.packages(), recursive=TRUE)',
    `packages <- unique(c(${basePackages.map((pkg) => `'${pkg}'`).join(',')}, direct, unlist(deps)))`,
    'packages <- packages[nzchar(vapply(packages, function(package) system.file(package=package), character(1)))]',
    'for (package in packages) cat(package, "\\t", normalizePath(system.file(package=package)), "\\n", sep="")',
  ].join('; '));
  const packageLocations = packageLines.split(/\r?\n/).filter(Boolean).map((line) => {
    const tab = line.indexOf('\t');
    if (tab < 1) throw new Error(`Could not parse package location: ${line}`);
    return { name: line.slice(0, tab), source: line.slice(tab + 1) };
  });
  console.log(`Bundled R packages: ${packageLocations.length}`);

  section('Copying minimal R home');
  fs.rmSync(portableDir, { recursive: true, force: true });
  const sourceLibrary = path.join(rHome, 'library');
  fs.cpSync(rHome, portableDir, {
    recursive: true,
    dereference: true,
    preserveTimestamps: true,
    filter: (source) => source !== sourceLibrary && !source.startsWith(`${sourceLibrary}${path.sep}`),
  });
  const portableLibrary = path.join(portableDir, 'library');
  fs.mkdirSync(portableLibrary, { recursive: true });
  for (const pkg of packageLocations) {
    fs.cpSync(pkg.source, path.join(portableLibrary, pkg.name), {
      recursive: true,
      dereference: true,
      preserveTimestamps: true,
    });
  }

  // X11/TclTk are not used by the headless API and would introduce external
  // /opt dependencies. grDevices itself remains available for base R startup.
  fs.rmSync(path.join(portableLibrary, 'grDevices', 'libs', 'cairo.so'), { force: true });
  fs.rmSync(path.join(portableDir, 'modules', 'R_X11.so'), { force: true });
  fs.rmSync(path.join(portableDir, 'modules', 'R_de.so'), { force: true });

  writeLaunchers();
  patchMachOFiles();

  const packageNames = packageLocations.map((pkg) => `'${pkg.name}'`).join(',');
  const notices = rEval([
    `packages <- c(${packageNames})`,
    'for (package in packages) { d <- packageDescription(package); cat(package, " ", d$Version, " | ", d$License, "\\n", sep="") }',
  ].join('; '));
  fs.writeFileSync(
    path.join(portableDir, 'THIRD_PARTY_NOTICES.txt'),
    `DDM-UI bundled R runtime\n\nR is distributed under GPL-2 | GPL-3. See COPYING in this directory.\n\nBundled packages:\n${notices}\n`,
  );

  verifyPortableRuntime();
  const size = execFileSync('du', ['-sh', portableDir], { encoding: 'utf8' }).trim().split(/\s+/)[0];
  console.log(`\nPortable R runtime ready: ${portableDir} (${size})`);
}

try {
  main();
} catch (error) {
  console.error(`\nPortable R setup failed: ${error.message}`);
  process.exit(1);
}
