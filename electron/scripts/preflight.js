#!/usr/bin/env node

const { spawnSync } = require('child_process');
const fs = require('fs');
const os = require('os');
const path = require('path');

const target = (process.argv[2] || '').toLowerCase();
if (!['mac', 'win'].includes(target)) {
  console.error('Usage: node scripts/preflight.js <mac|win>');
  process.exit(2);
}

const electronDir = path.resolve(__dirname, '..');
const root = path.resolve(electronDir, '..');
const apiDir = path.join(root, 'api');
const required = [
  path.join(electronDir, 'main.js'),
  path.join(root, 'frontend', 'dist', 'index.html'),
  path.join(apiDir, 'plumber.R'),
  path.join(apiDir, 'ddm_inspector.R'),
  path.join(apiDir, 'helpers.R'),
  path.join(apiDir, 'templates.R'),
  path.join(apiDir, 'self_tests.R'),
  path.join(apiDir, 'simulation.R.dtd-backup'),
];

const runtime = target === 'mac'
  ? path.join(electronDir, 'R-portable')
  : path.join(electronDir, 'R-portable-win');
const rscript = target === 'mac'
  ? path.join(runtime, 'bin', 'Rscript')
  : path.join(runtime, 'bin', 'Rscript.exe');

if (target === 'mac') {
  required.push(path.join(root, 'icon.icns'));
  required.push(path.join(runtime, 'THIRD_PARTY_NOTICES.txt'));
} else {
  required.push(path.join(root, 'icon.ico'));
}
required.push(rscript);

const missing = required.filter((file) => !fs.existsSync(file));
if (missing.length) {
  console.error(`Packaging preflight failed for ${target}:`);
  for (const file of missing) console.error(`  missing: ${path.relative(root, file)}`);
  console.error(target === 'mac'
    ? 'Run npm run setup:r and compile the frontend first.'
    : 'Run npm run setup:r:win in a supported Windows build environment first.');
  process.exit(1);
}

function fail(message) {
  console.error(`Packaging preflight failed for ${target}: ${message}`);
  process.exit(1);
}

if (target === 'mac') {
  try {
    fs.accessSync(rscript, fs.constants.X_OK);
  } catch {
    fail('the bundled Rscript launcher is not executable');
  }

  const executable = path.join(runtime, 'bin', 'exec', 'R');
  const architecture = spawnSync('file', [executable], { encoding: 'utf8' });
  if (architecture.status !== 0 || !architecture.stdout.includes('arm64')) {
    fail(`bundled R is not arm64: ${architecture.stdout || architecture.stderr}`);
  }

  const candidates = [];
  const visit = (directory) => {
    for (const entry of fs.readdirSync(directory, { withFileTypes: true })) {
      const full = path.join(directory, entry.name);
      if (entry.isDirectory()) visit(full);
      else if (entry.isFile()) {
        const extension = path.extname(full).toLowerCase();
        if (extension === '.so' || extension === '.dylib' || (fs.statSync(full).mode & 0o111)) {
          candidates.push(full);
        }
      }
    }
  };
  visit(runtime);

  const badLinks = [];
  for (const file of candidates) {
    const linked = spawnSync('otool', ['-L', file], { encoding: 'utf8' });
    if (linked.status !== 0) continue;
    for (const line of linked.stdout.split(/\r?\n/).slice(1)) {
      const dependency = line.trim().split(/\s+/)[0];
      if (/^\/Library\/Frameworks\/R\.framework\//.test(dependency) ||
          /^\/(?:opt|usr\/local)\//.test(dependency)) {
        badLinks.push(`${path.relative(runtime, file)} -> ${dependency}`);
      }
    }
  }
  if (badLinks.length) fail(`non-portable libraries remain:\n${badLinks.join('\n')}`);

  const library = path.join(runtime, 'library');
  const env = {
    PATH: '/usr/bin:/bin:/usr/sbin:/sbin',
    HOME: os.tmpdir(),
    LANG: process.env.LANG || 'en_US.UTF-8',
    R_HOME: runtime,
    R_LIBS: library,
    R_LIBS_USER: library,
  };
  const expression = [
    'stopifnot(normalizePath(R.home()) == normalizePath(Sys.getenv("R_HOME")))',
    'library(plumber)',
    'library(jsonlite)',
    `setwd(${JSON.stringify(apiDir)})`,
    'source("simulation.R.dtd-backup")',
    'source("ddm_inspector.R")',
    'source("helpers.R")',
    'source("templates.R")',
    'source("self_tests.R")',
    'result <- DDM.run_self_tests()',
    'stopifnot(isTRUE(result$success), result$summary$passed == 6)',
    'cat("PREFLIGHT_SELF_TESTS=6/6\\n")',
  ].join('; ');
  const check = spawnSync(rscript, ['-e', expression], {
    cwd: apiDir,
    env,
    encoding: 'utf8',
    timeout: 300000,
    maxBuffer: 32 * 1024 * 1024,
  });
  if (check.status !== 0 || !check.stdout.includes('PREFLIGHT_SELF_TESTS=6/6')) {
    fail(`portable R self-test failed:\n${check.stdout || ''}\n${check.stderr || ''}`);
  }
  console.log(check.stdout.trim());
}

console.log(`Packaging preflight passed for ${target}.`);
