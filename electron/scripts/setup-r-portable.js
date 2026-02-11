#!/usr/bin/env node
// ══════════════════════════════════════════════════════════════════════════════
// DDM-UI — Setup R Portable
// Creates a self-contained R installation with all required packages
// for bundling inside the Electron app.
//
// Usage: node scripts/setup-r-portable.js
//
// This copies your system R into electron/R-portable/ along with
// only the packages needed by the DDM-UI API.
// ══════════════════════════════════════════════════════════════════════════════

const { execSync, spawnSync } = require('child_process');
const fs = require('fs');
const path = require('path');

const R_PORTABLE_DIR = path.join(__dirname, '..', 'R-portable');
const REQUIRED_PACKAGES = [
  'plumber', 'jsonlite', 'dplyr', 'igraph', 'tidygraph', 'ggraph',
  'readr', 'stringr', 'tidyr', 'tools', 'visNetwork',
  // Dependencies that plumber needs
  'httpuv', 'webutils', 'swagger', 'crayon', 'promises', 'later',
  'Rcpp', 'R6', 'magrittr', 'rlang', 'cli', 'glue', 'lifecycle',
  'vctrs', 'pillar', 'tibble', 'tidyselect', 'generics', 'fansi',
  'utf8', 'pkgconfig', 'withr',
];

function run(cmd) {
  console.log(`  $ ${cmd}`);
  return execSync(cmd, { encoding: 'utf-8', stdio: 'pipe' }).trim();
}

function main() {
  console.log('═══════════════════════════════════════════════════════');
  console.log('  DDM-UI — R Portable Setup');
  console.log('═══════════════════════════════════════════════════════\n');

  // 1. Find system R
  const platform = process.platform;
  let rHome;
  try {
    rHome = run('R RHOME');
  } catch {
    console.error('❌ R is not installed or not in PATH');
    process.exit(1);
  }
  console.log(`✅ Found R at: ${rHome}`);

  // 2. Check required packages are installed
  console.log('\nChecking required packages...');
  const missingPkgs = [];
  for (const pkg of REQUIRED_PACKAGES) {
    try {
      run(`Rscript -e "library(${pkg})"`);
    } catch {
      missingPkgs.push(pkg);
    }
  }
  if (missingPkgs.length > 0) {
    console.log(`\n⚠️  Missing packages: ${missingPkgs.join(', ')}`);
    console.log('Installing missing packages...');
    const installCmd = `Rscript -e "install.packages(c(${missingPkgs.map(p => `'${p}'`).join(', ')}), repos='https://cran.r-project.org')"`;
    try {
      execSync(installCmd, { stdio: 'inherit' });
    } catch (err) {
      console.error('❌ Failed to install packages');
      process.exit(1);
    }
  }
  console.log('✅ All required packages available');

  // 3. Create R-portable directory
  if (fs.existsSync(R_PORTABLE_DIR)) {
    console.log(`\nRemoving existing R-portable...`);
    fs.rmSync(R_PORTABLE_DIR, { recursive: true });
  }
  fs.mkdirSync(R_PORTABLE_DIR, { recursive: true });

  // 4. Copy R installation
  console.log(`\nCopying R from ${rHome} to R-portable...`);

  if (platform === 'darwin') {
    // macOS: Copy R framework
    run(`cp -R "${rHome}/." "${R_PORTABLE_DIR}/"`);
  } else if (platform === 'win32') {
    // Windows: Copy R directory
    run(`xcopy "${rHome}" "${R_PORTABLE_DIR}" /E /I /Q`);
  } else {
    // Linux
    run(`cp -R "${rHome}/." "${R_PORTABLE_DIR}/"`);
  }

  // 5. Copy library packages
  console.log('\nCopying R packages...');
  const libPaths = run('Rscript -e "cat(.libPaths(), sep=\\"\\n\\")"').split('\n');
  console.log(`  R library paths: ${libPaths.join(', ')}`);

  const portableLib = path.join(R_PORTABLE_DIR, 'library');
  if (!fs.existsSync(portableLib)) {
    fs.mkdirSync(portableLib, { recursive: true });
  }

  // Get all required packages INCLUDING their dependencies
  const pkgList = REQUIRED_PACKAGES.map(p => `'${p}'`).join(',');
  const allDeps = run(
    `Rscript -e "options(repos='https://cran.r-project.org'); pkgs <- c(${pkgList}); deps <- tools::package_dependencies(pkgs, recursive=TRUE); all_pkgs <- unique(c(pkgs, unlist(deps))); cat(all_pkgs, sep='\\n')"`
  ).split('\n').filter(Boolean);

  console.log(`  Total packages to copy (including deps): ${allDeps.length}`);

  let copied = 0;
  for (const pkg of allDeps) {
    for (const libPath of libPaths) {
      const pkgDir = path.join(libPath, pkg);
      const destDir = path.join(portableLib, pkg);
      if (fs.existsSync(pkgDir) && !fs.existsSync(destDir)) {
        if (platform === 'win32') {
          spawnSync('xcopy', [pkgDir, destDir, '/E', '/I', '/Q'], { stdio: 'pipe' });
        } else {
          spawnSync('cp', ['-R', pkgDir, destDir], { stdio: 'pipe' });
        }
        copied++;
        break;
      }
    }
  }
  console.log(`  Copied ${copied} packages`);

  // 6. Verify
  console.log('\nVerifying R-portable...');
  const rscript = platform === 'win32'
    ? path.join(R_PORTABLE_DIR, 'bin', 'Rscript.exe')
    : path.join(R_PORTABLE_DIR, 'bin', 'Rscript');

  if (!fs.existsSync(rscript)) {
    console.error(`❌ Rscript not found at ${rscript}`);
    // List what's in bin/
    const binDir = path.join(R_PORTABLE_DIR, 'bin');
    if (fs.existsSync(binDir)) {
      console.log('  Contents of bin/:', fs.readdirSync(binDir).join(', '));
    }
    process.exit(1);
  }

  try {
    const env = { ...process.env, R_LIBS_USER: portableLib, R_LIBS: portableLib };
    const result = spawnSync(rscript, ['-e', 'library(plumber); cat("OK")'], {
      env, encoding: 'utf-8', stdio: 'pipe'
    });
    if (result.stdout.includes('OK')) {
      console.log('✅ R-portable verified — plumber loads correctly');
    } else {
      console.log('⚠️  Plumber load test output:', result.stdout, result.stderr);
    }
  } catch (err) {
    console.error('⚠️  Verification failed:', err.message);
  }

  // 7. Show size
  const size = run(`du -sh "${R_PORTABLE_DIR}"`).split('\t')[0];
  console.log(`\n📦 R-portable size: ${size}`);
  console.log('\n✅ R-portable setup complete!');
}

main();
