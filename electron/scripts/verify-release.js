#!/usr/bin/env node

const { execFileSync, spawnSync } = require('child_process');
const fs = require('fs');
const os = require('os');
const path = require('path');

const electronDir = path.resolve(__dirname, '..');
const releaseDir = path.resolve(electronDir, process.argv[2] || 'release');

function run(file, args, options = {}) {
  const result = spawnSync(file, args, {
    encoding: 'utf8',
    maxBuffer: 32 * 1024 * 1024,
    ...options,
  });
  if (result.status !== 0) {
    throw new Error(`${file} failed:\n${result.stdout || ''}\n${result.stderr || ''}`);
  }
  return (result.stdout || '').trim();
}

function findFirst(directory, predicate) {
  if (!fs.existsSync(directory)) return null;
  for (const entry of fs.readdirSync(directory, { withFileTypes: true })) {
    const full = path.join(directory, entry.name);
    if (predicate(full, entry)) return full;
    if (entry.isDirectory()) {
      const nested = findFirst(full, predicate);
      if (nested) return nested;
    }
  }
  return null;
}

function verifyApp(appPath, label) {
  console.log(`Verifying ${label}: ${appPath}`);
  run('codesign', ['--verify', '--deep', '--strict', '--verbose=2', appPath]);
  const resources = path.join(appPath, 'Contents', 'Resources');
  const runtime = path.join(resources, 'R-portable');
  const rscript = path.join(runtime, 'bin', 'Rscript');
  const api = path.join(resources, 'api');
  const library = path.join(runtime, 'library');
  const expression = [
    'stopifnot(normalizePath(R.home()) == normalizePath(Sys.getenv("R_HOME")))',
    'library(plumber)',
    `setwd(${JSON.stringify(api)})`,
    'source("simulation.R.dtd-backup")',
    'source("ddm_inspector.R")',
    'source("helpers.R")',
    'source("templates.R")',
    'source("self_tests.R")',
    'stopifnot(length(DDM.self_test_catalog()) == 6)',
    'cat("PACKAGED_RUNTIME_OK\\n")',
  ].join('; ');
  const output = run(rscript, ['-e', expression], {
    cwd: api,
    env: {
      PATH: '/usr/bin:/bin:/usr/sbin:/sbin',
      HOME: os.tmpdir(),
      LANG: process.env.LANG || 'en_US.UTF-8',
      R_HOME: runtime,
      R_LIBS: library,
      R_LIBS_USER: library,
    },
    timeout: 120000,
  });
  if (!output.includes('PACKAGED_RUNTIME_OK')) throw new Error(`${label} runtime marker missing.`);
  console.log(output);
}

try {
  const dmg = findFirst(releaseDir, (full, entry) => entry.isFile() && full.endsWith('.dmg'));
  const unpackedApp = findFirst(releaseDir, (full, entry) => entry.isDirectory() && full.endsWith('.app'));
  if (!dmg) throw new Error('No DMG was produced.');
  if (!unpackedApp) throw new Error('No unpacked .app was produced.');

  run('hdiutil', ['verify', dmg], { timeout: 300000 });
  verifyApp(unpackedApp, 'unpacked application');

  const mountPoint = fs.mkdtempSync(path.join(os.tmpdir(), 'ddm-ui-dmg-'));
  try {
    run('hdiutil', ['attach', '-nobrowse', '-readonly', '-mountpoint', mountPoint, dmg], { timeout: 300000 });
    const mountedApp = findFirst(mountPoint, (full, entry) => entry.isDirectory() && full.endsWith('.app'));
    if (!mountedApp) throw new Error('The mounted DMG does not contain an application.');
    verifyApp(mountedApp, 'application mounted from DMG');
  } finally {
    try { execFileSync('hdiutil', ['detach', mountPoint], { stdio: 'ignore' }); } catch { /* best effort */ }
    fs.rmSync(mountPoint, { recursive: true, force: true });
  }

  const size = fs.statSync(dmg).size / (1024 * 1024);
  console.log(`Release verification passed: ${dmg} (${size.toFixed(1)} MB)`);
} catch (error) {
  console.error(`Release verification failed: ${error.message}`);
  process.exit(1);
}
