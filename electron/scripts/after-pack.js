#!/usr/bin/env node
const { execFileSync } = require('child_process');
const fs = require('fs');
const path = require('path');

function log(message) {
  console.log(`[after-pack] ${message}`);
}

function run(command, args, options = {}) {
  execFileSync(command, args, {
    stdio: options.stdio || 'pipe',
    ...options,
  });
}

function chmodIfExists(filePath) {
  if (!fs.existsSync(filePath)) return;
  const stat = fs.statSync(filePath);
  fs.chmodSync(filePath, stat.mode | 0o755);
}

function chmodFilesIn(dirPath) {
  if (!fs.existsSync(dirPath)) return;

  for (const entry of fs.readdirSync(dirPath, { withFileTypes: true })) {
    const entryPath = path.join(dirPath, entry.name);
    if (entry.isDirectory()) {
      chmodFilesIn(entryPath);
    } else if (entry.isFile()) {
      chmodIfExists(entryPath);
    }
  }
}

function stripExtendedAttributes(appPath) {
  try {
    run('/usr/bin/xattr', ['-cr', appPath]);
    log('Removed extended attributes from app bundle');
  } catch {
    log('xattr cleanup skipped');
  }
}

function signAdHoc(appPath) {
  run('/usr/bin/codesign', ['--force', '--deep', '--sign', '-', appPath], {
    stdio: 'inherit',
  });
  run('/usr/bin/codesign', ['--verify', '--deep', '--strict', '--verbose=2', appPath], {
    stdio: 'inherit',
  });
  log('Applied ad-hoc macOS code signature');
}

module.exports = async function afterPack(context) {
  if (context.electronPlatformName !== 'darwin' || process.platform !== 'darwin') {
    return;
  }

  const macConfig = context.packager.config.mac || {};
  if (macConfig.identity !== null) {
    log('Skipping ad-hoc signing because a signing identity may be used');
    return;
  }

  const appName = `${context.packager.appInfo.productFilename}.app`;
  const appPath = path.join(context.appOutDir, appName);
  if (!fs.existsSync(appPath)) {
    throw new Error(`Expected app bundle not found: ${appPath}`);
  }

  const resourcesPath = path.join(appPath, 'Contents', 'Resources');
  chmodFilesIn(path.join(resourcesPath, 'R-portable', 'bin'));
  chmodIfExists(path.join(resourcesPath, 'R-portable', 'R'));
  chmodIfExists(path.join(resourcesPath, 'R-portable', 'Rscript'));
  chmodFilesIn(path.join(appPath, 'Contents', 'MacOS'));

  stripExtendedAttributes(appPath);
  signAdHoc(appPath);
};
