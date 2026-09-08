#!/usr/bin/env node

const { execFileSync } = require('child_process');
const fs = require('fs');
const path = require('path');

module.exports = async function afterPack(context) {
  if (context.electronPlatformName !== 'darwin') return;
  const appName = `${context.packager.appInfo.productFilename}.app`;
  const appPath = path.join(context.appOutDir, appName);
  if (!fs.existsSync(appPath)) throw new Error(`Packed app not found: ${appPath}`);

  // No Developer ID is configured in this workspace. Ad-hoc signing keeps the
  // Apple-Silicon bundle internally consistent and executable for local use.
  // Public distribution still requires Developer ID signing and notarization.
  execFileSync('codesign', ['--force', '--deep', '--sign', '-', appPath], { stdio: 'inherit' });
  execFileSync('codesign', ['--verify', '--deep', '--strict', '--verbose=2', appPath], { stdio: 'inherit' });
  console.log(`Ad-hoc signed and verified: ${appPath}`);
};
