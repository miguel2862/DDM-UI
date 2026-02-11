#!/usr/bin/env node
// ==============================================================================
// DDM-UI — Setup R Portable for WINDOWS (cross-platform)
// Downloads R for Windows from CRAN, extracts it with innoextract,
// and downloads Windows binary packages from CRAN.
//
// Prerequisites: innoextract (brew install innoextract)
// Usage: node scripts/setup-r-portable-win.js
// ==============================================================================

const { execSync } = require('child_process');
const fs = require('fs');
const path = require('path');
const https = require('https');

const R_VERSION = '4.4.2';
const R_INSTALLER_URL = `https://cran.r-project.org/bin/windows/base/old/${R_VERSION}/R-${R_VERSION}-win.exe`;
const R_INSTALLER_FILE = path.join(__dirname, `R-${R_VERSION}-win.exe`);
const R_PORTABLE_WIN = path.join(__dirname, '..', 'R-portable-win');
const EXTRACT_DIR = path.join(__dirname, 'r-win-extracted');
const TMP_PKG_DIR = path.join(__dirname, 'win-pkg-tmp');

const REQUIRED_PACKAGES = [
  'plumber', 'jsonlite', 'dplyr', 'igraph', 'tidygraph', 'ggraph',
  'readr', 'stringr', 'tidyr', 'visNetwork',
];

function run(cmd) {
  return execSync(cmd, { encoding: 'utf-8', stdio: 'pipe', maxBuffer: 50 * 1024 * 1024 }).trim();
}

function download(url, dest) {
  return new Promise((resolve, reject) => {
    console.log(`  Downloading ${url}`);
    const file = fs.createWriteStream(dest);
    const request = (url) => {
      https.get(url, (res) => {
        if (res.statusCode === 301 || res.statusCode === 302) {
          request(res.headers.location);
          return;
        }
        const total = parseInt(res.headers['content-length'] || '0');
        let downloaded = 0;
        res.on('data', (chunk) => {
          downloaded += chunk.length;
          file.write(chunk);
          if (total > 0) {
            const pct = ((downloaded / total) * 100).toFixed(1);
            process.stdout.write(`\r  Progress: ${pct}% (${(downloaded / 1024 / 1024).toFixed(1)} MB)`);
          }
        });
        res.on('end', () => {
          file.end();
          console.log('');
          resolve();
        });
      }).on('error', reject);
    };
    request(url);
  });
}

async function main() {
  console.log('===============================================================');
  console.log('  DDM-UI — R Portable for Windows (Cross-Platform Setup)');
  console.log('===============================================================\n');

  // 0. Check innoextract
  try { run('which innoextract'); } catch {
    console.error('innoextract not found. Install: brew install innoextract');
    process.exit(1);
  }

  // 1. Download R for Windows
  if (fs.existsSync(R_INSTALLER_FILE)) {
    console.log(`R installer already cached: ${path.basename(R_INSTALLER_FILE)}`);
  } else {
    console.log(`Downloading R ${R_VERSION} for Windows...`);
    await download(R_INSTALLER_URL, R_INSTALLER_FILE);
  }

  // 2. Extract with innoextract
  if (fs.existsSync(EXTRACT_DIR)) fs.rmSync(EXTRACT_DIR, { recursive: true });
  fs.mkdirSync(EXTRACT_DIR, { recursive: true });
  console.log('\nExtracting R installer...');
  run(`innoextract "${R_INSTALLER_FILE}" -d "${EXTRACT_DIR}"`);

  const appDir = path.join(EXTRACT_DIR, 'app');
  if (!fs.existsSync(appDir)) {
    console.error('Extraction failed. Contents:', fs.readdirSync(EXTRACT_DIR).join(', '));
    process.exit(1);
  }
  console.log('  R extracted successfully');

  // 3. Copy to R-portable-win
  if (fs.existsSync(R_PORTABLE_WIN)) fs.rmSync(R_PORTABLE_WIN, { recursive: true });
  console.log(`\nCopying to R-portable-win...`);
  run(`cp -R "${appDir}" "${R_PORTABLE_WIN}"`);

  // 4. Download ALL Windows binary packages in one R call
  console.log('\nDownloading Windows binary packages from CRAN...');
  if (fs.existsSync(TMP_PKG_DIR)) fs.rmSync(TMP_PKG_DIR, { recursive: true });
  fs.mkdirSync(TMP_PKG_DIR, { recursive: true });

  const winLibDir = path.join(R_PORTABLE_WIN, 'library');
  const pkgList = REQUIRED_PACKAGES.map(p => `'${p}'`).join(',');

  // Use a single R script to download all packages + dependencies
  const rScript = `
    options(repos = 'https://cran.r-project.org')

    # Get available Windows binary packages for R ${R_VERSION}
    r_ver <- '${R_VERSION.split('.').slice(0, 2).join('.')}'
    contrib_url <- paste0('https://cran.r-project.org/bin/windows/contrib/', r_ver)

    # Required packages
    pkgs <- c(${pkgList})

    # Resolve all dependencies
    db <- available.packages(contriburl = contrib_url, type = 'win.binary')
    deps <- tools::package_dependencies(pkgs, db = db, recursive = TRUE)
    all_pkgs <- unique(c(pkgs, unlist(deps)))

    # Filter out base packages that ship with R
    base_pkgs <- installed.packages(priority = c('base', 'recommended'))[,'Package']
    to_download <- setdiff(all_pkgs, base_pkgs)

    cat('Packages to download:', length(to_download), '\\n')
    cat(paste(to_download, collapse = ', '), '\\n\\n')

    # Download all at once
    dest <- '${TMP_PKG_DIR.replace(/'/g, "\\'")}'
    result <- download.packages(to_download, destdir = dest, contriburl = contrib_url, type = 'win.binary')

    cat('\\nDownloaded', nrow(result), 'packages\\n')

    # Print file paths
    for (i in seq_len(nrow(result))) {
      cat('FILE:', result[i, 2], '\\n')
    }
  `.replace(/\n/g, '\n');

  // Write R script to temp file and run it
  const rScriptFile = path.join(__dirname, 'download-win-pkgs.R');
  fs.writeFileSync(rScriptFile, rScript);

  console.log('  Running R to download packages...');
  const dlOutput = execSync(`Rscript "${rScriptFile}"`, {
    encoding: 'utf-8',
    stdio: ['pipe', 'pipe', 'pipe'],
    maxBuffer: 50 * 1024 * 1024,
  });

  console.log(dlOutput);

  // Extract each downloaded .zip into the library
  const zipFiles = dlOutput.split('\n')
    .filter(line => line.startsWith('FILE:'))
    .map(line => line.replace('FILE:', '').trim());

  console.log(`\nExtracting ${zipFiles.length} packages to library...`);
  let extracted = 0;
  for (const zipFile of zipFiles) {
    if (fs.existsSync(zipFile)) {
      try {
        run(`unzip -o -q "${zipFile}" -d "${winLibDir}"`);
        extracted++;
      } catch (err) {
        console.log(`  Warning: failed to extract ${path.basename(zipFile)}`);
      }
    }
  }
  console.log(`  Extracted ${extracted} packages`);

  // 5. Verify plumber exists in library
  const plumberDir = path.join(winLibDir, 'plumber');
  if (fs.existsSync(plumberDir)) {
    console.log('  plumber package found in library');
  } else {
    console.log('  WARNING: plumber not found in library!');
  }

  // 6. Clean up temp files
  console.log('\nCleaning up...');
  if (fs.existsSync(TMP_PKG_DIR)) fs.rmSync(TMP_PKG_DIR, { recursive: true });
  if (fs.existsSync(EXTRACT_DIR)) fs.rmSync(EXTRACT_DIR, { recursive: true });
  if (fs.existsSync(rScriptFile)) fs.unlinkSync(rScriptFile);

  // 7. Show size
  const size = run(`du -sh "${R_PORTABLE_WIN}"`).split('\t')[0];
  console.log(`\n R-portable-win size: ${size}`);
  console.log(' Setup complete!\n');
}

main().catch(err => {
  console.error('Fatal error:', err);
  process.exit(1);
});
