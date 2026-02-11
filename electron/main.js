// ══════════════════════════════════════════════════════════════════════════════
// DDM-UI — Electron Main Process
// Launches R Plumber API + serves React frontend in a desktop window
// ══════════════════════════════════════════════════════════════════════════════

const { app, BrowserWindow, dialog } = require('electron');
const path = require('path');
const { spawn, execSync } = require('child_process');
const http = require('http');
const fs = require('fs');

let mainWindow = null;
let rProcess = null;
let rProcessPid = null; // Keep PID separately so SIGKILL fallback works
const R_PORT = 8000;
const isDev = !app.isPackaged;

// ── Paths ────────────────────────────────────────────────────────────────────
function getResourcePath(...segments) {
  if (isDev) {
    return path.join(__dirname, '..', ...segments);
  }
  return path.join(process.resourcesPath, ...segments);
}

function getRPath() {
  if (isDev) {
    // In dev mode, use system R
    return process.platform === 'win32' ? 'Rscript.exe' : 'Rscript';
  }

  // In production, use bundled R-portable
  const rPortable = path.join(process.resourcesPath, 'R-portable');

  if (process.platform === 'win32') {
    return path.join(rPortable, 'bin', 'Rscript.exe');
  } else {
    // macOS
    return path.join(rPortable, 'bin', 'Rscript');
  }
}

function getApiPath() {
  if (isDev) {
    return path.join(__dirname, '..', 'api');
  }
  return path.join(process.resourcesPath, 'api');
}

function getFrontendPath() {
  if (isDev) {
    return null; // Use Vite dev server
  }
  return path.join(process.resourcesPath, 'frontend');
}

// ── Kill orphan R processes on our port ──────────────────────────────────────
function killOrphanR() {
  try {
    if (process.platform === 'win32') {
      // Find and kill any process listening on our port
      const result = execSync(`netstat -ano | findstr :${R_PORT} | findstr LISTENING`, { encoding: 'utf8', timeout: 5000 });
      const lines = result.trim().split('\n');
      for (const line of lines) {
        const pid = line.trim().split(/\s+/).pop();
        if (pid && pid !== '0') {
          console.log(`[DDM-UI] Killing orphan process on port ${R_PORT} (PID ${pid})`);
          try { execSync(`taskkill /pid ${pid} /f /t`, { timeout: 5000 }); } catch { /* ignore */ }
        }
      }
    } else {
      // macOS/Linux: find process on port and kill it
      const result = execSync(`lsof -ti :${R_PORT}`, { encoding: 'utf8', timeout: 5000 });
      const pids = result.trim().split('\n').filter(Boolean);
      for (const pid of pids) {
        console.log(`[DDM-UI] Killing orphan process on port ${R_PORT} (PID ${pid})`);
        try { execSync(`kill -9 ${pid}`, { timeout: 5000 }); } catch { /* ignore */ }
      }
    }
  } catch {
    // No process found on the port — good, nothing to clean up
  }
}

// ── Start R Plumber API ──────────────────────────────────────────────────────
function startR() {
  return new Promise((resolve, reject) => {
    // First, clean up any orphan R processes from a previous bad shutdown
    killOrphanR();

    const rscript = getRPath();
    const apiDir = getApiPath();
    const plumberFile = path.join(apiDir, 'plumber.R');

    console.log(`[DDM-UI] Starting R Plumber...`);
    console.log(`[DDM-UI]   Rscript: ${rscript}`);
    console.log(`[DDM-UI]   API dir: ${apiDir}`);
    console.log(`[DDM-UI]   Plumber: ${plumberFile}`);

    if (!fs.existsSync(plumberFile)) {
      reject(new Error(`plumber.R not found at ${plumberFile}`));
      return;
    }

    // Set R_LIBS_USER to the bundled packages directory
    const env = { ...process.env };
    if (!isDev) {
      const rLibs = path.join(process.resourcesPath, 'R-portable', 'library');
      env.R_LIBS_USER = rLibs;
      env.R_LIBS = rLibs;
    }

    // Suppress macOS Java dialog that appears when R loads packages referencing rJava.
    // NOAWT prevents AWT/Swing from initializing (which triggers the "install Java" popup).
    // JAVA_TOOL_OPTIONS silences the JVM startup banner.
    if (process.platform === 'darwin') {
      env.NOAWT = '1';
      env.JAVA_TOOL_OPTIONS = '-Djava.awt.headless=true';
      // Tell R not to look for Java at all
      env.R_JAVA_LD_LIBRARY_PATH = '';
    }

    const rCmd = `plumber::plumb('${plumberFile.replace(/\\/g, '/')}')$run(host='127.0.0.1', port=${R_PORT})`;

    // Use detached + process group so we can kill all child processes together
    rProcess = spawn(rscript, ['-e', rCmd], {
      cwd: apiDir,
      env,
      stdio: ['ignore', 'pipe', 'pipe'],
      detached: process.platform !== 'win32', // Create process group on macOS/Linux
    });

    rProcessPid = rProcess.pid;
    console.log(`[DDM-UI] R process started with PID ${rProcessPid}`);

    // Guard: only resolve/reject once
    let settled = false;
    let poll = null;

    const finish = (fn, arg) => {
      if (settled) return;
      settled = true;
      if (poll) clearInterval(poll);
      fn(arg);
    };

    rProcess.stdout.on('data', (data) => {
      console.log(`[R] ${data.toString().trim()}`);
    });

    rProcess.stderr.on('data', (data) => {
      const msg = data.toString().trim();
      console.log(`[R] ${msg}`);
      // Plumber logs "Running plumber API" or "Starting server" to stderr
      if (msg.includes('Running') || msg.includes('Starting') || msg.includes('port')) {
        finish(resolve);
      }
    });

    rProcess.on('error', (err) => {
      console.error(`[DDM-UI] Failed to start R:`, err);
      finish(reject, err);
    });

    rProcess.on('exit', (code) => {
      console.log(`[DDM-UI] R process exited with code ${code}`);
      rProcess = null;
      rProcessPid = null;
    });

    // Fallback: poll the health endpoint
    let attempts = 0;
    const maxAttempts = 30; // 15 seconds
    poll = setInterval(() => {
      attempts++;
      http.get(`http://127.0.0.1:${R_PORT}/api/health`, (res) => {
        if (res.statusCode === 200) {
          console.log(`[DDM-UI] R API is ready (attempt ${attempts})`);
          finish(resolve);
        }
      }).on('error', () => {
        if (attempts >= maxAttempts) {
          finish(reject, new Error('R API failed to start within 15 seconds'));
        }
      });
    }, 500);
  });
}

// ── Stop R ───────────────────────────────────────────────────────────────────
function stopR() {
  const pid = rProcessPid;
  if (!pid && !rProcess) return;

  console.log(`[DDM-UI] Stopping R process (PID ${pid})...`);

  if (process.platform === 'win32') {
    // Windows: taskkill with /t kills the entire process tree
    if (pid) {
      try { execSync(`taskkill /pid ${pid} /f /t`, { timeout: 5000 }); } catch { /* ignore */ }
    }
  } else {
    // macOS/Linux: kill the entire process group (negative PID)
    // This kills R + all child processes (httpuv, etc.) at once
    if (pid) {
      try {
        // Kill process group: negative PID kills all processes in the group
        process.kill(-pid, 'SIGTERM');
      } catch {
        // If process group kill fails, try killing just the main process
        try { process.kill(pid, 'SIGTERM'); } catch { /* ignore */ }
      }

      // Force kill after 3 seconds if still alive
      setTimeout(() => {
        try {
          process.kill(-pid, 0); // Check if still alive (signal 0 = test)
          console.log(`[DDM-UI] R still alive after SIGTERM, sending SIGKILL...`);
          try { process.kill(-pid, 'SIGKILL'); } catch { /* ignore */ }
        } catch {
          // Process is already dead — good
        }
      }, 3000);
    }
  }

  rProcess = null;
  rProcessPid = null;
}

// ── Create Window ────────────────────────────────────────────────────────────
function createWindow() {
  mainWindow = new BrowserWindow({
    width: 1400,
    height: 900,
    minWidth: 1024,
    minHeight: 700,
    title: 'DDM-UI',
    icon: path.join(__dirname, '..', 'icon.ico'),
    webPreferences: {
      nodeIntegration: false,
      contextIsolation: true,
    },
    show: false, // Show after content loads
    backgroundColor: '#0f172a', // Match app dark background
  });

  // Remove menu bar (optional, keeps it clean)
  mainWindow.setMenuBarVisibility(false);

  if (isDev) {
    // Dev mode: load from Vite dev server
    mainWindow.loadURL('http://localhost:5173');
    // mainWindow.webContents.openDevTools();
  } else {
    // Production: load built frontend, proxy API calls
    const frontendPath = getFrontendPath();
    const indexPath = path.join(frontendPath, 'index.html');
    mainWindow.loadFile(indexPath);
  }

  mainWindow.once('ready-to-show', () => {
    mainWindow.show();
  });

  mainWindow.on('closed', () => {
    mainWindow = null;
  });
}

// ── App Lifecycle ────────────────────────────────────────────────────────────
app.whenReady().then(async () => {
  try {
    // Show splash/loading while R starts
    createWindow();

    // Start R API
    await startR();
    console.log('[DDM-UI] R API started successfully');

    // In production, we need to intercept /api calls and proxy them to R
    if (!isDev) {
      const { session } = require('electron');
      session.defaultSession.webRequest.onBeforeRequest(
        { urls: ['file://*/api/*'] },
        (details, callback) => {
          const apiPath = details.url.replace(/^file:\/\/.*?\/api\//, '');
          callback({ redirectURL: `http://127.0.0.1:${R_PORT}/api/${apiPath}` });
        }
      );
    }
  } catch (err) {
    console.error('[DDM-UI] Startup error:', err);
    dialog.showErrorBox(
      'DDM-UI — Startup Error',
      `Failed to start the simulation engine.\n\n${err.message}\n\nPlease ensure R is installed or contact support.`
    );
    app.quit();
  }
});

app.on('window-all-closed', () => {
  stopR();
  app.quit();
});

app.on('before-quit', () => {
  stopR();
});

app.on('activate', () => {
  if (BrowserWindow.getAllWindows().length === 0) {
    createWindow();
  }
});
