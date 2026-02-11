// ══════════════════════════════════════════════════════════════════════════════
// DDM-UI — API Client with retry logic, timeout, and connection monitoring
// ══════════════════════════════════════════════════════════════════════════════

// In Electron production: file:// protocol can't use relative /api paths,
// so we connect directly to the R Plumber server.
// In dev (Vite): the proxy handles /api → localhost:8000
const isElectron = typeof window !== 'undefined' && window.location.protocol === 'file:';
const API_BASE = isElectron ? 'http://127.0.0.1:8000/api' : '/api';

const DEFAULT_TIMEOUT_MS = 10000;
const SIMULATION_TIMEOUT_MS = 300000; // 5 minutes for large simulations
const MAX_RETRIES = 3;
const RETRY_DELAY_MS = 1000;

// ── Connection state ────────────────────────────────────────────────────────
type ConnectionListener = (connected: boolean) => void;
let _connected = true;
const _listeners = new Set<ConnectionListener>();

export function onConnectionChange(fn: ConnectionListener): () => void {
  _listeners.add(fn);
  return () => _listeners.delete(fn);
}

function setConnected(val: boolean) {
  if (val !== _connected) {
    _connected = val;
    _listeners.forEach(fn => fn(val));
  }
}

export function isConnected(): boolean {
  return _connected;
}

// ── Fetch with timeout ──────────────────────────────────────────────────────
async function fetchWithTimeout(url: string, options: RequestInit = {}, timeoutMs: number): Promise<Response> {
  const controller = new AbortController();
  const timer = setTimeout(() => controller.abort(), timeoutMs);

  try {
    const res = await fetch(url, { ...options, signal: controller.signal });
    clearTimeout(timer);
    return res;
  } catch (err: any) {
    clearTimeout(timer);
    if (err.name === 'AbortError') {
      throw new Error(`Request timed out after ${Math.round(timeoutMs / 1000)}s`);
    }
    throw err;
  }
}

// ── Retry wrapper ───────────────────────────────────────────────────────────
async function fetchWithRetry<T>(
  url: string,
  options: RequestInit = {},
  { timeout = DEFAULT_TIMEOUT_MS, retries = MAX_RETRIES } = {}
): Promise<T> {
  let lastError: Error | null = null;

  for (let attempt = 0; attempt <= retries; attempt++) {
    try {
      const res = await fetchWithTimeout(`${API_BASE}${url}`, {
        headers: { 'Content-Type': 'application/json' },
        ...options,
      }, timeout);

      setConnected(true);

      if (!res.ok) {
        const body = await res.text().catch(() => '');
        throw new Error(`API error ${res.status}: ${body || res.statusText}`);
      }

      return await res.json();
    } catch (err: any) {
      lastError = err;

      // Never retry POST requests — simulations are not idempotent
      if (options.method === 'POST') {
        throw err;
      }

      // Don't retry on client errors (4xx)
      if (err.message.includes('API error 4')) {
        throw err;
      }

      if (attempt < retries) {
        await new Promise(r => setTimeout(r, RETRY_DELAY_MS * (attempt + 1)));
      }
    }
  }

  setConnected(false);
  throw lastError || new Error('Request failed after retries');
}

// ── Public API ──────────────────────────────────────────────────────────────

export async function getHealth() {
  return fetchWithRetry<{ status: string }>('/health', {}, { timeout: 5000, retries: 0 });
}

export async function getTemplates() {
  return fetchWithRetry<Record<string, {
    id: string;
    name: string;
    description: string;
    category: string;
    available: boolean;
    phases: string;
    npeCount: number;
    connectionCount: number;
  }>>('/templates');
}

export async function getTemplate(id: string) {
  return fetchWithRetry<any>(`/templates/${id}`);
}

export async function validateNetwork(npes: any[], connections: any[]) {
  return fetchWithRetry<{ valid: boolean; error?: string; message?: string }>('/validate', {
    method: 'POST',
    body: JSON.stringify({ npes, connections }),
  });
}

export async function runSimulation(params: {
  npes: any;
  connections: any;
  trials: Record<string, string[]>;
  contingencies: string[];
  hasITI: boolean[];
  numNetworks: number;
  threshold?: string;
  disc?: number;
  pupdate?: string;
}) {
  return fetchWithRetry<any>('/simulate', {
    method: 'POST',
    body: JSON.stringify(params),
  }, { timeout: SIMULATION_TIMEOUT_MS, retries: 1 });
}

// Run a single network simulation (for real-time progress tracking)
export async function runSimulationOne(params: {
  npes: any;
  connections: any;
  trials: Record<string, string[]>;
  contingencies: string[];
  hasITI: boolean[];
  threshold?: string;
  disc?: number;
  pupdate?: string;
}) {
  return fetchWithRetry<{ success: boolean; result?: any; error?: string }>('/simulate-one', {
    method: 'POST',
    body: JSON.stringify(params),
  }, { timeout: SIMULATION_TIMEOUT_MS, retries: 0 });
}
