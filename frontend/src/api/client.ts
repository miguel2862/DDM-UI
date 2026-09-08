// ══════════════════════════════════════════════════════════════════════════════
// DDM-UI — API Client with retry logic, timeout, and connection monitoring
// ══════════════════════════════════════════════════════════════════════════════

import type {
  Connection,
  NPE,
  SimulationMetadata,
  SimulationResult,
  TemplateData,
} from '../types/ddm';
import type { SimulationInspector } from '../types/inspector';

// In Electron production: file:// protocol can't use relative /api paths,
// so we connect directly to the R Plumber server.
// In dev (Vite): the proxy handles /api → localhost:8000
const isElectron = typeof window !== 'undefined' && window.location.protocol === 'file:';
const electronParams = isElectron ? new URLSearchParams(window.location.search) : null;
const electronPort = electronParams?.get('apiPort') || '8000';
const electronToken = electronParams?.get('apiToken') || '';
const API_BASE = isElectron ? `http://127.0.0.1:${electronPort}/api` : '/api';

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
  } catch (err: unknown) {
    clearTimeout(timer);
    if (err instanceof Error && err.name === 'AbortError') {
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
        headers: {
          'Content-Type': 'application/json',
          ...(electronToken ? { 'X-DDM-Token': electronToken } : {}),
        },
        ...options,
      }, timeout);

      setConnected(true);

      if (!res.ok) {
        const body = await res.text().catch(() => '');
        throw new Error(`API error ${res.status}: ${body || res.statusText}`);
      }

      return await res.json();
    } catch (err: unknown) {
      const error = err instanceof Error ? err : new Error(String(err));
      lastError = error;

      // Never retry POST requests — simulations are not idempotent
      if (options.method === 'POST') {
        throw error;
      }

      // Don't retry on client errors (4xx)
      if (error.message.includes('API error 4')) {
        throw error;
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
  return fetchWithRetry<Partial<TemplateData> & { error?: string }>(`/templates/${id}`);
}

export interface NPEColumnPayload {
  NPE: string[];
  Type: NPE['type'][];
  Layer: NPE['layer'][];
  Activation: number[];
  'Temporal.Summation': number[];
  'Activation.Decay': number[];
  mu: number[];
  sigma: number[];
  logisSigma: number[];
}

export interface ConnectionColumnPayload {
  PreSinapticNPE: string[];
  PostSinapticNPE: string[];
  Weight: number[];
  alpha: number[];
  beta: number[];
  alpha_prime: number[];
  beta_prime: number[];
}

/** Serialize the row-oriented editor state into the column-oriented R payload. */
export function serializeNetwork(npes: NPE[], connections: Connection[]): {
  npes: NPEColumnPayload;
  connections: ConnectionColumnPayload;
} {
  return {
    npes: {
      NPE: npes.map((n) => n.name),
      Type: npes.map((n) => n.type),
      Layer: npes.map((n) => n.layer),
      Activation: npes.map((n) => n.activation),
      'Temporal.Summation': npes.map((n) => n.temporalSummation),
      'Activation.Decay': npes.map((n) => n.activationDecay),
      mu: npes.map((n) => n.mu),
      sigma: npes.map((n) => n.sigma),
      logisSigma: npes.map((n) => n.logisSigma),
    },
    connections: {
      PreSinapticNPE: connections.map((c) => c.presynapticNPE),
      PostSinapticNPE: connections.map((c) => c.postsynapticNPE),
      Weight: connections.map((c) => c.weight),
      alpha: connections.map((c) => c.alpha),
      beta: connections.map((c) => c.beta),
      alpha_prime: connections.map((c) => c.alphaPrime),
      beta_prime: connections.map((c) => c.betaPrime),
    },
  };
}

export async function validateNetwork(
  npes: NPE[],
  connections: Connection[],
) {
  const network = serializeNetwork(npes, connections);
  return fetchWithRetry<{ valid: boolean; error?: string; message?: string }>('/validate', {
    method: 'POST',
    body: JSON.stringify({ ...network, model: 'dtd' }),
  });
}

interface SimulationParams {
  npes: unknown;
  connections: unknown;
  trials: Record<string, string[]>;
  contingencies: string[];
  hasITI: boolean[];
  threshold?: string;
  disc?: number;
  pupdate?: string;
  model?: 'dtd';
  inspector?: { enabled: boolean; maxTimesteps: number };
}

interface SimulationSuccessResponse {
  success: true;
  results: SimulationResult[][];
  metadata: SimulationMetadata;
}

interface SimulationErrorResponse {
  success: false;
  error?: string;
}

type SimulationResponse = SimulationSuccessResponse | SimulationErrorResponse;

interface SingleSimulationMetadata extends Partial<SimulationMetadata> {
  sampledParameters?: Record<string, number>;
}

type SingleSimulationResponse =
  | {
      success: true;
      result: SimulationResult[];
      metadata?: SingleSimulationMetadata;
      inspector?: SimulationInspector;
    }
  | SimulationErrorResponse;

export async function runSimulation(params: SimulationParams & { numNetworks: number }) {
  return fetchWithRetry<SimulationResponse>('/simulate', {
    method: 'POST',
    body: JSON.stringify(params),
  }, { timeout: SIMULATION_TIMEOUT_MS, retries: 1 });
}

// Run a single network simulation (for real-time progress tracking)
export async function runSimulationOne(params: SimulationParams) {
  return fetchWithRetry<SingleSimulationResponse>('/simulate-one', {
    method: 'POST',
    body: JSON.stringify(params),
  }, { timeout: SIMULATION_TIMEOUT_MS, retries: 0 });
}
