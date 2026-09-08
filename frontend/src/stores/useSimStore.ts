import { create } from 'zustand';
import { assertDDMPayload } from '../utils/ddmCompatibility';
import type { SimulationInspector } from '../types/inspector';
import type {
  NPE, Connection, Phase, SimStatus, AppMode, SimulationResult,
  SimulationMetadata, ThresholdPreset, ModelKind,
} from '../types/ddm';

// Maps each unified preset to its R-backend threshold type + NPE mu/sigma
const THRESHOLD_CONFIGS: Record<ThresholdPreset, { threshold: string; mu: number; sigma: number }> = {
  'gaussian_ddmui':       { threshold: 'gaussian', mu: 0.2, sigma: 0.15 },
  'gaussian_donahoe1993': { threshold: 'gaussian', mu: 0.0, sigma: 1.0 },
  'beta_ddmui':           { threshold: 'beta',     mu: 0.2, sigma: 0.15 },
};





type UnknownRecord = Record<string, unknown>;

const VALID_NPE_LAYERS: readonly NPE['layer'][] = [
  'US', 'PrimarySensory', 'AssociativeSensory', 'Hippocampal', 'AssociativeMotor',
  'PrimaryMotor', 'Dopaminergic',
];

function isRecord(value: unknown): value is UnknownRecord {
  return typeof value === 'object' && value !== null && !Array.isArray(value);
}

function asRecord(value: unknown): UnknownRecord {
  return isRecord(value) ? value : {};
}

function unboxScalar(value: unknown): unknown {
  return Array.isArray(value) ? value[0] : value;
}

function toStringValue(value: unknown, fallback = ''): string {
  const scalar = unboxScalar(value);
  return scalar === undefined || scalar === null ? fallback : String(scalar);
}

function toNumberValue(value: unknown, fallback: number): number {
  const parsed = Number(unboxScalar(value));
  return Number.isFinite(parsed) ? parsed : fallback;
}

function toBooleanValue(value: unknown, fallback = false): boolean {
  const scalar = unboxScalar(value);
  if (typeof scalar === 'boolean') return scalar;
  if (typeof scalar === 'number') return scalar !== 0;
  if (typeof scalar === 'string') {
    if (scalar.toLowerCase() === 'true') return true;
    if (scalar.toLowerCase() === 'false') return false;
  }
  return fallback;
}

function toStringArray(value: unknown, fallback: string[] = []): string[] {
  if (value === undefined || value === null) return [...fallback];
  const values = Array.isArray(value) ? value : [value];
  return values.map((entry) => toStringValue(entry));
}

function toBooleanArray(value: unknown): boolean[] {
  const values = Array.isArray(value) ? value : value === undefined ? [] : [value];
  return values.map((entry) => toBooleanValue(entry));
}

function columnsToRows(value: unknown): UnknownRecord[] {
  if (Array.isArray(value)) return value.filter(isRecord);
  if (!isRecord(value)) return [];

  const keys = Object.keys(value);
  if (keys.length === 0) return [];
  const firstValue = value[keys[0]];
  if (!Array.isArray(firstValue)) return [value];

  return Array.from({ length: firstValue.length }, (_, index) => {
    const row: UnknownRecord = {};
    for (const key of keys) {
      const column = value[key];
      row[key] = Array.isArray(column) ? column[index] : column;
    }
    return row;
  });
}

function parseNPEType(value: unknown): NPE['type'] {
  return toStringValue(value) === 'Inhibitory' ? 'Inhibitory' : 'Excitatory';
}

function parseNPELayer(value: unknown): NPE['layer'] {
  const layer = toStringValue(value);
  if (!VALID_NPE_LAYERS.includes(layer as NPE['layer'])) throw new Error(`Unsupported DDM layer: ${layer}`);
  return layer as NPE['layer'];
}

function parseNPEs(value: unknown): NPE[] {
  return columnsToRows(value).map((row) => ({
    name: toStringValue(row.NPE ?? row.name),
    type: parseNPEType(row.Type ?? row.type),
    layer: parseNPELayer(row.Layer ?? row.layer),
    activation: toNumberValue(row.Activation ?? row.activation, 0),
    temporalSummation: toNumberValue(row['Temporal.Summation'] ?? row.temporalSummation, 0),
    activationDecay: toNumberValue(row['Activation.Decay'] ?? row.activationDecay, 0),
    mu: toNumberValue(row.mu, 0.2),
    sigma: toNumberValue(row.sigma, 0.15),
    logisSigma: toNumberValue(row.logisSigma, 1),
  }));
}

function parseConnections(value: unknown): Connection[] {
  return columnsToRows(value).map((row) => ({
    presynapticNPE: toStringValue(row.PreSinapticNPE ?? row.presynapticNPE),
    postsynapticNPE: toStringValue(row.PostSinapticNPE ?? row.postsynapticNPE),
    weight: toNumberValue(row.Weight ?? row.weight, 0),
    alpha: toNumberValue(row.alpha, 0),
    beta: toNumberValue(row.beta, 0),
    alphaPrime: toNumberValue(row.alpha_prime ?? row.alphaPrime, 0),
    betaPrime: toNumberValue(row.beta_prime ?? row.betaPrime, 0),
  }));
}

function parseTrials(value: unknown): Record<string, string[]> {
  return Object.fromEntries(
    Object.entries(asRecord(value)).map(([name, timesteps]) => [name, toStringArray(timesteps)]),
  );
}

function parseLockedLayout(value: unknown): Record<string, { x: number; y: number }> | null {
  if (!isRecord(value)) return null;
  const entries = Object.entries(value).flatMap(([name, rawPosition]) => {
    if (!isRecord(rawPosition)) return [];
    const x = toNumberValue(rawPosition.x, Number.NaN);
    const y = toNumberValue(rawPosition.y, Number.NaN);
    return Number.isFinite(x) && Number.isFinite(y) ? [[name, { x, y }] as const] : [];
  });
  return entries.length > 0 ? Object.fromEntries(entries) : null;
}





interface SimStore {
  // Model configuration
  modelKind: ModelKind;

  npes: NPE[];
  connections: Connection[];
  trials: Record<string, string[]>;
  contingencies: string[];
  hasITI: boolean[];
  phases: Phase[];

  // Simulation state
  simulationResults: SimulationResult[][] | null;
  simulationMetadata: SimulationMetadata | null;
  simulationInspector: SimulationInspector | null;
  simStatus: SimStatus;
  simError: string | null;
  numNetworks: number;
  thresholdPreset: ThresholdPreset;
  disc: number;
  pupdate: string;

  // UI state
  appMode: AppMode;
  selectedTemplate: string | null;
  selectedNetwork: number;

  // Layout
  lockedLayout: Record<string, { x: number; y: number }> | null;

  // Playback
  playbackIndex: number;
  isPlaying: boolean;
  playbackSpeed: number;

  // Derived getter — R-backend threshold type
  getThresholdType: () => string;

  // Actions
  setNPEs: (npes: NPE[]) => void;
  setModelKind: (model: ModelKind) => void;


  addNPE: (npe: NPE) => void;
  removeNPE: (name: string) => void;
  setConnections: (connections: Connection[]) => void;
  addConnection: (connection: Connection) => void;
  removeConnection: (pre: string, post: string) => void;
  setTrials: (trials: Record<string, string[]>) => void;
  addTrial: (name: string, timesteps: string[]) => void;
  removeTrial: (name: string) => void;
  setContingencies: (contingencies: string[]) => void;
  setHasITI: (hasITI: boolean[]) => void;
  setPhases: (phases: Phase[]) => void;
  setSimResults: (results: SimulationResult[][], metadata: SimulationMetadata, inspector?: SimulationInspector | null) => void;
  setSimStatus: (status: SimStatus) => void;
  setSimError: (error: string | null) => void;
  setNumNetworks: (n: number) => void;
  setThresholdPreset: (preset: ThresholdPreset) => void;
  setDisc: (d: number) => void;
  setPupdate: (p: string) => void;
  setAppMode: (mode: AppMode) => void;
  setSelectedTemplate: (id: string | null) => void;
  setSelectedNetwork: (n: number) => void;
  setPlaybackIndex: (i: number) => void;
  setIsPlaying: (p: boolean) => void;
  setPlaybackSpeed: (s: number) => void;
  setLockedLayout: (layout: Record<string, { x: number; y: number }> | null) => void;
  loadTemplate: (data: unknown) => void;
  loadExperiment: (data: unknown) => void;
  reset: () => void;
}

const initialState = {
  modelKind: 'dtd' as ModelKind,

  npes: [],
  connections: [],
  trials: {},
  contingencies: [],
  hasITI: [],
  phases: [],
  simulationResults: null,
  simulationMetadata: null,
  simulationInspector: null,
  simStatus: 'idle' as SimStatus,
  simError: null,
  numNetworks: 5,
  thresholdPreset: 'gaussian_ddmui' as ThresholdPreset,
  disc: 0.0015,
  pupdate: 'async_random',
  appMode: 'beginner' as AppMode,
  selectedTemplate: null,
  selectedNetwork: 0,
  lockedLayout: null,
  playbackIndex: 0,
  isPlaying: false,
  playbackSpeed: 1,
};

export const useSimStore = create<SimStore>((set, get) => ({
  ...initialState,

  getThresholdType: () => THRESHOLD_CONFIGS[get().thresholdPreset].threshold,

  setModelKind: (modelKind) => set({ modelKind, simulationInspector: null }),


  setNPEs: (npes) => set({ npes, simulationInspector: null }),
  addNPE: (npe) => set((s) => ({ npes: [...s.npes, npe], simulationInspector: null })),
  removeNPE: (name) => set((s) => ({
    simulationInspector: null,
    npes: s.npes.filter((n) => n.name !== name),
    connections: s.connections.filter((c) => c.presynapticNPE !== name && c.postsynapticNPE !== name),
  })),
  setConnections: (connections) => set({ connections, simulationInspector: null }),
  addConnection: (connection) => set((s) => ({ connections: [...s.connections, connection], simulationInspector: null })),
  removeConnection: (pre, post) => set((s) => ({
    simulationInspector: null,
    connections: s.connections.filter((c) => !(c.presynapticNPE === pre && c.postsynapticNPE === post)),
  })),
  setTrials: (trials) => set({ trials, simulationInspector: null }),
  addTrial: (name, timesteps) => set((s) => ({ trials: { ...s.trials, [name]: timesteps }, simulationInspector: null })),
  removeTrial: (name) => set((s) => {
    const trials = { ...s.trials };
    delete trials[name];
    return { trials, simulationInspector: null };
  }),
  setContingencies: (contingencies) => set({ contingencies, simulationInspector: null }),
  setHasITI: (hasITI) => set({ hasITI, simulationInspector: null }),
  setPhases: (phases) => set({ phases, simulationInspector: null }),
  setSimResults: (results, metadata, inspector = null) => set({ simulationResults: results, simulationMetadata: metadata, simulationInspector: inspector, simStatus: 'complete' }),
  setSimStatus: (simStatus) => set({ simStatus }),
  setSimError: (simError) => set({ simError, simStatus: simError ? 'error' : 'idle' }),
  setNumNetworks: (numNetworks) => set({ numNetworks }),
  setThresholdPreset: (thresholdPreset) => set((s) => {
    const { mu, sigma } = THRESHOLD_CONFIGS[thresholdPreset];
    return {
      thresholdPreset,
      simulationInspector: null,
      npes: s.npes.map((npe) => ({ ...npe, mu, sigma })),
    };
  }),
  setDisc: (disc) => set({ disc, simulationInspector: null }),
  setPupdate: (pupdate) => set({ pupdate, simulationInspector: null }),
  setAppMode: (appMode) => set({ appMode }),
  setSelectedTemplate: (selectedTemplate) => set({ selectedTemplate }),
  setSelectedNetwork: (selectedNetwork) => set({ selectedNetwork }),
  setPlaybackIndex: (playbackIndex) => set({ playbackIndex }),
  setIsPlaying: (isPlaying) => set({ isPlaying }),
  setPlaybackSpeed: (playbackSpeed) => set({ playbackSpeed }),
  setLockedLayout: (lockedLayout) => set({ lockedLayout }),

  loadTemplate: (payload) => {
    // R plumber can return dataframes in two formats:
    // 1. Array of row-objects: [{ NPE: "US", Type: "Excitatory", ... }, ...]
    // 2. Column-arrays (named list): { NPE: ["US", "D"], Type: ["Excitatory", "Excitatory"], ... }
    // We handle both formats.
    assertDDMPayload(payload);
    const data = asRecord(payload);
    const npes = parseNPEs(data.npes);
    const connections = parseConnections(data.connections);
    const templateId = toStringValue(data.id) || null;
    const modelKind: ModelKind = 'dtd';
    const contingencies = toStringArray(data.contingencies);
    const parsedHasITI = toBooleanArray(data.hasITI);
    const hasITI = parsedHasITI.length > 0 ? parsedHasITI : contingencies.map(() => false);

    // Auto-detect threshold preset from NPE mu/sigma values
    const firstMu = npes.length > 0 ? npes[0].mu : 0.2;
    const firstSigma = npes.length > 0 ? npes[0].sigma : 0.15;
    const detectedPreset: ThresholdPreset =
      (firstMu === 0 && firstSigma === 1) ? 'gaussian_donahoe1993' : 'gaussian_ddmui';

    set({
      npes,
      modelKind,

      connections,
      trials: parseTrials(data.trials),
      contingencies,
      hasITI,
      selectedTemplate: templateId,
      simulationInspector: null,
      thresholdPreset: detectedPreset,
      simulationResults: null,
      simulationMetadata: null,
      simStatus: 'idle',
      simError: null,
    });
  },

  loadExperiment: (payload) => {
    assertDDMPayload(payload);
    const data = asRecord(payload);
    // Backward compatibility: convert old separate threshold+thresholdPreset to unified format
    const rawPreset = toStringValue(data.thresholdPreset, 'gaussian_ddmui');
    let preset: ThresholdPreset;
    if (rawPreset === 'ddm-ui' || rawPreset === 'donahoe1993') {
      // Old format: 'ddm-ui' / 'donahoe1993' were separate from threshold type
      const oldThreshold = toStringValue(data.threshold, 'gaussian');
      if (rawPreset === 'donahoe1993') {
        preset = 'gaussian_donahoe1993';
      } else if (oldThreshold === 'beta') {
        preset = 'beta_ddmui';
      } else {
        preset = 'gaussian_ddmui';
      }
    } else {
      // New unified format — validate or default
      preset = (['gaussian_ddmui', 'gaussian_donahoe1993', 'beta_ddmui'] as ThresholdPreset[]).includes(rawPreset as ThresholdPreset)
        ? (rawPreset as ThresholdPreset)
        : 'gaussian_ddmui';
    }

    // Load from a saved experiment JSON file (already in frontend format)
    set({
      npes: parseNPEs(data.npes),
      modelKind: 'dtd',

      connections: parseConnections(data.connections),
      trials: parseTrials(data.trials),
      contingencies: toStringArray(data.contingencies),
      hasITI: toBooleanArray(data.hasITI),
      numNetworks: Math.max(1, Math.floor(toNumberValue(data.numNetworks, 5))),
      thresholdPreset: preset,
      disc: toNumberValue(data.disc, 0.0015),
      pupdate: toStringValue(data.pupdate, 'async_random'),
      lockedLayout: parseLockedLayout(data.lockedLayout),
      selectedTemplate: null,
      simulationInspector: null,
      simulationResults: null,
      simulationMetadata: null,
      simStatus: 'idle',
      simError: null,
    });
  },

  reset: () => set(initialState),
}));
