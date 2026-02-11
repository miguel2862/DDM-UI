import { create } from 'zustand';
import type { NPE, Connection, Phase, SimStatus, AppMode, SimulationResult, SimulationMetadata, ThresholdPreset } from '../types/ddm';

// Maps each unified preset to its R-backend threshold type + NPE mu/sigma
const THRESHOLD_CONFIGS: Record<ThresholdPreset, { threshold: string; mu: number; sigma: number }> = {
  'gaussian_ddmui':       { threshold: 'gaussian', mu: 0.2, sigma: 0.15 },
  'gaussian_donahoe1993': { threshold: 'gaussian', mu: 0.0, sigma: 1.0 },
  'beta_ddmui':           { threshold: 'beta',     mu: 0.2, sigma: 0.15 },
};

interface SimStore {
  // Model configuration
  npes: NPE[];
  connections: Connection[];
  trials: Record<string, string[]>;
  contingencies: string[];
  hasITI: boolean[];
  phases: Phase[];

  // Simulation state
  simulationResults: SimulationResult[][] | null;
  simulationMetadata: SimulationMetadata | null;
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
  setSimResults: (results: SimulationResult[][], metadata: SimulationMetadata) => void;
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
  loadTemplate: (data: any) => void;
  loadExperiment: (data: any) => void;
  reset: () => void;
}

const initialState = {
  npes: [],
  connections: [],
  trials: {},
  contingencies: [],
  hasITI: [],
  phases: [],
  simulationResults: null,
  simulationMetadata: null,
  simStatus: 'idle' as SimStatus,
  simError: null,
  numNetworks: 5,
  thresholdPreset: 'gaussian_ddmui' as ThresholdPreset,
  disc: 0.001,
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

  setNPEs: (npes) => set({ npes }),
  addNPE: (npe) => set((s) => ({ npes: [...s.npes, npe] })),
  removeNPE: (name) => set((s) => ({
    npes: s.npes.filter((n) => n.name !== name),
    connections: s.connections.filter((c) => c.presynapticNPE !== name && c.postsynapticNPE !== name),
  })),
  setConnections: (connections) => set({ connections }),
  addConnection: (connection) => set((s) => ({ connections: [...s.connections, connection] })),
  removeConnection: (pre, post) => set((s) => ({
    connections: s.connections.filter((c) => !(c.presynapticNPE === pre && c.postsynapticNPE === post)),
  })),
  setTrials: (trials) => set({ trials }),
  addTrial: (name, timesteps) => set((s) => ({ trials: { ...s.trials, [name]: timesteps } })),
  removeTrial: (name) => set((s) => {
    const { [name]: _, ...rest } = s.trials;
    return { trials: rest };
  }),
  setContingencies: (contingencies) => set({ contingencies }),
  setHasITI: (hasITI) => set({ hasITI }),
  setPhases: (phases) => set({ phases }),
  setSimResults: (results, metadata) => set({ simulationResults: results, simulationMetadata: metadata, simStatus: 'complete' }),
  setSimStatus: (simStatus) => set({ simStatus }),
  setSimError: (simError) => set({ simError, simStatus: simError ? 'error' : 'idle' }),
  setNumNetworks: (numNetworks) => set({ numNetworks }),
  setThresholdPreset: (thresholdPreset) => set((s) => {
    const { mu, sigma } = THRESHOLD_CONFIGS[thresholdPreset];
    return {
      thresholdPreset,
      npes: s.npes.map((npe) => ({ ...npe, mu, sigma })),
    };
  }),
  setDisc: (disc) => set({ disc }),
  setPupdate: (pupdate) => set({ pupdate }),
  setAppMode: (appMode) => set({ appMode }),
  setSelectedTemplate: (selectedTemplate) => set({ selectedTemplate }),
  setSelectedNetwork: (selectedNetwork) => set({ selectedNetwork }),
  setPlaybackIndex: (playbackIndex) => set({ playbackIndex }),
  setIsPlaying: (isPlaying) => set({ isPlaying }),
  setPlaybackSpeed: (playbackSpeed) => set({ playbackSpeed }),
  setLockedLayout: (lockedLayout) => set({ lockedLayout }),

  loadTemplate: (data: any) => {
    // R plumber can return dataframes in two formats:
    // 1. Array of row-objects: [{ NPE: "US", Type: "Excitatory", ... }, ...]
    // 2. Column-arrays (named list): { NPE: ["US", "D"], Type: ["Excitatory", "Excitatory"], ... }
    // We handle both formats.

    const unbox = (v: any) => Array.isArray(v) ? v[0] : v;

    // Helper: convert column-arrays format to row-objects format
    const columnsToRows = (obj: any): any[] => {
      if (Array.isArray(obj)) return obj; // Already row-objects
      const keys = Object.keys(obj);
      if (keys.length === 0) return [];
      const firstVal = obj[keys[0]];
      if (!Array.isArray(firstVal)) return [obj]; // Single row wrapped in object
      const len = firstVal.length;
      const rows: any[] = [];
      for (let i = 0; i < len; i++) {
        const row: any = {};
        for (const key of keys) {
          row[key] = Array.isArray(obj[key]) ? obj[key][i] : obj[key];
        }
        rows.push(row);
      }
      return rows;
    };

    const rawNpes = columnsToRows(data.npes || []);
    const npes: NPE[] = rawNpes.map((row: any) => ({
      name: row.NPE,
      type: row.Type as NPE['type'],
      layer: row.Layer as NPE['layer'],
      activation: Number(row.Activation),
      temporalSummation: Number(row['Temporal.Summation']),
      activationDecay: Number(row['Activation.Decay']),
      mu: Number(row.mu),
      sigma: Number(row.sigma),
      logisSigma: Number(row.logisSigma),
    }));

    const rawConns = columnsToRows(data.connections || []);
    const connections: Connection[] = rawConns.map((row: any) => ({
      presynapticNPE: row.PreSinapticNPE,
      postsynapticNPE: row.PostSinapticNPE,
      weight: Number(row.Weight),
      alpha: Number(row.alpha),
      beta: Number(row.beta),
      alphaPrime: Number(row.alpha_prime),
      betaPrime: Number(row.beta_prime),
    }));

    const templateId = unbox(data.id);
    const contingencies = data.contingencies || [];
    const hasITI = data.hasITI || contingencies.map(() => false);

    // Auto-detect threshold preset from NPE mu/sigma values
    const firstMu = npes.length > 0 ? npes[0].mu : 0.2;
    const firstSigma = npes.length > 0 ? npes[0].sigma : 0.15;
    const detectedPreset: ThresholdPreset =
      (firstMu === 0 && firstSigma === 1) ? 'gaussian_donahoe1993' : 'gaussian_ddmui';

    set({
      npes,
      connections,
      trials: data.trials,
      contingencies,
      hasITI,
      selectedTemplate: templateId,
      thresholdPreset: detectedPreset,
      simulationResults: null,
      simulationMetadata: null,
      simStatus: 'idle',
      simError: null,
    });
  },

  loadExperiment: (data: any) => {
    // Backward compatibility: convert old separate threshold+thresholdPreset to unified format
    const rawPreset: string = data.thresholdPreset || 'gaussian_ddmui';
    let preset: ThresholdPreset;
    if (rawPreset === 'ddm-ui' || rawPreset === 'donahoe1993') {
      // Old format: 'ddm-ui' / 'donahoe1993' were separate from threshold type
      const oldThreshold = data.threshold || 'gaussian';
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
      npes: data.npes || [],
      connections: data.connections || [],
      trials: data.trials || {},
      contingencies: data.contingencies || [],
      hasITI: data.hasITI || [],
      numNetworks: data.numNetworks || 5,
      thresholdPreset: preset,
      disc: data.disc || 0.001,
      pupdate: data.pupdate || 'async_random',
      lockedLayout: data.lockedLayout || null,
      selectedTemplate: null,
      simulationResults: null,
      simulationMetadata: null,
      simStatus: 'idle',
      simError: null,
    });
  },

  reset: () => set(initialState),
}));
