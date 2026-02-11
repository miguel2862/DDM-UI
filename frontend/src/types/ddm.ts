export type NPELayer =
  | 'US'
  | 'PrimarySensory'
  | 'AssociativeSensory'
  | 'Hippocampal'
  | 'AssociativeMotor'
  | 'PrimaryMotor'
  | 'Dopaminergic';

export type NPEType = 'Excitatory' | 'Inhibitory';

export interface NPE {
  name: string;
  type: NPEType;
  layer: NPELayer;
  activation: number;
  temporalSummation: number;
  activationDecay: number;
  mu: number;
  sigma: number;
  logisSigma: number;
}

export interface Connection {
  presynapticNPE: string;
  postsynapticNPE: string;
  weight: number;
  alpha: number;
  beta: number;
  alphaPrime: number;
  betaPrime: number;
}

export interface Trial {
  name: string;
  timesteps: string[];
}

export interface Phase {
  name: string;
  trialOrder: 'Random' | 'In bulk' | 'Alternated';
  trialTypes: string[];
  trialCounts: number[];
  hasITI: boolean;
  minITI?: number;
  maxITI?: number;
  itiTrialName?: string;
}

export interface SimulationResult {
  Phase: string;
  Trial: number;
  TimeStep: number;
  [key: string]: string | number;
}

export interface SimulationMetadata {
  numNetworks: number;
  phases: string[];
  units: string[];
  connections: string[];
  totalTrials: number;
  totalTimesteps: number;
  duration: number;
  disc?: number;
}

export interface SimulationResponse {
  success: boolean;
  results?: SimulationResult[][];
  metadata?: SimulationMetadata;
  error?: string;
}

export interface TemplateMetadata {
  id: string;
  name: string;
  description: string;
  category: string;
  available: boolean;
  phases: string;
  npeCount: number;
  connectionCount: number;
}

export interface TemplateData {
  id: string;
  name: string;
  available: boolean;
  npes: NPE[];
  connections: Connection[];
  trials: Record<string, string[]>;
  contingencies: string[];
  hasITI: boolean[];
}

export type ThresholdPreset = 'gaussian_ddmui' | 'gaussian_donahoe1993' | 'beta_ddmui';
export type AppMode = 'beginner' | 'advanced';
export type SimStatus = 'idle' | 'running' | 'complete' | 'error';
export interface PlaybackState {
  isPlaying: boolean;
  speed: number;
  currentIndex: number;
  maxIndex: number;
}
