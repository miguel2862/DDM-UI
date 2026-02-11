// ══════════════════════════════════════════════════════════════════════════════
// DDM-UI — Type Definitions (direct mapping from R S4 classes)
// ══════════════════════════════════════════════════════════════════════════════

export type NPEType = 'Excitatory' | 'Inhibitory';

export type Layer =
  | 'US'
  | 'PrimarySensory'
  | 'AssociativeSensory'
  | 'Hippocampal'
  | 'Dopaminergic'
  | 'AssociativeMotor'
  | 'PrimaryMotor';

/** Matches R S4 class "Connection" */
export interface Connection {
  weight: number;
  Name: string;
  alpha: number;
  beta: number;
  alpha_prime: number;
  beta_prime: number;
  preSinapticNPE: string;
  p: number;
}

/** Matches R S4 class "NPE" */
export interface NPE {
  Activation: number;
  PreviousActivation: number;
  ActivationDecay: number;
  ExcitatoryInput: number;
  InhibitoryInput: number;
  PreviousExcitatoryInput: number;
  TemporalSummation: number;
  Name: string;
  Type: NPEType;
  Layer: Layer;
  Threshold: number;
  mu: number;
  sigma: number;
  logisSigma: number;
  InputConnections: Record<string, Connection>;
  r: [number, number];
}

/** NPE row from the input dataframe (9 columns) */
export interface NPERow {
  NPE: string;
  Type: NPEType;
  Layer: Layer;
  Activation: number;
  'Temporal.Summation': number;
  'Activation.Decay': number;
  mu: number;
  sigma: number;
  logisSigma: number;
}

/** Connection row from the input dataframe (7 columns) */
export interface ConnectionRow {
  PreSinapticNPE: string;
  PostSinapticNPE: string;
  Weight: number;
  alpha: number;
  beta: number;
  alpha_prime: number;
  beta_prime: number;
}

/** A single row in the TimeSteps table */
export interface TimeStepRow {
  Phase: string;
  Trial: number;
  TimeStep: number;
  /** Pairs of [npeName, activationValue, ...] followed by learningActive boolean */
  stimuli: Array<{ npe: string; activation: number }>;
  learningActive: boolean;
}

/** Options for the simulation */
export interface SimulationOptions {
  threshold?: 'gaussian' | 'beta';
  disc?: number;
  pupdate?: 'async_random' | 'async_sequential' | 'sync_random' | 'sync_sequential';
  saveData?: {
    Elements?: string[];
    TimeSteps?: number[];
  };
  /** Fixed seed for reproducibility (optional) */
  seed?: number;
}

/** One row of simulation output */
export interface SimulationOutputRow {
  Phase: string;
  Trial: number;
  TimeStep: number;
  [key: string]: number | string; // NPE activations + Connection weights + dVTA + dH
}
