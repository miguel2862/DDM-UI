/** Recorded operands from the original R engine, not a frontend re-simulation. */
export interface InspectorUnit {
  name: string;
  layer: string;
  order: number;
  branch: 'external' | 'unconditional' | 'suprathreshold' | 'subthreshold' | 'inhibited';
  previousActivation: number;
  activation: number;
  previousExcitatoryInput: number;
  excInput: number | null;
  inhInput: number | null;
  logisticExc: number | null;
  logisticInh: number | null;
  threshold: number | null;
  mu: number;
  sigma: number;
  logisSigma: number;
  temporalSummation: number;
  activationDecay: number;
}

export interface InspectorConnection {
  name: string;
  pre: string;
  post: string;
  order: number;
  preType: string;
  preActivation: number;
  postActivation: number;
  weightBefore: number;
  weightUnclipped: number;
  weightAfter: number;
  deltaWeight: number;
  branch: 'potentiation' | 'decrement' | 'fixedUS' | 'learningOff';
  signalKind: 'dD' | 'dH' | null;
  signal: number | null;
  disc: number;
  capacity: number | null;
  proportion: number | null;
  alpha: number | null;
  beta: number | null;
  inputDenominator: number | null;
}

export interface InspectorStep {
  rowIndex: number;
  phase: string;
  trial: number;
  timestep: number;
  resetApplied: boolean;
  learningEnabled: boolean;
  activationOrder: string[];
  learningOrder: string[];
  dD: number | null;
  dH: number | null;
  previousDH: number;
  units: InspectorUnit[];
  connections: InspectorConnection[];
}

export interface SimulationInspector {
  schemaVersion: 1;
  engine: string;
  engineHash: string;
  recordedTimesteps: number;
  totalTimesteps: number;
  truncated: boolean;
  steps: InspectorStep[];
}
