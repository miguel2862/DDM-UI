// ══════════════════════════════════════════════════════════════════════════════
// DDM-UI — Templates — Port from api/templates.R
// ══════════════════════════════════════════════════════════════════════════════

import type { NPERow, ConnectionRow } from './types.js';

export interface TemplateData {
  id: string;
  name: string;
  npes: NPERow[];
  connections: ConnectionRow[];
  trials: Record<string, string[]>;
  contingencies: string[];
  hasITI: boolean[];
}

// ── Helper to create NPE rows ─────────────────────────────────────────────────
function makeNPEs(
  specs: Array<[string, string, string]>,
  defaults: { activation?: number; ts?: number; decay?: number; mu?: number; sigma?: number; logisSigma?: number } = {}
): NPERow[] {
  const {
    activation = 0, ts = 0.1, decay = 0.1,
    mu = 0.2, sigma = 0.15, logisSigma = 0.1
  } = defaults;
  return specs.map(([NPE, Type, Layer]) => ({
    NPE, Type: Type as NPERow['Type'], Layer: Layer as NPERow['Layer'],
    Activation: activation, 'Temporal.Summation': ts, 'Activation.Decay': decay,
    mu, sigma, logisSigma,
  }));
}

// ── Helper to create Connection rows ──────────────────────────────────────────
function makeConns(
  specs: Array<[string, string, number]>,
  defaults: { alpha?: number; beta?: number; alpha_prime?: number; beta_prime?: number } = {}
): ConnectionRow[] {
  const { alpha = 0.5, beta = 0.12, alpha_prime = 0.5, beta_prime = 0.12 } = defaults;
  return specs.map(([PreSinapticNPE, PostSinapticNPE, Weight]) => ({
    PreSinapticNPE, PostSinapticNPE, Weight, alpha, beta, alpha_prime, beta_prime,
  }));
}

// ═══════════════════════════════════════════════════════════════════════════════
// Templates
// ═══════════════════════════════════════════════════════════════════════════════

export function getExtinctionTemplate(): TemplateData {
  return {
    id: 'extinction',
    name: 'Extinction',
    npes: makeNPEs([
      ['US', 'Excitatory', 'US'],
      ['D', 'Excitatory', 'Dopaminergic'],
      ['S1', 'Excitatory', 'PrimarySensory'],
      ['S..1', 'Excitatory', 'AssociativeSensory'],
      ['H1', 'Excitatory', 'Hippocampal'],
      ['M..1', 'Excitatory', 'AssociativeMotor'],
      ['M.1', 'Excitatory', 'PrimaryMotor'],
    ]),
    connections: makeConns([
      ['S1', 'S..1', 0.1],
      ['S..1', 'H1', 0.1],
      ['S..1', 'M..1', 0.1],
      ['M..1', 'D', 0.1],
      ['M..1', 'M.1', 0.1],
      ['US', 'D', 1.0],
    ]),
    trials: {
      Training: [
        'US,0.00,S1,1.00,True',
        'US,0.00,S1,1.00,True',
        'US,0.00,S1,1.00,True',
        'US,0.00,S1,1.00,True',
        'US,1.00,S1,1.00,True',
      ],
      Extinction: [
        'US,0.00,S1,1.00,True',
        'US,0.00,S1,1.00,True',
        'US,0.00,S1,1.00,True',
        'US,0.00,S1,1.00,True',
        'US,0.00,S1,1.00,True',
      ],
    },
    contingencies: [
      'training, Random, Training, 100, False',
      'extinction, Random, Extinction, 100, False',
    ],
    hasITI: [false, false],
  };
}

export function getAcquisitionTemplate(): TemplateData {
  return {
    id: 'acquisition',
    name: 'Acquisition',
    npes: makeNPEs([
      ['US', 'Excitatory', 'US'],
      ['D', 'Excitatory', 'Dopaminergic'],
      ['S1', 'Excitatory', 'PrimarySensory'],
      ['S..1', 'Excitatory', 'AssociativeSensory'],
      ['H1', 'Excitatory', 'Hippocampal'],
      ['M..1', 'Excitatory', 'AssociativeMotor'],
      ['M.1', 'Excitatory', 'PrimaryMotor'],
    ]),
    connections: makeConns([
      ['S1', 'S..1', 0.1],
      ['S..1', 'H1', 0.1],
      ['S..1', 'M..1', 0.1],
      ['M..1', 'D', 0.1],
      ['M..1', 'M.1', 0.1],
      ['US', 'D', 1.0],
    ]),
    trials: {
      Training: [
        'US,0.00,S1,1.00,True',
        'US,0.00,S1,1.00,True',
        'US,0.00,S1,1.00,True',
        'US,0.00,S1,1.00,True',
        'US,1.00,S1,1.00,True',
      ],
    },
    contingencies: [
      'training, Random, Training, 100, False',
    ],
    hasITI: [false],
  };
}

export function getSpontaneousRecoveryTemplate(): TemplateData {
  return {
    id: 'spontaneous_recovery',
    name: 'Spontaneous Recovery',
    npes: makeNPEs([
      ['US', 'Excitatory', 'US'],
      ['D', 'Excitatory', 'Dopaminergic'],
      ['S1', 'Excitatory', 'PrimarySensory'],
      ['S..1', 'Excitatory', 'AssociativeSensory'],
      ['H1', 'Excitatory', 'Hippocampal'],
      ['M..1', 'Excitatory', 'AssociativeMotor'],
      ['M.1', 'Excitatory', 'PrimaryMotor'],
    ]),
    connections: makeConns([
      ['S1', 'S..1', 0.1],
      ['S..1', 'H1', 0.1],
      ['S..1', 'M..1', 0.1],
      ['M..1', 'D', 0.1],
      ['M..1', 'M.1', 0.1],
      ['US', 'D', 1.0],
    ]),
    trials: {
      Training: [
        'US,0.00,S1,1.00,True',
        'US,0.00,S1,1.00,True',
        'US,0.00,S1,1.00,True',
        'US,0.00,S1,1.00,True',
        'US,1.00,S1,1.00,True',
      ],
      Extinction: [
        'US,0.00,S1,1.00,True',
        'US,0.00,S1,1.00,True',
        'US,0.00,S1,1.00,True',
        'US,0.00,S1,1.00,True',
        'US,0.00,S1,1.00,True',
      ],
      Rest: ['US,0.00,S1,0.00,True'],
      Test: [
        'US,0.00,S1,1.00,False',
        'US,0.00,S1,1.00,False',
        'US,0.00,S1,1.00,False',
        'US,0.00,S1,1.00,False',
        'US,0.00,S1,1.00,False',
      ],
    },
    contingencies: [
      'training, Random, Training, 100, False',
      'extinction, Random, Extinction, 100, False',
      'rest, Random, Rest, 50, False',
      'test, In bulk, Test, 25, False',
    ],
    hasITI: [false, false, false, false],
  };
}

export function getLatentInhibitionTemplate(): TemplateData {
  return {
    id: 'latent_inhibition',
    name: 'Latent Inhibition',
    npes: makeNPEs([
      ['US', 'Excitatory', 'US'],
      ['D', 'Excitatory', 'Dopaminergic'],
      ['S1', 'Excitatory', 'PrimarySensory'],
      ['S..1', 'Excitatory', 'AssociativeSensory'],
      ['H1', 'Excitatory', 'Hippocampal'],
      ['M..1', 'Excitatory', 'AssociativeMotor'],
      ['M.1', 'Excitatory', 'PrimaryMotor'],
    ]),
    connections: makeConns([
      ['S1', 'S..1', 0.1],
      ['S..1', 'H1', 0.1],
      ['S..1', 'M..1', 0.1],
      ['M..1', 'D', 0.1],
      ['M..1', 'M.1', 0.1],
      ['US', 'D', 1.0],
    ]),
    trials: {
      PreExposure: [
        'US,0.00,S1,1.00,True',
        'US,0.00,S1,1.00,True',
        'US,0.00,S1,1.00,True',
        'US,0.00,S1,1.00,True',
        'US,0.00,S1,1.00,True',
      ],
      Training: [
        'US,0.00,S1,1.00,True',
        'US,0.00,S1,1.00,True',
        'US,0.00,S1,1.00,True',
        'US,0.00,S1,1.00,True',
        'US,1.00,S1,1.00,True',
      ],
    },
    contingencies: [
      'pre-exposure, Random, PreExposure, 100, False',
      'training, Random, Training, 100, False',
    ],
    hasITI: [false, false],
  };
}

export function getBlockingTemplate(): TemplateData {
  const npes = makeNPEs([
    ['US', 'Excitatory', 'US'],
    ['D', 'Excitatory', 'Dopaminergic'],
    ['A', 'Excitatory', 'PrimarySensory'],
    ['C', 'Excitatory', 'PrimarySensory'],
    ['X', 'Excitatory', 'PrimarySensory'],
    ['S..1', 'Excitatory', 'AssociativeSensory'],
    ['S..2', 'Excitatory', 'AssociativeSensory'],
    ['H1', 'Excitatory', 'Hippocampal'],
    ['H2', 'Excitatory', 'Hippocampal'],
    ['M..1', 'Excitatory', 'AssociativeMotor'],
    ['M..2', 'Excitatory', 'AssociativeMotor'],
    ['M.1', 'Excitatory', 'PrimaryMotor'],
    ['M.2', 'Excitatory', 'PrimaryMotor'],
  ]);
  const connections = makeConns([
    ['A', 'S..1', 0.2], ['C', 'S..1', 0.2], ['C', 'S..2', 0.2], ['X', 'S..2', 0.2],
    ['S..1', 'H1', 0.2], ['S..1', 'M..1', 0.2], ['S..2', 'H2', 0.2], ['S..2', 'M..2', 0.2],
    ['S..1', 'M..2', 0.2], ['S..2', 'M..1', 0.2],
    ['M..1', 'D', 0.2], ['M..2', 'D', 0.2], ['M..1', 'M.1', 0.2], ['M..2', 'M.2', 0.2],
    ['US', 'D', 1.0], ['US', 'M.1', 1.0], ['US', 'M.2', 1.0],
  ]);
  return {
    id: 'burgos_donahoe_blocking',
    name: 'Blocking (Burgos & Donahoe, 2016)',
    npes,
    connections,
    trials: {
      'A+': [
        'US,0.00,A,1.00,C,0.90,X,0.00,True',
        'US,0.00,A,1.00,C,0.90,X,0.00,True',
        'US,0.00,A,1.00,C,0.90,X,0.00,True',
        'US,0.00,A,1.00,C,0.90,X,0.00,True',
        'US,1.00,A,1.00,C,0.90,X,0.00,True',
      ],
      'AX+': [
        'US,0.00,A,1.00,C,0.90,X,1.00,True',
        'US,0.00,A,1.00,C,0.90,X,1.00,True',
        'US,0.00,A,1.00,C,0.90,X,1.00,True',
        'US,0.00,A,1.00,C,0.90,X,1.00,True',
        'US,1.00,A,1.00,C,0.90,X,1.00,True',
      ],
      'X TST': [
        'US,0.00,A,0.00,C,0.90,X,1.00,False',
        'US,0.00,A,0.00,C,0.90,X,1.00,False',
        'US,0.00,A,0.00,C,0.90,X,1.00,False',
        'US,0.00,A,0.00,C,0.90,X,1.00,False',
        'US,0.00,A,0.00,C,0.90,X,1.00,False',
      ],
    },
    contingencies: [
      'Entrenamiento, Random, A+, 100, False',
      'Bloqueo, Random, AX+, 100, False',
      'Test, In bulk, X TST, 25, False',
    ],
    hasITI: [false, false, false],
  };
}

/** All templates by ID */
export const TEMPLATES: Record<string, () => TemplateData> = {
  extinction: getExtinctionTemplate,
  acquisition: getAcquisitionTemplate,
  spontaneous_recovery: getSpontaneousRecoveryTemplate,
  latent_inhibition: getLatentInhibitionTemplate,
  burgos_donahoe_blocking: getBlockingTemplate,
};
