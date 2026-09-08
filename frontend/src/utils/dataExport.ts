import { assertDDMPayload } from './ddmCompatibility';
import type {
  NPE,
  Connection,
  ModelKind,
  SimulationResult,
} from '../types/ddm';

// Download simulation results as CSV for a single network
export function downloadResultsCSV(results: SimulationResult[], filename = 'ddm-results.csv') {
  if (!results || results.length === 0) return;

  // Get all column names from first row
  const columns = Object.keys(results[0]);

  // Build CSV string
  const header = columns.join(',');
  const rows = results.map(row =>
    columns.map(col => {
      const val = row[col];
      if (typeof val === 'string') return `"${val}"`;
      if (typeof val === 'number') return val.toFixed(6);
      return String(val);
    }).join(',')
  );

  const csv = [header, ...rows].join('\n');
  const blob = new Blob([csv], { type: 'text/csv;charset=utf-8;' });
  const url = URL.createObjectURL(blob);
  const link = document.createElement('a');
  link.href = url;
  link.download = filename;
  link.click();
  URL.revokeObjectURL(url);
}

// Download ALL networks combined into a single CSV (with a Network column)
export function downloadAllNetworksCSV(allResults: SimulationResult[][], filename = 'ddm-results-all-networks.csv') {
  if (!allResults || allResults.length === 0) return;

  const firstRow = allResults[0]?.[0];
  if (!firstRow) return;

  const baseColumns = Object.keys(firstRow);
  const header = ['Network', ...baseColumns].join(',');

  const allRows: string[] = [];
  allResults.forEach((networkData, netIdx) => {
    networkData.forEach(row => {
      const values = baseColumns.map(col => {
        const val = row[col];
        if (typeof val === 'string') return `"${val}"`;
        if (typeof val === 'number') return val.toFixed(6);
        return String(val);
      });
      allRows.push([netIdx + 1, ...values].join(','));
    });
  });

  const csv = [header, ...allRows].join('\n');
  const blob = new Blob([csv], { type: 'text/csv;charset=utf-8;' });
  const url = URL.createObjectURL(blob);
  const link = document.createElement('a');
  link.href = url;
  link.download = filename;
  link.click();
  URL.revokeObjectURL(url);
}

// Download architecture/experiment as unified JSON format
// Both Network "Export" and Simulation "Save Experiment" now produce the same format
export function downloadArchitectureJSON(
  npes: NPE[],
  connections: Connection[],
  trials: Record<string, string[]>,
  contingencies: string[],
  hasITI: boolean[],
  filename = 'ddm-experiment.json',
  simParams?: {
    numNetworks?: number;
    thresholdPreset?: string;
    disc?: number;
    pupdate?: string;
    modelKind?: ModelKind;
  },
  lockedLayout?: Record<string, { x: number; y: number }> | null
) {
  const data: Record<string, unknown> = {
    _type: 'ddm-ui-experiment',
    _version: '3.0',
    exportedAt: new Date().toISOString(),
    npes,
    connections,
    trials,
    contingencies,
    hasITI,
  };

  // Include simulation parameters if provided
  if (simParams) {
    if (simParams.numNetworks !== undefined) data.numNetworks = simParams.numNetworks;
    if (simParams.thresholdPreset !== undefined) data.thresholdPreset = simParams.thresholdPreset;
    if (simParams.disc !== undefined) data.disc = simParams.disc;
    if (simParams.pupdate !== undefined) data.pupdate = simParams.pupdate;
    if (simParams.modelKind !== undefined) data.modelKind = simParams.modelKind;
  }

  // Include locked layout if available
  if (lockedLayout) {
    data.lockedLayout = lockedLayout;
  }

  const json = JSON.stringify(data, null, 2);
  const blob = new Blob([json], { type: 'application/json' });
  const url = URL.createObjectURL(blob);
  const link = document.createElement('a');
  link.href = url;
  link.download = filename;
  link.click();
  URL.revokeObjectURL(url);
}

// Valid NPE layers and types for validation
const VALID_LAYERS = [
  'US', 'PrimarySensory', 'AssociativeSensory', 'Hippocampal', 'AssociativeMotor', 'PrimaryMotor', 'Dopaminergic',
];
const VALID_NPE_TYPES = ['Excitatory', 'Inhibitory'];

// Deep validation for imported experiment files
function isRecord(value: unknown): value is Record<string, unknown> {
  return typeof value === 'object' && value !== null && !Array.isArray(value);
}

function validateExperimentData(data: unknown): string | null {
  if (!isRecord(data)) return 'Experiment file must contain a JSON object';
  try { assertDDMPayload(data); } catch (error) { return error instanceof Error ? error.message : 'Not a DDM experiment'; }

  // Validate NPEs array
  if (!Array.isArray(data.npes) || data.npes.length === 0) {
    return 'File must contain at least one NPE';
  }
  for (let i = 0; i < data.npes.length; i++) {
    const npe = data.npes[i];
    if (!isRecord(npe)) return `NPE at index ${i} must be an object`;
    if (!npe.name || typeof npe.name !== 'string') return `NPE at index ${i} is missing a valid name`;
    if (typeof npe.type !== 'string' || !VALID_NPE_TYPES.includes(npe.type)) {
      return `NPE "${npe.name}" has invalid type "${String(npe.type)}"`;
    }
    if (typeof npe.layer !== 'string' || !VALID_LAYERS.includes(npe.layer)) {
      return `NPE "${npe.name}" has invalid layer "${String(npe.layer)}"`;
    }
    if (typeof npe.mu !== 'number' || typeof npe.sigma !== 'number') return `NPE "${npe.name}" is missing mu/sigma values`;
  }

  // Check for duplicate NPE names
  const npeNames = new Set<string>();
  for (const npe of data.npes) {
    if (!isRecord(npe) || typeof npe.name !== 'string') continue;
    if (npeNames.has(npe.name)) return `Duplicate NPE name: "${npe.name}"`;
    npeNames.add(npe.name);
  }

  // Validate connections array
  if (!Array.isArray(data.connections)) return 'Connections must be an array';
  for (let i = 0; i < data.connections.length; i++) {
    const conn = data.connections[i];
    if (!isRecord(conn)) return `Connection at index ${i} must be an object`;
    if (typeof conn.presynapticNPE !== 'string' || typeof conn.postsynapticNPE !== 'string') {
      return `Connection at index ${i} is missing source or target`;
    }
    if (!npeNames.has(conn.presynapticNPE)) return `Connection references unknown NPE "${conn.presynapticNPE}"`;
    if (!npeNames.has(conn.postsynapticNPE)) return `Connection references unknown NPE "${conn.postsynapticNPE}"`;
    if (typeof conn.weight !== 'number') return `Connection ${conn.presynapticNPE}→${conn.postsynapticNPE} has invalid weight`;
  }

  // Validate trials (optional but if present must be valid)
  if (data.trials && isRecord(data.trials)) {
    for (const [name, timesteps] of Object.entries(data.trials)) {
      if (!Array.isArray(timesteps)) return `Trial "${name}" has invalid timesteps (expected array)`;
    }
  }

  // Validate contingencies (optional)
  if (data.contingencies && !Array.isArray(data.contingencies)) {
    return 'Contingencies must be an array of phase definition strings';
  }

  return null; // All valid
}

// Parse uploaded JSON — accepts both old architecture-only and new unified experiment format
export function parseArchitectureJSON(jsonString: string): {
  npes: NPE[];
  connections: Connection[];
  trials: Record<string, string[]>;
  contingencies: string[];
  hasITI: boolean[];
  // Optional simulation parameters (present in unified format)
  numNetworks?: number;
  thresholdPreset?: string;
  disc?: number;
  pupdate?: string;
  modelKind?: ModelKind;
  lockedLayout?: Record<string, { x: number; y: number }>;
  validationError?: string;
} | null {
  try {
    const data: unknown = JSON.parse(jsonString);
    if (!isRecord(data)) {
      throw new Error('Invalid experiment/architecture file: expected a JSON object');
    }
    if (!data.npes || !data.connections) {
      throw new Error('Invalid experiment/architecture file: missing npes or connections');
    }

    // Deep validation
    const validationError = validateExperimentData(data);
    if (validationError) {
      return {
        npes: [],
        connections: [],
        trials: {},
        contingencies: [],
        hasITI: [],
        validationError,
      };
    }

    return {
      npes: data.npes as NPE[],
      connections: data.connections as Connection[],
      trials: isRecord(data.trials) ? data.trials as Record<string, string[]> : {},
      contingencies: Array.isArray(data.contingencies) ? data.contingencies as string[] : [],
      hasITI: Array.isArray(data.hasITI) ? data.hasITI as boolean[] : [],
      numNetworks: typeof data.numNetworks === 'number' ? data.numNetworks : undefined,
      thresholdPreset: typeof data.thresholdPreset === 'string' ? data.thresholdPreset : undefined,
      disc: typeof data.disc === 'number' ? data.disc : undefined,
      pupdate: typeof data.pupdate === 'string' ? data.pupdate : undefined,
      modelKind: 'dtd',
      lockedLayout: isRecord(data.lockedLayout)
        ? data.lockedLayout as Record<string, { x: number; y: number }>
        : undefined,
    };
  } catch (err) {
    console.error('Failed to parse file:', err);
    return null;
  }
}
