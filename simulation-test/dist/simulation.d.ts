import type { NPERow, ConnectionRow, SimulationOptions, SimulationOutputRow } from './types.js';
export declare function simulateDBP(npesData: NPERow[], connectionsData: ConnectionRow[], timeStepsRaw: Record<string, unknown>[], hasITI: boolean[], options?: SimulationOptions): SimulationOutputRow[];
