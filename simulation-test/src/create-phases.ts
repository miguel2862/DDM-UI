// ══════════════════════════════════════════════════════════════════════════════
// DDM-UI — Create.Phases() — Line-by-line port from R to TypeScript
// Original: api/simulation.R (Create.Phases function, lines 404-499)
// ══════════════════════════════════════════════════════════════════════════════

import { RNG } from './random.js';

/**
 * Creates the TimeSteps table from phase definitions and trial types.
 * Each row represents one timestep with stimulus activations and learning flag.
 *
 * @param phases   Array of phase definition strings, e.g.:
 *                 "training, Random, Training, 100, False"
 *                 "Entrenamiento, Random, SS/LL, 100-100, True, 30, 30, IEEn"
 * @param trials   Map of trial type name → array of timestep strings, e.g.:
 *                 { Training: ["US,0.00,S1,1.00,True", ...] }
 * @param seed     Optional seed for reproducibility (affects Random order & ITI sampling)
 */
export function createPhases(
  phases: string[],
  trials: Record<string, string[]>,
  seed?: number
): Record<string, unknown>[] {
  const rng = new RNG(seed);
  const timesteps: Record<string, unknown>[] = [];

  // Determine number of columns from first trial's first timestep
  const firstTrialKey = Object.keys(trials)[0];
  const firstRow = trials[firstTrialKey][0].split(',');
  const nCols = 3 + firstRow.length; // Phase, Trial, TimeStep + stimulus columns

  for (let i = 0; i < phases.length; i++) {
    let currentTs = 1;
    let currentTrial = 1;

    const currentPhase = phases[i].split(',').map(s => s.trim());
    const phaseName = currentPhase[0];
    const trialOrder = currentPhase[1].toLowerCase();
    const trialTypes = currentPhase[2].split('/').map(s => s.trim());
    const trialNumbers = currentPhase[3].split('-').map(s => parseInt(s.trim()));
    const hasIti = currentPhase[4].toLowerCase() === 'true';

    let minITI = 0;
    let maxITI = 0;
    let itiTimestep = '';

    if (hasIti) {
      minITI = parseInt(currentPhase[5].trim());
      maxITI = parseInt(currentPhase[6].trim());
      itiTimestep = currentPhase[7].trim();
    }

    // Helper to add one timestep row
    const addRow = (trialName: string, tsIdx: number) => {
      const values = trials[trialName][tsIdx].split(',').map(s => s.trim());
      const row: Record<string, unknown> = {};
      row['0'] = phaseName;
      row['1'] = currentTrial;
      row['2'] = currentTs;
      for (let c = 0; c < values.length; c++) {
        // Try to parse as number, otherwise keep as string
        const num = Number(values[c]);
        row[String(c + 3)] = isNaN(num) ? values[c] : num;
      }
      timesteps.push(row);
      currentTs++;
    };

    // Helper to add ITI timesteps
    const addITI = () => {
      if (!hasIti) {
        currentTs = 1;
        return;
      }
      if (!(itiTimestep in trials)) throw new Error('ITI time steps not in trials');
      const currentITI = minITI === maxITI ? minITI : rng.randInt(minITI, maxITI);
      for (let ts = 0; ts < currentITI; ts++) {
        const values = trials[itiTimestep][0].split(',').map(s => s.trim());
        const row: Record<string, unknown> = {};
        row['0'] = phaseName;
        row['1'] = currentTrial;
        row['2'] = currentTs;
        for (let c = 0; c < values.length; c++) {
          const num = Number(values[c]);
          row[String(c + 3)] = isNaN(num) ? values[c] : num;
        }
        timesteps.push(row);
        currentTs++;
      }
    };

    if (trialOrder === 'in bulk') {
      // R lines 426-446: each trial type in sequence
      for (let j = 0; j < trialTypes.length; j++) {
        for (let k = 0; k < trialNumbers[j]; k++) {
          addITI();
          if (!hasIti) currentTs = 1;
          for (let ts = 0; ts < trials[trialTypes[j]].length; ts++) {
            addRow(trialTypes[j], ts);
          }
          currentTrial++;
          currentTs = 1;
        }
      }
    } else if (trialOrder === 'alternated') {
      // R lines 447-467: alternate between trial types
      for (let j = 0; j < trialNumbers[0]; j++) {
        for (let k = 0; k < trialTypes.length; k++) {
          addITI();
          if (!hasIti) currentTs = 1;
          for (let ts = 0; ts < trials[trialTypes[k]].length; ts++) {
            addRow(trialTypes[k], ts);
          }
          currentTrial++;
          currentTs = 1;
        }
      }
    } else if (trialOrder === 'random') {
      // R lines 468-491: random order
      const trialSequence: number[] = [];
      for (let n = 0; n < trialNumbers.length; n++) {
        for (let r = 0; r < trialNumbers[n]; r++) {
          trialSequence.push(n);
        }
      }
      rng.shuffle(trialSequence);

      for (const j of trialSequence) {
        addITI();
        if (!hasIti) currentTs = 1;
        for (let ts = 0; ts < trials[trialTypes[j]].length; ts++) {
          addRow(trialTypes[j], ts);
        }
        currentTrial++;
        currentTs = 1;
      }
    } else {
      throw new Error(`Unknown trial order: "${trialOrder}"`);
    }
  }

  // No duplicate keys — simulation.ts reads values by position using Object.values()

  return timesteps;
}
