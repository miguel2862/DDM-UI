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
export declare function createPhases(phases: string[], trials: Record<string, string[]>, seed?: number): Record<string, unknown>[];
