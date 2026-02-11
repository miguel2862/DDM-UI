// ══════════════════════════════════════════════════════════════════════════════
// DDM-UI — Run Extinction Template and display results
// Quick smoke test to verify the TS simulation produces reasonable output
// ══════════════════════════════════════════════════════════════════════════════
import { simulateDBP } from './simulation.js';
import { createPhases } from './create-phases.js';
import { getExtinctionTemplate } from './templates.js';
const template = getExtinctionTemplate();
console.log('═══════════════════════════════════════════════════════');
console.log('  DDM Simulation Engine — TypeScript Port');
console.log('  Template: Extinction (100 training + 100 extinction)');
console.log('═══════════════════════════════════════════════════════\n');
// 1. Create timesteps
console.log('Creating timesteps...');
const timeSteps = createPhases(template.contingencies, template.trials);
console.log(`  Total timestep rows: ${timeSteps.length}`);
console.log(`  First row:`, JSON.stringify(timeSteps[0]));
console.log(`  Last row:`, JSON.stringify(timeSteps[timeSteps.length - 1]));
// 2. Run simulation
console.log('\nRunning simulation...');
const start = performance.now();
const results = simulateDBP(template.npes, template.connections, timeSteps, template.hasITI, { threshold: 'gaussian', disc: 0.001, pupdate: 'async_random' });
const elapsed = performance.now() - start;
console.log(`  Done in ${elapsed.toFixed(1)}ms`);
console.log(`  Output rows: ${results.length}`);
// 3. Show summary
const phases = [...new Set(results.map(r => r.Phase))];
console.log(`  Phases: ${phases.join(', ')}`);
// 4. Show M.1 activation (the conditioned response) at key points
console.log('\n── M.1 (Conditioned Response) ──');
console.log('Training phase (trials 1, 25, 50, 75, 100):');
const trainingRows = results.filter(r => r.Phase === 'training');
for (const trialNum of [1, 25, 50, 75, 100]) {
    const trialRows = trainingRows.filter(r => Number(r.Trial) === trialNum);
    if (trialRows.length > 0) {
        const lastTs = trialRows[trialRows.length - 1];
        console.log(`  Trial ${trialNum}: M.1 = ${Number(lastTs['M.1'] || 0).toFixed(6)}`);
    }
}
console.log('\nExtinction phase (trials 1, 25, 50, 75, 100):');
const extinctionRows = results.filter(r => r.Phase === 'extinction');
for (const trialNum of [1, 25, 50, 75, 100]) {
    const trialRows = extinctionRows.filter(r => Number(r.Trial) === trialNum);
    if (trialRows.length > 0) {
        const lastTs = trialRows[trialRows.length - 1];
        console.log(`  Trial ${trialNum}: M.1 = ${Number(lastTs['M.1'] || 0).toFixed(6)}`);
    }
}
// 5. Show connection weights at end
console.log('\n── Connection Weights (end of training) ──');
const lastTraining = trainingRows[trainingRows.length - 1];
for (const key of Object.keys(lastTraining)) {
    if (key.includes('-')) {
        console.log(`  ${key} = ${Number(lastTraining[key] || 0).toFixed(6)}`);
    }
}
console.log('\n── Connection Weights (end of extinction) ──');
const lastExtinction = extinctionRows[extinctionRows.length - 1];
for (const key of Object.keys(lastExtinction)) {
    if (key.includes('-')) {
        console.log(`  ${key} = ${Number(lastExtinction[key] || 0).toFixed(6)}`);
    }
}
// 6. Show learning signals
console.log('\n── Learning Signals (last 5 training timesteps) ──');
const last5Training = trainingRows.slice(-5);
for (const row of last5Training) {
    console.log(`  ts=${row.TimeStep} dVTA=${Number(row.dVTA || 0).toFixed(6)} dH=${Number(row.dH || 0).toFixed(6)}`);
}
console.log('\n✅ Simulation completed successfully!');
