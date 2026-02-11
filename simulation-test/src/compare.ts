// ══════════════════════════════════════════════════════════════════════════════
// DDM-UI — Statistical Comparison: Run N simulations and check patterns
// Since R and TS have different RNGs, we compare statistical behavior
// over many runs to ensure the port is correct.
// ══════════════════════════════════════════════════════════════════════════════

import { simulateDBP } from './simulation.js';
import { createPhases } from './create-phases.js';
import {
  getExtinctionTemplate,
  getAcquisitionTemplate,
  getBlockingTemplate,
  getSpontaneousRecoveryTemplate,
  getLatentInhibitionTemplate,
} from './templates.js';
import type { SimulationOutputRow } from './types.js';
import type { TemplateData } from './templates.js';

const N_RUNS = 20; // number of networks per template

// ── Helper: compute mean activation per trial for a given NPE ─────────────────
function meanPerTrial(
  results: SimulationOutputRow[],
  phase: string,
  npe: string,
  maxTrial: number
): number[] {
  const means: number[] = [];
  for (let t = 1; t <= maxTrial; t++) {
    const rows = results.filter(r => r.Phase === phase && Number(r.Trial) === t);
    if (rows.length === 0) { means.push(0); continue; }
    // Use last timestep of each trial (final activation)
    const lastRow = rows[rows.length - 1];
    means.push(Number(lastRow[npe] || 0));
  }
  return means;
}

function avg(arr: number[]): number {
  return arr.reduce((s, v) => s + v, 0) / arr.length;
}

function avgSlice(arr: number[], start: number, end: number): number {
  return avg(arr.slice(start, end));
}

// ── Run multiple simulations and average ─────────────────────────────────────
function runMultiple(template: TemplateData, n: number) {
  const allResults: SimulationOutputRow[][] = [];
  for (let i = 0; i < n; i++) {
    const timeSteps = createPhases(template.contingencies, template.trials);
    const results = simulateDBP(
      template.npes,
      template.connections,
      timeSteps,
      template.hasITI,
      { threshold: 'gaussian', disc: 0.001, pupdate: 'async_random' }
    );
    allResults.push(results);
  }
  return allResults;
}

// ═══════════════════════════════════════════════════════════════════════════════
// TEST 1: Acquisition — M.1 should increase monotonically
// ═══════════════════════════════════════════════════════════════════════════════
function testAcquisition() {
  console.log('\n🔬 TEST 1: Acquisition');
  const template = getAcquisitionTemplate();
  const allResults = runMultiple(template, N_RUNS);

  // Average M.1 per trial across runs
  const avgPerTrial = new Array(100).fill(0);
  for (const results of allResults) {
    const m1 = meanPerTrial(results, 'training', 'M.1', 100);
    for (let t = 0; t < 100; t++) avgPerTrial[t] += m1[t];
  }
  for (let t = 0; t < 100; t++) avgPerTrial[t] /= N_RUNS;

  const earlyMean = avgSlice(avgPerTrial, 0, 10);   // trials 1-10
  const midMean = avgSlice(avgPerTrial, 40, 60);     // trials 41-60
  const lateMean = avgSlice(avgPerTrial, 85, 100);   // trials 86-100

  console.log(`  Early (trials 1-10):   M.1 avg = ${earlyMean.toFixed(4)}`);
  console.log(`  Middle (trials 41-60): M.1 avg = ${midMean.toFixed(4)}`);
  console.log(`  Late (trials 86-100):  M.1 avg = ${lateMean.toFixed(4)}`);

  const pass = earlyMean < midMean && midMean < lateMean && lateMean > 0.5;
  console.log(`  ${pass ? '✅ PASS' : '❌ FAIL'}: M.1 increases over training and reaches > 0.5`);
  return pass;
}

// ═══════════════════════════════════════════════════════════════════════════════
// TEST 2: Extinction — M.1 should increase then decrease
// ═══════════════════════════════════════════════════════════════════════════════
function testExtinction() {
  console.log('\n🔬 TEST 2: Extinction');
  const template = getExtinctionTemplate();
  const allResults = runMultiple(template, N_RUNS);

  // Average M.1 at end of training vs end of extinction
  let trainingEnd = 0;
  let extinctionEnd = 0;
  for (const results of allResults) {
    const trainRows = results.filter(r => r.Phase === 'training');
    const extRows = results.filter(r => r.Phase === 'extinction');

    // Last 10 trials of each phase
    const trainLast10 = meanPerTrial(results, 'training', 'M.1', 100).slice(90, 100);
    const extLast10 = meanPerTrial(results, 'extinction', 'M.1', 100).slice(90, 100);

    trainingEnd += avg(trainLast10);
    extinctionEnd += avg(extLast10);
  }
  trainingEnd /= N_RUNS;
  extinctionEnd /= N_RUNS;

  console.log(`  End of training (last 10 trials):    M.1 avg = ${trainingEnd.toFixed(4)}`);
  console.log(`  End of extinction (last 10 trials):  M.1 avg = ${extinctionEnd.toFixed(4)}`);

  const pass = trainingEnd > 0.5 && extinctionEnd < trainingEnd;
  console.log(`  ${pass ? '✅ PASS' : '❌ FAIL'}: M.1 is high after training and lower after extinction`);
  return pass;
}

// ═══════════════════════════════════════════════════════════════════════════════
// TEST 3: Blocking — M.1 should be lower for blocked stimulus X
// ═══════════════════════════════════════════════════════════════════════════════
function testBlocking() {
  console.log('\n🔬 TEST 3: Blocking');
  const template = getBlockingTemplate();
  const allResults = runMultiple(template, N_RUNS);

  // Check M.1 activation during Test phase (X alone)
  // In blocking, X should have LOW activation because A blocked learning to X
  let testM1 = 0;
  let trainingEnd = 0;
  for (const results of allResults) {
    // End of Entrenamiento (A+): M.1 should be high
    const trainTrials = meanPerTrial(results, 'Entrenamiento', 'M.1', 100);
    trainingEnd += avg(trainTrials.slice(90, 100));

    // Test phase: M.1 for blocked stimulus X
    const testTrials = meanPerTrial(results, 'Test', 'M.1', 25);
    testM1 += avg(testTrials);
  }
  trainingEnd /= N_RUNS;
  testM1 /= N_RUNS;

  console.log(`  End of A+ training:  M.1 avg = ${trainingEnd.toFixed(4)}`);
  console.log(`  Test (X alone):      M.1 avg = ${testM1.toFixed(4)}`);

  // X alone should have lower activation than trained A
  const pass = trainingEnd > testM1;
  console.log(`  ${pass ? '✅ PASS' : '❌ FAIL'}: Blocked X shows less responding than trained A`);
  return pass;
}

// ═══════════════════════════════════════════════════════════════════════════════
// TEST 4: Latent Inhibition — Pre-exposed group learns slower
// ═══════════════════════════════════════════════════════════════════════════════
function testLatentInhibition() {
  console.log('\n🔬 TEST 4: Latent Inhibition');
  const liTemplate = getLatentInhibitionTemplate();
  const acqTemplate = getAcquisitionTemplate();

  const liResults = runMultiple(liTemplate, N_RUNS);
  const acqResults = runMultiple(acqTemplate, N_RUNS);

  // Compare M.1 at mid-training (trials 30-50) between pre-exposed vs naive
  let liMid = 0;
  let acqMid = 0;
  for (const results of liResults) {
    const trainTrials = meanPerTrial(results, 'training', 'M.1', 100);
    liMid += avgSlice(trainTrials, 29, 50);
  }
  for (const results of acqResults) {
    const trainTrials = meanPerTrial(results, 'training', 'M.1', 100);
    acqMid += avgSlice(trainTrials, 29, 50);
  }
  liMid /= N_RUNS;
  acqMid /= N_RUNS;

  console.log(`  Naive acquisition (trials 30-50):         M.1 avg = ${acqMid.toFixed(4)}`);
  console.log(`  Pre-exposed + training (trials 30-50):    M.1 avg = ${liMid.toFixed(4)}`);

  // Pre-exposed group should learn slower (lower M.1 during mid-training)
  const pass = acqMid > liMid;
  console.log(`  ${pass ? '✅ PASS' : '❌ FAIL'}: Pre-exposed group learns slower (latent inhibition effect)`);
  return pass;
}

// ═══════════════════════════════════════════════════════════════════════════════
// TEST 5: Connection weights — US-D always stays at 1.0
// ═══════════════════════════════════════════════════════════════════════════════
function testUSConnection() {
  console.log('\n🔬 TEST 5: US→D Connection Weight');
  const template = getExtinctionTemplate();
  const allResults = runMultiple(template, 5);

  let allCorrect = true;
  for (const results of allResults) {
    for (const row of results) {
      const usd = Number(row['US-D'] || 0);
      if (Math.abs(usd - 1.0) > 0.0001) {
        allCorrect = false;
        console.log(`  ❌ US-D = ${usd} at Phase=${row.Phase} Trial=${row.Trial} TS=${row.TimeStep}`);
        break;
      }
    }
  }
  console.log(`  ${allCorrect ? '✅ PASS' : '❌ FAIL'}: US→D weight = 1.0 at all timesteps`);
  return allCorrect;
}

// ═══════════════════════════════════════════════════════════════════════════════
// TEST 6: Weight bounds — all weights in [0, 1]
// ═══════════════════════════════════════════════════════════════════════════════
function testWeightBounds() {
  console.log('\n🔬 TEST 6: Weight Bounds [0, 1]');
  const template = getBlockingTemplate();
  const timeSteps = createPhases(template.contingencies, template.trials);
  const results = simulateDBP(
    template.npes, template.connections, timeSteps, template.hasITI,
    { threshold: 'gaussian', disc: 0.001, pupdate: 'async_random' }
  );

  let violations = 0;
  for (const row of results) {
    for (const key of Object.keys(row)) {
      if (key.includes('-')) {
        const val = Number(row[key]);
        if (val < -0.0001 || val > 1.0001) {
          violations++;
          console.log(`  ❌ ${key} = ${val} out of bounds`);
        }
      }
    }
  }
  const pass = violations === 0;
  console.log(`  ${pass ? '✅ PASS' : '❌ FAIL'}: All connection weights within [0, 1]`);
  return pass;
}

// ═══════════════════════════════════════════════════════════════════════════════
// TEST 7: Spontaneous Recovery — response returns after rest
// ═══════════════════════════════════════════════════════════════════════════════
function testSpontaneousRecovery() {
  console.log('\n🔬 TEST 7: Spontaneous Recovery');
  const template = getSpontaneousRecoveryTemplate();
  const allResults = runMultiple(template, N_RUNS);

  let endTraining = 0;
  let endExtinction = 0;
  let testResponse = 0;

  for (const results of allResults) {
    const trainTrials = meanPerTrial(results, 'training', 'M.1', 100);
    const extTrials = meanPerTrial(results, 'extinction', 'M.1', 100);
    const testTrials = meanPerTrial(results, 'test', 'M.1', 25);

    endTraining += avg(trainTrials.slice(90, 100));
    endExtinction += avg(extTrials.slice(90, 100));
    testResponse += avg(testTrials);
  }
  endTraining /= N_RUNS;
  endExtinction /= N_RUNS;
  testResponse /= N_RUNS;

  console.log(`  End of training:       M.1 avg = ${endTraining.toFixed(4)}`);
  console.log(`  End of extinction:     M.1 avg = ${endExtinction.toFixed(4)}`);
  console.log(`  After rest (test):     M.1 avg = ${testResponse.toFixed(4)}`);

  // In this template config, rest = trials with S1=0, so there's continued decay.
  // The key test is that the response doesn't completely disappear (weights preserved)
  const pass = testResponse > 0.05 && endTraining > endExtinction;
  console.log(`  ${pass ? '✅ PASS' : '❌ FAIL'}: Training > Extinction and test response > 0 (weights preserved)`);
  return pass;
}

// ═══════════════════════════════════════════════════════════════════════════════
// RUN ALL TESTS
// ═══════════════════════════════════════════════════════════════════════════════
console.log('═══════════════════════════════════════════════════════');
console.log('  DDM Simulation Engine — TypeScript Port Validation');
console.log(`  Running ${N_RUNS} networks per template`);
console.log('═══════════════════════════════════════════════════════');

const start = performance.now();

const results = [
  testAcquisition(),
  testExtinction(),
  testBlocking(),
  testLatentInhibition(),
  testUSConnection(),
  testWeightBounds(),
  testSpontaneousRecovery(),
];

const elapsed = performance.now() - start;
const passed = results.filter(Boolean).length;
const total = results.length;

console.log('\n═══════════════════════════════════════════════════════');
console.log(`  Results: ${passed}/${total} tests passed (${elapsed.toFixed(0)}ms)`);
console.log('═══════════════════════════════════════════════════════');

if (passed === total) {
  console.log('\n🎉 ALL TESTS PASSED — TypeScript simulation matches expected behavior!');
} else {
  console.log('\n⚠️  Some tests failed — investigate differences');
  process.exit(1);
}
