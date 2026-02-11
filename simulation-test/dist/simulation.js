// ══════════════════════════════════════════════════════════════════════════════
// DDM-UI — Simulate.DBP() — Line-by-line port from R to TypeScript
// Original: api/simulation.R (Simulate.DBP function, lines 4-401)
// ══════════════════════════════════════════════════════════════════════════════
import { RNG } from './random.js';
// ── Helper: Logistic function (R line 85-87) ────────────────────────────────
function L(x, sigma) {
    return 1 / (1 + Math.exp((-x + 0.5) / sigma));
}
// ── Helper: Estimate Beta distribution parameters (R line 89-95) ─────────────
function estBetaParams(mu, sigma) {
    const variance = sigma * sigma;
    const a = mu * (mu - mu * mu - variance) / variance;
    const b = (mu * (1 - mu) / variance - 1) * (1 - mu);
    if (a <= 0 || b <= 0) {
        throw new Error('This combination of mean and standard deviation for the threshold results in invalid parameters for a beta distribution');
    }
    return { alpha: a, beta: b };
}
// ═══════════════════════════════════════════════════════════════════════════════
// Simulate.DBP — Main simulation function
// ═══════════════════════════════════════════════════════════════════════════════
export function simulateDBP(npesData, connectionsData, timeStepsRaw, hasITI, options = {}) {
    const { threshold = 'gaussian', disc = 0.001, pupdate = 'async_random', seed, } = options;
    const rng = new RNG(seed);
    // ── Validate inputs (R lines 47-68) ─────────────────────────────────────────
    if (npesData.length === 0)
        throw new Error('NPEs array is empty');
    if (connectionsData.length === 0)
        throw new Error('Connections array is empty');
    // ── Initialize network (R lines 140-172) ────────────────────────────────────
    const network = {};
    const networkOrder = []; // preserve insertion order
    for (const row of npesData) {
        const name = String(row.NPE);
        networkOrder.push(name);
        network[name] = {
            Name: name,
            Type: row.Type,
            Layer: row.Layer,
            Activation: Number(row.Activation),
            TemporalSummation: Number(row['Temporal.Summation']),
            ActivationDecay: Number(row['Activation.Decay']),
            mu: Number(row.mu),
            sigma: Number(row.sigma),
            logisSigma: Number(row.logisSigma),
            PreviousExcitatoryInput: 0,
            ExcitatoryInput: 0,
            InhibitoryInput: 0,
            PreviousActivation: 0,
            Threshold: 0,
            InputConnections: {},
            r: [0, 0],
        };
    }
    for (const row of connectionsData) {
        const pre = String(row.PreSinapticNPE);
        const post = String(row.PostSinapticNPE);
        const connName = `${pre}-${post}`;
        network[post].InputConnections[connName] = {
            Name: connName,
            weight: Number(row.Weight),
            alpha: Number(row.alpha),
            beta: Number(row.beta),
            alpha_prime: Number(row.alpha_prime),
            beta_prime: Number(row.beta_prime),
            preSinapticNPE: pre,
            p: 0,
        };
    }
    const nLen = networkOrder.length;
    // ── Determine input layers (R line 174) ─────────────────────────────────────
    const isInput = {};
    for (const name of networkOrder) {
        isInput[name] = network[name].Layer === 'PrimarySensory' || network[name].Layer === 'US';
    }
    // ── Build saveData (R lines 176-192) ────────────────────────────────────────
    const saveElements = new Set();
    if (options.saveData?.Elements) {
        for (const el of options.saveData.Elements)
            saveElements.add(el);
    }
    else {
        // Default: all non-input NPEs + all connections
        for (const name of networkOrder) {
            if (network[name].Layer !== 'PrimarySensory' && network[name].Layer !== 'US') {
                saveElements.add(name);
            }
            for (const connName of Object.keys(network[name].InputConnections)) {
                saveElements.add(connName);
            }
        }
    }
    // Always include input NPEs (R line 192)
    for (const name of networkOrder) {
        if (isInput[name])
            saveElements.add(name);
    }
    const timeSteps = [];
    for (const row of timeStepsRaw) {
        const values = Object.values(row);
        const phase = String(values[0]);
        const trial = Number(values[1]);
        const ts = Number(values[2]);
        const stimuli = [];
        // Last column is LearningActive (R line 228)
        const lastVal = values[values.length - 1];
        // Pairs from column 4 to second-to-last, stepping by 2 (R: seq(4, ncol(TimeSteps)-1, 2))
        // The last value is the learning flag, pairs are everything between col 3 and the last
        for (let c = 3; c < values.length - 1; c += 2) {
            stimuli.push({
                npe: String(values[c]),
                activation: Number(values[c + 1]),
            });
        }
        const learningActive = lastVal === true || lastVal === 'True' || lastVal === 'TRUE' || lastVal === 1 || lastVal === '1';
        timeSteps.push({ Phase: phase, Trial: trial, TimeStep: ts, stimuli, learningActive });
    }
    // ── Determine which timesteps to save ───────────────────────────────────────
    const saveTimeSteps = new Set();
    if (options.saveData?.TimeSteps) {
        for (const ts of options.saveData.TimeSteps)
            saveTimeSteps.add(ts);
    }
    else {
        for (const row of timeSteps)
            saveTimeSteps.add(row.TimeStep);
    }
    // ── Allocate output data (R lines 194-195) ──────────────────────────────────
    const dataSim = [];
    const saveElementsList = Array.from(saveElements);
    // ── ComputeInputs (R lines 70-83) ──────────────────────────────────────────
    function computeInputs(npe) {
        let exc = 0;
        let inh = 0;
        for (const connName in npe.InputConnections) {
            const conn = npe.InputConnections[connName];
            const pre = network[conn.preSinapticNPE];
            if (pre.Type === 'Excitatory') {
                exc += pre.Activation * conn.weight;
            }
            else {
                inh += pre.Activation * conn.weight;
            }
        }
        return [exc, inh];
    }
    // ── dVTA (R lines 97-108) ──────────────────────────────────────────────────
    function dVTA() {
        let d = 0;
        let n = 0;
        for (const name of networkOrder) {
            if (network[name].Layer === 'Dopaminergic') {
                d += network[name].Activation - network[name].PreviousActivation;
                n++;
            }
        }
        return n === 0 ? 0 : d / n;
    }
    // ── dCA1 (R lines 110-122) ─────────────────────────────────────────────────
    function dCA1(dVTAVal, previousdCA1) {
        let d = 0;
        let n = 0;
        for (const name of networkOrder) {
            if (network[name].Layer === 'Hippocampal') {
                d += Math.abs(network[name].Activation - network[name].PreviousActivation);
                n++;
            }
        }
        const dH = n === 0 ? 0 : d / n;
        return dH + dVTAVal * (1 - previousdCA1);
    }
    // ── Compute.r (R lines 124-138) ────────────────────────────────────────────
    function computeR(npeName) {
        let sumWeightsExc = 0;
        let sumWeightsInh = 0;
        for (const connName in network[npeName].InputConnections) {
            const conn = network[npeName].InputConnections[connName];
            const pre = network[conn.preSinapticNPE];
            if (pre.Type === 'Excitatory') {
                sumWeightsExc += pre.Layer === 'US' ? 0 : conn.weight;
            }
            else {
                sumWeightsInh += conn.weight;
            }
        }
        return [1 - sumWeightsExc, 1 - sumWeightsInh];
    }
    // ── Main simulation loop (R lines 197-396) ─────────────────────────────────
    let learningRuleIsActive = false;
    let previousdCA1 = 0;
    let currentPhase = 0; // 0-indexed (R uses 1-indexed)
    for (let ts = 0; ts < timeSteps.length; ts++) {
        const row = timeSteps[ts];
        const shouldSave = saveTimeSteps.has(row.TimeStep);
        // Initialize output row if saving (R lines 204-209)
        let outputRow = null;
        if (shouldSave) {
            outputRow = {
                Phase: row.Phase,
                Trial: row.Trial,
                TimeStep: row.TimeStep,
            };
        }
        // Phase transition (R lines 211-213)
        if (ts > 0 && timeSteps[ts - 1].Phase !== row.Phase) {
            currentPhase++;
        }
        // Reset activations at trial start if no ITI (R lines 215-220)
        if (row.TimeStep === 1 && !hasITI[currentPhase]) {
            for (const name of networkOrder) {
                network[name].Activation = L(0, network[name].logisSigma);
                network[name].ExcitatoryInput = L(0, network[name].logisSigma);
            }
        }
        // Zero out input layer activations (R lines 222-226)
        for (const name of networkOrder) {
            if (network[name].Layer === 'US' || network[name].Layer === 'PrimarySensory') {
                network[name].Activation = 0;
            }
        }
        // Set learning rule active flag (R line 228)
        learningRuleIsActive = row.learningActive;
        // Apply stimuli (R lines 230-232)
        for (const stim of row.stimuli) {
            network[stim.npe].Activation = stim.activation;
        }
        // Scramble NPE order for activation computation (R line 234)
        const indices = Array.from({ length: nLen }, (_, i) => i);
        const scrambledNPEs = rng.shuffled(indices);
        // ── Activation phase (R lines 236-294) ───────────────────────────────────
        for (const i of scrambledNPEs) {
            const npeName = networkOrder[i];
            const npe = network[npeName];
            npe.PreviousActivation = npe.Activation;
            npe.PreviousExcitatoryInput = npe.ExcitatoryInput;
            // Skip input layers (R lines 241-244)
            if (npe.Layer === 'US' || npe.Layer === 'PrimarySensory') {
                if (shouldSave && saveElements.has(npeName)) {
                    outputRow[npeName] = npe.Activation;
                }
                continue;
            }
            // Check unconditional activation for Dopaminergic/PrimaryMotor (R lines 246-258)
            let npeIsUnconditionallyActivated = false;
            let usActivation = 0;
            if ((npe.Layer === 'Dopaminergic' || npe.Layer === 'PrimaryMotor') &&
                Object.keys(npe.InputConnections).length > 0) {
                for (const connName in npe.InputConnections) {
                    const conn = npe.InputConnections[connName];
                    const pre = network[conn.preSinapticNPE];
                    if (pre.Layer === 'US' && pre.Activation > 0) {
                        npeIsUnconditionallyActivated = true;
                        usActivation = pre.Activation;
                        break;
                    }
                }
            }
            if (npeIsUnconditionallyActivated) {
                // R line 261
                npe.Activation = usActivation;
            }
            else {
                // Compute threshold (R lines 263-268)
                if (threshold === 'gaussian') {
                    npe.Threshold = rng.rnorm(npe.mu, npe.sigma);
                }
                else {
                    const p = estBetaParams(npe.mu, npe.sigma);
                    npe.Threshold = rng.rbeta(p.alpha, p.beta);
                }
                // Compute inputs (R lines 270-272)
                const inputs = computeInputs(npe);
                npe.ExcitatoryInput = inputs[0];
                npe.InhibitoryInput = inputs[1];
                const pEpsp = L(npe.ExcitatoryInput, npe.logisSigma);
                const pIpsp = L(npe.InhibitoryInput, npe.logisSigma);
                // Activation rule (R lines 277-287)
                // Note: PrimarySensory/US already skipped via continue above
                if (pEpsp > pIpsp) {
                    if (pEpsp >= npe.Threshold) {
                        // R line 280
                        npe.Activation =
                            pEpsp +
                                npe.TemporalSummation * L(npe.PreviousExcitatoryInput, npe.logisSigma) * (1 - pEpsp) -
                                pIpsp;
                    }
                    else {
                        // R line 282
                        npe.Activation =
                            L(npe.PreviousExcitatoryInput, npe.logisSigma) -
                                npe.ActivationDecay * L(npe.PreviousExcitatoryInput, npe.logisSigma);
                    }
                }
                else {
                    npe.Activation = 0;
                }
            }
            // Save activation (R lines 290-293)
            if (shouldSave && saveElements.has(npeName)) {
                outputRow[npeName] = npe.Activation;
            }
        }
        // Initialize learning signals (R lines 296-300)
        if (shouldSave) {
            outputRow['dVTA'] = 0;
            outputRow['dH'] = 0;
        }
        // ── Learning phase (R lines 302-395) ─────────────────────────────────────
        if (learningRuleIsActive) {
            const dD = dVTA();
            const dH = dCA1(dD, previousdCA1);
            previousdCA1 = dH;
            // Save learning signals (R lines 308-311)
            if (shouldSave) {
                outputRow['dVTA'] = dD;
                outputRow['dH'] = dH;
            }
            // PUpdate procedure: order NPEs (R lines 313-318)
            let learningNPEOrder;
            if (pupdate === 'async_random' || pupdate === 'sync_random') {
                learningNPEOrder = rng.shuffled(indices);
            }
            else {
                learningNPEOrder = [...indices];
            }
            for (const i of learningNPEOrder) {
                const npeName = networkOrder[i];
                const npe = network[npeName];
                // Recompute inputs (R lines 322-324)
                const inputs = computeInputs(npe);
                npe.ExcitatoryInput = inputs[0];
                npe.InhibitoryInput = inputs[1];
                // Skip input layers (R lines 326-328)
                if (npe.Layer === 'US' || npe.Layer === 'PrimarySensory') {
                    continue;
                }
                // Compute r (R line 330)
                npe.r = computeR(npeName);
                // PUpdate procedure: order connections (R lines 332-337)
                const connNames = Object.keys(npe.InputConnections);
                let connectionOrder;
                if (pupdate === 'async_random' || pupdate === 'sync_random') {
                    connectionOrder = rng.shuffled(Array.from({ length: connNames.length }, (_, k) => k));
                }
                else {
                    connectionOrder = Array.from({ length: connNames.length }, (_, k) => k);
                }
                // Select discrepancy signal (R lines 339-343)
                let d;
                if (npe.Layer === 'AssociativeSensory' || npe.Layer === 'Hippocampal') {
                    d = dH;
                }
                else {
                    d = dD;
                }
                // Connection weight update (R lines 345-379)
                for (const j of connectionOrder) {
                    const connName = connNames[j];
                    const conn = npe.InputConnections[connName];
                    const pre = network[conn.preSinapticNPE];
                    if (pre.Layer !== 'US') {
                        if (d >= disc) {
                            // Reinforcement (R lines 350-360)
                            if (pre.Type === 'Excitatory') {
                                conn.p =
                                    npe.ExcitatoryInput === 0
                                        ? 0
                                        : (pre.Activation * conn.weight) / npe.ExcitatoryInput;
                                conn.weight =
                                    conn.weight +
                                        conn.alpha * npe.r[0] * npe.Activation * d * conn.p;
                            }
                            else {
                                conn.p =
                                    npe.InhibitoryInput === 0
                                        ? 0
                                        : (pre.Activation * conn.weight) / npe.InhibitoryInput;
                                conn.weight =
                                    conn.weight +
                                        conn.alpha_prime * npe.r[1] * npe.Activation * d * conn.p;
                            }
                        }
                        else {
                            // Extinction (R lines 362-370)
                            const decayRate = pre.Type === 'Excitatory' ? conn.beta : conn.beta_prime;
                            conn.weight =
                                conn.weight - decayRate * conn.weight * pre.Activation * npe.Activation;
                        }
                    }
                    // Clamp weight to [0, 1] (R line 373)
                    conn.weight = Math.min(Math.max(conn.weight, 0), 1);
                    // Save connection weight (R lines 375-378)
                    if (shouldSave && saveElements.has(connName)) {
                        outputRow[connName] = conn.weight;
                    }
                }
            }
        }
        else {
            // No learning — just save connection weights (R lines 382-394)
            for (const i of scrambledNPEs) {
                const npeName = networkOrder[i];
                const npe = network[npeName];
                const connNames = Object.keys(npe.InputConnections);
                let connectionOrder;
                if (pupdate === 'async_random' || pupdate === 'sync_random') {
                    connectionOrder = rng.shuffled(Array.from({ length: connNames.length }, (_, k) => k));
                }
                else {
                    connectionOrder = Array.from({ length: connNames.length }, (_, k) => k);
                }
                for (const j of connectionOrder) {
                    const connName = connNames[j];
                    const conn = npe.InputConnections[connName];
                    if (shouldSave && saveElements.has(connName)) {
                        outputRow[connName] = conn.weight;
                    }
                }
            }
        }
        // Push output row (R equivalent of data.sim[t,] assignment)
        if (shouldSave && outputRow) {
            dataSim.push(outputRow);
        }
    }
    return dataSim;
}
