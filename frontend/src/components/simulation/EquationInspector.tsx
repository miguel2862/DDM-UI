import { useMemo, useState } from 'react';
import { ChevronLeft, ChevronRight, Download } from 'lucide-react';
import type { InspectorUnit, SimulationInspector } from '../../types/inspector';
import { getDisplayName } from '../../utils/displayNames';
import './equationInspector.css';

interface Props {
  trace: SimulationInspector;
  rowIndex: number;
  onRowChange: (row: number) => void;
  selectedUnit?: string;
  onSelectUnit: (name: string) => void;
  language: string;
}

// Precision is only a display choice; neither operands nor results are rounded in the trace.
const number = (value: number | null | undefined) => value == null ? '—'
  : Number.isFinite(value) ? Number(value.toPrecision(6)).toString() : String(value);
const logistic = (x: number, sigma: number) => 1 / (1 + Math.exp((0.5 - x) / sigma));

function ActivationCurve({ unit, es }: { unit: InspectorUnit; es: boolean }) {
  if (unit.excInput == null || unit.inhInput == null || unit.logisSigma <= 0) return null;
  const low = Math.min(0, unit.excInput, unit.inhInput);
  const high = Math.max(1, unit.excInput, unit.inhInput);
  const x = (v: number) => 44 + (v - low) / (high - low) * 330;
  const y = (v: number) => 151 - v * 123;
  const curve = Array.from({ length: 121 }, (_, i) => {
    const input = low + (high - low) * i / 120;
    return `${i ? 'L' : 'M'}${x(input)},${y(logistic(input, unit.logisSigma))}`;
  }).join(' ');
  const thresholdVisible = unit.threshold != null && unit.threshold >= 0 && unit.threshold <= 1;
  return <figure className="equation-chart">
    <svg viewBox="0 0 410 198" role="img" aria-label={es ? 'Función logística y entradas realmente utilizadas' : 'Logistic function and actual recorded inputs'}>
      {[0, 0.5, 1].map(v => <g key={v}>
        <line x1="44" x2="374" y1={y(v)} y2={y(v)} stroke="#E2E8F0" />
        <text x="34" y={y(v) + 4} textAnchor="end">{v}</text>
      </g>)}
      <path d={curve} fill="none" stroke="#45535D" strokeWidth="2.3" />
      {thresholdVisible && <g>
        <line x1="44" x2="374" y1={y(unit.threshold!)} y2={y(unit.threshold!)} stroke="#B88817" strokeDasharray="5 4" />
        <text x="383" y={y(unit.threshold!) + 4} fill="#8A6510">θ</text>
      </g>}
      {[[unit.excInput, unit.logisticExc, '#0072B2'], [unit.inhInput, unit.logisticInh, '#8E679C']].map(([input, output, color], i) => output != null && <g key={i}>
        <line x1={x(input as number)} x2={x(input as number)} y1={151} y2={y(output as number)} stroke={color as string} strokeDasharray="2 3" />
        <circle cx={x(input as number)} cy={y(output as number)} r={i === 0 ? 5 : 3.3} fill={color as string} stroke="white" strokeWidth="1" />
      </g>)}
      <text x="44" y="171" textAnchor="middle">{number(low)}</text>
      <text x="374" y="171" textAnchor="middle">{number(high)}</text>
      <text x="209" y="191" textAnchor="middle">{es ? 'Entrada ponderada, x' : 'Weighted input, x'}</text>
      <text x="44" y="15">L(x)</text>
    </svg>
    <figcaption><span className="equation-exc">● E = L(e)</span><span className="equation-inh">● I = L(i)</span><span>θ = {number(unit.threshold)}{!thresholdVisible && (es ? ' (fuera de 0–1)' : ' (outside 0–1)')}</span></figcaption>
    <p className="equation-note">{es ? 'Curva analítica de la función, no una serie temporal. Los puntos son las entradas registradas en este timestep.' : 'Analytical function curve, not a time series. Points are the inputs recorded at this timestep.'}</p>
  </figure>;
}

function WeightCurve({ points, rowIndex, es }: { points: { row: number; value: number }[]; rowIndex: number; es: boolean }) {
  if (!points.length) return null;
  const maxRow = Math.max(1, points[points.length - 1].row);
  const low = Math.min(...points.map(p => p.value));
  const high = Math.max(...points.map(p => p.value));
  const pad = Math.max(0.005, (high - low) * 0.12);
  const bottom = Math.max(0, low - pad), top = Math.min(1, high + pad);
  const span = Math.max(0.000001, top - bottom);
  const x = (row: number) => 60 + row / maxRow * 314;
  const y = (value: number) => 151 - (value - bottom) / span * 123;
  const selected = points.find(p => p.row === rowIndex);
  return <figure className="equation-chart">
    <svg viewBox="0 0 410 198" role="img" aria-label={es ? 'Pesos registrados al final de cada timestep' : 'Weights recorded at the end of each timestep'}>
      {[bottom, (bottom + top) / 2, top].map((v, i) => <g key={i}>
        <line x1="60" x2="374" y1={y(v)} y2={y(v)} stroke="#E2E8F0" />
        <text x="52" y={y(v) + 4} textAnchor="end">{Number(v.toPrecision(3))}</text>
      </g>)}
      <path d={points.map((p, i) => `${i ? 'L' : 'M'}${x(p.row)},${y(p.value)}`).join(' ')} fill="none" stroke="#008D79" strokeWidth="2.3" />
      {selected && <g>
        <line x1={x(selected.row)} x2={x(selected.row)} y1="25" y2="151" stroke="#008D79" strokeDasharray="3 4" />
        <circle cx={x(selected.row)} cy={y(selected.value)} r="4.5" fill="#008D79" stroke="white" strokeWidth="1.5" />
      </g>}
      <text x="60" y="171" textAnchor="middle">1</text>
      <text x="374" y="171" textAnchor="middle">{maxRow + 1}</text>
      <text x="217" y="191" textAnchor="middle">{es ? 'Timestep acumulado (registro)' : 'Cumulative timestep (record)'}</text>
      <text x="60" y="15">w</text>
    </svg>
    <p className="equation-note">{es ? 'Peso al terminar cada timestep. El eje vertical se ajusta para hacer visibles cambios pequeños; no representa ensayos ni segundos.' : 'Weight at the end of each timestep. The vertical axis adapts to reveal small changes; the horizontal axis is neither trials nor seconds.'}</p>
  </figure>;
}

function ActivationEquation({ unit, es }: { unit: InspectorUnit; es: boolean }) {
  const previous = logistic(unit.previousExcitatoryInput, unit.logisSigma);
  if (unit.branch === 'external' || unit.branch === 'unconditional') return <div className="equation-formula">
    <p>{es ? 'Activación asignada por el protocolo o la entrada incondicional.' : 'Activation assigned by the protocol or unconditional input.'}</p>
    <strong>a = {number(unit.activation)}</strong>
    <p className="equation-note">{es ? 'Aquí no se evaluaron un umbral aleatorio ni la ecuación de activación ordinaria.' : 'No random threshold or ordinary activation equation was evaluated here.'}</p>
  </div>;
  return <div className="equation-formula">
    <p>L(x) = 1 / [1 + exp((0.5 − x) / σL)]</p>
    {unit.branch === 'suprathreshold' ? <>
      <p className="equation-branch">E &gt; I; E ≥ θ</p>
      <p>a = E + τ · L(e {es ? 'anterior' : 'previous'}) · (1 − E) − I</p>
      <p className="equation-substitution">{number(unit.logisticExc)} + {number(unit.temporalSummation)} × {number(previous)} × (1 − {number(unit.logisticExc)}) − {number(unit.logisticInh)}</p>
    </> : unit.branch === 'subthreshold' ? <>
      <p className="equation-branch">E &gt; I; E &lt; θ</p>
      <p>a = (1 − κ) · L(e {es ? 'anterior' : 'previous'})</p>
      <p className="equation-substitution">(1 − {number(unit.activationDecay)}) × {number(previous)}</p>
    </> : <><p className="equation-branch">E ≤ I</p><p>a = 0</p></>}
    <strong>a = {number(unit.activation)}</strong>
  </div>;
}

export function EquationInspector({ trace, rowIndex, onRowChange, selectedUnit, onSelectUnit, language }: Props) {
  const es = language === 'es';
  const [connectionName, setConnectionName] = useState('');
  const step = trace.steps.find(s => s.rowIndex === rowIndex);
  // Defaults must not jump when the engine's random traversal order changes.
  const firstUnits = trace.steps[0]?.units ?? [];
  const defaultUnit = firstUnits.find(u => u.layer === 'AssociativeSensory') ?? firstUnits.find(u => !['external', 'unconditional'].includes(u.branch)) ?? firstUnits[0];
  const unit = step?.units.find(u => u.name === (selectedUnit ?? defaultUnit?.name));
  const stableConnections = [...(step?.connections ?? [])].sort((a, b) => a.name.localeCompare(b.name));
  const connection = stableConnections.find(c => c.name === connectionName) ?? stableConnections.find(c => c.post === unit?.name && c.branch !== 'fixedUS') ?? stableConnections[0];
  const selectedConnectionName = connection?.name;
  const trajectories = useMemo(() => {
    const result = new Map<string, { row: number; value: number }[]>();
    for (const s of trace.steps) for (const c of s.connections) {
      const series = result.get(c.name) ?? [];
      series.push({ row: s.rowIndex, value: c.weightAfter });
      result.set(c.name, series);
    }
    return result;
  }, [trace]);
  const points = trajectories.get(selectedConnectionName ?? '') ?? [];
  const trialStarts = useMemo(() => trace.steps.filter((s, i, all) => i === 0 || s.phase !== all[i - 1].phase || s.trial !== all[i - 1].trial), [trace]);
  const currentTrialIndex = step ? trialStarts.findIndex(s => s.phase === step.phase && s.trial === step.trial) : -1;
  const previousTrial = trialStarts[currentTrialIndex - 1];
  const nextTrial = trialStarts[currentTrialIndex + 1];
  const selectUnit = (name: string) => { onSelectUnit(name); setConnectionName(''); };
  const downloadTrace = () => {
    const url = URL.createObjectURL(new Blob([JSON.stringify(trace)], { type: 'application/json' }));
    const a = document.createElement('a');
    a.href = url; a.download = 'DDM-inspector-network-1.json'; a.click();
    window.setTimeout(() => URL.revokeObjectURL(url), 1000);
  };

  return <section className="equation-inspector" aria-labelledby="equation-inspector-title">
    <header className="equation-heading">
      <div><h2 id="equation-inspector-title">{es ? 'Dentro de la actualización' : 'Inside the update'}</h2>
        <p>{es ? 'Red 1 · Operaciones registradas por el motor DDM' : 'Network 1 · Operations recorded by the DDM engine'}</p></div>
      <button type="button" onClick={downloadTrace}><Download size={15} />{es ? 'Guardar registro' : 'Save trace'}</button>
    </header>
    {trace.truncated && <p className="equation-notice">{es ? `Registro limitado a ${trace.recordedTimesteps} de ${trace.totalTimesteps} timesteps. La simulación sí se ejecutó completa.` : `Trace limited to ${trace.recordedTimesteps} of ${trace.totalTimesteps} timesteps. The full simulation still ran.`}</p>}
    {!step ? <div className="equation-empty"><p>{es ? 'Este timestep no está incluido en el registro de ecuaciones.' : 'This timestep is not included in the equation trace.'}</p>
      {trace.steps.length > 0 && <button type="button" onClick={() => onRowChange(trace.steps[trace.steps.length - 1].rowIndex)}>{es ? 'Ir al último timestep registrado' : 'Go to the last recorded timestep'}</button>}</div> : <>
      <div className="equation-navigation">
        <div className="equation-trial-buttons">
          <button type="button" disabled={!previousTrial} onClick={() => previousTrial && onRowChange(previousTrial.rowIndex)}><ChevronLeft size={15} />{es ? 'Ensayo anterior' : 'Previous trial'}</button>
          <button type="button" disabled={!nextTrial} onClick={() => nextTrial && onRowChange(nextTrial.rowIndex)}>{es ? 'Ensayo siguiente' : 'Next trial'}<ChevronRight size={15} /></button>
        </div>
        <p>{step.phase} · {es ? 'Ensayo' : 'Trial'} {step.trial} · timestep {step.timestep}</p>
      </div>
      <p className="equation-note equation-context">{step.resetApplied
        ? (es ? 'Inicio de ensayo sin ITI: se reiniciaron activaciones y entradas según el motor; los pesos se conservaron.' : 'Trial start without ITI: activations and inputs were reset by the engine; weights were retained.')
        : (es ? 'Sin reinicio al inicio de este timestep.' : 'No reset at the start of this timestep.')} {es ? 'La red superior muestra el estado final; los operandos de abajo corresponden al orden real de actualización.' : 'The network above shows the final state; the operands below follow the actual update order.'}</p>
      <div className="equation-columns">
        <article>
          <label htmlFor="inspector-unit">{es ? 'Activación de la unidad' : 'Unit activation'}</label>
          <select id="inspector-unit" value={unit?.name ?? ''} onChange={e => selectUnit(e.target.value)}>{[...step.units].sort((a, b) => a.name.localeCompare(b.name)).map(u => <option key={u.name} value={u.name}>{getDisplayName(u.name)} · {u.layer}</option>)}</select>
          {unit && <>
            <dl className="equation-values">
              <div><dt>{es ? 'Orden de activación' : 'Activation order'}</dt><dd>{unit.order} / {step.units.length}</dd></div>
              <div><dt>{es ? 'a anterior → actual' : 'Previous → current a'}</dt><dd>{number(unit.previousActivation)} → {number(unit.activation)}</dd></div>
              <div><dt>e / i</dt><dd>{number(unit.excInput)} / {number(unit.inhInput)}</dd></div>
              <div><dt>E / I</dt><dd>{number(unit.logisticExc)} / {number(unit.logisticInh)}</dd></div>
              <div><dt>θ / σL</dt><dd>{number(unit.threshold)} / {number(unit.logisSigma)}</dd></div>
              <div><dt>τ / κ</dt><dd>{number(unit.temporalSummation)} / {number(unit.activationDecay)}</dd></div>
            </dl>
            <ActivationEquation unit={unit} es={es} />
            <ActivationCurve unit={unit} es={es} />
          </>}
        </article>
        <article>
          <label htmlFor="inspector-connection">{es ? 'Aprendizaje de la conexión' : 'Connection learning'}</label>
          <select id="inspector-connection" value={connection?.name ?? ''} onChange={e => setConnectionName(e.target.value)} disabled={!step.connections.length}>{stableConnections.map(c => <option key={c.name} value={c.name}>{getDisplayName(c.pre)} → {getDisplayName(c.post)}</option>)}</select>
          {!connection ? <p className="equation-empty">{es ? 'No hay conexiones registradas en este timestep.' : 'No connections recorded at this timestep.'}</p> : <>
            <dl className="equation-values">
              <div><dt>{es ? 'Orden de conexión' : 'Connection order'}</dt><dd>{connection.order} / {step.connections.length}</dd></div>
              <div><dt>w {es ? 'antes → después' : 'before → after'}</dt><dd>{number(connection.weightBefore)} → {number(connection.weightAfter)}</dd></div>
              <div><dt>a pre / a post</dt><dd>{number(connection.preActivation)} / {number(connection.postActivation)}</dd></div>
              <div><dt>{connection.signalKind ?? (es ? 'Señal no evaluada' : 'Signal not evaluated')}</dt><dd>{number(connection.signal)}</dd></div>
              <div><dt>r / p</dt><dd>{number(connection.capacity)} / {number(connection.proportion)}</dd></div>
              <div><dt>{es ? 'Discrepancia' : 'Discrepancy'} / {es ? 'denominador' : 'denominator'}</dt><dd>{number(connection.disc)} / {number(connection.inputDenominator)}</dd></div>
            </dl>
            <div className="equation-formula">
              {connection.branch === 'potentiation' ? <>
                <p className="equation-branch">d ≥ {number(connection.disc)} · {es ? 'Incremento' : 'Potentiation'}</p>
                <p>w̃ = w + α{connection.preType === 'Inhibitory' ? '′' : ''} · r · a post · d · p</p>
                <p className="equation-substitution">{number(connection.weightBefore)} + {number(connection.alpha)} × {number(connection.capacity)} × {number(connection.postActivation)} × {number(connection.signal)} × {number(connection.proportion)}</p>
              </> : connection.branch === 'decrement' ? <>
                <p className="equation-branch">d &lt; {number(connection.disc)} · {es ? 'Decremento' : 'Decrement'}</p>
                <p>w̃ = w − β{connection.preType === 'Inhibitory' ? '′' : ''} · w · a pre · a post</p>
                <p className="equation-substitution">{number(connection.weightBefore)} − {number(connection.beta)} × {number(connection.weightBefore)} × {number(connection.preActivation)} × {number(connection.postActivation)}</p>
              </> : <p>{connection.branch === 'fixedUS' ? (es ? 'Conexión desde S*: peso fijo, sin regla de aprendizaje.' : 'Connection from S*: fixed weight, no learning rule.') : (es ? 'Aprendizaje desactivado: el peso se conserva.' : 'Learning disabled: the weight is retained.')}</p>}
              {connection.branch !== 'learningOff' && <p>w̃ = {number(connection.weightUnclipped)}; w {es ? 'nuevo' : 'new'} = min(1, max(0, w̃))</p>}
              <strong>w = {number(connection.weightAfter)} <span>· Δw = {number(connection.deltaWeight)}</span></strong>
            </div>
            <WeightCurve points={points} rowIndex={rowIndex} es={es} />
          </>}
        </article>
      </div>
      <details className="equation-details">
        <summary>{es ? 'Señales moduladoras y orden de actualización' : 'Modulatory signals and update order'}</summary>
        <p>dD = {number(step.dD)} · dH = {number(step.dH)} · dH {es ? 'anterior' : 'previous'} = {number(step.previousDH)}</p>
        {step.learningEnabled ? <><p>dD = mean(Δa D); dH = mean(|Δa H|) + dD · (1 − dH {es ? 'anterior' : 'previous'})</p><p className="equation-note">{es ? 'r es la capacidad restante de la clase de entrada correspondiente; p = a pre · w / entrada total de esa clase. El motor conserva sus reglas y el momento en que calcula cada operando.' : 'r is the remaining capacity of the relevant input class; p = a pre · w / total input of that class. The engine retains its rules and the timing of each operand.'}</p></> : <p>{es ? 'Aprendizaje desactivado: dD y dH no se calcularon en este timestep.' : 'Learning disabled: dD and dH were not calculated at this timestep.'}</p>}
        <p><b>{es ? 'Activación' : 'Activation'}:</b> {step.activationOrder.map(getDisplayName).join(' → ')}</p>
        <p><b>{es ? 'Aprendizaje (unidades)' : 'Learning (units)'}:</b> {step.learningOrder.length ? step.learningOrder.map(getDisplayName).join(' → ') : '—'}</p>
      </details>
      <footer className="equation-note">{es ? 'Seis cifras significativas en pantalla; precisión completa en el registro JSON. Navegar no ejecuta de nuevo el modelo ni sortea nuevos valores.' : 'Six significant digits on screen; full precision in the JSON trace. Navigation neither reruns the model nor draws new random values.'}</footer>
    </>}
  </section>;
}
