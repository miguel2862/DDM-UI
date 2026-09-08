import { useState, useCallback, useMemo, useRef } from 'react';
import {
  LineChart, Line, XAxis, YAxis, CartesianGrid, Tooltip as RechartsTooltip, Legend, ResponsiveContainer,
} from 'recharts';
import { Sliders, Play, Loader2, TrendingUp, Info, Square, Download } from 'lucide-react';
import { Card } from '../components/ui/Card';
import { Badge } from '../components/ui/Badge';
import { Tooltip } from '../components/ui/Tooltip';
import { useSimStore } from '../stores/useSimStore';
import { runSimulation } from '../api/client';
import { useI18n } from '../i18n';
import type { SimulationResult } from '../types/ddm';

const SWEEP_COLORS = [
  '#06b6d4', '#14b8a6', '#f59e0b', '#8b5cf6', '#ef4444',
  '#ec4899', '#10b981', '#3b82f6', '#f97316', '#a855f7',
];

interface SweepResult {
  paramValue: number;
  meanActivation: number;
  medianActivation: number;
}



export function ParameterSweep() {
  const { npes, connections, trials, contingencies, hasITI, getThresholdType, disc, modelKind } = useSimStore();
  const { t } = useI18n();

  const [sweepTarget, setSweepTarget] = useState<'connection' | 'npe'>(('connection'));
  const [selectedElement, setSelectedElement] = useState((''));
  const [selectedParam, setSelectedParam] = useState(('weight'));
  const [rangeMin, setRangeMin] = useState(0);
  const [rangeMax, setRangeMax] = useState(1);
  const [steps, setSteps] = useState(10);
  const [numNetworks, setNumNetworks] = useState(3);
  const [outputUnit, setOutputUnit] = useState((''));
  const [isRunning, setIsRunning] = useState(false);
  const [progress, setProgress] = useState(0);
  const [results, setResults] = useState<SweepResult[]>([]);
  const cancelRef = useRef(false);

  const canRun = npes.length >= 2 && connections.length >= 1 &&
    Object.keys(trials).length >= 1 && contingencies.length >= 1 &&
    selectedElement && outputUnit;

  const connectionParams = ['weight', 'alpha', 'beta', 'alpha_prime', 'beta_prime'];
  const npeParams = ['mu', 'sigma', 'temporalSummation', 'activationDecay', 'logisSigma'];


  const motorUnits = useMemo(() => (npes.filter(n => n.layer !== 'PrimarySensory' && n.layer !== 'US').map(n => n.name)),
    [npes]
  );

  const handleRunSweep = useCallback(async () => {
    if (!canRun) return;
    setIsRunning(true);
    setResults([]);
    setProgress(0);
    cancelRef.current = false;

    const stepSize = (rangeMax - rangeMin) / steps;
    const sweepResults: SweepResult[] = [];

    for (let i = 0; i <= steps; i++) {
      if (cancelRef.current) break;
      const paramValue = rangeMin + i * stepSize;
      setProgress(Math.round((i / steps) * 100));

      let modifiedConns = connections.map(c => ({ ...c }));
      let modifiedNpes = npes.map(n => ({ ...n }));


      if (sweepTarget === 'connection') {
        modifiedConns = modifiedConns.map(c => {
          const connName = `${c.presynapticNPE}-${c.postsynapticNPE}`;
          if (connName === selectedElement) {
            const paramKey = selectedParam === 'alpha_prime' ? 'alphaPrime' : selectedParam === 'beta_prime' ? 'betaPrime' : selectedParam;
            return { ...c, [paramKey]: paramValue };
          }
          return c;
        });
      } else {
        modifiedNpes = modifiedNpes.map(n => {
          if (n.name === selectedElement) {
            return { ...n, [selectedParam]: paramValue };
          }
          return n;
        });
      }

      try {
        const npeData: Record<string, string[] | number[]> = {
          NPE: modifiedNpes.map(n => n.name),
          Type: modifiedNpes.map(n => n.type),
          Layer: modifiedNpes.map(n => n.layer),
          Activation: modifiedNpes.map(n => n.activation),
          'Temporal.Summation': modifiedNpes.map(n => n.temporalSummation),
          'Activation.Decay': modifiedNpes.map(n => n.activationDecay),
          mu: modifiedNpes.map(n => n.mu),
          sigma: modifiedNpes.map(n => n.sigma),
          logisSigma: modifiedNpes.map(n => n.logisSigma),
        };
        const connData: Record<string, string[] | number[]> = {
          PreSinapticNPE: modifiedConns.map(c => c.presynapticNPE),
          PostSinapticNPE: modifiedConns.map(c => c.postsynapticNPE),
          Weight: modifiedConns.map(c => c.weight),
          alpha: modifiedConns.map(c => c.alpha),
          beta: modifiedConns.map(c => c.beta),
          alpha_prime: modifiedConns.map(c => c.alphaPrime),
          beta_prime: modifiedConns.map(c => c.betaPrime),
        };

        const result = await runSimulation({
          npes: npeData,
          connections: connData,
          trials,
          contingencies,
          hasITI,
          numNetworks,
          threshold: getThresholdType(),
          disc,
          model: modelKind,
        });

        if (result.success) {
          const allNetMeans: number[] = [];
          for (const netData of result.results) {
            const lastPhase = netData[netData.length - 1]?.Phase;
            const lastPhaseData = netData.filter((r: SimulationResult) =>
              r.Phase === lastPhase && (true)
            );
            const vals = lastPhaseData.map((r: SimulationResult) => typeof r[outputUnit] === 'number' ? Number(r[outputUnit]) : 0);
            const mean = vals.length > 0 ? vals.reduce((a: number, b: number) => a + b, 0) / vals.length : 0;
            allNetMeans.push(mean);
          }
          const overallMean = allNetMeans.reduce((a, b) => a + b, 0) / allNetMeans.length;
          const sorted = [...allNetMeans].sort((a, b) => a - b);
          const overallMedian = sorted[Math.floor(sorted.length / 2)];

          sweepResults.push({
            paramValue: Number(paramValue.toFixed(4)),
            meanActivation: Number(overallMean.toFixed(6)),
            medianActivation: Number(overallMedian.toFixed(6)),
          });
        }
      } catch (err) {
        console.error(`Sweep step ${i} failed:`, err);
      }
    }

    setResults(sweepResults);
    setProgress(100);
    setIsRunning(false);
  }, [canRun, sweepTarget, selectedElement, selectedParam, rangeMin, rangeMax, steps,
      npes, connections, trials, contingencies, hasITI, numNetworks, getThresholdType, disc, outputUnit,
      modelKind]);

  return (
    <div className="max-w-7xl mx-auto space-y-6">
      <div className="flex items-center gap-3">
        <div className="w-10 h-10 rounded-xl bg-amber-500/10 flex items-center justify-center">
          <Sliders size={20} className="text-amber-600" />
        </div>
        <div>
          <h1 className="text-2xl font-bold text-slate-800">{t.sweep.pageTitle}</h1>
          <p className="text-sm text-slate-500">{t.sweep.pageSubtitle}</p>
        </div>
      </div>

      {/* Explanation card */}
      <Card className="p-4 border-amber-200 bg-amber-50/50">
        <div className="flex gap-3">
          <Info size={18} className="text-amber-600 flex-shrink-0 mt-0.5" />
          <div className="text-xs text-slate-600 leading-relaxed space-y-1">
            <p className="font-bold text-slate-700">{t.sweep.howItWorks}</p>
            <p dangerouslySetInnerHTML={{ __html: t.sweep.howDescription1 }} />
            <p dangerouslySetInnerHTML={{ __html: t.sweep.howDescription2 }} />
          </div>
        </div>
      </Card>

      <div className="grid grid-cols-1 lg:grid-cols-3 gap-6">
        <Card className="p-5 space-y-4">
          <h3 className="text-sm font-bold text-slate-700 flex items-center gap-2">
            <Sliders size={14} className="text-amber-600" /> {t.sweep.sweepConfig}
          </h3>

          <div>
            <label className="flex items-center gap-1.5 text-xs font-semibold text-slate-500 mb-1">
              {t.sweep.targetType}
              <Tooltip content={t.sweep.targetTypeTooltip} />
            </label>
            <div className="flex gap-2">
              {((['connection', 'npe'] as const)).map(typ => (
                <button key={typ} onClick={() => {
                  setSweepTarget(typ);
                  setSelectedElement('');
                  setSelectedParam(typ === 'connection' ? 'weight' : 'mu');
                }}
                  className={`flex-1 py-1.5 rounded-lg text-xs font-bold border transition-all ${
                    sweepTarget === typ ? 'bg-amber-50 text-amber-600 border-amber-200' : 'text-slate-400 border-slate-200'
                  }`}>
                  {typ === 'connection' ? t.sweep.connection : t.sweep.unitNPE}
                </button>
              ))}
            </div>
          </div>

          <div>
            <label className="flex items-center gap-1.5 text-xs font-semibold text-slate-500 mb-1">
              {t.sweep.element}
              <Tooltip content={t.sweep.elementTooltip} />
            </label>
            <select value={selectedElement} onChange={(e) => setSelectedElement(e.target.value)}
              className="w-full px-3 py-2 rounded-lg bg-white border border-slate-200 text-slate-800 text-sm focus:border-amber-500/50 focus:outline-none">
              <option value="">{t.sweep.select}</option>
              {(sweepTarget === 'connection'
                ? connections.map(c => {
                    const name = `${c.presynapticNPE}-${c.postsynapticNPE}`;
                    return <option key={name} value={name}>{name}</option>;
                  })
                : npes.map(n => <option key={n.name} value={n.name}>{n.name}</option>))
              }
            </select>
          </div>

          <div>
            <label className="flex items-center gap-1.5 text-xs font-semibold text-slate-500 mb-1">
              {t.sweep.parameter}
              <Tooltip content={t.sweep.parameterTooltip} />
            </label>
            <select value={selectedParam} onChange={(e) => setSelectedParam(e.target.value)}
              className="w-full px-3 py-2 rounded-lg bg-white border border-slate-200 text-slate-800 text-sm focus:border-amber-500/50 focus:outline-none">
              {((sweepTarget === 'connection' ? connectionParams : npeParams)).map(p => (
                <option key={p} value={p}>{p}</option>
              ))}
            </select>
          </div>

          <div className="grid grid-cols-2 gap-2">
            <div>
              <label className="flex items-center gap-1 text-[10px] font-semibold text-slate-400 mb-0.5">{t.sweep.min} <Tooltip content={t.sweep.minTooltip} iconSize={10} /></label>
              <input type="number" step="0.01" value={rangeMin}
                onChange={(e) => setRangeMin(parseFloat(e.target.value) || 0)}
                className="w-full px-2 py-1.5 rounded-md bg-white border border-slate-200 text-slate-700 text-xs focus:border-amber-500/50 focus:outline-none" />
            </div>
            <div>
              <label className="flex items-center gap-1 text-[10px] font-semibold text-slate-400 mb-0.5">{t.sweep.max} <Tooltip content={t.sweep.maxTooltip} iconSize={10} /></label>
              <input type="number" step="0.01" value={rangeMax}
                onChange={(e) => setRangeMax(parseFloat(e.target.value) || 0)}
                className="w-full px-2 py-1.5 rounded-md bg-white border border-slate-200 text-slate-700 text-xs focus:border-amber-500/50 focus:outline-none" />
            </div>
          </div>

          <div className="grid grid-cols-2 gap-2">
            <div>
              <label className="flex items-center gap-1 text-[10px] font-semibold text-slate-400 mb-0.5">{t.sweep.steps} <Tooltip content={t.sweep.stepsTooltip} iconSize={10} /></label>
              <input type="number" min="2" max="50" value={steps}
                onChange={(e) => setSteps(parseInt(e.target.value) || 5)}
                className="w-full px-2 py-1.5 rounded-md bg-white border border-slate-200 text-slate-700 text-xs focus:border-amber-500/50 focus:outline-none" />
            </div>
            <div>
              <label className="flex items-center gap-1 text-[10px] font-semibold text-slate-400 mb-0.5">{t.sweep.networksPerStep} <Tooltip content={t.sweep.networksPerStepTooltip} iconSize={10} /></label>
              <input type="number" min="1" max="20" value={numNetworks}
                onChange={(e) => setNumNetworks(parseInt(e.target.value) || 1)}
                className="w-full px-2 py-1.5 rounded-md bg-white border border-slate-200 text-slate-700 text-xs focus:border-amber-500/50 focus:outline-none" />
            </div>
          </div>

          <div>
            <label className="flex items-center gap-1.5 text-xs font-semibold text-slate-500 mb-1">
              {t.sweep.outputUnit}
              <Tooltip content={t.sweep.outputUnitTooltip} />
            </label>
            <select value={outputUnit} onChange={(e) => setOutputUnit(e.target.value)}
              className="w-full px-3 py-2 rounded-lg bg-white border border-slate-200 text-slate-800 text-sm focus:border-amber-500/50 focus:outline-none">
              <option value="">{t.sweep.select}</option>
              {motorUnits.map(u => <option key={u} value={u}>{u}</option>)}
            </select>
          </div>

          <button onClick={handleRunSweep} disabled={!canRun || isRunning}
            className="w-full py-2.5 rounded-xl bg-gradient-to-r from-amber-500 to-amber-600 text-white font-bold text-sm disabled:opacity-40 disabled:cursor-not-allowed hover:shadow-lg hover:shadow-amber-500/20 transition-shadow flex items-center justify-center gap-2">
            {isRunning ? <><Loader2 size={14} className="animate-spin" /> {t.sweep.running} ({progress}%)...</> : <><Play size={14} /> {t.sweep.runSweep}</>}
          </button>
          {isRunning && (
            <button onClick={() => { cancelRef.current = true; }}
              className="w-full py-2.5 rounded-xl bg-slate-100 text-slate-600 font-bold text-sm border border-slate-200 hover:bg-slate-200 transition-colors flex items-center justify-center gap-2 mt-2">
              <Square size={14} /> Cancel
            </button>
          )}
          {results.length > 0 && !isRunning && (
            <button
              onClick={() => {
                const csv = ['paramValue,meanActivation,medianActivation',
                  ...results.map(r => `${r.paramValue},${r.meanActivation},${r.medianActivation}`)
                ].join('\n');
                const blob = new Blob([csv], { type: 'text/csv' });
                const url = URL.createObjectURL(blob);
                const a = document.createElement('a');
                a.href = url;
                a.download = `sweep-${selectedElement}-${selectedParam}.csv`;
                a.click();
                URL.revokeObjectURL(url);
              }}
              className="w-full py-2.5 rounded-xl bg-white text-slate-600 font-bold text-sm border border-slate-200 hover:bg-slate-50 transition-colors flex items-center justify-center gap-2 mt-2">
              <Download size={14} /> Export CSV
            </button>
          )}
        </Card>

        <div className="lg:col-span-2">
          <Card className="p-6">
            <h3 className="text-lg font-bold text-slate-800 mb-4 flex items-center gap-2">
              <TrendingUp size={18} className="text-amber-600" />
              {t.sweep.sensitivityAnalysis}
              {results.length > 0 && <Badge variant="info">{results.length} {t.sweep.points}</Badge>}
            </h3>
            <div className="h-[450px]">
              {results.length > 0 ? (
                <ResponsiveContainer
                  width="100%"
                  height="100%"
                  minWidth={0}
                  minHeight={450}
                  initialDimension={{ width: 800, height: 450 }}
                >
                  <LineChart data={results}>
                    <CartesianGrid strokeDasharray="3 3" stroke="#e2e8f0" />
                    <XAxis dataKey="paramValue" stroke="#64748b" fontSize={11}
                      label={{ value: selectedParam, position: 'bottom', offset: -5, fill: '#64748b', fontSize: 12 }} />
                    <YAxis stroke="#64748b" fontSize={11}
                      label={{ value: `${outputUnit} ${t.sweep.activationLabel}`, angle: -90, position: 'insideLeft', fill: '#64748b', fontSize: 12 }} />
                    <RechartsTooltip contentStyle={{ backgroundColor: '#fff', border: '1px solid #e2e8f0', borderRadius: 12, fontSize: 12 }} />
                    <Legend wrapperStyle={{ fontSize: '12px' }} />
                    <Line type="monotone" dataKey="meanActivation" stroke={SWEEP_COLORS[0]} strokeWidth={2.5} name={t.sweep.meanLabel} dot={{ r: 3, fill: SWEEP_COLORS[0] }} />
                    <Line type="monotone" dataKey="medianActivation" stroke={SWEEP_COLORS[1]} strokeWidth={2} strokeDasharray="5 5" name={t.sweep.medianLabel} dot={{ r: 3, fill: SWEEP_COLORS[1] }} />
                  </LineChart>
                </ResponsiveContainer>
              ) : (
                <div className="h-full flex items-center justify-center">
                  <div className="text-center space-y-3">
                    <div className="w-16 h-16 rounded-2xl bg-slate-50 border border-slate-200 mx-auto flex items-center justify-center">
                      <Sliders size={28} className="text-slate-300" />
                    </div>
                    <p className="text-sm text-slate-400 max-w-sm">
                      {isRunning ? t.sweep.runningSweep : t.sweep.emptyState}
                    </p>
                  </div>
                </div>
              )}
            </div>
          </Card>
        </div>
      </div>
    </div>
  );
}
