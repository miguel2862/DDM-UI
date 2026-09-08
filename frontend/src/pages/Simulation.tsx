import { useState, useCallback, useEffect, useRef, useMemo } from 'react';
import { useNavigate } from 'react-router-dom';
import { motion, AnimatePresence } from 'framer-motion';
import {
  Play, Pause, SkipForward, RotateCcw, CheckCircle2, AlertCircle,
  Loader2, Brain, Activity, Gauge, ChevronRight, ChevronLeft, Download, Upload, ChevronDown,
} from 'lucide-react';
import { Card } from '../components/ui/Card';
import { StatCard } from '../components/ui/StatCard';
import { Badge } from '../components/ui/Badge';
import { Tooltip } from '../components/ui/Tooltip';
import { useSimStore } from '../stores/useSimStore';
import { useConfirm } from '../components/ui/ConfirmDialog';
import { useToast } from '../components/ui/Toast';
import { runSimulationOne, serializeNetwork, validateNetwork } from '../api/client';
import { downloadArchitectureJSON } from '../utils/dataExport';
import { useI18n, type Translations } from '../i18n';
import type { SimulationMetadata, SimulationResult, ThresholdPreset } from '../types/ddm';
import type { SimulationInspector } from '../types/inspector';
import { PublicationNetwork } from '../components/network/PublicationNetwork';
import { computePublicationLayout } from '../utils/publicationNetwork';
import { EquationInspector } from '../components/simulation/EquationInspector';

const pupdateOptions = [
  { value: 'async_random', labelKey: 'asyncRandom' as const, tooltipKey: 'asyncRandomTooltip' as const },
  { value: 'async_sequential', labelKey: 'asyncSequential' as const, tooltipKey: 'asyncSequentialTooltip' as const },
  { value: 'sync_random', labelKey: 'syncRandom' as const, tooltipKey: 'syncRandomTooltip' as const },
  { value: 'sync_sequential', labelKey: 'syncSequential' as const, tooltipKey: 'syncSequentialTooltip' as const },
];

function PUpdateDropdown({ pupdate, setPupdate, t }: { pupdate: string; setPupdate: (v: string) => void; t: Translations }) {
  const [open, setOpen] = useState(false);
  const [hoveredIdx, setHoveredIdx] = useState<number | null>(null);
  const dropdownRef = useRef<HTMLDivElement>(null);

  // Close on outside click
  useEffect(() => {
    const handler = (e: MouseEvent) => {
      if (dropdownRef.current && !dropdownRef.current.contains(e.target as HTMLElement)) {
        setOpen(false);
      }
    };
    document.addEventListener('mousedown', handler);
    return () => document.removeEventListener('mousedown', handler);
  }, []);

  const selectedOpt = pupdateOptions.find(o => o.value === pupdate) || pupdateOptions[0];

  return (
    <div ref={dropdownRef}>
      <label className="flex items-center gap-1.5 text-xs font-semibold text-slate-500 mb-1">
        {t.sim.updateProcedure}
        <Tooltip content={t.sim.updateProcedureTooltip} />
      </label>
      <div className="relative">
        <button
          type="button"
          onClick={() => setOpen(!open)}
          className="w-full flex items-center justify-between px-3 py-1.5 rounded-lg bg-slate-100 border border-slate-200 text-slate-800 text-sm focus:border-cyan-500/50 focus:outline-none hover:bg-slate-50 transition-colors"
        >
          <span>{t.sim[selectedOpt.labelKey]}</span>
          <ChevronDown size={14} className={`text-slate-400 transition-transform ${open ? 'rotate-180' : ''}`} />
        </button>
        <AnimatePresence>
          {open && (
            <motion.div
              initial={{ opacity: 0, y: -4 }}
              animate={{ opacity: 1, y: 0 }}
              exit={{ opacity: 0, y: -4 }}
              transition={{ duration: 0.15 }}
              className="absolute z-50 top-full left-0 right-0 mt-1 rounded-xl bg-white border border-slate-200 shadow-lg overflow-hidden"
            >
              {pupdateOptions.map((opt, idx) => (
                <button
                  key={opt.value}
                  type="button"
                  onClick={() => { setPupdate(opt.value); setOpen(false); }}
                  onMouseEnter={() => setHoveredIdx(idx)}
                  onMouseLeave={() => setHoveredIdx(null)}
                  className={`w-full text-left px-3 py-2 text-sm transition-colors ${
                    pupdate === opt.value
                      ? 'bg-cyan-50 text-cyan-700 font-semibold'
                      : 'text-slate-600 hover:bg-slate-50'
                  }`}
                >
                  {t.sim[opt.labelKey]}
                </button>
              ))}
              {/* Tooltip panel for hovered option */}
              <AnimatePresence>
                {hoveredIdx !== null && (
                  <motion.div
                    key={hoveredIdx}
                    initial={{ opacity: 0, height: 0 }}
                    animate={{ opacity: 1, height: 'auto' }}
                    exit={{ opacity: 0, height: 0 }}
                    transition={{ duration: 0.15 }}
                    className="border-t border-slate-100 bg-slate-50 px-3 py-2"
                  >
                    <p className="text-[11px] text-slate-500 leading-relaxed">
                      {t.sim[pupdateOptions[hoveredIdx].tooltipKey]}
                    </p>
                  </motion.div>
                )}
              </AnimatePresence>
            </motion.div>
          )}
        </AnimatePresence>
      </div>
    </div>
  );
}

export function Simulation() {
  const navigate = useNavigate();
  const {
    npes, connections, trials, contingencies, hasITI, modelKind,
    numNetworks, thresholdPreset, disc, pupdate, appMode, simStatus, simError, simulationResults, simulationMetadata,
    getThresholdType,
    setNumNetworks, setThresholdPreset, setDisc, setPupdate, setSimStatus, setSimResults, setSimError,
    playbackIndex, isPlaying, playbackSpeed,
    setPlaybackIndex, setIsPlaying, setPlaybackSpeed,
    lockedLayout, simulationInspector,
  } = useSimStore();
  const { t, language } = useI18n();
  const { confirm } = useConfirm();
  const toast = useToast();

  const [progress, setProgress] = useState(0);
  const [saveLoadMsg, setSaveLoadMsg] = useState<string | null>(null);
  const [isCancelling, setIsCancelling] = useState(false);
  const [simStartTime, setSimStartTime] = useState<number | null>(null);
  const [completedElapsedMs, setCompletedElapsedMs] = useState(0);
  const [etaNow, setEtaNow] = useState(0);
  const [inspectorEnabled, setInspectorEnabled] = useState(false);
  const [inspectorLimit, setInspectorLimit] = useState(2000);
  const [selectedInspectorUnit, setSelectedInspectorUnit] = useState<string>();
  const playbackRef = useRef<number | null>(null);
  const fileInputRef = useRef<HTMLInputElement>(null);
  const cancelledRef = useRef(false);

  const canRun = npes.length >= 2 && connections.length >= 1 && Object.keys(trials).length >= 1 && contingencies.length >= 1;

  const handleSaveExperiment = useCallback(() => {
    downloadArchitectureJSON(
      npes, connections, trials, contingencies, hasITI,
      `${modelKind}-experiment-${new Date().toISOString().slice(0, 10)}.json`,
      { numNetworks, thresholdPreset, disc, pupdate, modelKind },
      lockedLayout
    );
    setSaveLoadMsg(t.sim.experimentSaved);
    setTimeout(() => setSaveLoadMsg(null), 3000);
  }, [npes, connections, trials, contingencies, hasITI, numNetworks, thresholdPreset, disc, pupdate, modelKind, lockedLayout, t]);

  const handleLoadExperiment = useCallback((e: React.ChangeEvent<HTMLInputElement>) => {
    const file = e.target.files?.[0];
    if (!file) return;
    const reader = new FileReader();
    reader.onload = (ev) => {
      try {
        const data = JSON.parse(ev.target?.result as string);
        if (data._type !== 'ddm-ui-experiment' || !data.npes || !data.connections) {
          setSaveLoadMsg(t.sim.invalidFile);
          toast.error(t.sim.invalidFile);
          setTimeout(() => setSaveLoadMsg(null), 4000);
          return;
        }
        // Load all experiment data into the store
        const store = useSimStore.getState();
        store.loadExperiment(data);
        if (data.numNetworks) setNumNetworks(data.numNetworks);
        if (data.thresholdPreset) setThresholdPreset(data.thresholdPreset);
        if (data.disc) setDisc(data.disc);
        if (data.pupdate) setPupdate(data.pupdate);
        // Restore locked layout if present
        if (data.lockedLayout) {
          store.setLockedLayout(data.lockedLayout);
        }
        setSaveLoadMsg(t.sim.experimentLoaded);
        setTimeout(() => setSaveLoadMsg(null), 3000);
      } catch {
        setSaveLoadMsg(t.sim.invalidFile);
        toast.error(t.sim.invalidFile);
        setTimeout(() => setSaveLoadMsg(null), 4000);
      }
    };
    reader.readAsText(file);
    // Reset file input so same file can be loaded again
    e.target.value = '';
  }, [t, toast, setNumNetworks, setThresholdPreset, setDisc, setPupdate]);

  const [completedNetworks, setCompletedNetworks] = useState(0);

  // Cancel simulation handler
  const handleCancelSimulation = useCallback(() => {
    cancelledRef.current = true;
    setIsCancelling(true);
  }, []);

  const handleRunSimulation = useCallback(async () => {
    if (!canRun) return;

    // ── Frontend validation ──
    const hasUSD = connections.some(c => c.presynapticNPE === 'US' && c.postsynapticNPE === 'D');
    if (!hasUSD) {
      toast.warning(t.sim.validationNoUSD);
    }

    const connectedNPEs = new Set<string>();
    connections.forEach(c => { connectedNPEs.add(c.presynapticNPE); connectedNPEs.add(c.postsynapticNPE); });
    const disconnected = npes.filter(n => !connectedNPEs.has(n.name));
    if (disconnected.length > 0) {
      toast.warning(t.sim.validationDisconnected.replace('{units}', disconnected.map(n => n.name).join(', ')));
    }

    try {
      const validation = await validateNetwork(npes, connections);
      if (!validation.valid) {
        const message = validation.error || t.sim.simulationFailed;
        setSimError(message);
        toast.error(message);
        return;
      }
    } catch (error: unknown) {
      const message = error instanceof Error ? error.message : t.sim.errorNetwork;
      setSimError(message);
      toast.error(message);
      return;
    }

    setSimStatus('running');
    setProgress(0);
    setCompletedNetworks(0);
    cancelledRef.current = false;
    setIsCancelling(false);
    const startMs = Date.now();
    setSimStartTime(startMs);
    setCompletedElapsedMs(0);
    setEtaNow(startMs);

    // Build R-compatible data (shared across all network calls)
    const { npes: npeData, connections: connData } = serializeNetwork(npes, connections);

    const baseParams = {
      npes: npeData,
      connections: connData,
      trials,
      contingencies,
      hasITI,
      threshold: getThresholdType(),
      disc,
      pupdate,
      model: modelKind,
    };

    try {
      // Run networks one at a time for real progress tracking
      const results: SimulationResult[][] = [];
      type NetworkMetadata = Partial<SimulationMetadata> & {
        sampledParameters?: Record<string, number>;
      };
      const networkMetadata: NetworkMetadata[] = [];
      let firstInspector: SimulationInspector | null = null;

      for (let i = 0; i < numNetworks; i++) {
        // Check if user cancelled
        if (cancelledRef.current) {
          setSimStatus('idle');
          setIsCancelling(false);
          setSimStartTime(null);
          toast.info(t.toast.simulationCancelled);
          return;
        }

        const res = await runSimulationOne({
          ...baseParams,
          inspector: inspectorEnabled && i === 0
            ? { enabled: true, maxTimesteps: inspectorLimit } : undefined,
        });

        if (cancelledRef.current) {
          setSimStatus('idle');
          setIsCancelling(false);
          setSimStartTime(null);
          toast.info(t.toast.simulationCancelled);
          return;
        }

        if (!res.success) {
          setSimError(res.error || t.sim.simulationFailed);
          setSimStartTime(null);
          return;
        }

        results.push(res.result);
        if (i === 0 && res.inspector) firstInspector = res.inspector;
        if (res.metadata) networkMetadata.push(res.metadata);
        const done = i + 1;
        setCompletedNetworks(done);
        const now = Date.now();
        setCompletedElapsedMs(now - startMs);
        setEtaNow(now);
        // Real progress: percentage of completed networks
        setProgress(Math.round((done / numNetworks) * 100));
      }

      if (results.length === 0) {
        throw new Error(t.sim.simulationFailed);
      }

      const duration = ((Date.now() - startMs) / 1000).toFixed(2);
      setSimStartTime(null);

      const allCols = Object.keys(results[0][0] || {});
      const metaCols = ['Phase', 'Trial', 'TimeStep'];
      const signalCols = ['dVTA', 'dH'];
      const dataCols = allCols.filter(c => !metaCols.includes(c) && !signalCols.includes(c));
      const unitCols = dataCols.filter(c => !c.includes('-'));
      const connectionCols = dataCols.filter(c => c.includes('-'));


      setSimResults(results, {
        numNetworks,
        model: 'DTD',
        phases: Array.from(new Set(results[0].map((r) => String(r.Phase)))),
        units: unitCols,
        connections: connectionCols,
        totalTrials: Math.max(...results[0].map((r) => Number(r.Trial))),
        totalTimesteps: results[0].length,
        duration: parseFloat(duration),
        disc,
        signals: signalCols,
        networkParameters: networkMetadata.map(meta => meta.sampledParameters || {}),
      }, firstInspector);
      setPlaybackIndex(0);
      setIsPlaying(false);
      if (inspectorEnabled && !firstInspector) {
        toast.info(language === 'es' ? 'Este servidor no incluye el inspector. Los resultados se conservaron; ejecuta la versión actualizada para registrar las ecuaciones.' : 'This server does not include the inspector. Results were preserved; use the updated version to record equations.');
      }
    } catch (err: unknown) {
      setSimStartTime(null);
      if (cancelledRef.current) {
        setSimStatus('idle');
        setIsCancelling(false);
        toast.info(t.toast.simulationCancelled);
        return;
      }
      let errorMsg = err instanceof Error ? err.message : t.sim.errorFallback;
      if (errorMsg.includes('timed out')) {
        errorMsg = t.sim.errorTimeout.replace('{networks}', String(numNetworks)).replace('{phases}', String(contingencies.length));
      } else if (errorMsg.includes('Failed to fetch') || errorMsg.includes('NetworkError')) {
        errorMsg = t.sim.errorNetwork;
      } else if (errorMsg.includes('500')) {
        errorMsg = t.sim.errorInternal;
      }
      setSimError(errorMsg);
    }
  }, [canRun, npes, connections, trials, contingencies, hasITI, numNetworks, getThresholdType, disc, pupdate, modelKind, setSimStatus, setSimResults, setSimError, setPlaybackIndex, setIsPlaying, inspectorEnabled, inspectorLimit, language, t, toast]);

  // Live countdown ETA — ticks every second so it never "freezes"
  useEffect(() => {
    if (!simStartTime || completedNetworks === 0 || completedNetworks >= numNetworks) return;
    const id = window.setInterval(() => setEtaNow(Date.now()), 1000);
    return () => clearInterval(id);
  }, [simStartTime, completedNetworks, numNetworks]);

  const etaText = useMemo(() => {
    if (!simStartTime || completedNetworks === 0 || completedNetworks >= numNetworks || completedElapsedMs <= 0) {
      return '';
    }
    const perNetwork = completedElapsedMs / completedNetworks;
    const remaining = Math.max(0, perNetwork * numNetworks - (etaNow - simStartTime));
    if (remaining < 500) return '';
    const secs = Math.round(remaining / 1000);
    const time = secs < 60
      ? `${secs}s`
      : `${Math.floor(secs / 60)}m ${secs % 60}s`;
    return t.sim.etaRemaining.replace('{time}', time);
  }, [simStartTime, completedNetworks, numNetworks, completedElapsedMs, etaNow, t]);

  // Playback logic
  const currentData = useMemo(() => {
    if (!simulationResults || simulationResults.length === 0) return null;
    return simulationResults[0];
  }, [simulationResults]);

  const maxIndex = currentData ? currentData.length - 1 : 0;

  const publicationPositions = useMemo(() => {
    const automatic = computePublicationLayout(npes, connections);
    if (!lockedLayout) return automatic;
    return Object.fromEntries(npes.map(npe => {
      const point = lockedLayout[npe.name];
      return [npe.name, point ? { x: point.x + 30, y: point.y + 30 } : automatic[npe.name]];
    }));
  }, [npes, connections, lockedLayout]);

  const inspectRow = useCallback((index: number) => {
    setIsPlaying(false);
    setPlaybackIndex(Math.max(0, Math.min(maxIndex, index)));
  }, [maxIndex, setIsPlaying, setPlaybackIndex]);

  useEffect(() => {
    if (isPlaying && currentData) {
      playbackRef.current = window.setInterval(() => {
        const current = useSimStore.getState().playbackIndex;
        if (current >= maxIndex) {
          setIsPlaying(false);
          setPlaybackIndex(maxIndex);
        } else {
          setPlaybackIndex(current + 1);
        }
      }, 200 / playbackSpeed);
    }
    return () => {
      if (playbackRef.current) clearInterval(playbackRef.current);
    };
  }, [isPlaying, playbackSpeed, maxIndex, currentData, setPlaybackIndex, setIsPlaying]);

  // Build animated network nodes




  const currentRow = currentData && playbackIndex < currentData.length ? currentData[playbackIndex] : null;

  return (
    <div className="max-w-7xl mx-auto space-y-6">
      <div className="flex items-center gap-3">
        <div className="w-10 h-10 rounded-xl bg-emerald-500/10 flex items-center justify-center">
          <Activity size={20} className="text-emerald-600" />
        </div>
        <div>
          <h1 className="text-2xl font-bold text-slate-800">{t.sim.pageTitle}</h1>
          <p className="text-sm text-slate-500">{t.sim.pageSubtitle}</p>
        </div>
      </div>

      {/* Config + Controls */}
      <div className="grid grid-cols-1 md:grid-cols-4 gap-4">
        <StatCard icon={Brain} label={t.sim.networks} value={numNetworks} color="cyan" />
        <StatCard
          icon={Gauge}
          label={(t.sim.thresholdLabel)}
          value={(t.sim[`preset_${thresholdPreset}` as keyof typeof t.sim] || thresholdPreset)}
          color="teal"
        />
        <StatCard icon={Activity} label={t.sim.status} value={simStatus} color={simStatus === 'complete' ? 'emerald' : simStatus === 'error' ? 'rose' : 'amber'} />
        {appMode === 'advanced' && (
          <StatCard icon={Gauge} label={t.sim.discCriterion} value={disc} color="violet" />
        )}
      </div>

      {/* Simulation Parameters */}
      <Card className="p-4">

        <div className={`grid grid-cols-1 ${appMode === 'advanced' ? 'md:grid-cols-4' : 'md:grid-cols-2'} gap-4`}>
          <div>
            <label className="block text-xs font-semibold text-slate-500 mb-1">{t.sim.networks}</label>
            <input
              type="number"
              min="1"
              max="100"
              value={numNetworks}
              onChange={(e) => setNumNetworks(parseInt(e.target.value) || 1)}
              className="w-full px-3 py-1.5 rounded-lg bg-slate-100 border border-slate-200 text-slate-800 text-sm focus:border-cyan-500/50 focus:outline-none"
            />
          </div>
          {<div>
            <label className="flex items-center gap-1.5 text-xs font-semibold text-slate-500 mb-1">
              {t.sim.thresholdPreset}
              <Tooltip content={t.sim.thresholdPresetTooltip} />
            </label>
            <select
              value={thresholdPreset}
              onChange={async (e) => {
                const newPreset = e.target.value as ThresholdPreset;
                if (npes.length > 0 && newPreset !== thresholdPreset) {
                  // Only warn if μ/σ values will actually change
                  // gaussian_ddmui and beta_ddmui share μ=0.2, σ=0.15 — switching between them doesn't change NPE params
                  const sameParams = new Set(['gaussian_ddmui', 'beta_ddmui']);
                  const paramsWillChange = !(sameParams.has(thresholdPreset) && sameParams.has(newPreset));
                  if (paramsWillChange) {
                    const ok = await confirm({
                      title: t.confirm.thresholdChangeTitle,
                      message: t.confirm.thresholdChangeMessage,
                      confirmLabel: t.confirm.confirm,
                      cancelLabel: t.confirm.cancel,
                      variant: 'warning',
                    });
                    if (!ok) return;
                  }
                }
                setThresholdPreset(newPreset);
              }}
              className="w-full px-3 py-1.5 rounded-lg bg-slate-100 border border-slate-200 text-slate-800 text-sm focus:border-cyan-500/50 focus:outline-none"
            >
              <option value="gaussian_ddmui">{t.sim.preset_gaussian_ddmui}</option>
              <option value="gaussian_donahoe1993">{t.sim.preset_gaussian_donahoe1993}</option>
              <option value="beta_ddmui">{t.sim.preset_beta_ddmui}</option>
            </select>
            <p className="text-[10px] text-slate-400 mt-0.5">
              {thresholdPreset === 'gaussian_ddmui' && 'Gaussian: θ ~ N(μ=0.2, σ=0.15)'}
              {thresholdPreset === 'gaussian_donahoe1993' && 'Gaussian: θ ~ N(μ=0.0, σ=1.0)'}
              {thresholdPreset === 'beta_ddmui' && 'Beta: θ ~ Beta(μ=0.2, σ=0.15)'}
            </p>
          </div>}
          {appMode === 'advanced' && (
            <>
              <PUpdateDropdown pupdate={pupdate} setPupdate={setPupdate} t={t} />
              <div>
                <label className="flex items-center gap-1.5 text-xs font-semibold text-slate-500 mb-1">
                  {t.sim.discrepancyCriterion}: {disc}
                  <Tooltip content={t.sim.discrepancyTooltip} />
                </label>
                <input
                  type="range"
                  min="0.0001"
                  max="0.1"
                  step="0.0001"
                  value={disc}
                  onChange={(e) => setDisc(parseFloat(e.target.value))}
                  className="w-full accent-cyan-500"
                />
                <div className="flex justify-between text-[10px] text-slate-400 mt-0.5">
                  <span>0.0001</span>
                  <span>0.1</span>
                </div>
              </div>
            </>
          )}
        </div>

      </Card>

      {/* Save / Load Experiment */}
      <div className="flex items-center gap-3">
        <button
          onClick={handleSaveExperiment}
          disabled={!canRun}
          className="flex items-center gap-2 px-4 py-2.5 rounded-xl bg-white text-slate-600 font-semibold text-sm border border-slate-200 hover:bg-slate-50 hover:border-slate-300 transition-all disabled:opacity-30 disabled:cursor-not-allowed"
          title={t.sim.saveExperimentTooltip}
        >
          <Download size={16} />
          {t.sim.saveExperiment}
        </button>
        <button
          onClick={() => fileInputRef.current?.click()}
          className="flex items-center gap-2 px-4 py-2.5 rounded-xl bg-white text-slate-600 font-semibold text-sm border border-slate-200 hover:bg-slate-50 hover:border-slate-300 transition-all"
          title={t.sim.loadExperimentTooltip}
        >
          <Upload size={16} />
          {t.sim.loadExperiment}
        </button>
        <input
          ref={fileInputRef}
          type="file"
          accept=".json"
          className="hidden"
          onChange={handleLoadExperiment}
        />
        <AnimatePresence>
          {saveLoadMsg && (
            <motion.span
              initial={{ opacity: 0, x: -10 }}
              animate={{ opacity: 1, x: 0 }}
              exit={{ opacity: 0 }}
              className={`text-xs font-semibold ${
                saveLoadMsg === t.sim.invalidFile ? 'text-rose-500' : 'text-emerald-500'
              }`}
            >
              {saveLoadMsg}
            </motion.span>
          )}
        </AnimatePresence>
      </div>

      {(
        <section className="rounded-xl border border-slate-200 bg-white p-4" aria-label={language === 'es' ? 'Registro didáctico' : 'Teaching trace'}>
          <div className="flex flex-wrap items-start justify-between gap-4">
            <div className="flex-1 min-w-56">
              <label className="flex items-center gap-3 text-sm font-semibold text-slate-800">
                <input type="checkbox" checked={inspectorEnabled} onChange={event => setInspectorEnabled(event.target.checked)}
                  disabled={simStatus === 'running'} className="h-4 w-4 accent-cyan-700" />
                {language === 'es' ? 'Registrar ecuaciones paso a paso' : 'Record equations step by step'}
              </label>
              <p className="mt-2 text-sm leading-relaxed text-slate-600 max-w-prose">
                {language === 'es' ? 'Registra los valores reales de la primera red: umbrales, activaciones, señales difusas y cambios de peso. Después podrás avanzar y retroceder por timesteps y ensayos sin ejecutar de nuevo el modelo.' : 'Record the first network’s actual thresholds, activations, diffuse signals and weight changes. Afterwards, move through timesteps and trials without running the model again.'}
              </p>
            </div>
            <label className="flex flex-col gap-1.5 text-xs font-medium text-slate-600">
              {language === 'es' ? 'Límite del registro detallado' : 'Detailed trace limit'}
              <select value={inspectorLimit} onChange={event => setInspectorLimit(Number(event.target.value))}
                disabled={!inspectorEnabled || simStatus === 'running'}
                className="rounded-lg border border-slate-300 px-3 py-2 text-sm text-slate-800 bg-white disabled:opacity-50 focus-visible:outline-2 focus-visible:outline-cyan-600">
                {[2000, 5000, 10000].map(limit => <option key={limit} value={limit}>{limit.toLocaleString()} timesteps</option>)}
              </select>
              <span>{language === 'es' ? 'El límite no acorta la simulación.' : 'The limit does not shorten the simulation.'}</span>
            </label>
          </div>
        </section>
      )}

      {/* Run button / Progress */}
      <AnimatePresence mode="wait">
        {simStatus === 'idle' && (
          <motion.div key="run" initial={{ opacity: 0 }} animate={{ opacity: 1 }} exit={{ opacity: 0 }}>
            <button
              onClick={handleRunSimulation}
              disabled={!canRun}
              className="w-full py-4 rounded-2xl bg-gradient-to-r from-cyan-500 to-teal-500 text-white font-bold text-lg disabled:opacity-30 disabled:cursor-not-allowed hover:shadow-xl hover:shadow-cyan-500/20 transition-all flex items-center justify-center gap-3"
            >
              <Play size={24} />
              {t.sim.runSimulation}
            </button>
            {!canRun && (
              <p className="text-xs text-amber-400 mt-2 text-center">
                {t.sim.completeBeforeRunning}
              </p>
            )}
          </motion.div>
        )}

        {simStatus === 'running' && (
          <motion.div key="progress" initial={{ opacity: 0 }} animate={{ opacity: 1 }} exit={{ opacity: 0 }}>
            <Card className="p-6">
              <div className="flex items-center gap-4 mb-4">
                <Loader2 size={24} className="text-cyan-600 animate-spin" />
                <div className="flex-1">
                  <p className="text-sm font-bold text-slate-700">{t.sim.simulating}</p>
                  <p className="text-xs text-slate-500">
                    {t.results.networkLabel} {completedNetworks + (completedNetworks < numNetworks ? 1 : 0)}/{numNetworks}
                    {completedNetworks > 0 && completedNetworks < numNetworks && (
                      <span className="text-emerald-500 ml-2">✓ {completedNetworks} {t.sim.completed}</span>
                    )}
                    {etaText && (
                      <span className="text-cyan-500 ml-2">{etaText}</span>
                    )}
                    {completedNetworks === 0 && (
                      <span className="text-slate-400 ml-2">{t.sim.etaCalculating}</span>
                    )}
                  </p>
                </div>
                <span className="text-lg font-bold text-cyan-600">{Math.round(progress)}%</span>
                <button
                  onClick={handleCancelSimulation}
                  disabled={isCancelling}
                  className="px-3 py-1.5 rounded-lg text-xs font-bold text-rose-600 border border-rose-200 bg-rose-50 hover:bg-rose-100 transition-colors disabled:opacity-50"
                >
                  {isCancelling ? t.sim.cancelling : t.sim.cancelSimulation}
                </button>
              </div>
              <div className="w-full h-3 bg-slate-100 rounded-full overflow-hidden">
                <motion.div
                  className="h-full bg-gradient-to-r from-cyan-500 to-teal-500 rounded-full"
                  animate={{ width: `${progress}%` }}
                  transition={{ duration: 0.3 }}
                />
              </div>
              {numNetworks > 1 && (
                <div className="flex gap-1 mt-3">
                  {Array.from({ length: numNetworks }, (_, i) => (
                    <div
                      key={i}
                      className={`h-1.5 flex-1 rounded-full transition-colors duration-300 ${
                        i < completedNetworks
                          ? 'bg-emerald-400'
                          : i === completedNetworks
                            ? 'bg-cyan-400 animate-pulse'
                            : 'bg-slate-200'
                      }`}
                    />
                  ))}
                </div>
              )}

              {/* Neural animation during simulation */}
              <div className="mt-6 flex justify-center">
                <div className="relative w-[280px] h-[160px]">
                  <svg viewBox="0 0 280 160" className="w-full h-full" overflow="visible">
                    <defs>
                      <filter id="sim_glow" x="-50%" y="-50%" width="200%" height="200%">
                        <feGaussianBlur in="SourceGraphic" stdDeviation="3" />
                      </filter>
                    </defs>

                    {/* Connection lines — animated pulses */}
                    {/* Layer 1 → Layer 2 */}
                    <motion.line x1="40" y1="50" x2="100" y2="35" stroke="#06b6d4" strokeOpacity="0.4"
                      animate={{ strokeWidth: [0.5, 2.5, 0.5], strokeOpacity: [0.2, 0.6, 0.2] }}
                      transition={{ duration: 1.5, repeat: Infinity, ease: 'easeInOut', delay: 0 }} />
                    <motion.line x1="40" y1="50" x2="100" y2="80" stroke="#06b6d4" strokeOpacity="0.4"
                      animate={{ strokeWidth: [0.5, 2, 0.5], strokeOpacity: [0.2, 0.5, 0.2] }}
                      transition={{ duration: 1.8, repeat: Infinity, ease: 'easeInOut', delay: 0.3 }} />
                    <motion.line x1="40" y1="110" x2="100" y2="80" stroke="#06b6d4" strokeOpacity="0.4"
                      animate={{ strokeWidth: [0.5, 2, 0.5], strokeOpacity: [0.15, 0.5, 0.15] }}
                      transition={{ duration: 1.6, repeat: Infinity, ease: 'easeInOut', delay: 0.5 }} />
                    <motion.line x1="40" y1="110" x2="100" y2="125" stroke="#06b6d4" strokeOpacity="0.4"
                      animate={{ strokeWidth: [0.5, 2, 0.5], strokeOpacity: [0.2, 0.5, 0.2] }}
                      transition={{ duration: 1.7, repeat: Infinity, ease: 'easeInOut', delay: 0.2 }} />

                    {/* Layer 2 → Layer 3 */}
                    <motion.line x1="100" y1="35" x2="175" y2="55" stroke="#14b8a6" strokeOpacity="0.4"
                      animate={{ strokeWidth: [0.5, 2.5, 0.5], strokeOpacity: [0.2, 0.6, 0.2] }}
                      transition={{ duration: 1.4, repeat: Infinity, ease: 'easeInOut', delay: 0.8 }} />
                    <motion.line x1="100" y1="80" x2="175" y2="55" stroke="#14b8a6" strokeOpacity="0.4"
                      animate={{ strokeWidth: [0.5, 2, 0.5], strokeOpacity: [0.2, 0.5, 0.2] }}
                      transition={{ duration: 1.6, repeat: Infinity, ease: 'easeInOut', delay: 1.0 }} />
                    <motion.line x1="100" y1="80" x2="175" y2="105" stroke="#14b8a6" strokeOpacity="0.4"
                      animate={{ strokeWidth: [0.5, 2, 0.5], strokeOpacity: [0.15, 0.5, 0.15] }}
                      transition={{ duration: 1.5, repeat: Infinity, ease: 'easeInOut', delay: 0.7 }} />
                    <motion.line x1="100" y1="125" x2="175" y2="105" stroke="#14b8a6" strokeOpacity="0.4"
                      animate={{ strokeWidth: [0.5, 2.5, 0.5], strokeOpacity: [0.2, 0.6, 0.2] }}
                      transition={{ duration: 1.3, repeat: Infinity, ease: 'easeInOut', delay: 1.2 }} />

                    {/* Layer 3 → Output */}
                    <motion.line x1="175" y1="55" x2="240" y2="80" stroke="#8b5cf6" strokeOpacity="0.4"
                      animate={{ strokeWidth: [0.5, 3, 0.5], strokeOpacity: [0.2, 0.7, 0.2] }}
                      transition={{ duration: 1.2, repeat: Infinity, ease: 'easeInOut', delay: 1.5 }} />
                    <motion.line x1="175" y1="105" x2="240" y2="80" stroke="#8b5cf6" strokeOpacity="0.4"
                      animate={{ strokeWidth: [0.5, 3, 0.5], strokeOpacity: [0.2, 0.7, 0.2] }}
                      transition={{ duration: 1.4, repeat: Infinity, ease: 'easeInOut', delay: 1.7 }} />

                    {/* Nodes — Layer 1 (input) */}
                    <motion.circle cx="40" cy="50" r="8" fill="#06b6d4" filter="url(#sim_glow)"
                      animate={{ fillOpacity: [0.3, 0.8, 0.3], r: [7, 9, 7] }}
                      transition={{ duration: 1.2, repeat: Infinity, ease: 'easeInOut', delay: 0 }} />
                    <motion.circle cx="40" cy="110" r="8" fill="#06b6d4" filter="url(#sim_glow)"
                      animate={{ fillOpacity: [0.3, 0.8, 0.3], r: [7, 9, 7] }}
                      transition={{ duration: 1.4, repeat: Infinity, ease: 'easeInOut', delay: 0.4 }} />

                    {/* Nodes — Layer 2 (hidden) */}
                    <motion.circle cx="100" cy="35" r="7" fill="#14b8a6" filter="url(#sim_glow)"
                      animate={{ fillOpacity: [0.25, 0.75, 0.25], r: [6, 8, 6] }}
                      transition={{ duration: 1.3, repeat: Infinity, ease: 'easeInOut', delay: 0.6 }} />
                    <motion.circle cx="100" cy="80" r="7" fill="#14b8a6" filter="url(#sim_glow)"
                      animate={{ fillOpacity: [0.25, 0.75, 0.25], r: [6, 8, 6] }}
                      transition={{ duration: 1.5, repeat: Infinity, ease: 'easeInOut', delay: 0.9 }} />
                    <motion.circle cx="100" cy="125" r="7" fill="#14b8a6" filter="url(#sim_glow)"
                      animate={{ fillOpacity: [0.25, 0.7, 0.25], r: [6, 8, 6] }}
                      transition={{ duration: 1.1, repeat: Infinity, ease: 'easeInOut', delay: 0.3 }} />

                    {/* Nodes — Layer 3 (association) */}
                    <motion.circle cx="175" cy="55" r="7" fill="#8b5cf6" filter="url(#sim_glow)"
                      animate={{ fillOpacity: [0.2, 0.7, 0.2], r: [6, 8, 6] }}
                      transition={{ duration: 1.6, repeat: Infinity, ease: 'easeInOut', delay: 1.2 }} />
                    <motion.circle cx="175" cy="105" r="7" fill="#8b5cf6" filter="url(#sim_glow)"
                      animate={{ fillOpacity: [0.2, 0.7, 0.2], r: [6, 8, 6] }}
                      transition={{ duration: 1.3, repeat: Infinity, ease: 'easeInOut', delay: 1.0 }} />

                    {/* Output node — larger, prominent pulse */}
                    <motion.circle cx="240" cy="80" r="10" fill="#f59e0b" filter="url(#sim_glow)"
                      animate={{ fillOpacity: [0.3, 0.9, 0.3], r: [9, 12, 9] }}
                      transition={{ duration: 1.0, repeat: Infinity, ease: 'easeInOut', delay: 1.8 }} />

                    {/* Travelling signal dots */}
                    <motion.circle r="2.5" fill="#06b6d4"
                      animate={{ cx: [40, 100], cy: [50, 35], opacity: [0.8, 0] }}
                      transition={{ duration: 1.0, repeat: Infinity, ease: 'easeIn', delay: 0 }} />
                    <motion.circle r="2.5" fill="#14b8a6"
                      animate={{ cx: [100, 175], cy: [80, 55], opacity: [0.8, 0] }}
                      transition={{ duration: 1.0, repeat: Infinity, ease: 'easeIn', delay: 0.8 }} />
                    <motion.circle r="3" fill="#8b5cf6"
                      animate={{ cx: [175, 240], cy: [55, 80], opacity: [0.9, 0] }}
                      transition={{ duration: 0.8, repeat: Infinity, ease: 'easeIn', delay: 1.5 }} />
                  </svg>
                </div>
              </div>
            </Card>
          </motion.div>
        )}

        {simStatus === 'error' && (
          <motion.div key="error" initial={{ opacity: 0 }} animate={{ opacity: 1 }} exit={{ opacity: 0 }}>
            <Card className="p-6 border-rose-500/30">
              <div className="flex items-center gap-3 mb-2">
                <AlertCircle size={20} className="text-rose-600" />
                <p className="text-sm font-bold text-rose-600">{t.sim.simulationFailed}</p>
              </div>
              <p className="text-xs text-slate-500 mb-4">
                {simError || t.sim.apiNotRunning}
              </p>
              <button
                onClick={() => { setSimStatus('idle'); setSimError(null); }}
                className="px-4 py-2 rounded-lg bg-white text-slate-600 text-sm font-bold border border-slate-200 hover:bg-slate-100 transition-colors"
              >
                {t.sim.tryAgain}
              </button>
            </Card>
          </motion.div>
        )}
      </AnimatePresence>

      {/* Playback view */}
      {simStatus === 'complete' && currentData && (
        <motion.div
          initial={{ opacity: 0, y: 20 }}
          animate={{ opacity: 1, y: 0 }}
          transition={{ duration: 0.5 }}
          className="space-y-4"
        >
          {/* Success banner */}
          <Card className="p-4 border-emerald-200 bg-emerald-50">
            <div className="flex items-center justify-between">
              <div className="flex items-center gap-3">
                <CheckCircle2 size={20} className="text-emerald-600" />
                <div>
                  <p className="text-sm font-bold text-emerald-600">{t.sim.simulationComplete}</p>
                  <p className="text-xs text-slate-500">
                    {simulationMetadata?.numNetworks} {t.sim.networksSuffix} {t.sim.durationConnector} {simulationMetadata?.duration}s
                  </p>
                </div>
              </div>
              <div className="flex items-center gap-2">
                <button
                  onClick={() => { setSimStatus('idle'); setProgress(0); }}
                  className="px-4 py-2 rounded-xl bg-white text-slate-600 font-bold text-sm border border-slate-200 hover:bg-slate-100 transition-colors flex items-center gap-1"
                >
                  <RotateCcw size={14} /> {t.sim.runAgain}
                </button>
                <button
                  onClick={() => navigate('/results')}
                  className="px-4 py-2 rounded-xl bg-emerald-50 text-emerald-600 font-bold text-sm border border-emerald-200 hover:bg-emerald-200 transition-colors flex items-center gap-1"
                >
                  {t.sim.viewResults} <ChevronRight size={14} />
                </button>
              </div>
            </div>
          </Card>

          {/* Network Playback */}
          <Card className="p-0 overflow-hidden" glow="cyan">
            {/* Playback controls */}
            <div className="flex flex-wrap items-center gap-3 p-3 border-b border-slate-200 bg-slate-50">
              <button
                onClick={() => setIsPlaying(!isPlaying)}
                aria-label={isPlaying ? (language === 'es' ? 'Pausar reproducción' : 'Pause playback') : (language === 'es' ? 'Reproducir' : 'Play')}
                className="w-10 h-10 rounded-xl bg-cyan-50 text-cyan-600 flex items-center justify-center hover:bg-cyan-200 transition-colors"
              >
                {isPlaying ? <Pause size={18} /> : <Play size={18} />}
              </button>
              <button
                onClick={() => { setPlaybackIndex(0); setIsPlaying(false); }}
                aria-label={language === 'es' ? 'Volver al primer timestep' : 'Go to first timestep'}
                className="w-10 h-10 rounded-xl bg-slate-100 text-slate-500 flex items-center justify-center hover:bg-slate-100 transition-colors"
              >
                <RotateCcw size={16} />
              </button>
              <button
                onClick={() => inspectRow(playbackIndex - 1)} disabled={playbackIndex === 0}
                aria-label={language === 'es' ? 'Timestep anterior' : 'Previous timestep'}
                className="w-10 h-10 rounded-xl bg-slate-100 text-slate-600 flex items-center justify-center hover:bg-slate-200 disabled:opacity-40"
              >
                <ChevronLeft size={16} />
              </button>
              <button
                onClick={() => inspectRow(playbackIndex + 1)} disabled={playbackIndex >= maxIndex}
                aria-label={language === 'es' ? 'Siguiente timestep' : 'Next timestep'}
                className="w-10 h-10 rounded-xl bg-slate-100 text-slate-500 flex items-center justify-center hover:bg-slate-100 transition-colors"
              >
                <SkipForward size={16} />
              </button>

              {/* Speed control */}
              <div className="flex items-center gap-1 ml-2">
                {[1, 2, 5, 10].map(speed => (
                  <button
                    key={speed}
                    onClick={() => setPlaybackSpeed(speed)}
                    className={`px-2 py-1 rounded-md text-xs font-bold transition-all ${
                      playbackSpeed === speed
                        ? 'bg-cyan-100 text-cyan-600'
                        : 'text-slate-400 hover:text-slate-600'
                    }`}
                  >
                    {speed}x
                  </button>
                ))}
              </div>

              {/* Progress slider */}
              <div className="flex-1 min-w-40 mx-2">
                <input
                  type="range"
                  aria-label={language === 'es' ? 'Timestep de reproducción' : 'Playback timestep'}
                  min="0"
                  max={maxIndex}
                  value={playbackIndex}
                  onChange={(e) => { setPlaybackIndex(parseInt(e.target.value)); setIsPlaying(false); }}
                  className="w-full accent-cyan-500"
                />
              </div>

              {/* Current info */}
              {currentRow && (
                <div className="flex items-center gap-2 text-xs">
                  <Badge variant="info">{String(currentRow.Phase)}</Badge>
                  <span className="text-slate-500">T{String(currentRow.Trial)}</span>
                  <span className="text-slate-400">t{String(currentRow.TimeStep)}</span>
                </div>
              )}
            </div>

            {/* Animated network */}
            <div className="h-[450px]">
              {((
                <PublicationNetwork npes={npes} connections={connections} positions={publicationPositions}
                  activations={Object.fromEntries(npes.map(npe => [npe.name, Number(currentRow?.[npe.name] ?? npe.activation)]))}
                  weights={Object.fromEntries(connections.map(connection => {
                    const key = `${connection.presynapticNPE}-${connection.postsynapticNPE}`;
                    return [key, Number(currentRow?.[key] ?? connection.weight)];
                  }))}
                  showValues selectedUnit={selectedInspectorUnit} onSelectUnit={setSelectedInspectorUnit} language={language} />
              ))}
            </div>
            {<p className="px-4 py-2 text-xs text-slate-600 border-t border-slate-200">
              {language === 'es' ? 'Red 1 · Estado al final del timestep. Selecciona una unidad para inspeccionar su cálculo.' : 'Network 1 · State at the end of the timestep. Select a unit to inspect its calculation.'}
            </p>}
          </Card>
          {simulationInspector && (
            <EquationInspector trace={simulationInspector} rowIndex={playbackIndex} onRowChange={inspectRow}
              selectedUnit={selectedInspectorUnit} onSelectUnit={setSelectedInspectorUnit} language={language} />
          )}
          {!simulationInspector && <p className="text-sm text-slate-600">
            {language === 'es' ? 'Esta ejecución no tiene registro de ecuaciones. Activa «Registrar ecuaciones paso a paso» y ejecuta de nuevo para inspeccionarlas.' : 'This run has no equation trace. Enable “Record equations step by step” and run again to inspect them.'}
          </p>}
        </motion.div>
      )}
    </div>
  );
}
