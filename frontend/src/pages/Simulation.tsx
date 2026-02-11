import { useState, useCallback, useEffect, useRef, useMemo } from 'react';
import { useNavigate } from 'react-router-dom';
import { motion, AnimatePresence } from 'framer-motion';
import {
  ReactFlow,
  Background,
  type Node,
  type Edge,
  MarkerType,
} from '@xyflow/react';
import '@xyflow/react/dist/style.css';
import {
  Play, Pause, SkipForward, RotateCcw, CheckCircle2, AlertCircle,
  Loader2, Brain, Activity, Gauge, ChevronRight, Download, Upload, ChevronDown,
} from 'lucide-react';
import { Card } from '../components/ui/Card';
import { StatCard } from '../components/ui/StatCard';
import { Badge } from '../components/ui/Badge';
import { Tooltip } from '../components/ui/Tooltip';
import { useSimStore } from '../stores/useSimStore';
import { useConfirm } from '../components/ui/ConfirmDialog';
import { useToast } from '../components/ui/Toast';
import { getDisplayName } from '../utils/displayNames';
import { runSimulationOne } from '../api/client';
import { downloadArchitectureJSON } from '../utils/dataExport';
import { useI18n } from '../i18n';

const layerPositions: Record<string, { x: number; y: number }> = {
  US: { x: 50, y: 350 },
  PrimarySensory: { x: 50, y: 100 },
  AssociativeSensory: { x: 250, y: 100 },
  Hippocampal: { x: 250, y: 320 },
  Dopaminergic: { x: 450, y: 350 },
  AssociativeMotor: { x: 450, y: 150 },
  PrimaryMotor: { x: 650, y: 150 },
};

function getNodeShape(layer: string, type: string): string {
  if (layer === 'US') return '8px';
  if (layer === 'PrimarySensory' || layer === 'PrimaryMotor') return '4px';
  if (type === 'Inhibitory') return '2px';
  if (layer === 'Dopaminergic') return '12px';
  return '50%';
}

function activationColor(val: number): string {
  if (val < 0.3) return '#3b82f6';
  if (val < 0.6) return '#eab308';
  return '#ef4444';
}

function activationGlow(val: number): string {
  const color = activationColor(val);
  return `0 0 ${Math.round(val * 30)}px ${color}88`;
}

const pupdateOptions = [
  { value: 'async_random', labelKey: 'asyncRandom' as const, tooltipKey: 'asyncRandomTooltip' as const },
  { value: 'async_sequential', labelKey: 'asyncSequential' as const, tooltipKey: 'asyncSequentialTooltip' as const },
  { value: 'sync_random', labelKey: 'syncRandom' as const, tooltipKey: 'syncRandomTooltip' as const },
  { value: 'sync_sequential', labelKey: 'syncSequential' as const, tooltipKey: 'syncSequentialTooltip' as const },
];

function PUpdateDropdown({ pupdate, setPupdate, t }: { pupdate: string; setPupdate: (v: string) => void; t: any }) {
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
    npes, connections, trials, contingencies, hasITI,
    numNetworks, thresholdPreset, disc, pupdate, appMode, simStatus, simError, simulationResults, simulationMetadata,
    getThresholdType,
    setNumNetworks, setThresholdPreset, setDisc, setPupdate, setSimStatus, setSimResults, setSimError,
    playbackIndex, isPlaying, playbackSpeed,
    setPlaybackIndex, setIsPlaying, setPlaybackSpeed,
    lockedLayout,
  } = useSimStore();
  const { t } = useI18n();
  const { confirm } = useConfirm();
  const toast = useToast();

  const [progress, setProgress] = useState(0);
  const [saveLoadMsg, setSaveLoadMsg] = useState<string | null>(null);
  const [isCancelling, setIsCancelling] = useState(false);
  const [simStartTime, setSimStartTime] = useState<number | null>(null);
  const playbackRef = useRef<number | null>(null);
  const fileInputRef = useRef<HTMLInputElement>(null);
  const progressIntervalRef = useRef<number | null>(null);
  const cancelledRef = useRef(false);

  // Clean up progress interval on unmount (prevents memory leak if user navigates away during simulation)
  useEffect(() => {
    return () => {
      if (progressIntervalRef.current) clearInterval(progressIntervalRef.current);
    };
  }, []);

  const canRun = npes.length >= 2 && connections.length >= 1 && Object.keys(trials).length >= 1 && contingencies.length >= 1;

  const handleSaveExperiment = useCallback(() => {
    downloadArchitectureJSON(
      npes, connections, trials, contingencies, hasITI,
      `ddm-experiment-${new Date().toISOString().slice(0, 10)}.json`,
      { numNetworks, thresholdPreset, disc, pupdate },
      lockedLayout
    );
    setSaveLoadMsg(t.sim.experimentSaved);
    setTimeout(() => setSaveLoadMsg(null), 3000);
  }, [npes, connections, trials, contingencies, hasITI, numNetworks, thresholdPreset, disc, pupdate, lockedLayout, t]);

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

    setSimStatus('running');
    setProgress(0);
    setCompletedNetworks(0);
    cancelledRef.current = false;
    setIsCancelling(false);
    const startMs = Date.now();
    setSimStartTime(startMs);

    // Clear any leftover interval
    if (progressIntervalRef.current) clearInterval(progressIntervalRef.current);

    // Build R-compatible data (shared across all network calls)
    const npeData: any = {
      NPE: npes.map(n => n.name),
      Type: npes.map(n => n.type),
      Layer: npes.map(n => n.layer),
      Activation: npes.map(n => n.activation),
      'Temporal.Summation': npes.map(n => n.temporalSummation),
      'Activation.Decay': npes.map(n => n.activationDecay),
      mu: npes.map(n => n.mu),
      sigma: npes.map(n => n.sigma),
      logisSigma: npes.map(n => n.logisSigma),
    };

    const connData: any = {
      PreSinapticNPE: connections.map(c => c.presynapticNPE),
      PostSinapticNPE: connections.map(c => c.postsynapticNPE),
      Weight: connections.map(c => c.weight),
      alpha: connections.map(c => c.alpha),
      beta: connections.map(c => c.beta),
      alpha_prime: connections.map(c => c.alphaPrime),
      beta_prime: connections.map(c => c.betaPrime),
    };

    const baseParams = {
      npes: npeData,
      connections: connData,
      trials,
      contingencies,
      hasITI,
      threshold: getThresholdType(),
      disc,
      pupdate,
    };

    try {
      // Run networks one at a time for real progress tracking
      const results: any[] = [];

      for (let i = 0; i < numNetworks; i++) {
        // Check if user cancelled
        if (cancelledRef.current) {
          setSimStatus('idle');
          setIsCancelling(false);
          setSimStartTime(null);
          toast.info(t.toast.simulationCancelled);
          return;
        }

        const res = await runSimulationOne(baseParams);

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
        const done = i + 1;
        setCompletedNetworks(done);
        // Real progress: percentage of completed networks
        setProgress(Math.round((done / numNetworks) * 100));
      }

      const duration = ((Date.now() - startMs) / 1000).toFixed(2);
      setSimStartTime(null);

      // Identify unit columns vs connection columns vs signal columns
      const allCols = Object.keys(results[0][0] || {});
      const metaCols = ['Phase', 'Trial', 'TimeStep'];
      const signalCols = ['dVTA', 'dH'];
      const dataCols = allCols.filter(c => !metaCols.includes(c) && !signalCols.includes(c));
      const unitCols = dataCols.filter(c => !c.includes('-'));
      const connectionCols = dataCols.filter(c => c.includes('-'));

      setSimResults(results, {
        numNetworks,
        phases: Array.from(new Set(results[0].map((r: any) => String(r.Phase)))),
        units: unitCols,
        connections: connectionCols,
        totalTrials: Math.max(...results[0].map((r: any) => Number(r.Trial))),
        totalTimesteps: results[0].length,
        duration: parseFloat(duration),
        disc,
      });
      setPlaybackIndex(0);
    } catch (err: any) {
      setSimStartTime(null);
      if (cancelledRef.current) {
        setSimStatus('idle');
        setIsCancelling(false);
        toast.info(t.toast.simulationCancelled);
        return;
      }
      let errorMsg = err.message || t.sim.errorFallback;
      if (errorMsg.includes('timed out')) {
        errorMsg = t.sim.errorTimeout.replace('{networks}', String(numNetworks)).replace('{phases}', String(contingencies.length));
      } else if (errorMsg.includes('Failed to fetch') || errorMsg.includes('NetworkError')) {
        errorMsg = t.sim.errorNetwork;
      } else if (errorMsg.includes('500')) {
        errorMsg = t.sim.errorInternal;
      }
      setSimError(errorMsg);
    }
  }, [canRun, npes, connections, trials, contingencies, hasITI, numNetworks, thresholdPreset, getThresholdType, disc, pupdate, setSimStatus, setSimResults, setSimError, setPlaybackIndex, t, toast]);

  // Live countdown ETA — ticks every second so it never "freezes"
  const [etaText, setEtaText] = useState('');
  useEffect(() => {
    if (!simStartTime || completedNetworks === 0 || completedNetworks >= numNetworks) {
      setEtaText('');
      return;
    }
    // Compute per-network average once when completedNetworks changes
    const elapsed = Date.now() - simStartTime;
    const perNetwork = elapsed / completedNetworks;

    const tick = () => {
      const now = Date.now();
      const totalElapsed = now - simStartTime;
      const estimated = perNetwork * numNetworks;
      const remaining = Math.max(0, estimated - totalElapsed);
      if (remaining < 500) { setEtaText(''); return; }
      const secs = Math.round(remaining / 1000);
      if (secs < 60) {
        setEtaText(t.sim.etaRemaining.replace('{time}', `${secs}s`));
      } else {
        const mins = Math.floor(secs / 60);
        const remSecs = secs % 60;
        setEtaText(t.sim.etaRemaining.replace('{time}', `${mins}m ${remSecs}s`));
      }
    };
    tick(); // immediate update
    const id = window.setInterval(tick, 1000);
    return () => clearInterval(id);
  }, [simStartTime, completedNetworks, numNetworks, t]);

  // Playback logic
  const currentData = useMemo(() => {
    if (!simulationResults || simulationResults.length === 0) return null;
    return simulationResults[0];
  }, [simulationResults]);

  const maxIndex = currentData ? currentData.length - 1 : 0;

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
  const animatedNodes: Node[] = useMemo(() => {
    // Pre-count NPEs per layer for proper spacing (fallback when no locked layout)
    const layerTotals: Record<string, number> = {};
    npes.forEach(npe => { layerTotals[npe.layer] = (layerTotals[npe.layer] || 0) + 1; });
    const layerIdx: Record<string, number> = {};

    return npes.map((npe) => {
      const idx = layerIdx[npe.layer] || 0;
      layerIdx[npe.layer] = idx + 1;
      const total = layerTotals[npe.layer] || 1;

      // Use locked layout from NetworkBuilder if available, otherwise fallback to default
      let pos: { x: number; y: number };
      if (lockedLayout && lockedLayout[npe.name]) {
        pos = lockedLayout[npe.name];
      } else {
        const basePos = layerPositions[npe.layer] || { x: 300, y: 300 };
        const yOffset = total > 1 ? (idx - (total - 1) / 2) * 80 : 0;
        pos = { x: basePos.x, y: basePos.y + yOffset };
      }

      let activation = 0;
      if (currentData && playbackIndex < currentData.length) {
        const row = currentData[playbackIndex];
        const val = row[npe.name];
        if (typeof val === 'number') activation = val;
      }

      const color = activationColor(activation);

      return {
        id: npe.name,
        position: pos,
        data: {
          label: (
            <div className="text-center">
              <div className="text-[10px] font-bold italic">{getDisplayName(npe.name)}</div>
              <div className="text-[9px] opacity-70">{activation.toFixed(2)}</div>
            </div>
          ),
        },
        style: {
          background: color,
          color: '#fff',
          border: `3px solid ${color}`,
          borderRadius: getNodeShape(npe.layer, npe.type),
          width: 65,
          height: 65,
          display: 'flex',
          alignItems: 'center',
          justifyContent: 'center',
          fontSize: '11px',
          fontWeight: 700,
          boxShadow: activationGlow(activation),
          transition: 'all 0.15s ease',
        },
      };
    });
  }, [npes, currentData, playbackIndex, lockedLayout]);

  const animatedEdges: Edge[] = useMemo(() => {
    return connections.map((conn) => {
      let weight = conn.weight;
      if (currentData && playbackIndex < currentData.length) {
        const row = currentData[playbackIndex];
        const connName = `${conn.presynapticNPE}-${conn.postsynapticNPE}`;
        const val = row[connName];
        if (typeof val === 'number') weight = val;
      }

      return {
        id: `${conn.presynapticNPE}-${conn.postsynapticNPE}`,
        source: conn.presynapticNPE,
        target: conn.postsynapticNPE,
        animated: isPlaying,
        label: weight.toFixed(2),
        labelStyle: { fill: '#64748b', fontSize: 9, fontWeight: 700 },
        labelBgStyle: { fill: '#ffffff', fillOpacity: 0.85 },
        style: {
          stroke: weight >= 0.999 ? '#ef4444' : activationColor(weight),
          strokeWidth: Math.max(1.5, weight * 6),
          transition: 'all 0.15s ease',
        },
        markerEnd: { type: MarkerType.ArrowClosed, color: weight >= 0.999 ? '#ef4444' : '#64748b' },
      };
    });
  }, [connections, currentData, playbackIndex, isPlaying]);

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
        <StatCard icon={Gauge} label={t.sim.thresholdLabel} value={t.sim[`preset_${thresholdPreset}` as keyof typeof t.sim] || thresholdPreset} color="teal" />
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
          <div>
            <label className="flex items-center gap-1.5 text-xs font-semibold text-slate-500 mb-1">
              {t.sim.thresholdPreset}
              <Tooltip content={t.sim.thresholdPresetTooltip} />
            </label>
            <select
              value={thresholdPreset}
              onChange={async (e) => {
                const newPreset = e.target.value as any;
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
          </div>
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
                    {simulationMetadata?.numNetworks} {t.sim.networksSuffix} in {simulationMetadata?.duration}s
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
            <div className="flex items-center gap-3 p-3 border-b border-slate-200 bg-slate-50">
              <button
                onClick={() => setIsPlaying(!isPlaying)}
                className="w-10 h-10 rounded-xl bg-cyan-50 text-cyan-600 flex items-center justify-center hover:bg-cyan-200 transition-colors"
              >
                {isPlaying ? <Pause size={18} /> : <Play size={18} />}
              </button>
              <button
                onClick={() => { setPlaybackIndex(0); setIsPlaying(false); }}
                className="w-10 h-10 rounded-xl bg-slate-100 text-slate-500 flex items-center justify-center hover:bg-slate-100 transition-colors"
              >
                <RotateCcw size={16} />
              </button>
              <button
                onClick={() => setPlaybackIndex(Math.min(playbackIndex + 10, maxIndex))}
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
              <div className="flex-1 mx-4">
                <input
                  type="range"
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
              <ReactFlow
                nodes={animatedNodes}
                edges={animatedEdges}
                fitView
                proOptions={{ hideAttribution: true }}
                nodesDraggable={false}
                nodesConnectable={false}
                elementsSelectable={false}
              >
                <Background color="#e2e8f0" gap={20} />
              </ReactFlow>
            </div>
          </Card>
        </motion.div>
      )}
    </div>
  );
}
