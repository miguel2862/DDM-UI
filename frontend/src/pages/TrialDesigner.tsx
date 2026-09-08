import { useState, useCallback, useEffect, useRef } from 'react';
import { motion, AnimatePresence } from 'framer-motion';
import {
  ListChecks, Trash2, Clock, Shuffle, ArrowRight, Layers, Pencil, ArrowUp, ArrowDown,
  Timer, ToggleLeft, ToggleRight, Eraser, CheckSquare,
} from 'lucide-react';
import { Card } from '../components/ui/Card';
import { Badge } from '../components/ui/Badge';
import { Tooltip } from '../components/ui/Tooltip';
import { PageTransition } from '../components/ui/Skeleton';
import { useSimStore } from '../stores/useSimStore';
import { useConfirm } from '../components/ui/ConfirmDialog';
import { useI18n } from '../i18n';

export function TrialDesigner() {
  const {
    npes, trials, contingencies, addTrial, removeTrial, setContingencies, setHasITI,
  } = useSimStore();
  const { t } = useI18n();
  const { confirm } = useConfirm();

  // Trial form mode: 'trial' = normal, 'iti' = ITI builder
  const [trialMode, setTrialMode] = useState<'trial' | 'iti'>('trial');

  // Trial form
  const [trialName, setTrialName] = useState('');
  const [numTimesteps, setNumTimesteps] = useState(5);
  const [timestepData, setTimestepData] = useState<{ stimuli: Record<string, string>; learning: boolean }[]>([]);

  // ITI form (single timestep)
  const [itiName, setItiName] = useState('ITI');
  const [itiStimuli, setItiStimuli] = useState<Record<string, string>>({});
  const [itiLearning, setItiLearning] = useState(false);

  // Trial edit state
  const [editingTrial, setEditingTrial] = useState<string | null>(null);

  // Contingency form
  const [phaseName, setPhaseName] = useState('');
  const [trialOrder, setTrialOrder] = useState('Random');
  const [selectedTrialTypes, setSelectedTrialTypes] = useState<string[]>([]);
  const [trialCounts, setTrialCounts] = useState<string>('100');
  // "Reset Activations" = original R logic: true means HasITI=false, false means HasITI=true
  const [resetActivations, setResetActivations] = useState(true);
  const [itiMinTrials, setItiMinTrials] = useState(30);
  const [itiMaxTrials, setItiMaxTrials] = useState(30);
  const [itiTrialName, setItiTrialName] = useState('');

  // Phase edit state
  const [editingPhaseIdx, setEditingPhaseIdx] = useState<number | null>(null);

  const isOutcomeLayer = (layer: string) => layer === 'US';
  const primarySensoryNPEs = npes.filter(n =>
    n.layer === 'PrimarySensory' || isOutcomeLayer(n.layer)
  );

  // Find ITI trials (single timestep, typically used as ITI)
  const itiTrialNames = Object.entries(trials)
    .filter(([, ts]) => ts.length === 1)
    .map(([name]) => name);

  // Presentation order label mapping (values stay as-is for R backend, display is translated)
  const orderLabels: Record<string, string> = {
    'Random': t.trial.random,
    'In bulk': t.trial.inBulk,
    'Alternated': t.trial.alternated,
  };

  // Auto-initialize timestep grid when NPEs are available and no data exists yet
  // Using a ref to track if we've auto-inited so it only happens once
  const autoInitDone = useRef(false);
  useEffect(() => {
    if (primarySensoryNPEs.length > 0 && timestepData.length === 0 && editingTrial === null && trialMode === 'trial' && !autoInitDone.current) {
      autoInitDone.current = true;
      const steps = [];
      for (let i = 0; i < numTimesteps; i++) {
        const stimuli: Record<string, string> = {};
        primarySensoryNPEs.forEach(npe => {
          stimuli[npe.name] = i === numTimesteps - 1 && isOutcomeLayer(npe.layer) ? '1.00' : isOutcomeLayer(npe.layer) ? '0.00' : '1.00';
        });
        steps.push({ stimuli, learning: true });
      }
      setTimestepData(steps);
    }
    // Reset auto-init flag when NPEs change (e.g., new template loaded)
    if (primarySensoryNPEs.length === 0) {
      autoInitDone.current = false;
    }
  }, [primarySensoryNPEs, numTimesteps, timestepData.length, editingTrial, trialMode]);

  // Auto-initialize ITI stimuli when NPEs change
  useEffect(() => {
    if (primarySensoryNPEs.length > 0 && Object.keys(itiStimuli).length === 0) {
      const stim: Record<string, string> = {};
      primarySensoryNPEs.forEach(npe => {
        stim[npe.name] = '0.00';
      });
      setItiStimuli(stim);
    }
  }, [primarySensoryNPEs.length]); // eslint-disable-line react-hooks/exhaustive-deps

  // Auto-set ITI trial name when an ITI trial is available
  useEffect(() => {
    if (!resetActivations && !itiTrialName && itiTrialNames.length > 0) {
      setItiTrialName(itiTrialNames[0]);
    }
  }, [resetActivations, itiTrialName, itiTrialNames]);

  const initTimesteps = useCallback((count: number) => {
    const steps = [];
    for (let i = 0; i < count; i++) {
      const stimuli: Record<string, string> = {};
      primarySensoryNPEs.forEach(npe => {
        stimuli[npe.name] = i === count - 1 && isOutcomeLayer(npe.layer) ? '1.00' : isOutcomeLayer(npe.layer) ? '0.00' : '1.00';
      });
      steps.push({ stimuli, learning: true });
    }
    setTimestepData(steps);
    setNumTimesteps(count);
  }, [primarySensoryNPEs]);

  // --- Bulk actions ---
  const handleToggleAllLearning = useCallback((enabled: boolean) => {
    setTimestepData(prev => prev.map(ts => ({ ...ts, learning: enabled })));
  }, []);

  const handleFillAllStimuli = useCallback(() => {
    setTimestepData(prev => prev.map((ts, i) => {
      const stimuli: Record<string, string> = {};
      primarySensoryNPEs.forEach(npe => {
        stimuli[npe.name] = i === prev.length - 1 && isOutcomeLayer(npe.layer) ? '1.00' : isOutcomeLayer(npe.layer) ? '0.00' : '1.00';
      });
      return { ...ts, stimuli };
    }));
  }, [primarySensoryNPEs]);

  const handleClearAllStimuli = useCallback(() => {
    setTimestepData(prev => prev.map(ts => {
      const stimuli: Record<string, string> = {};
      primarySensoryNPEs.forEach(npe => {
        stimuli[npe.name] = '0.00';
      });
      return { ...ts, stimuli };
    }));
  }, [primarySensoryNPEs]);

  const handleAddTrial = useCallback(() => {
    if (!trialName.trim() || timestepData.length === 0) return;
    const timestepStrings = timestepData.map(ts => {
      const parts: string[] = [];
      Object.entries(ts.stimuli).forEach(([name, val]) => {
        parts.push(name, val);
      });
      parts.push(ts.learning ? 'True' : 'False');
      return parts.join(',');
    });
    addTrial(trialName.trim(), timestepStrings);
    setTrialName('');
    setTimestepData([]);
    setEditingTrial(null);
  }, [trialName, timestepData, addTrial]);

  const handleAddITI = useCallback(() => {
    if (!itiName.trim()) return;
    // ITI is a single-timestep trial
    const parts: string[] = [];
    Object.entries(itiStimuli).forEach(([name, val]) => {
      parts.push(name, val);
    });
    parts.push(itiLearning ? 'True' : 'False');
    addTrial(itiName.trim(), [parts.join(',')]);
    // Reset
    setItiName('ITI');
    const stim: Record<string, string> = {};
    primarySensoryNPEs.forEach(npe => { stim[npe.name] = '0.00'; });
    setItiStimuli(stim);
    setItiLearning(false);
  }, [itiName, itiStimuli, itiLearning, primarySensoryNPEs, addTrial]);

  const handleEditTrial = useCallback((name: string, timesteps: string[]) => {
    setTrialMode('trial');
    setEditingTrial(name);
    setTrialName(name);
    const parsed = timesteps.map(tsStr => {
      const parts = tsStr.split(',');
      const stimuli: Record<string, string> = {};
      for (let i = 0; i < parts.length - 1; i += 2) {
        stimuli[parts[i]] = parts[i + 1];
      }
      const learning = parts[parts.length - 1] === 'True';
      return { stimuli, learning };
    });
    setTimestepData(parsed);
    setNumTimesteps(parsed.length);
  }, []);

  const handleCancelEditTrial = useCallback(() => {
    setEditingTrial(null);
    setTrialName('');
    setTimestepData([]);
  }, []);

  const handleAddContingency = useCallback(() => {
    if (!phaseName.trim() || selectedTrialTypes.length === 0) return;
    const trialTypesStr = selectedTrialTypes.join('/');
    // resetActivations=true → hasITI=false, resetActivations=false → hasITI=true
    const hasITIVal = !resetActivations;
    let spec = `${phaseName.trim()}, ${trialOrder}, ${trialTypesStr}, ${trialCounts}, ${hasITIVal ? 'True' : 'False'}`;
    if (hasITIVal) {
      spec += `, ${itiMinTrials}, ${itiMaxTrials}, ${itiTrialName}`;
    }
    if (editingPhaseIdx !== null) {
      const updated = [...contingencies];
      updated[editingPhaseIdx] = spec;
      setContingencies(updated);
      const store = useSimStore.getState();
      const updatedHasITI = [...(store.hasITI || contingencies.map(() => false))];
      updatedHasITI[editingPhaseIdx] = hasITIVal;
      setHasITI(updatedHasITI);
      setEditingPhaseIdx(null);
    } else {
      setContingencies([...contingencies, spec]);
      const store = useSimStore.getState();
      setHasITI([...(store.hasITI || contingencies.map(() => false)), hasITIVal]);
    }

    setPhaseName('');
    setSelectedTrialTypes([]);
    setTrialCounts('100');
    setTrialOrder('Random');
    setResetActivations(true);
    setItiMinTrials(30);
    setItiMaxTrials(30);
    setItiTrialName('');
  }, [phaseName, trialOrder, selectedTrialTypes, trialCounts, resetActivations, itiMinTrials, itiMaxTrials, itiTrialName,
      contingencies, setContingencies, setHasITI, editingPhaseIdx]);

  const handleEditPhase = useCallback((spec: string, idx: number) => {
    const parts = spec.split(',').map(s => s.trim());
    setPhaseName(parts[0] || '');
    setTrialOrder(parts[1] || 'Random');
    setSelectedTrialTypes(parts[2] ? parts[2].split('/') : []);
    setTrialCounts(parts[3] || '100');
    const hasItiVal = parts[4]?.toLowerCase() === 'true';
    setResetActivations(!hasItiVal); // Invert: hasITI=true → resetActivations=false
    if (hasItiVal && parts.length >= 8) {
      setItiMinTrials(parseInt(parts[5]) || 30);
      setItiMaxTrials(parseInt(parts[6]) || 30);
      setItiTrialName(parts[7] || '');
    } else {
      setItiMinTrials(30);
      setItiMaxTrials(30);
      setItiTrialName('');
    }

    setEditingPhaseIdx(idx);
  }, []);

  const handleCancelEditPhase = useCallback(() => {
    setEditingPhaseIdx(null);
    setPhaseName('');
    setSelectedTrialTypes([]);
    setTrialCounts('100');
    setTrialOrder('Random');
    setResetActivations(true);
    setItiMinTrials(30);
    setItiMaxTrials(30);
    setItiTrialName('');
  }, []);

  const handleMovePhase = useCallback((idx: number, direction: 'up' | 'down') => {
    const newIdx = direction === 'up' ? idx - 1 : idx + 1;
    if (newIdx < 0 || newIdx >= contingencies.length) return;
    const updated = [...contingencies];
    [updated[idx], updated[newIdx]] = [updated[newIdx], updated[idx]];
    setContingencies(updated);
    // Keep hasITI array in sync
    const store = useSimStore.getState();
    const updatedITI = [...(store.hasITI || contingencies.map(() => false))];
    [updatedITI[idx], updatedITI[newIdx]] = [updatedITI[newIdx], updatedITI[idx]];
    setHasITI(updatedITI);
    if (editingPhaseIdx === idx) {
      setEditingPhaseIdx(newIdx);
    } else if (editingPhaseIdx === newIdx) {
      setEditingPhaseIdx(idx);
    }
  }, [contingencies, setContingencies, setHasITI, editingPhaseIdx]);

  const updateTimestepStimulus = (tsIndex: number, npeName: string, value: string) => {
    setTimestepData(prev => {
      const copy = [...prev];
      copy[tsIndex] = { ...copy[tsIndex], stimuli: { ...copy[tsIndex].stimuli, [npeName]: value } };
      return copy;
    });
  };

  return (
    <PageTransition>
    <div className="max-w-6xl mx-auto space-y-6">
      <div className="flex items-center gap-3">
        <div className="w-10 h-10 rounded-xl bg-amber-500/10 flex items-center justify-center">
          <ListChecks size={20} className="text-amber-600" />
        </div>
        <div>
          <h1 className="text-2xl font-bold text-slate-800">{t.trial.pageTitle}</h1>
          <p className="text-sm text-slate-500">{t.trial.pageSubtitle}</p>
        </div>
      </div>

      <div className="grid grid-cols-1 lg:grid-cols-2 gap-6">
        {/* Trial Builder */}
        <div className="space-y-4">
          {/* Mode toggle: Trial / ITI */}
          <div className="flex gap-1 p-1 rounded-xl bg-slate-100 border border-slate-200">
            <button
              onClick={() => { setTrialMode('trial'); setEditingTrial(null); }}
              className={`flex-1 flex items-center justify-center gap-2 py-2 rounded-lg text-sm font-bold transition-all ${
                trialMode === 'trial' ? 'bg-amber-50 text-amber-600 shadow-sm' : 'text-slate-500 hover:text-slate-400'
              }`}
            >
              <Clock size={14} /> {t.trial.addTrial}
            </button>
            <button
              onClick={() => { setTrialMode('iti'); setEditingTrial(null); }}
              className={`flex-1 flex items-center justify-center gap-2 py-2 rounded-lg text-sm font-bold transition-all ${
                trialMode === 'iti' ? 'bg-violet-50 text-violet-600 shadow-sm' : 'text-slate-500 hover:text-slate-400'
              }`}
            >
              <Timer size={14} /> {t.trial.addITI}
            </button>
          </div>

          <AnimatePresence mode="wait">
            {trialMode === 'trial' ? (
              <motion.div key="trial-form" initial={{ opacity: 0, x: -10 }} animate={{ opacity: 1, x: 0 }} exit={{ opacity: 0, x: 10 }}>
                <Card className="p-5 space-y-4">
                  <h3 className="text-sm font-bold text-slate-700 flex items-center gap-2">
                    <Clock size={14} className="text-amber-600" /> {t.trial.createTrial}
                  </h3>

                  <div>
                    <label className="flex items-center gap-1.5 text-xs font-semibold text-slate-500 mb-1">
                      {t.trial.trialName}
                      <Tooltip content={t.trial.trialNameTooltip} />
                    </label>
                    <input
                      type="text"
                      value={trialName}
                      onChange={(e) => setTrialName(e.target.value)}
                      disabled={editingTrial !== null}
                      placeholder={t.trial.trialNamePlaceholder}
                      className={`w-full px-3 py-2 rounded-lg bg-white border border-slate-200 text-slate-800 text-sm focus:border-cyan-500/50 focus:outline-none ${editingTrial !== null ? 'opacity-60 cursor-not-allowed' : ''}`}
                    />
                  </div>

                  <div>
                    <label className="flex items-center gap-1.5 text-xs font-semibold text-slate-500 mb-1">
                      {t.trial.timesteps}: {numTimesteps}
                      <Tooltip content={t.trial.timestepsTooltip} />
                    </label>
                    <div className="flex items-center gap-2">
                      <input
                        type="range"
                        min="1"
                        max="10"
                        value={numTimesteps}
                        onChange={(e) => initTimesteps(parseInt(e.target.value))}
                        className="flex-1 accent-amber-500"
                      />
                      <button
                        onClick={() => initTimesteps(numTimesteps)}
                        className="px-3 py-1.5 rounded-lg bg-amber-50 text-amber-600 text-xs font-bold border border-amber-200 hover:bg-amber-500/25 transition-colors"
                      >
                        {t.trial.generate}
                      </button>
                    </div>
                  </div>

                  {/* Timestep grid */}
                  {timestepData.length > 0 && (
                    <div className="space-y-3">
                      <div className="text-xs font-semibold text-slate-400">{t.trial.stimulusActivations}</div>
                      {/* Bulk action buttons — single compact row */}
                      <div className="flex gap-1.5">
                        <button
                          onClick={handleFillAllStimuli}
                          title={t.trial.fillTitle}
                          className="flex items-center gap-1 py-1 px-2 rounded-md text-[10px] font-bold border border-amber-200 bg-amber-50 text-amber-600 hover:bg-amber-100 transition-colors"
                        >
                          <CheckSquare size={11} />
                          {t.trial.fill}
                        </button>
                        <button
                          onClick={handleClearAllStimuli}
                          title={t.trial.clearTitle}
                          className="flex items-center gap-1 py-1 px-2 rounded-md text-[10px] font-bold border border-slate-200 bg-slate-50 text-slate-500 hover:bg-slate-100 transition-colors"
                        >
                          <Eraser size={11} />
                          {t.trial.clear}
                        </button>
                        <button
                          onClick={() => handleToggleAllLearning(true)}
                          title={t.trial.learnOnTitle}
                          className="flex items-center gap-1 py-1 px-2 rounded-md text-[10px] font-bold border border-emerald-200 bg-emerald-50 text-emerald-600 hover:bg-emerald-100 transition-colors"
                        >
                          <ToggleRight size={11} />
                          {t.trial.learnOn}
                        </button>
                        <button
                          onClick={() => handleToggleAllLearning(false)}
                          title={t.trial.learnOffTitle}
                          className="flex items-center gap-1 py-1 px-2 rounded-md text-[10px] font-bold border border-rose-200 bg-rose-50 text-rose-500 hover:bg-rose-100 transition-colors"
                        >
                          <ToggleLeft size={11} />
                          {t.trial.learnOff}
                        </button>
                      </div>
                      <div className="overflow-x-auto">
                        <table className="w-full text-xs">
                          <thead>
                            <tr className="text-slate-400">
                              <th className="text-left py-1 pr-2">{t.trial.tHeader}</th>
                              {primarySensoryNPEs.map(npe => (
                                <th key={npe.name} className="text-center py-1 px-1">{npe.name}</th>
                              ))}
                              <th className="text-center py-1 px-1">{t.trial.learnHeader}</th>
                            </tr>
                          </thead>
                          <tbody>
                            {timestepData.map((ts, i) => (
                              <tr key={i} className="border-t border-slate-100">
                                <td className="py-1.5 pr-2 font-bold text-slate-500">{i + 1}</td>
                                {primarySensoryNPEs.map(npe => (
                                  <td key={npe.name} className="py-1.5 px-1">
                                    <input
                                      type="number"
                                      step="0.1"
                                      min="0"
                                      max="1"
                                      value={ts.stimuli[npe.name] || '0'}
                                      onChange={(e) => updateTimestepStimulus(i, npe.name, e.target.value)}
                                      className="w-14 px-1.5 py-1 rounded bg-white border border-slate-200 text-slate-700 text-center text-xs focus:border-cyan-500/50 focus:outline-none"
                                    />
                                  </td>
                                ))}
                                <td className="py-1.5 px-1 text-center">
                                  <input
                                    type="checkbox"
                                    checked={ts.learning}
                                    onChange={(e) => {
                                      setTimestepData(prev => {
                                        const copy = [...prev];
                                        copy[i] = { ...copy[i], learning: e.target.checked };
                                        return copy;
                                      });
                                    }}
                                    className="accent-cyan-500"
                                  />
                                </td>
                              </tr>
                            ))}
                          </tbody>
                        </table>
                      </div>
                    </div>
                  )}

                  <button
                    onClick={handleAddTrial}
                    disabled={!trialName.trim() || timestepData.length === 0}
                    className="w-full py-2.5 rounded-xl bg-gradient-to-r from-amber-500 to-amber-600 text-white font-bold text-sm disabled:opacity-40 disabled:cursor-not-allowed hover:shadow-lg hover:shadow-amber-500/20 transition-shadow"
                  >
                    {editingTrial ? t.trial.updateTrial : t.trial.addTrial}
                  </button>
                  {editingTrial && (
                    <button
                      onClick={handleCancelEditTrial}
                      className="w-full py-2 rounded-xl text-slate-500 text-sm font-bold hover:bg-slate-100 transition-colors"
                    >
                      {t.trial.cancelEdit}
                    </button>
                  )}
                </Card>
              </motion.div>
            ) : (
              <motion.div key="iti-form" initial={{ opacity: 0, x: 10 }} animate={{ opacity: 1, x: 0 }} exit={{ opacity: 0, x: -10 }}>
                <Card className="p-5 space-y-4">
                  <h3 className="text-sm font-bold text-slate-700 flex items-center gap-2">
                    <Timer size={14} className="text-violet-600" /> {t.trial.configureITI}
                  </h3>
                  <p className="text-[11px] text-slate-400 leading-relaxed">
                    {t.trial.itiDescription}
                  </p>

                  <div>
                    <label className="flex items-center gap-1.5 text-xs font-semibold text-slate-500 mb-1">
                      {t.trial.itiName}
                      <Tooltip content={t.trial.itiNameTooltip} />
                    </label>
                    <input
                      type="text"
                      value={itiName}
                      onChange={(e) => setItiName(e.target.value)}
                      placeholder={t.trial.itiNamePlaceholder}
                      className="w-full px-3 py-2 rounded-lg bg-white border border-slate-200 text-slate-800 text-sm focus:border-violet-500/50 focus:outline-none"
                    />
                  </div>

                  {/* Single timestep ITI values */}
                  <div className="space-y-2">
                    <div className="text-xs font-semibold text-slate-400">{t.trial.npeActivations}</div>
                    <div className="grid grid-cols-2 gap-2">
                      {primarySensoryNPEs.map(npe => (
                        <div key={npe.name} className="flex items-center gap-2">
                          <label className="text-[10px] font-bold text-slate-500 w-10">{npe.name}</label>
                          <input
                            type="number"
                            step="0.1"
                            min="0"
                            max="1"
                            value={itiStimuli[npe.name] || '0'}
                            onChange={(e) => setItiStimuli(prev => ({ ...prev, [npe.name]: e.target.value }))}
                            className="flex-1 px-2 py-1.5 rounded-md bg-white border border-slate-200 text-slate-700 text-xs text-center focus:border-violet-500/50 focus:outline-none"
                          />
                        </div>
                      ))}
                    </div>
                    <div className="flex items-center gap-2 pt-1">
                      <input
                        type="checkbox"
                        checked={itiLearning}
                        onChange={(e) => setItiLearning(e.target.checked)}
                        className="accent-violet-500"
                      />
                      <label className="text-xs font-semibold text-slate-500">{t.trial.enableLearning}</label>
                    </div>
                  </div>

                  <button
                    onClick={handleAddITI}
                    disabled={!itiName.trim()}
                    className="w-full py-2.5 rounded-xl bg-gradient-to-r from-violet-500 to-violet-600 text-white font-bold text-sm disabled:opacity-40 disabled:cursor-not-allowed hover:shadow-lg hover:shadow-violet-500/20 transition-shadow"
                  >
                    {t.trial.addITI}
                  </button>
                </Card>
              </motion.div>
            )}
          </AnimatePresence>

          {/* Trial list */}
          <div className="space-y-2">
            <h4 className="text-xs font-bold text-slate-400 uppercase tracking-wider">{t.trial.definedTrials} ({Object.keys(trials).length})</h4>
            {Object.entries(trials).map(([name, timesteps]) => {
              const isItiTrial = timesteps.length === 1;
              return (
                <motion.div
                  key={name}
                  initial={{ opacity: 0, x: -10 }}
                  animate={{ opacity: 1, x: 0 }}
                  className="flex items-center justify-between p-3 rounded-xl bg-slate-50 border border-slate-200"
                >
                  <div className="flex items-center gap-3">
                    <div className={`w-8 h-8 rounded-lg flex items-center justify-center ${isItiTrial ? 'bg-violet-50' : 'bg-amber-50'}`}>
                      {isItiTrial ? <Timer size={14} className="text-violet-600" /> : <Clock size={14} className="text-amber-600" />}
                    </div>
                    <div>
                      <span className="text-sm font-bold text-slate-700">{name}</span>
                      <div className="flex gap-1.5 mt-0.5">
                        <Badge variant={isItiTrial ? 'muted' : 'warning'}>
                          {isItiTrial ? t.trial.timestepITI : `${timesteps.length} ${t.trial.timestepLabel}`}
                        </Badge>
                      </div>
                    </div>
                  </div>
                  <div className="flex items-center gap-1.5">
                    <button
                      onClick={() => handleEditTrial(name, timesteps)}
                      className="flex items-center gap-1 px-2.5 py-1.5 rounded-lg text-[11px] font-bold text-slate-500 border border-slate-200 hover:text-amber-600 hover:border-amber-200 hover:bg-amber-50 transition-colors"
                    >
                      <Pencil size={12} /> {t.trial.edit}
                    </button>
                    <button
                      onClick={async () => {
                        const isUsed = contingencies.some(c => c.includes(name));
                        if (isUsed) {
                          const ok = await confirm({
                            title: t.confirm.deleteTrialTitle,
                            message: t.confirm.deleteTrialMessage,
                            confirmLabel: t.common.delete,
                            cancelLabel: t.confirm.cancel,
                            variant: 'warning',
                          });
                          if (!ok) return;
                        }
                        removeTrial(name);
                      }}
                      className="flex items-center gap-1 px-2.5 py-1.5 rounded-lg text-[11px] font-bold text-slate-500 border border-slate-200 hover:text-rose-600 hover:border-rose-200 hover:bg-rose-50 transition-colors"
                    >
                      <Trash2 size={12} /> {t.trial.delete}
                    </button>
                  </div>
                </motion.div>
              );
            })}
          </div>
        </div>

        {/* Contingency Builder */}
        <div className="space-y-4">
          <Card className="p-5 space-y-4">
            <h3 className="text-sm font-bold text-slate-700 flex items-center gap-2">
              <Layers size={14} className="text-violet-600" /> {t.trial.configurePhase}
            </h3>

            <div>
              <label className="flex items-center gap-1.5 text-xs font-semibold text-slate-500 mb-1">
                {t.trial.phaseName}
                <Tooltip content={t.trial.phaseNameTooltip} />
              </label>
              <input
                type="text"
                value={phaseName}
                onChange={(e) => setPhaseName(e.target.value)}
                placeholder={t.trial.phaseNamePlaceholder}
                className="w-full px-3 py-2 rounded-lg bg-white border border-slate-200 text-slate-800 text-sm focus:border-cyan-500/50 focus:outline-none"
              />
            </div>



            <div>
              <label className="flex items-center gap-1.5 text-xs font-semibold text-slate-500 mb-1">
                {t.trial.presentationOrder}
                <Tooltip content={t.trial.presentationOrderTooltip} />
              </label>
              <div className="grid grid-cols-3 gap-2">
                {['Random', 'In bulk', 'Alternated'].map(order => (
                  <button
                    key={order}
                    onClick={() => setTrialOrder(order)}
                    className={`py-2 rounded-lg text-xs font-bold border transition-all ${
                      trialOrder === order
                        ? 'bg-violet-50 text-violet-600 border-violet-200'
                        : 'text-slate-500 border-slate-200 hover:bg-slate-100'
                    }`}
                  >
                    {order === 'Random' && <Shuffle size={12} className="inline mr-1" />}
                    {orderLabels[order] || order}
                  </button>
                ))}
              </div>
            </div>

            <div>
              <label className="block text-xs font-semibold text-slate-500 mb-1">{t.trial.trialTypesLabel}</label>
              <div className="flex flex-wrap gap-2">
                {Object.keys(trials).map(name => (
                  <button
                    key={name}
                    onClick={() => {
                      setSelectedTrialTypes(prev =>
                        prev.includes(name) ? prev.filter(t => t !== name) : [...prev, name]
                      );
                    }}
                    className={`px-3 py-1.5 rounded-lg text-xs font-bold border transition-all ${
                      selectedTrialTypes.includes(name)
                        ? 'bg-cyan-50 text-cyan-600 border-cyan-200'
                        : 'text-slate-500 border-slate-200 hover:bg-slate-100'
                    }`}
                  >
                    {name}
                  </button>
                ))}
                {Object.keys(trials).length === 0 && (
                  <span className="text-xs text-slate-400 italic">{t.trial.noTrialsDefined}</span>
                )}
              </div>
            </div>

            <div>
              <label className="flex items-center gap-1.5 text-xs font-semibold text-slate-500 mb-1">
                {t.trial.trialCounts}
                <Tooltip content={t.trial.trialCountsTooltip} />
              </label>
              <input
                type="text"
                value={trialCounts}
                onChange={(e) => setTrialCounts(e.target.value)}
                placeholder={t.trial.trialCountsPlaceholder}
                className="w-full px-3 py-2 rounded-lg bg-white border border-slate-200 text-slate-800 text-sm focus:border-cyan-500/50 focus:outline-none"
              />
            </div>

            {/* Reset Activations / ITI Configuration */}
            <div className="space-y-2 pt-2 border-t border-slate-200">
              <div className="flex items-center gap-2">
                <input
                  type="checkbox"
                  checked={resetActivations}
                  onChange={(e) => setResetActivations(e.target.checked)}
                  className="accent-violet-500"
                />
                <label className="flex items-center gap-1.5 text-xs font-semibold text-slate-500">
                  {t.trial.resetActivations}
                  <Tooltip content={t.trial.resetActivationsTooltip} />
                </label>
              </div>
              {!resetActivations && (
                <div className="grid grid-cols-3 gap-2 pl-5">
                  <div>
                    <label className="block text-[10px] font-semibold text-slate-400 mb-0.5">{t.trial.minITI}</label>
                    <input
                      type="number"
                      min="1"
                      value={itiMinTrials}
                      onChange={(e) => setItiMinTrials(parseInt(e.target.value) || 1)}
                      className="w-full px-2 py-1.5 rounded-md bg-white border border-slate-200 text-slate-700 text-xs focus:border-violet-500/50 focus:outline-none"
                    />
                  </div>
                  <div>
                    <label className="block text-[10px] font-semibold text-slate-400 mb-0.5">{t.trial.maxITI}</label>
                    <input
                      type="number"
                      min="1"
                      value={itiMaxTrials}
                      onChange={(e) => setItiMaxTrials(parseInt(e.target.value) || 1)}
                      className="w-full px-2 py-1.5 rounded-md bg-white border border-slate-200 text-slate-700 text-xs focus:border-violet-500/50 focus:outline-none"
                    />
                  </div>
                  <div>
                    <label className="block text-[10px] font-semibold text-slate-400 mb-0.5">{t.trial.itiTrial}</label>
                    <select
                      value={itiTrialName}
                      onChange={(e) => setItiTrialName(e.target.value)}
                      className="w-full px-2 py-1.5 rounded-md bg-white border border-slate-200 text-slate-700 text-xs focus:border-violet-500/50 focus:outline-none"
                    >
                      <option value="">Select...</option>
                      {Object.keys(trials).map(name => (
                        <option key={name} value={name}>{name}</option>
                      ))}
                    </select>
                  </div>
                  {itiTrialNames.length === 0 && (
                    <div className="col-span-3">
                      <p className="text-[10px] text-amber-500 font-semibold">
                        {t.trial.noITITrialDefined}
                      </p>
                    </div>
                  )}
                </div>
              )}
            </div>

            <button
              onClick={handleAddContingency}
              disabled={!phaseName.trim() || selectedTrialTypes.length === 0 || (!resetActivations && !itiTrialName)}
              className="w-full py-2.5 rounded-xl bg-gradient-to-r from-violet-500 to-violet-600 text-white font-bold text-sm disabled:opacity-40 disabled:cursor-not-allowed hover:shadow-lg hover:shadow-violet-500/20 transition-shadow"
            >
              {editingPhaseIdx !== null ? t.trial.updatePhase : t.trial.addPhase}
            </button>
            {!resetActivations && !itiTrialName && (
              <p className="text-[10px] text-amber-500 text-center">{t.trial.selectITITrial}</p>
            )}
            {editingPhaseIdx !== null && (
              <button
                onClick={handleCancelEditPhase}
                className="w-full py-2 rounded-xl text-slate-500 text-sm font-bold hover:bg-slate-100 transition-colors"
              >
                {t.trial.cancelEdit}
              </button>
            )}
          </Card>

          {/* Contingency list */}
          <div className="space-y-2">
            <h4 className="text-xs font-bold text-slate-400 uppercase tracking-wider">{t.trial.phasesLabel} ({contingencies.length})</h4>
            {contingencies.map((spec, i) => {
              const parts = spec.split(',').map(s => s.trim());
              const hasIti = parts[4]?.toLowerCase() === 'true';
              return (
                <motion.div
                  key={i}
                  initial={{ opacity: 0, x: 10 }}
                  animate={{ opacity: 1, x: 0 }}
                  className="p-4 rounded-xl bg-slate-50 border border-slate-200"
                >
                  <div className="flex items-center justify-between mb-2">
                    <div className="flex items-center gap-2">
                      <span className="w-6 h-6 rounded-full bg-violet-50 text-violet-600 text-xs font-bold flex items-center justify-center">
                        {i + 1}
                      </span>
                      <span className="text-sm font-bold text-slate-700">{parts[0]}</span>
                    </div>
                    <div className="flex items-center gap-1.5">
                      <div className="flex items-center gap-0.5 mr-1">
                        <button
                          onClick={() => handleMovePhase(i, 'up')}
                          disabled={i === 0}
                          title={t.trial.moveUp}
                          className="p-1.5 rounded-md text-slate-400 hover:text-violet-600 hover:bg-violet-50 transition-colors disabled:opacity-20 disabled:cursor-not-allowed"
                        >
                          <ArrowUp size={14} />
                        </button>
                        <button
                          onClick={() => handleMovePhase(i, 'down')}
                          disabled={i === contingencies.length - 1}
                          title={t.trial.moveDown}
                          className="p-1.5 rounded-md text-slate-400 hover:text-violet-600 hover:bg-violet-50 transition-colors disabled:opacity-20 disabled:cursor-not-allowed"
                        >
                          <ArrowDown size={14} />
                        </button>
                      </div>
                      <button
                        onClick={() => handleEditPhase(spec, i)}
                        className="flex items-center gap-1 px-2.5 py-1.5 rounded-lg text-[11px] font-bold text-slate-500 border border-slate-200 hover:text-amber-600 hover:border-amber-200 hover:bg-amber-50 transition-colors"
                      >
                        <Pencil size={12} /> {t.trial.edit}
                      </button>
                      <button
                        onClick={() => {

                          setContingencies(contingencies.filter((_, j) => j !== i));
                          // Keep hasITI array in sync with contingencies
                          const store = useSimStore.getState();
                          setHasITI((store.hasITI || []).filter((_, j) => j !== i));
                          if (editingPhaseIdx === i) setEditingPhaseIdx(null);
                          else if (editingPhaseIdx !== null && editingPhaseIdx > i) setEditingPhaseIdx(editingPhaseIdx - 1);
                        }}
                        className="flex items-center gap-1 px-2.5 py-1.5 rounded-lg text-[11px] font-bold text-slate-500 border border-slate-200 hover:text-rose-600 hover:border-rose-200 hover:bg-rose-50 transition-colors"
                      >
                        <Trash2 size={12} /> {t.trial.delete}
                      </button>
                    </div>
                  </div>
                  <div className="flex items-center gap-2 text-xs text-slate-500">
                    <Badge variant="info">{parts[1]}</Badge>
                    <ArrowRight size={12} />
                    <span className="font-mono">{parts[2]}</span>
                    <span>({parts[3]} {t.trial.trials})</span>
                    {hasIti ? (
                      <Badge variant="warning">ITI: {parts[5]}-{parts[6]} ({parts[7]})</Badge>
                    ) : (
                      <Badge variant="muted">{t.trial.reset}</Badge>
                    )}

                  </div>
                </motion.div>
              );
            })}
          </div>
        </div>
      </div>
    </div>
    </PageTransition>
  );
}
