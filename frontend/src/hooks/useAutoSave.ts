import { useEffect } from 'react';
import { useSimStore } from '../stores/useSimStore';

const AUTOSAVE_KEY = 'ddm-ui-autosave';

/**
 * Checks if a previous session exists in localStorage.
 * Returns the saved data (or null) without restoring it.
 */
export function getSavedSession(): any | null {
  try {
    const saved = localStorage.getItem(AUTOSAVE_KEY);
    if (!saved) return null;
    const data = JSON.parse(saved);
    if (!data || !data.npes || data.npes.length === 0) return null;
    return data;
  } catch {
    return null;
  }
}

/**
 * Restores a saved session into the Zustand store.
 */
export function restoreSession(data: any) {
  const store = useSimStore.getState();
  store.setNPEs(data.npes || []);
  store.setConnections(data.connections || []);
  store.setTrials(data.trials || {});
  store.setContingencies(data.contingencies || []);
  store.setHasITI(data.hasITI || []);
  store.setPhases(data.phases || []);
  if (data.numNetworks) store.setNumNetworks(data.numNetworks);
  if (data.thresholdPreset) store.setThresholdPreset(data.thresholdPreset);
  if (data.disc) store.setDisc(data.disc);
  if (data.pupdate) store.setPupdate(data.pupdate);
  if (data.appMode) store.setAppMode(data.appMode);
  // Restore simulation results if available
  if (data.simulationResults && data.simulationMetadata) {
    store.setSimResults(data.simulationResults, data.simulationMetadata);
  }
}

/**
 * Clears any saved session data from localStorage.
 */
export function clearAutoSave() {
  try {
    localStorage.removeItem(AUTOSAVE_KEY);
  } catch { /* ignore */ }
}

/**
 * Saves current state to localStorage.
 * Only saves if there's actual content (at least 1 NPE).
 */
function saveCurrentState() {
  try {
    const state = useSimStore.getState();
    if (state.npes.length === 0) return;

    const payload: Record<string, any> = {
      npes: state.npes,
      connections: state.connections,
      trials: state.trials,
      contingencies: state.contingencies,
      hasITI: state.hasITI,
      phases: state.phases,
      numNetworks: state.numNetworks,
      thresholdPreset: state.thresholdPreset,
      disc: state.disc,
      pupdate: state.pupdate,
      appMode: state.appMode,
      _timestamp: Date.now(),
    };

    // Include simulation results if available
    if (state.simulationResults && state.simulationMetadata) {
      payload.simulationResults = state.simulationResults;
      payload.simulationMetadata = state.simulationMetadata;
    }

    localStorage.setItem(AUTOSAVE_KEY, JSON.stringify(payload));
  } catch (e) {
    console.warn('[DDM-UI] Save failed:', e);
  }
}

/**
 * Hook that auto-saves the session when the user leaves/closes the page.
 * Does NOT auto-restore — use getSavedSession() + restoreSession() with
 * an explicit user prompt for that.
 */
export function useAutoSave() {
  useEffect(() => {
    const handleVisibility = () => {
      if (document.visibilityState === 'hidden') saveCurrentState();
    };
    const handleBeforeUnload = () => saveCurrentState();

    document.addEventListener('visibilitychange', handleVisibility);
    window.addEventListener('beforeunload', handleBeforeUnload);

    return () => {
      document.removeEventListener('visibilitychange', handleVisibility);
      window.removeEventListener('beforeunload', handleBeforeUnload);
    };
  }, []);
}
