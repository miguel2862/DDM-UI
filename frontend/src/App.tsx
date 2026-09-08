import { useState, useCallback } from 'react';
import { HashRouter, Routes, Route, Navigate } from 'react-router-dom';
import { motion, AnimatePresence } from 'framer-motion';
import { RotateCcw, Plus } from 'lucide-react';
import { Layout } from './components/layout/Layout';
import { Dashboard } from './pages/Dashboard';
import { NetworkBuilder } from './pages/NetworkBuilder';
import { TrialDesigner } from './pages/TrialDesigner';
import { Simulation } from './pages/Simulation';
import { Results } from './pages/Results';
import { ParameterSweep } from './pages/ParameterSweep';
import { Help } from './pages/Help';
import { SplashScreen } from './components/SplashScreen';
import { ErrorBoundary } from './components/ui/ErrorBoundary';
import { ToastProvider, ConnectionMonitor } from './components/ui/Toast';
import { ConfirmDialogProvider } from './components/ui/ConfirmDialog';
import { useAutoSave, getSavedSession, restoreSession, clearAutoSave } from './hooks/useAutoSave';
import { usePlaybackKeys } from './hooks/usePlaybackKeys';
import { useI18n } from './i18n';

type SavedSession = NonNullable<ReturnType<typeof getSavedSession>>;

function AppContent() {
  // Auto-save model state on page close/navigate away
  useAutoSave();

  // Global keyboard shortcuts for playback
  usePlaybackKeys();

  return (
    <HashRouter>
      <Routes>
        <Route element={<Layout />}>
          <Route path="/" element={<Dashboard />} />
          <Route path="/network" element={<NetworkBuilder />} />
          <Route path="/trials" element={<TrialDesigner />} />
          <Route path="/simulation" element={<Simulation />} />
          <Route path="/results" element={<Results />} />
          <Route path="/sweep" element={<ParameterSweep />} />
          <Route path="/help" element={<Help />} />
          <Route path="*" element={<Navigate to="/" replace />} />
        </Route>
      </Routes>
    </HashRouter>
  );
}

/** Dialog asking user if they want to restore a previous session */
function RestoreSessionDialog({ savedData, onRestore, onNewSession }: {
  savedData: SavedSession;
  onRestore: () => void;
  onNewSession: () => void;
}) {
  const { t } = useI18n();
  const npeCount = savedData?.npes?.length || 0;
  const connCount = savedData?.connections?.length || 0;
  const phaseCount = savedData?.contingencies?.length || 0;

  return (
    <motion.div
      initial={{ opacity: 0 }}
      animate={{ opacity: 1 }}
      exit={{ opacity: 0 }}
      className="fixed inset-0 z-[90] flex items-center justify-center bg-black/50 backdrop-blur-sm"
    >
      <motion.div
        initial={{ opacity: 0, scale: 0.95, y: 20 }}
        animate={{ opacity: 1, scale: 1, y: 0 }}
        exit={{ opacity: 0, scale: 0.95, y: 20 }}
        transition={{ duration: 0.25 }}
        className="max-w-md w-full mx-4 rounded-2xl bg-white border border-slate-200 shadow-2xl overflow-hidden"
      >
        <div className="p-6 space-y-4">
          <div className="flex items-center gap-3">
            <div className="w-10 h-10 rounded-xl bg-cyan-50 flex items-center justify-center">
              <RotateCcw size={20} className="text-cyan-600" />
            </div>
            <div>
              <h2 className="text-base font-bold text-slate-800">{t.session.restoreTitle}</h2>
              <p className="text-xs text-slate-500">{t.session.restoreDescription}</p>
            </div>
          </div>

          {/* Summary of saved data */}
          <div className="p-3 rounded-xl bg-slate-50 border border-slate-100">
            <div className="flex gap-4 text-xs">
              <div className="text-center">
                <p className="text-lg font-bold text-cyan-600">{npeCount}</p>
                <p className="text-slate-500">NPEs</p>
              </div>
              <div className="text-center">
                <p className="text-lg font-bold text-teal-600">{connCount}</p>
                <p className="text-slate-500">{t.dashboard.connections}</p>
              </div>
              <div className="text-center">
                <p className="text-lg font-bold text-violet-600">{phaseCount}</p>
                <p className="text-slate-500">{t.dashboard.phases}</p>
              </div>
            </div>
          </div>

          {/* Buttons */}
          <div className="flex gap-3">
            <button
              onClick={onRestore}
              className="flex-1 flex items-center justify-center gap-2 px-4 py-2.5 rounded-xl bg-gradient-to-r from-cyan-500 to-teal-500 text-white font-bold text-sm hover:shadow-lg hover:shadow-cyan-500/20 transition-all"
            >
              <RotateCcw size={14} />
              {t.session.restore}
            </button>
            <button
              onClick={onNewSession}
              className="flex-1 flex items-center justify-center gap-2 px-4 py-2.5 rounded-xl bg-white text-slate-600 font-bold text-sm border border-slate-200 hover:bg-slate-50 transition-all"
            >
              <Plus size={14} />
              {t.session.newSession}
            </button>
          </div>
        </div>
      </motion.div>
    </motion.div>
  );
}

function App() {
  const [splashDone, setSplashDone] = useState(false);
  const [savedData, setSavedData] = useState<ReturnType<typeof getSavedSession>>(
    () => getSavedSession(),
  );
  const [sessionDecided, setSessionDecided] = useState(() => savedData === null);

  const handleRestore = useCallback(() => {
    if (savedData) restoreSession(savedData);
    setSavedData(null);
    setSessionDecided(true);
  }, [savedData]);

  const handleNewSession = useCallback(() => {
    clearAutoSave();
    setSavedData(null);
    setSessionDecided(true);
  }, []);

  if (!splashDone) {
    return <SplashScreen onComplete={() => setSplashDone(true)} />;
  }

  return (
    <ErrorBoundary>
      <ToastProvider>
        <ConfirmDialogProvider>
          <ConnectionMonitor />
          {sessionDecided && <AppContent />}
          <AnimatePresence>
            {savedData && !sessionDecided && (
              <RestoreSessionDialog
                savedData={savedData}
                onRestore={handleRestore}
                onNewSession={handleNewSession}
              />
            )}
          </AnimatePresence>
        </ConfirmDialogProvider>
      </ToastProvider>
    </ErrorBoundary>
  );
}

export default App;
