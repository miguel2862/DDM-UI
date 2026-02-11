import { useState, useEffect, useRef } from 'react';
import { motion, AnimatePresence } from 'framer-motion';
import { History, X } from 'lucide-react';
import { useI18n } from '../i18n';
import { VERSION_HISTORY } from '../data/versionHistory';
import { getHealth } from '../api/client';

const MIN_SPLASH_MS = 3500;
const MAX_SPLASH_MS = 30000;
const POLL_INTERVAL_MS = 800;
const VISUAL_PROGRESS_INTERVAL = 100; // Update visual progress every 100ms

export function SplashScreen({ onComplete }: { onComplete: () => void }) {
  const [show, setShow] = useState(true);
  const [showHistory, setShowHistory] = useState(false);
  const [apiReady, setApiReady] = useState(false);
  const [progress, setProgress] = useState(0);
  const { t } = useI18n();
  const startTime = useRef(Date.now());
  const apiReadyTime = useRef<number | null>(null);

  // Smooth visual progress: always animate from 0→100 over MIN_SPLASH_MS
  useEffect(() => {
    const interval = setInterval(() => {
      const elapsed = Date.now() - startTime.current;
      const fraction = elapsed / MIN_SPLASH_MS;

      if (apiReadyTime.current) {
        // API is ready — accelerate to 100%
        const sinceReady = Date.now() - apiReadyTime.current;
        const currentProgress = Math.min(90, fraction * 90);
        const rampUp = Math.min(10, sinceReady / 30); // Fill remaining 10% over ~300ms
        setProgress(Math.min(100, currentProgress + rampUp));
      } else {
        // Still waiting for API — fill up to 85% over MIN_SPLASH_MS
        setProgress(Math.min(85, fraction * 85));
      }
    }, VISUAL_PROGRESS_INTERVAL);

    return () => clearInterval(interval);
  }, []);

  // Poll R API health until it responds
  useEffect(() => {
    let cancelled = false;

    const poll = async () => {
      if (cancelled) return;
      try {
        await getHealth();
        if (!cancelled) {
          apiReadyTime.current = Date.now();
          setApiReady(true);
        }
      } catch {
        // R not ready yet — keep polling
        if (!cancelled) {
          setTimeout(poll, POLL_INTERVAL_MS);
        }
      }
    };
    poll();
    return () => { cancelled = true; };
  }, []);

  // Once API is ready AND minimum time passed, close splash
  useEffect(() => {
    if (!apiReady) return;
    if (showHistory) return;

    const elapsed = Date.now() - startTime.current;
    const remaining = Math.max(0, MIN_SPLASH_MS - elapsed);

    const timer = setTimeout(() => {
      setProgress(100);
      // Small delay after reaching 100% so user sees the full bar
      setTimeout(() => {
        setShow(false);
        setTimeout(onComplete, 500);
      }, 300);
    }, remaining);

    return () => clearTimeout(timer);
  }, [apiReady, onComplete, showHistory]);

  // Safety timeout — don't wait forever
  useEffect(() => {
    const safety = setTimeout(() => {
      if (!apiReady) {
        setShow(false);
        setTimeout(onComplete, 500);
      }
    }, MAX_SPLASH_MS);
    return () => clearTimeout(safety);
  }, [apiReady, onComplete]);

  const handleCloseHistory = () => {
    setShowHistory(false);
    if (apiReady) {
      setTimeout(() => {
        setShow(false);
        setTimeout(onComplete, 500);
      }, 300);
    }
  };

  return (
    <AnimatePresence>
      {show && (
        <motion.div
          initial={{ opacity: 0 }}
          animate={{ opacity: 1 }}
          exit={{ opacity: 0 }}
          transition={{ duration: 0.5 }}
          className="fixed inset-0 z-[100] flex items-center justify-center bg-gradient-to-b from-slate-900 via-slate-800 to-slate-900"
        >
          {!showHistory ? (
            <div className="text-center space-y-8 px-6">
              {/* Animated network icon */}
              <motion.div
                initial={{ scale: 0.5, opacity: 0 }}
                animate={{ scale: 1, opacity: 1 }}
                transition={{ delay: 0.2, duration: 0.6, ease: 'easeOut' }}
                className="mx-auto w-24 h-24 rounded-2xl bg-gradient-to-br from-cyan-500 to-teal-500 flex items-center justify-center shadow-2xl shadow-cyan-500/30"
              >
                <svg viewBox="0 0 24 24" fill="none" className="w-12 h-12">
                  <circle cx="4" cy="7" r="2.5" fill="white" opacity="0.9"/>
                  <circle cx="4" cy="17" r="2.5" fill="white" opacity="0.9"/>
                  <circle cx="12" cy="6" r="2" fill="white" opacity="0.8"/>
                  <circle cx="12" cy="18" r="2" fill="white" opacity="0.8"/>
                  <circle cx="20" cy="12" r="2.5" fill="white" opacity="0.9"/>
                  <line x1="6.5" y1="7" x2="10" y2="6" stroke="white" strokeWidth="0.8" opacity="0.5"/>
                  <line x1="6.5" y1="17" x2="10" y2="18" stroke="white" strokeWidth="0.8" opacity="0.5"/>
                  <line x1="6.5" y1="7" x2="10" y2="18" stroke="white" strokeWidth="0.5" opacity="0.3"/>
                  <line x1="6.5" y1="17" x2="10" y2="6" stroke="white" strokeWidth="0.5" opacity="0.3"/>
                  <line x1="14" y1="6" x2="17.5" y2="12" stroke="white" strokeWidth="0.8" opacity="0.5"/>
                  <line x1="14" y1="18" x2="17.5" y2="12" stroke="white" strokeWidth="0.8" opacity="0.5"/>
                </svg>
              </motion.div>

              {/* Title */}
              <motion.div
                initial={{ y: 20, opacity: 0 }}
                animate={{ y: 0, opacity: 1 }}
                transition={{ delay: 0.5, duration: 0.5 }}
              >
                <h1 className="text-5xl font-extrabold text-white tracking-tight">
                  DDM<span className="text-cyan-400">-UI</span>
                </h1>
                <p className="text-xl text-cyan-300/80 font-medium mt-2">
                  {t.app.subtitle}
                </p>
                <motion.p
                  initial={{ opacity: 0 }}
                  animate={{ opacity: 1 }}
                  transition={{ delay: 1.0, duration: 0.5 }}
                  className="text-xs text-slate-500 mt-1 italic"
                >
                  {t.app.subtitleFormer}
                </motion.p>
              </motion.div>

              {/* Tagline */}
              <motion.div
                initial={{ y: 20, opacity: 0 }}
                animate={{ y: 0, opacity: 1 }}
                transition={{ delay: 0.8, duration: 0.5 }}
                className="space-y-2"
              >
                <p className="text-base text-slate-400">
                  {t.app.tagline}
                </p>
                <p className="text-sm text-slate-500">
                  {t.app.taglineLong}
                </p>
              </motion.div>

              {/* Credits */}
              <motion.div
                initial={{ opacity: 0 }}
                animate={{ opacity: 1 }}
                transition={{ delay: 1.2, duration: 0.5 }}
                className="space-y-1.5 pt-2"
              >
                <p className="text-sm text-slate-300">
                  {t.app.credits}
                </p>
                <p className="text-sm text-slate-400">
                  {t.app.university}
                </p>
                <p className="text-xs text-slate-400 mt-2">
                  {t.app.lastModified}
                </p>
              </motion.div>

              {/* Version history button */}
              <motion.div
                initial={{ opacity: 0 }}
                animate={{ opacity: 1 }}
                transition={{ delay: 1.6, duration: 0.4 }}
              >
                <button
                  onClick={(e) => { e.stopPropagation(); setShowHistory(true); }}
                  className="inline-flex items-center gap-2 px-5 py-2 rounded-lg text-xs font-semibold text-cyan-400/80 border border-cyan-500/20 hover:bg-cyan-500/10 hover:text-cyan-300 transition-colors"
                >
                  <History size={14} />
                  {t.app.version} — {t.splash.versionHistory}
                </button>
              </motion.div>

              {/* Loading bar — smoothly fills from 0 to reflect actual R startup */}
              <motion.div
                initial={{ opacity: 0 }}
                animate={{ opacity: 1 }}
                transition={{ delay: 0.3, duration: 0.3 }}
                className="pt-4"
              >
                <div className="w-40 h-1.5 mx-auto bg-slate-700 rounded-full overflow-hidden">
                  <motion.div
                    className="h-full bg-gradient-to-r from-cyan-500 to-teal-500 rounded-full"
                    initial={{ width: '0%' }}
                    animate={{ width: `${progress}%` }}
                    transition={{ duration: 0.6, ease: 'easeOut' }}
                  />
                </div>
                {!apiReady && (
                  <p className="text-xs text-slate-600 mt-3">
                    {t.splash.loading}
                  </p>
                )}
              </motion.div>
            </div>
          ) : (
            /* Version History Panel */
            <motion.div
              initial={{ opacity: 0, y: 20 }}
              animate={{ opacity: 1, y: 0 }}
              transition={{ duration: 0.3 }}
              className="max-w-lg w-full mx-4 max-h-[80vh] overflow-y-auto rounded-2xl bg-slate-800/95 border border-slate-700 shadow-2xl"
            >
              <div className="sticky top-0 flex items-center justify-between p-4 border-b border-slate-700 bg-slate-800/95 backdrop-blur-sm rounded-t-2xl">
                <h2 className="text-base font-bold text-white flex items-center gap-2">
                  <History size={18} className="text-cyan-400" />
                  {t.splash.versionHistory}
                </h2>
                <button
                  onClick={handleCloseHistory}
                  className="w-8 h-8 rounded-lg flex items-center justify-center text-slate-400 hover:text-white hover:bg-slate-700 transition-colors"
                >
                  <X size={16} />
                </button>
              </div>
              <div className="p-4 space-y-5">
                {VERSION_HISTORY.map((entry, i) => (
                  <div key={i} className="space-y-2">
                    <div className="flex items-center gap-2">
                      <span className="px-2.5 py-0.5 rounded-md bg-cyan-500/10 text-cyan-400 text-xs font-bold border border-cyan-500/20">
                        v{entry.version}
                      </span>
                      <span className="text-xs text-slate-500 font-medium">{entry.date}</span>
                    </div>
                    <ul className="space-y-1 pl-3">
                      {entry.changes.map((change, j) => (
                        <li key={j} className="text-xs text-slate-400 leading-relaxed flex items-start gap-1.5">
                          <span className="w-1 h-1 rounded-full bg-slate-600 mt-1.5 flex-shrink-0" />
                          {change}
                        </li>
                      ))}
                    </ul>
                    {i < VERSION_HISTORY.length - 1 && (
                      <div className="border-b border-slate-700/50 pt-1" />
                    )}
                  </div>
                ))}
              </div>
            </motion.div>
          )}
        </motion.div>
      )}
    </AnimatePresence>
  );
}
