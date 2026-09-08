import { useState, useCallback, useEffect, createContext, useContext, type ReactNode } from 'react';
import { motion, AnimatePresence } from 'framer-motion';
import { AlertTriangle, X } from 'lucide-react';

// ── Types ───────────────────────────────────────────────────────────────────
interface ConfirmOptions {
  title: string;
  message: string;
  confirmLabel?: string;
  cancelLabel?: string;
  variant?: 'danger' | 'warning' | 'info';
}

interface ConfirmContextValue {
  confirm: (options: ConfirmOptions) => Promise<boolean>;
}

// ── Context ─────────────────────────────────────────────────────────────────
const ConfirmContext = createContext<ConfirmContextValue | null>(null);

export function ConfirmDialogProvider({ children }: { children: ReactNode }) {
  const [state, setState] = useState<{
    options: ConfirmOptions;
    resolve: (value: boolean) => void;
  } | null>(null);

  const confirm = useCallback((options: ConfirmOptions): Promise<boolean> => {
    return new Promise(resolve => {
      setState({ options, resolve });
    });
  }, []);

  const handleConfirm = useCallback(() => {
    state?.resolve(true);
    setState(null);
  }, [state]);

  const handleCancel = useCallback(() => {
    state?.resolve(false);
    setState(null);
  }, [state]);

  // Close on Escape key
  useEffect(() => {
    if (!state) return;
    const handler = (e: KeyboardEvent) => {
      if (e.key === 'Escape') handleCancel();
    };
    window.addEventListener('keydown', handler);
    return () => window.removeEventListener('keydown', handler);
  }, [state, handleCancel]);

  return (
    <ConfirmContext.Provider value={{ confirm }}>
      {children}
      <AnimatePresence>
        {state && (
          <motion.div
            initial={{ opacity: 0 }}
            animate={{ opacity: 1 }}
            exit={{ opacity: 0 }}
            transition={{ duration: 0.15 }}
            className="fixed inset-0 z-[150] flex items-center justify-center bg-black/50 backdrop-blur-sm"
            onClick={handleCancel}
          >
            <motion.div
              initial={{ opacity: 0, scale: 0.95, y: 10 }}
              animate={{ opacity: 1, scale: 1, y: 0 }}
              exit={{ opacity: 0, scale: 0.95, y: 10 }}
              transition={{ duration: 0.2 }}
              className="w-full max-w-sm mx-4 rounded-2xl bg-white border border-slate-200 shadow-2xl"
              onClick={e => e.stopPropagation()}
            >
              <div className="flex items-start gap-3 p-5">
                <div className={`w-10 h-10 rounded-xl flex items-center justify-center flex-shrink-0 ${
                  state.options.variant === 'danger' ? 'bg-rose-50' :
                  state.options.variant === 'warning' ? 'bg-amber-50' : 'bg-cyan-50'
                }`}>
                  <AlertTriangle size={20} className={
                    state.options.variant === 'danger' ? 'text-rose-500' :
                    state.options.variant === 'warning' ? 'text-amber-500' : 'text-cyan-500'
                  } />
                </div>
                <div className="flex-1 min-w-0">
                  <h3 className="text-sm font-bold text-slate-800 mb-1">
                    {state.options.title}
                  </h3>
                  <p className="text-xs text-slate-500 leading-relaxed">
                    {state.options.message}
                  </p>
                </div>
                <button
                  onClick={handleCancel}
                  className="flex-shrink-0 text-slate-400 hover:text-slate-600 transition-colors"
                >
                  <X size={16} />
                </button>
              </div>

              <div className="flex items-center gap-2 px-5 pb-5">
                <button
                  onClick={handleCancel}
                  className="flex-1 px-4 py-2.5 rounded-xl bg-slate-100 text-slate-600 font-semibold text-sm hover:bg-slate-200 transition-colors"
                >
                  {state.options.cancelLabel || 'Cancel'}
                </button>
                <button
                  onClick={handleConfirm}
                  className={`flex-1 px-4 py-2.5 rounded-xl font-semibold text-sm transition-colors ${
                    state.options.variant === 'danger'
                      ? 'bg-rose-500 text-white hover:bg-rose-400'
                      : state.options.variant === 'warning'
                        ? 'bg-amber-500 text-white hover:bg-amber-400'
                        : 'bg-cyan-500 text-white hover:bg-cyan-400'
                  }`}
                >
                  {state.options.confirmLabel || 'Confirm'}
                </button>
              </div>
            </motion.div>
          </motion.div>
        )}
      </AnimatePresence>
    </ConfirmContext.Provider>
  );
}

// The provider and its hook intentionally form one colocated context API.
// eslint-disable-next-line react-refresh/only-export-components
export function useConfirm(): ConfirmContextValue {
  const ctx = useContext(ConfirmContext);
  if (!ctx) throw new Error('useConfirm must be used within ConfirmDialogProvider');
  return ctx;
}
