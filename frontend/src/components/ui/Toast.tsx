import { useState, useEffect, useCallback, useRef, createContext, useContext, type ReactNode } from 'react';
import { motion, AnimatePresence } from 'framer-motion';
import { CheckCircle2, AlertCircle, Info, X } from 'lucide-react';
import { useI18n } from '../../i18n';

// ── Types ───────────────────────────────────────────────────────────────────
type ToastType = 'success' | 'error' | 'info' | 'warning';

interface Toast {
  id: number;
  type: ToastType;
  message: string;
  duration?: number;
}

interface ToastContextValue {
  addToast: (type: ToastType, message: string, duration?: number) => void;
  success: (message: string) => void;
  error: (message: string) => void;
  info: (message: string) => void;
  warning: (message: string) => void;
}

// ── Context ─────────────────────────────────────────────────────────────────
const ToastContext = createContext<ToastContextValue | null>(null);

let _nextId = 0;

export function ToastProvider({ children }: { children: ReactNode }) {
  const [toasts, setToasts] = useState<Toast[]>([]);

  const addToast = useCallback((type: ToastType, message: string, duration = 4000) => {
    const id = ++_nextId;
    setToasts(prev => [...prev, { id, type, message, duration }]);
  }, []);

  const removeToast = useCallback((id: number) => {
    setToasts(prev => prev.filter(t => t.id !== id));
  }, []);

  const success = useCallback((msg: string) => addToast('success', msg), [addToast]);
  const error = useCallback((msg: string) => addToast('error', msg, 6000), [addToast]);
  const info = useCallback((msg: string) => addToast('info', msg), [addToast]);
  const warning = useCallback((msg: string) => addToast('warning', msg, 5000), [addToast]);

  return (
    <ToastContext.Provider value={{ addToast, success, error, info, warning }}>
      {children}
      <div className="fixed bottom-4 right-4 z-[200] space-y-2 max-w-sm">
        <AnimatePresence>
          {toasts.map(toast => (
            <ToastItem key={toast.id} toast={toast} onDismiss={() => removeToast(toast.id)} />
          ))}
        </AnimatePresence>
      </div>
    </ToastContext.Provider>
  );
}

export function useToast(): ToastContextValue {
  const ctx = useContext(ToastContext);
  if (!ctx) throw new Error('useToast must be used within ToastProvider');
  return ctx;
}

// ── Toast Item ──────────────────────────────────────────────────────────────
const icons: Record<ToastType, typeof CheckCircle2> = {
  success: CheckCircle2,
  error: AlertCircle,
  info: Info,
  warning: AlertCircle,
};

const colors: Record<ToastType, string> = {
  success: 'border-emerald-200 bg-white/95',
  error: 'border-rose-200 bg-white/95',
  info: 'border-cyan-200 bg-white/95',
  warning: 'border-amber-200 bg-white/95',
};

const iconColors: Record<ToastType, string> = {
  success: 'text-emerald-500',
  error: 'text-rose-500',
  info: 'text-cyan-500',
  warning: 'text-amber-500',
};

function ToastItem({ toast, onDismiss }: { toast: Toast; onDismiss: () => void }) {
  const Icon = icons[toast.type];

  useEffect(() => {
    if (toast.duration && toast.duration > 0) {
      const timer = setTimeout(onDismiss, toast.duration);
      return () => clearTimeout(timer);
    }
  }, [toast.duration, onDismiss]);

  return (
    <motion.div
      initial={{ opacity: 0, x: 80, scale: 0.95 }}
      animate={{ opacity: 1, x: 0, scale: 1 }}
      exit={{ opacity: 0, x: 80, scale: 0.95 }}
      transition={{ duration: 0.25, ease: 'easeOut' }}
      className={`flex items-start gap-3 p-3 rounded-xl border backdrop-blur-lg shadow-2xl ${colors[toast.type]}`}
    >
      <Icon size={18} className={`flex-shrink-0 mt-0.5 ${iconColors[toast.type]}`} />
      <p className="text-sm text-slate-700 flex-1 leading-snug">{toast.message}</p>
      <button
        onClick={onDismiss}
        className="flex-shrink-0 text-slate-400 hover:text-slate-600 transition-colors"
      >
        <X size={14} />
      </button>
    </motion.div>
  );
}

// ── Connection Status Monitor (optional auto-toast for connection changes) ──
export function ConnectionMonitor() {
  const toast = useToast();
  const { t } = useI18n();
  const cancelledRef = useRef(false);
  const unsubRef = useRef<(() => void) | undefined>(undefined);

  useEffect(() => {
    cancelledRef.current = false;

    (async () => {
      const { onConnectionChange } = await import('../../api/client');
      if (cancelledRef.current) return; // component unmounted during import
      unsubRef.current = onConnectionChange((connected) => {
        if (connected) {
          toast.success(t.toast.connectionRestored);
        } else {
          toast.warning(t.toast.connectionLost);
        }
      });
    })();

    return () => {
      cancelledRef.current = true;
      unsubRef.current?.();
    };
  }, [toast, t]);

  return null;
}
