import { useCallback, useState, useEffect } from 'react';
import { useNavigate } from 'react-router-dom';
import { motion } from 'framer-motion';
import { Brain, Network, Zap, GitBranch, Activity, Loader2 } from 'lucide-react';
import { StatCard } from '../components/ui/StatCard';
import { PhenomenonGallery } from '../components/gallery/PhenomenonGallery';
import { PageTransition } from '../components/ui/Skeleton';
import { useSimStore } from '../stores/useSimStore';
import { useConfirm } from '../components/ui/ConfirmDialog';
import { useI18n } from '../i18n';
import { getTemplate } from '../api/client';

export function Dashboard() {
  const navigate = useNavigate();
  const { npes, connections, trials, contingencies, simStatus, simulationResults, loadTemplate } = useSimStore();
  const { t } = useI18n();
  const { confirm } = useConfirm();
  const [loadingTemplate, setLoadingTemplate] = useState<string | null>(null);
  const [templateError, setTemplateError] = useState<string | null>(null);
  const [nameRevealed, setNameRevealed] = useState(false);

  // Trigger name transition animation after hero appears
  useEffect(() => {
    const timer = setTimeout(() => setNameRevealed(true), 2200);
    return () => clearTimeout(timer);
  }, []);

  const handleLoadPhenomenon = useCallback(async (id: string) => {
    // If there's existing data, ask for confirmation
    const hasData = npes.length > 0 || Object.keys(trials).length > 0 || contingencies.length > 0 || simulationResults !== null;
    if (hasData) {
      const ok = await confirm({
        title: t.confirm.loadTemplateTitle,
        message: t.confirm.loadTemplateMessage,
        confirmLabel: t.confirm.confirm,
        cancelLabel: t.confirm.cancel,
        variant: 'warning',
      });
      if (!ok) return;
    }

    setLoadingTemplate(id);
    setTemplateError(null);
    try {
      const data = await getTemplate(id);
      if (data.error) {
        setTemplateError(data.error);
        setLoadingTemplate(null);
        return;
      }
      loadTemplate(data);
      setLoadingTemplate(null);
      navigate('/network');
    } catch (err: unknown) {
      setTemplateError(err instanceof Error ? err.message : 'Failed to connect to the R API.');
      setLoadingTemplate(null);
    }
  }, [loadTemplate, navigate, npes.length, trials, contingencies.length, simulationResults, confirm, t]);

  return (
    <PageTransition>
    <div className="max-w-7xl mx-auto space-y-8">
      {/* Hero */}
      <motion.div
        initial={{ opacity: 0, y: -20 }}
        animate={{ opacity: 1, y: 0 }}
        transition={{ duration: 0.6 }}
        className="relative overflow-hidden rounded-3xl border border-slate-200 p-8 md:p-12"
      >
        {/* Background gradient */}
        <div className="absolute inset-0 bg-gradient-to-b from-cyan-100/60 via-white to-white" />
        <div className="absolute top-0 right-0 w-96 h-96 bg-cyan-500/5 rounded-full blur-3xl" />
        <div className="absolute bottom-0 left-0 w-96 h-96 bg-teal-500/5 rounded-full blur-3xl" />

        {/* Floating DiffDiscM architecture */}
        <div className="absolute right-6 top-1/2 -translate-y-1/2 hidden lg:block">
          <div className="relative w-[300px] h-[230px]">
            <svg className="absolute inset-0 w-full h-full" viewBox="0 0 300 230" overflow="visible">
              <defs>
                {/* Soft blur for organic cloud edges */}
                <filter id="disc_blur" x="-40%" y="-40%" width="180%" height="180%">
                  <feGaussianBlur in="SourceGraphic" stdDeviation="12" />
                </filter>
              </defs>

              {/* ═══ Diffuse Discrepancy Signals — the hallmark of the model ═══ */}

              {/* Signal 1: S'' + H — hippocampal discrepancy cloud */}
              {/* S'' at (110,45), H at (148,118) → center ~(128, 80) */}
              <motion.ellipse cx="128" cy="82" rx="48" ry="54"
                fill="#94a3b8" filter="url(#disc_blur)"
                animate={{ opacity: [0.18, 0.30, 0.18], rx: [46, 50, 46], ry: [52, 56, 52] }}
                transition={{ duration: 5, repeat: Infinity, ease: 'easeInOut' }}
              />

              {/* Signal 2: M'' + D — dopaminergic discrepancy cloud (modulates Signal 1) */}
              {/* M'' at (200,45), D at (200,130) → center ~(200, 88) — slightly denser */}
              <motion.ellipse cx="200" cy="88" rx="44" ry="58"
                fill="#8494a7" filter="url(#disc_blur)"
                animate={{ opacity: [0.22, 0.35, 0.22], rx: [42, 46, 42], ry: [56, 60, 56] }}
                transition={{ duration: 5.5, repeat: Infinity, ease: 'easeInOut', delay: 0.5 }}
              />

              {/* ═══ Connections — animated strokeWidth (learning) ═══ */}
              {/* S' → S'' */}
              <motion.line x1="52" y1="45" x2="98" y2="45" stroke="#06b6d4" strokeOpacity="0.45"
                animate={{ strokeWidth: [0.8, 1.8, 0.8] }}
                transition={{ duration: 3, repeat: Infinity, ease: 'easeInOut', delay: 0 }} />
              {/* S'' → M'' */}
              <motion.line x1="122" y1="45" x2="188" y2="45" stroke="#06b6d4" strokeOpacity="0.45"
                animate={{ strokeWidth: [0.8, 2, 0.8] }}
                transition={{ duration: 3.5, repeat: Infinity, ease: 'easeInOut', delay: 0.5 }} />
              {/* M'' → M' */}
              <motion.line x1="212" y1="45" x2="255" y2="45" stroke="#06b6d4" strokeOpacity="0.45"
                animate={{ strokeWidth: [0.8, 1.8, 0.8] }}
                transition={{ duration: 3, repeat: Infinity, ease: 'easeInOut', delay: 1.0 }} />
              {/* S'' → H (inclined to the right) */}
              <motion.line x1="115" y1="56" x2="142" y2="108" stroke="#06b6d4" strokeOpacity="0.35"
                animate={{ strokeWidth: [0.6, 1.5, 0.6] }}
                transition={{ duration: 4, repeat: Infinity, ease: 'easeInOut', delay: 0.3 }} />
              {/* M'' → D */}
              <motion.line x1="200" y1="56" x2="200" y2="116" stroke="#06b6d4" strokeOpacity="0.35"
                animate={{ strokeWidth: [0.6, 1.5, 0.6] }}
                transition={{ duration: 4, repeat: Infinity, ease: 'easeInOut', delay: 0.8 }} />

              {/* US → D : curved connection from bottom-left sweeping up to D (fixed strong) */}
              <motion.path
                d="M 55 192 C 90 210, 170 195, 195 135"
                fill="none" stroke="#06b6d4" strokeOpacity="0.35"
                strokeLinecap="round"
                animate={{ strokeWidth: [1.5, 2.8, 1.5] }}
                transition={{ duration: 4, repeat: Infinity, ease: 'easeInOut', delay: 1.2 }}
              />

              {/* ═══ Nodes — all cyan, subtle opacity pulse ═══ */}
              {/* S' — square */}
              <motion.rect x="30" y="33" width="22" height="22" rx="3"
                fill="#06b6d4" stroke="#06b6d4" strokeWidth="1"
                animate={{ fillOpacity: [0.35, 0.6, 0.35], strokeOpacity: [0.25, 0.5, 0.25] }}
                transition={{ duration: 2.5, repeat: Infinity, ease: 'easeInOut', delay: 0 }}
              />
              {/* S'' — circle */}
              <motion.circle cx="110" cy="45" r="12"
                fill="#06b6d4" stroke="#06b6d4" strokeWidth="1"
                animate={{ fillOpacity: [0.35, 0.6, 0.35], strokeOpacity: [0.25, 0.5, 0.25] }}
                transition={{ duration: 2.5, repeat: Infinity, ease: 'easeInOut', delay: 0.3 }}
              />
              {/* M'' — circle */}
              <motion.circle cx="200" cy="45" r="12"
                fill="#06b6d4" stroke="#06b6d4" strokeWidth="1"
                animate={{ fillOpacity: [0.35, 0.6, 0.35], strokeOpacity: [0.25, 0.5, 0.25] }}
                transition={{ duration: 2.5, repeat: Infinity, ease: 'easeInOut', delay: 0.6 }}
              />
              {/* M' — square */}
              <motion.rect x="254" y="33" width="22" height="22" rx="3"
                fill="#06b6d4" stroke="#06b6d4" strokeWidth="1"
                animate={{ fillOpacity: [0.35, 0.6, 0.35], strokeOpacity: [0.25, 0.5, 0.25] }}
                transition={{ duration: 2.5, repeat: Infinity, ease: 'easeInOut', delay: 0.9 }}
              />
              {/* H — circle (positioned inclined to the right of S'') */}
              <motion.circle cx="148" cy="118" r="10"
                fill="#06b6d4" stroke="#06b6d4" strokeWidth="1"
                animate={{ fillOpacity: [0.3, 0.55, 0.3], strokeOpacity: [0.2, 0.45, 0.2] }}
                transition={{ duration: 2.5, repeat: Infinity, ease: 'easeInOut', delay: 0.5 }}
              />
              {/* D — circle, larger */}
              <motion.circle cx="200" cy="130" r="14"
                fill="#06b6d4" stroke="#06b6d4" strokeWidth="1"
                animate={{ fillOpacity: [0.3, 0.55, 0.3], strokeOpacity: [0.2, 0.45, 0.2] }}
                transition={{ duration: 2.5, repeat: Infinity, ease: 'easeInOut', delay: 0.8 }}
              />
              {/* US — hexagon */}
              <motion.polygon
                points="45,182 55,176 65,182 65,194 55,200 45,194"
                fill="#06b6d4" stroke="#06b6d4" strokeWidth="1"
                animate={{ fillOpacity: [0.3, 0.55, 0.3], strokeOpacity: [0.2, 0.45, 0.2] }}
                transition={{ duration: 2.5, repeat: Infinity, ease: 'easeInOut', delay: 1.1 }}
              />
            </svg>
          </div>
        </div>

        <div className="relative z-10 max-w-2xl">
          <h1 className="text-4xl md:text-5xl font-extrabold mb-4">
            {t.dashboard.heroTitle1 && (
              <span className="block text-slate-800">{t.dashboard.heroTitle1}</span>
            )}
            <span className="relative inline-block">
              {/* Old name — fades out mysteriously with blur */}
              <motion.span
                className="text-slate-800"
                initial={{ opacity: 1, filter: 'blur(0px)' }}
                animate={nameRevealed ? { opacity: 0, filter: 'blur(8px)' } : {}}
                transition={{ duration: 0.8, ease: 'easeIn' }}
                style={nameRevealed ? { position: 'absolute', left: 0, top: 0, whiteSpace: 'nowrap', pointerEvents: 'none' } : { whiteSpace: 'nowrap' }}
              >
                {t.dashboard.heroTitleOld}
              </motion.span>
              {/* New name — materializes from blur */}
              <motion.span
                className="gradient-text"
                initial={{ opacity: 0, filter: 'blur(8px)' }}
                animate={nameRevealed ? { opacity: 1, filter: 'blur(0px)' } : {}}
                transition={{ duration: 0.8, ease: 'easeOut', delay: 0.5 }}
                style={!nameRevealed ? { position: 'absolute', left: 0, top: 0, whiteSpace: 'nowrap', pointerEvents: 'none' } : { whiteSpace: 'nowrap' }}
              >
                {t.dashboard.heroTitleNew}
              </motion.span>
            </span>
          </h1>
          <p className="text-slate-500 text-lg leading-relaxed mb-6 max-w-xl">
            {t.dashboard.heroDescription}
          </p>
          <div className="flex flex-wrap gap-3">
            <motion.button
              whileHover={{ scale: loadingTemplate ? 1 : 1.02 }}
              whileTap={{ scale: loadingTemplate ? 1 : 0.98 }}
              onClick={() => handleLoadPhenomenon('extinction')}
              disabled={!!loadingTemplate}
              className="px-6 py-3 rounded-xl bg-gradient-to-r from-cyan-500 to-teal-500 text-white font-bold text-sm shadow-lg shadow-cyan-500/20 hover:shadow-cyan-500/30 transition-shadow disabled:opacity-70 disabled:cursor-not-allowed flex items-center gap-2"
            >
              {loadingTemplate && <Loader2 size={16} className="animate-spin" />}
              {loadingTemplate ? t.gallery.loading : t.dashboard.quickStart}
            </motion.button>

            <motion.button
              whileHover={{ scale: 1.02 }}
              whileTap={{ scale: 0.98 }}
              onClick={() => navigate('/network')}
              className="px-6 py-3 rounded-xl border border-slate-300 text-slate-600 font-bold text-sm hover:bg-slate-50 transition-colors"
            >
              {t.dashboard.buildScratch}
            </motion.button>
          </div>
        </div>
      </motion.div>

      {/* Stats */}
      <div className="grid grid-cols-2 md:grid-cols-4 gap-4">
        <StatCard icon={Brain} label={t.dashboard.npes} value={npes.length} color="cyan" delay={0.1} />
        <StatCard icon={GitBranch} label={t.dashboard.connections} value={connections.length} color="teal" delay={0.15} />
        <StatCard icon={Zap} label={t.dashboard.trialTypes} value={Object.keys(trials).length} color="amber" delay={0.2} />
        <StatCard icon={Activity} label={t.dashboard.phases} value={contingencies.length} color="violet" delay={0.25} />
      </div>

      {/* Workflow steps — each step fades in sequentially */}
      <div className="flex items-center gap-2 flex-wrap">
        {[
          { label: t.dashboard.architecture, done: npes.length > 0, icon: Network },
          { label: t.dashboard.trialsStep, done: Object.keys(trials).length > 0, icon: Zap },
          { label: t.dashboard.contingencies, done: contingencies.length > 0, icon: GitBranch },
          { label: t.dashboard.simulateStep, done: simStatus === 'complete', icon: Activity },
          { label: t.dashboard.resultsStep, done: simulationResults !== null && simulationResults.length > 0, icon: Activity },
        ].map((step, i) => (
          <motion.div
            key={i}
            initial={{ opacity: 0, x: -10 }}
            animate={{ opacity: 1, x: 0 }}
            transition={{ delay: 0.4 + i * 0.15, duration: 0.4, ease: 'easeOut' }}
            className="flex items-center gap-2 flex-shrink-0"
          >
            <div className={`flex items-center gap-2 px-4 py-2 rounded-xl border text-xs font-bold
              ${step.done
                ? 'border-emerald-200 bg-emerald-50 text-emerald-600'
                : 'border-slate-200 bg-white text-slate-400'
              }`}
            >
              <step.icon size={14} />
              {step.label}
            </div>
            {i < 4 && (
              <motion.div
                initial={{ opacity: 0, scaleX: 0 }}
                animate={{ opacity: 1, scaleX: 1 }}
                transition={{ delay: 0.55 + i * 0.15, duration: 0.3 }}
                className="w-6 h-px bg-slate-200 flex-shrink-0 origin-left"
              />
            )}
          </motion.div>
        ))}
      </div>

      {/* Template loading error */}
      {templateError && (
        <motion.div
          initial={{ opacity: 0, y: -10 }}
          animate={{ opacity: 1, y: 0 }}
          className="p-4 rounded-xl bg-rose-50 border border-rose-200 text-sm text-rose-600 font-medium flex items-center justify-between"
        >
          <span>⚠ {templateError}</span>
          <button onClick={() => setTemplateError(null)} className="text-rose-400 hover:text-rose-600 text-xs font-bold ml-4">✕</button>
        </motion.div>
      )}

      {/* Phenomenon Gallery */}
      <PhenomenonGallery onLoad={handleLoadPhenomenon} loadingId={loadingTemplate} />
    </div>
    </PageTransition>
  );
}
