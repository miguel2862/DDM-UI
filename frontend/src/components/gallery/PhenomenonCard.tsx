import { useState, useRef, useCallback } from 'react';
import { motion, AnimatePresence } from 'framer-motion';
import {
  Flame, ShieldOff, Eclipse, ShieldAlert, RotateCcw, Timer, Lock,
  Scale, Layers, ArrowRightLeft, Loader2, Zap, EyeOff,
} from 'lucide-react';
import { Badge } from '../ui/Badge';
import { useI18n } from '../../i18n';
import type { Phenomenon } from '../../data/phenomena';

const iconMap: Record<string, React.ComponentType<any>> = {
  Flame, ShieldOff, Eclipse, ShieldAlert, RotateCcw, Timer, Scale, Layers, ArrowRightLeft, Zap, EyeOff,
};

const colorStyles: Record<string, { gradient: string; iconBg: string; iconText: string; border: string; btnBg: string }> = {
  cyan: {
    gradient: 'from-cyan-100 to-cyan-50',
    iconBg: 'bg-cyan-100',
    iconText: 'text-cyan-600',
    border: 'hover:border-cyan-400',
    btnBg: 'bg-cyan-50 hover:bg-cyan-100 text-cyan-700 border-cyan-200',
  },
  violet: {
    gradient: 'from-violet-100 to-violet-50',
    iconBg: 'bg-violet-100',
    iconText: 'text-violet-600',
    border: 'hover:border-violet-400',
    btnBg: 'bg-violet-50 hover:bg-violet-100 text-violet-700 border-violet-200',
  },
  amber: {
    gradient: 'from-amber-100 to-amber-50',
    iconBg: 'bg-amber-100',
    iconText: 'text-amber-600',
    border: 'hover:border-amber-400',
    btnBg: 'bg-amber-50 hover:bg-amber-100 text-amber-700 border-amber-200',
  },
  rose: {
    gradient: 'from-rose-100 to-rose-50',
    iconBg: 'bg-rose-100',
    iconText: 'text-rose-600',
    border: 'hover:border-rose-400',
    btnBg: 'bg-rose-50 hover:bg-rose-100 text-rose-700 border-rose-200',
  },
  teal: {
    gradient: 'from-teal-100 to-teal-50',
    iconBg: 'bg-teal-100',
    iconText: 'text-teal-600',
    border: 'hover:border-teal-400',
    btnBg: 'bg-teal-50 hover:bg-teal-100 text-teal-700 border-teal-200',
  },
  emerald: {
    gradient: 'from-emerald-100 to-emerald-50',
    iconBg: 'bg-emerald-100',
    iconText: 'text-emerald-600',
    border: 'hover:border-emerald-400',
    btnBg: 'bg-emerald-50 hover:bg-emerald-100 text-emerald-700 border-emerald-200',
  },
};

// Map phenomenon IDs to translation keys
const phenomenonTranslationKeys: Record<string, string> = {
  extinction: 'extinction',
  alcala_2017: 'alcala_2017',
  burgos_donahoe_blocking: 'burgos_donahoe_blocking',
  burgos_donahoe_successive: 'burgos_donahoe_successive',
  acquisition: 'acquisition',
  latent_inhibition: 'latent_inhibition',
};

interface PhenomenonCardProps {
  phenomenon: Phenomenon;
  onLoad: (id: string) => void;
  delay?: number;
  isLoading?: boolean;
}

export function PhenomenonCard({ phenomenon, onLoad, delay = 0, isLoading = false }: PhenomenonCardProps) {
  const { t } = useI18n();
  const Icon = iconMap[phenomenon.icon] || Flame;
  const style = colorStyles[phenomenon.color] || colorStyles.cyan;
  const [showTooltip, setShowTooltip] = useState(false);
  const hoverTimerRef = useRef<ReturnType<typeof setTimeout> | null>(null);

  // Get translated name/description if available, fallback to data file
  const translationKey = phenomenonTranslationKeys[phenomenon.id] as keyof typeof t.phenomena;
  const translated = translationKey ? t.phenomena[translationKey] : null;
  const displayName = translated ? translated.name : phenomenon.name;
  const displayDescription = translated ? translated.description : phenomenon.description;

  const handleMouseEnter = useCallback(() => {
    hoverTimerRef.current = setTimeout(() => {
      setShowTooltip(true);
    }, 2000);
  }, []);

  const handleMouseLeave = useCallback(() => {
    if (hoverTimerRef.current) {
      clearTimeout(hoverTimerRef.current);
      hoverTimerRef.current = null;
    }
    setShowTooltip(false);
  }, []);

  return (
    <motion.div
      initial={{ opacity: 0, y: 20 }}
      animate={{ opacity: 1, y: 0 }}
      transition={{ duration: 0.4, delay }}
      whileHover={{ y: -4, transition: { duration: 0.2 } }}
      className={`
        group relative rounded-2xl border border-slate-200 bg-white shadow-sm
        overflow-hidden transition-all duration-300 ${style.border}
        ${!phenomenon.available ? 'opacity-60' : 'cursor-pointer hover:shadow-md'}
      `}
      onClick={() => phenomenon.available && onLoad(phenomenon.id)}
      onMouseEnter={handleMouseEnter}
      onMouseLeave={handleMouseLeave}
    >
      {/* Gradient header */}
      <div className={`h-2 bg-gradient-to-r ${style.gradient}`} />

      <div className="p-5">
        {/* Icon + Badge row */}
        <div className="flex items-start justify-between mb-4">
          <div className={`w-12 h-12 rounded-xl ${style.iconBg} flex items-center justify-center`}>
            <Icon size={24} className={style.iconText} />
          </div>
          {phenomenon.available ? (
            <Badge variant="success">{t.gallery.ready}</Badge>
          ) : (
            <Badge variant="muted">{t.gallery.comingSoon}</Badge>
          )}
        </div>

        {/* Title */}
        <h3 className="text-lg font-bold text-slate-800 mb-2">{displayName}</h3>

        {/* Description */}
        <p className="text-sm text-slate-500 leading-relaxed mb-4 line-clamp-3">
          {displayDescription}
        </p>

        {/* Meta */}
        <div className="flex items-center gap-4 text-xs text-slate-400">
          <span>{phenomenon.npeCount} {t.gallery.npesLabel}</span>
          <span>{phenomenon.connectionCount} {t.gallery.connectionsLabel}</span>
        </div>

        {/* Phases */}
        <div className="mt-3 px-3 py-2 rounded-lg bg-slate-50 text-xs text-slate-500 font-mono border border-slate-100">
          {phenomenon.phases}
        </div>

        {/* Load overlay */}
        {phenomenon.available && (
          <div className="mt-4 flex justify-center">
            <motion.button
              whileHover={{ scale: isLoading ? 1 : 1.02 }}
              whileTap={{ scale: isLoading ? 1 : 0.98 }}
              disabled={isLoading}
              className={`w-full py-2.5 rounded-xl font-semibold text-sm transition-all border ${style.btnBg} ${isLoading ? 'opacity-60 cursor-wait' : ''} flex items-center justify-center gap-2`}
            >
              {isLoading && <Loader2 size={14} className="animate-spin" />}
              {isLoading ? t.gallery.loading : t.gallery.loadTemplate}
            </motion.button>
          </div>
        )}

        {!phenomenon.available && (
          <div className="mt-4 flex justify-center items-center gap-2 text-slate-400 text-xs">
            <Lock size={12} />
            <span>{t.gallery.templateNotAvailable}</span>
          </div>
        )}
      </div>

      {/* Full description tooltip on hover (after 2 seconds) */}
      <AnimatePresence>
        {showTooltip && (
          <motion.div
            initial={{ opacity: 0, y: 5 }}
            animate={{ opacity: 1, y: 0 }}
            exit={{ opacity: 0, y: 5 }}
            transition={{ duration: 0.2 }}
            className="absolute inset-x-0 bottom-0 z-10 p-4 bg-slate-800/95 backdrop-blur-sm rounded-b-2xl border-t border-slate-700"
          >
            <p className="text-xs text-slate-200 leading-relaxed">
              {displayDescription}
            </p>
          </motion.div>
        )}
      </AnimatePresence>
    </motion.div>
  );
}
