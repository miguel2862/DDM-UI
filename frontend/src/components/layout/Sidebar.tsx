import { NavLink } from 'react-router-dom';
import { motion, AnimatePresence } from 'framer-motion';
import {
  Home,
  Network,
  ListChecks,
  Play,
  BarChart3,
  ChevronLeft,
  ChevronRight,
  HelpCircle,
  Sliders,
  Globe,
  History,
  X,
  GraduationCap,
  Lock,
} from 'lucide-react';
import { useState } from 'react';
import { useSimStore } from '../../stores/useSimStore';
import { useToast } from '../ui/Toast';
import { useI18n } from '../../i18n';
import { VERSION_HISTORY } from '../../data/versionHistory';

export function Sidebar() {
  const [collapsed, setCollapsed] = useState(false);
  const [showHistory, setShowHistory] = useState(false);
  const { appMode, setAppMode, simStatus } = useSimStore();
  const { t, language, setLanguage } = useI18n();
  const toast = useToast();
  const isSimRunning = simStatus === 'running';

  const navItems = [
    { to: '/', icon: Home, label: t.nav.dashboard, expertOnly: false },
    { to: '/network', icon: Network, label: t.nav.network, expertOnly: false },
    { to: '/trials', icon: ListChecks, label: t.nav.trials, expertOnly: false },
    { to: '/simulation', icon: Play, label: t.nav.simulate, expertOnly: false },
    { to: '/results', icon: BarChart3, label: t.nav.results, expertOnly: false },
    { to: '/sweep', icon: Sliders, label: t.nav.paramSweep, expertOnly: true },
    { to: '/help', icon: HelpCircle, label: t.nav.help, expertOnly: false },
  ];

  return (
    <>
      <motion.aside
        initial={false}
        animate={{ width: collapsed ? 72 : 240 }}
        transition={{ duration: 0.2, ease: 'easeInOut' }}
        className="fixed left-0 top-0 bottom-0 z-50 flex flex-col bg-gradient-to-b from-slate-900 via-slate-900 to-slate-950 shadow-xl"
      >
        {/* Logo */}
        <div className="flex items-center gap-3 px-4 py-5 border-b border-white/[0.06]">
          <div className="w-10 h-10 rounded-xl bg-gradient-to-br from-cyan-500 to-teal-500 flex items-center justify-center flex-shrink-0 shadow-lg shadow-cyan-500/20">
            <svg viewBox="0 0 24 24" fill="none" className="w-5 h-5">
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
          </div>
          {!collapsed && (
            <motion.div
              initial={{ opacity: 0 }}
              animate={{ opacity: 1 }}
              exit={{ opacity: 0 }}
            >
              <h1 className="text-sm font-bold text-white leading-tight">{t.app.title}</h1>
              <p className="text-[10px] text-slate-400 font-medium">{t.app.subtitle}</p>
            </motion.div>
          )}
        </div>

        {/* Navigation */}
        <nav className="flex-1 py-4 px-2 space-y-1">
          {navItems.filter(item => !item.expertOnly || appMode === 'advanced').map(({ to, icon: Icon, label }) => (
            <NavLink
              key={to}
              to={to}
              onClick={(e) => {
                if (isSimRunning && to !== '/simulation') {
                  e.preventDefault();
                  toast.warning(t.toast.navigationBlocked);
                }
              }}
              className={({ isActive }) =>
                `flex items-center gap-3 px-3 py-2.5 rounded-xl text-sm font-medium transition-all duration-200 group ${
                  isActive
                    ? 'bg-cyan-500/15 text-cyan-400 shadow-[inset_0_0_0_1px_rgba(34,211,238,0.15)]'
                    : isSimRunning && to !== '/simulation'
                      ? 'text-slate-600 cursor-not-allowed opacity-40'
                      : 'text-slate-400 hover:text-slate-200 hover:bg-white/[0.06]'
                }`
              }
            >
              {isSimRunning && to !== '/simulation' ? (
                <Lock size={16} className="flex-shrink-0 text-slate-600" />
              ) : (
                <Icon size={20} className="flex-shrink-0" />
              )}
              {!collapsed && (
                <motion.span
                  initial={{ opacity: 0 }}
                  animate={{ opacity: 1 }}
                  exit={{ opacity: 0 }}
                >
                  {label}
                </motion.span>
              )}
            </NavLink>
          ))}
        </nav>

        {/* Language Toggle - compact pill style */}
        {!collapsed && (
          <div className="px-3 pb-2">
            <label className="flex items-center gap-1.5 text-[10px] font-bold text-slate-500 uppercase tracking-wider mb-1.5">
              <Globe size={10} />
              {t.common.language}
            </label>
            <div className="flex gap-1 p-0.5 rounded-lg bg-white/[0.06] border border-white/[0.06]">
              <button
                onClick={() => setLanguage('en')}
                className={`flex-1 py-1.5 rounded-md text-xs font-bold transition-all ${
                  language === 'en' ? 'bg-white/10 text-cyan-400 shadow-sm' : 'text-slate-500 hover:text-slate-300'
                }`}
              >
                EN
              </button>
              <button
                onClick={() => setLanguage('es')}
                className={`flex-1 py-1.5 rounded-md text-xs font-bold transition-all ${
                  language === 'es' ? 'bg-white/10 text-cyan-400 shadow-sm' : 'text-slate-500 hover:text-slate-300'
                }`}
              >
                ES
              </button>
            </div>
          </div>
        )}

        {/* Mode Toggle - switch with animated indicator */}
        {!collapsed && (
          <div className="px-3 pb-3">
            <label className="flex items-center gap-1.5 text-[10px] font-bold text-slate-500 uppercase tracking-wider mb-1.5">
              <GraduationCap size={10} />
              {t.nav.mode}
            </label>
            <div className="relative rounded-xl bg-white/[0.04] border border-white/[0.06] p-1">
              {/* Animated slider background */}
              <motion.div
                className={`absolute top-1 bottom-1 w-[calc(50%-4px)] rounded-lg ${
                  appMode === 'beginner'
                    ? 'bg-gradient-to-r from-emerald-500/20 to-emerald-400/10 border border-emerald-500/25'
                    : 'bg-gradient-to-r from-violet-500/20 to-violet-400/10 border border-violet-500/25'
                }`}
                animate={{ x: appMode === 'beginner' ? 0 : '100%' }}
                transition={{ type: 'spring', stiffness: 400, damping: 30 }}
                style={{ marginLeft: appMode === 'advanced' ? 4 : 0 }}
              />
              <div className="relative flex">
                <button
                  onClick={() => setAppMode('beginner')}
                  className={`flex-1 flex items-center justify-center gap-1.5 py-2.5 rounded-lg text-xs font-bold transition-colors z-10 ${
                    appMode === 'beginner' ? 'text-emerald-400' : 'text-slate-500 hover:text-slate-300'
                  }`}
                >
                  {t.nav.beginner}
                </button>
                <button
                  onClick={() => setAppMode('advanced')}
                  className={`flex-1 flex items-center justify-center gap-1.5 py-2.5 rounded-lg text-xs font-bold transition-colors z-10 ${
                    appMode === 'advanced' ? 'text-violet-400' : 'text-slate-500 hover:text-slate-300'
                  }`}
                >
                  {t.nav.expert}
                </button>
              </div>
            </div>
          </div>
        )}

        {/* Version badge + History button */}
        <div className={`px-3 pb-2 ${collapsed ? 'text-center' : ''}`}>
          <button
            onClick={() => setShowHistory(true)}
            className="inline-flex items-center gap-1.5 px-2 py-1 rounded-md bg-cyan-500/10 text-cyan-400 text-[10px] font-bold border border-cyan-500/15 hover:bg-cyan-500/20 hover:border-cyan-500/25 transition-colors cursor-pointer"
            title={t.splash.versionHistory}
          >
            <History size={10} />
            {collapsed ? t.app.version : `${t.app.title} ${t.app.version}`}
          </button>
        </div>

        {/* Collapse toggle */}
        <button
          onClick={() => setCollapsed(!collapsed)}
          className="mx-2 mb-4 p-2 rounded-lg text-slate-500 hover:text-slate-300 hover:bg-white/[0.06] transition-colors"
        >
          {collapsed ? <ChevronRight size={18} /> : <ChevronLeft size={18} />}
        </button>
      </motion.aside>

      {/* Version History Modal */}
      <AnimatePresence>
        {showHistory && (
          <motion.div
            initial={{ opacity: 0 }}
            animate={{ opacity: 1 }}
            exit={{ opacity: 0 }}
            transition={{ duration: 0.2 }}
            className="fixed inset-0 z-[60] flex items-center justify-center bg-black/40 backdrop-blur-sm"
            onClick={() => setShowHistory(false)}
          >
            <motion.div
              initial={{ opacity: 0, scale: 0.95, y: 20 }}
              animate={{ opacity: 1, scale: 1, y: 0 }}
              exit={{ opacity: 0, scale: 0.95, y: 20 }}
              transition={{ duration: 0.25 }}
              className="max-w-lg w-full mx-4 max-h-[80vh] overflow-y-auto rounded-2xl bg-white border border-slate-200 shadow-2xl"
              onClick={(e) => e.stopPropagation()}
            >
              <div className="sticky top-0 flex items-center justify-between p-4 border-b border-slate-200 bg-white/95 backdrop-blur-sm rounded-t-2xl">
                <h2 className="text-sm font-bold text-slate-800 flex items-center gap-2">
                  <History size={16} className="text-cyan-600" />
                  {t.splash.versionHistory}
                </h2>
                <button
                  onClick={() => setShowHistory(false)}
                  className="w-8 h-8 rounded-lg flex items-center justify-center text-slate-400 hover:text-slate-700 hover:bg-slate-100 transition-colors"
                >
                  <X size={16} />
                </button>
              </div>
              <div className="p-4 space-y-5">
                {VERSION_HISTORY.map((entry, i) => (
                  <div key={i} className="space-y-2">
                    <div className="flex items-center gap-2">
                      <span className="px-2 py-0.5 rounded-md bg-cyan-50 text-cyan-600 text-[11px] font-bold border border-cyan-100">
                        v{entry.version}
                      </span>
                      <span className="text-[11px] text-slate-500 font-medium">{entry.date}</span>
                    </div>
                    <ul className="space-y-1 pl-3">
                      {entry.changes.map((change, j) => (
                        <li key={j} className="text-[11px] text-slate-500 leading-relaxed flex items-start gap-1.5">
                          <span className="w-1 h-1 rounded-full bg-slate-300 mt-1.5 flex-shrink-0" />
                          {change}
                        </li>
                      ))}
                    </ul>
                    {i < VERSION_HISTORY.length - 1 && (
                      <div className="border-b border-slate-100 pt-1" />
                    )}
                  </div>
                ))}
              </div>
            </motion.div>
          </motion.div>
        )}
      </AnimatePresence>
    </>
  );
}
