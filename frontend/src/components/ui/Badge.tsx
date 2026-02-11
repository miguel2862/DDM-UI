interface BadgeProps {
  children: React.ReactNode;
  variant: 'success' | 'warning' | 'info' | 'muted';
}

const variants = {
  success: 'bg-emerald-50 text-emerald-700 border-emerald-200',
  warning: 'bg-amber-50 text-amber-700 border-amber-200',
  info: 'bg-cyan-50 text-cyan-700 border-cyan-200',
  muted: 'bg-slate-100 text-slate-500 border-slate-200',
};

export function Badge({ children, variant }: BadgeProps) {
  return (
    <span className={`inline-flex items-center px-2.5 py-0.5 rounded-full text-[11px] font-bold border uppercase tracking-wider ${variants[variant]}`}>
      {children}
    </span>
  );
}
