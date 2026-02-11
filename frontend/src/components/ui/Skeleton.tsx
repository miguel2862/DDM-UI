import { motion } from 'framer-motion';

// ── Base Skeleton ───────────────────────────────────────────────────────────
export function Skeleton({ className = '', style }: { className?: string; style?: React.CSSProperties }) {
  return (
    <div className={`animate-pulse rounded-lg bg-slate-200/60 ${className}`} style={style} />
  );
}

// ── Card Skeleton ───────────────────────────────────────────────────────────
export function CardSkeleton({ className = '' }: { className?: string }) {
  return (
    <div className={`rounded-2xl border border-slate-200/60 bg-white p-5 ${className}`}>
      <Skeleton className="h-4 w-1/3 mb-4" />
      <Skeleton className="h-3 w-2/3 mb-2" />
      <Skeleton className="h-3 w-1/2" />
    </div>
  );
}

// ── Stat Cards Skeleton ─────────────────────────────────────────────────────
export function StatCardsSkeleton({ count = 4 }: { count?: number }) {
  return (
    <div className="grid grid-cols-1 md:grid-cols-4 gap-4">
      {Array.from({ length: count }).map((_, i) => (
        <div key={i} className="rounded-2xl border border-slate-200/60 bg-white p-4">
          <div className="flex items-center gap-3">
            <Skeleton className="w-10 h-10 rounded-xl" />
            <div className="flex-1">
              <Skeleton className="h-3 w-16 mb-2" />
              <Skeleton className="h-6 w-10" />
            </div>
          </div>
        </div>
      ))}
    </div>
  );
}

// ── Chart Skeleton ──────────────────────────────────────────────────────────
export function ChartSkeleton() {
  return (
    <div className="rounded-2xl border border-slate-200/60 bg-white p-5">
      <Skeleton className="h-4 w-40 mb-6" />
      <div className="flex items-end gap-2 h-48">
        {[40, 65, 35, 80, 55, 70, 45, 60, 75, 50, 85, 40].map((h, i) => (
          <Skeleton key={i} className="flex-1" style={{ height: `${h}%` }} />
        ))}
      </div>
    </div>
  );
}

// ── Page transition wrapper ─────────────────────────────────────────────────
export function PageTransition({ children }: { children: React.ReactNode }) {
  return (
    <motion.div
      initial={{ opacity: 0, y: 12 }}
      animate={{ opacity: 1, y: 0 }}
      transition={{ duration: 0.3, ease: 'easeOut' }}
    >
      {children}
    </motion.div>
  );
}
