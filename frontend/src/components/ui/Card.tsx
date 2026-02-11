import { type ReactNode } from 'react';
import { motion } from 'framer-motion';

interface CardProps {
  children: ReactNode;
  className?: string;
  hover?: boolean;
  glow?: 'cyan' | 'teal' | null;
}

export function Card({ children, className = '', hover = false, glow = null }: CardProps) {
  return (
    <motion.div
      initial={{ opacity: 0, y: 12 }}
      animate={{ opacity: 1, y: 0 }}
      transition={{ duration: 0.3 }}
      className={`
        rounded-2xl border border-slate-200 bg-white/80 backdrop-blur-sm shadow-sm
        ${hover ? 'transition-all duration-300 hover:border-cyan-400/40 hover:shadow-md cursor-pointer' : ''}
        ${glow === 'cyan' ? 'glow-cyan' : glow === 'teal' ? 'glow-teal' : ''}
        ${className}
      `}
    >
      {children}
    </motion.div>
  );
}
