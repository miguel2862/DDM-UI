import { useId, useState, useRef, type ReactNode } from 'react';
import { motion, AnimatePresence } from 'framer-motion';
import { HelpCircle } from 'lucide-react';

interface TooltipProps {
  content: string;
  children?: ReactNode;
  iconSize?: number;
}

export function Tooltip({ content, children, iconSize = 14 }: TooltipProps) {
  const [show, setShow] = useState(false);
  const tooltipRef = useRef<HTMLDivElement>(null);
  const triggerRef = useRef<HTMLSpanElement>(null);
  const [align, setAlign] = useState<'center' | 'left' | 'right'>('center');
  const tooltipId = useId();

  const showTooltip = () => {
    if (triggerRef.current) {
      const rect = triggerRef.current.getBoundingClientRect();
      const tooltipWidth = 256; // w-64
      const halfWidth = tooltipWidth / 2;

      if (rect.left < halfWidth + 8) {
        setAlign('left');
      } else if (window.innerWidth - rect.right < halfWidth + 8) {
        setAlign('right');
      } else {
        setAlign('center');
      }
    }
    setShow(true);
  };

  const alignClass =
    align === 'left'
      ? 'left-0'
      : align === 'right'
        ? 'right-0'
        : 'left-1/2 -translate-x-1/2';

  const arrowClass =
    align === 'left'
      ? 'left-4'
      : align === 'right'
        ? 'right-4'
        : 'left-1/2 -translate-x-1/2';

  return (
    <span className="relative inline-flex items-center" ref={triggerRef}>
      <span
        onMouseEnter={showTooltip}
        onMouseLeave={() => setShow(false)}
        onFocus={showTooltip}
        onBlur={() => setShow(false)}
        tabIndex={0}
        aria-describedby={show ? tooltipId : undefined}
        className="cursor-help inline-flex items-center"
      >
        {children || <HelpCircle size={iconSize} className="text-slate-400 hover:text-cyan-500 transition-colors" />}
      </span>
      <AnimatePresence>
        {show && (
          <motion.div
            id={tooltipId}
            role="tooltip"
            ref={tooltipRef}
            initial={{ opacity: 0, y: 4 }}
            animate={{ opacity: 1, y: 0 }}
            exit={{ opacity: 0, y: 4 }}
            transition={{ duration: 0.15 }}
            className={`absolute z-50 bottom-full ${alignClass} mb-2 px-3 py-2 rounded-xl bg-white border border-slate-200 shadow-lg text-xs text-slate-600 leading-relaxed w-64 pointer-events-none`}
          >
            {content}
            <div className={`absolute top-full ${arrowClass} -mt-1 w-2 h-2 bg-white border-r border-b border-slate-200 rotate-45`} />
          </motion.div>
        )}
      </AnimatePresence>
    </span>
  );
}
