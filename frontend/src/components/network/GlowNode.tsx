/**
 * GlowNode — Custom React Flow node for Modern & Transduction layouts.
 *
 * Features:
 * - Dark-themed circular node with radial gradient fill
 * - SVG glow filter (feGaussianBlur) per-node
 * - Pulsing outer ring animation via Framer Motion
 * - Glassmorphism inner ring with backdrop blur
 * - Animated "breathing" scale + opacity
 * - Layer-specific accent color inherited via data.color
 */
import { memo } from 'react';
import { Handle, Position, type Node, type NodeProps } from '@xyflow/react';
import { motion } from 'framer-motion';

interface GlowNodeData extends Record<string, unknown> {
  label: string;
  color: string;
  size: number;
  hasSource: boolean;
  hasTarget: boolean;
  layerIndex?: number;
}

type GlowNodeModel = Node<GlowNodeData, 'glowNode'>;

function GlowNodeComponent({ data, id }: NodeProps<GlowNodeModel>) {
  const { label, color, size, hasSource, hasTarget, layerIndex } = data;
  const r = size / 2;
  const svgSize = size + 30; // extra room for glow bleed
  const cx = svgSize / 2;
  const cy = svgSize / 2;
  const filterId = `glow_${id.replace(/[^a-zA-Z0-9]/g, '_')}`;
  const gradId = `grad_${id.replace(/[^a-zA-Z0-9]/g, '_')}`;
  const ringGradId = `ring_${id.replace(/[^a-zA-Z0-9]/g, '_')}`;
  const delay = (layerIndex || 0) * 0.35;

  return (
    <div
      style={{
        width: svgSize,
        height: svgSize,
        position: 'relative',
        display: 'flex',
        alignItems: 'center',
        justifyContent: 'center',
      }}
    >
      {/* Handles */}
      {hasTarget && (
        <Handle
          type="target"
          position={Position.Left}
          style={{
            background: 'transparent',
            border: 'none',
            width: 8,
            height: 8,
            left: -2,
          }}
        />
      )}
      {hasSource && (
        <Handle
          type="source"
          position={Position.Right}
          style={{
            background: 'transparent',
            border: 'none',
            width: 8,
            height: 8,
            right: -2,
          }}
        />
      )}

      {/* SVG glow layer */}
      <svg
        width={svgSize}
        height={svgSize}
        viewBox={`0 0 ${svgSize} ${svgSize}`}
        style={{ position: 'absolute', top: 0, left: 0, pointerEvents: 'none' }}
      >
        <defs>
          <filter id={filterId} x="-80%" y="-80%" width="260%" height="260%">
            <feGaussianBlur in="SourceGraphic" stdDeviation={r * 0.35} />
          </filter>
          <radialGradient id={gradId} cx="35%" cy="35%">
            <stop offset="0%" stopColor={color} stopOpacity="0.9" />
            <stop offset="70%" stopColor={color} stopOpacity="0.5" />
            <stop offset="100%" stopColor={color} stopOpacity="0.15" />
          </radialGradient>
          <radialGradient id={ringGradId} cx="50%" cy="50%">
            <stop offset="60%" stopColor={color} stopOpacity="0" />
            <stop offset="85%" stopColor={color} stopOpacity="0.25" />
            <stop offset="100%" stopColor={color} stopOpacity="0" />
          </radialGradient>
        </defs>

        {/* Outer glow */}
        <motion.circle
          cx={cx}
          cy={cy}
          r={r + 6}
          fill={color}
          filter={`url(#${filterId})`}
          initial={{ opacity: 0.15 }}
          animate={{ opacity: [0.1, 0.35, 0.1] }}
          transition={{
            duration: 2.5,
            repeat: Infinity,
            ease: 'easeInOut',
            delay,
          }}
        />

        {/* Pulsing ring */}
        <motion.circle
          cx={cx}
          cy={cy}
          r={r + 4}
          fill="none"
          stroke={color}
          strokeWidth={1.5}
          strokeOpacity={0.4}
          initial={{ r: r + 2 }}
          animate={{
            r: [r + 2, r + 8, r + 2],
            strokeOpacity: [0.35, 0.08, 0.35],
          }}
          transition={{
            duration: 3,
            repeat: Infinity,
            ease: 'easeInOut',
            delay: delay + 0.5,
          }}
        />

        {/* Core circle with gradient */}
        <circle
          cx={cx}
          cy={cy}
          r={r}
          fill={`url(#${gradId})`}
          stroke={color}
          strokeWidth={1.5}
          strokeOpacity={0.6}
        />

        {/* Inner highlight (specular) */}
        <circle
          cx={cx - r * 0.2}
          cy={cy - r * 0.2}
          r={r * 0.4}
          fill="white"
          opacity={0.08}
        />
      </svg>

      {/* Label overlay */}
      <div
        style={{
          position: 'absolute',
          inset: 0,
          display: 'flex',
          alignItems: 'center',
          justifyContent: 'center',
          pointerEvents: 'none',
          zIndex: 2,
        }}
      >
        <span
          style={{
            color: '#ffffff',
            fontSize: size > 60 ? 14 : size > 48 ? 12 : 10,
            fontWeight: 700,
            fontStyle: 'italic',
            textShadow: `0 0 8px ${color}, 0 1px 2px rgba(0,0,0,0.5)`,
            letterSpacing: '0.02em',
          }}
        >
          {label}
        </span>
      </div>
    </div>
  );
}

export const GlowNode = memo(GlowNodeComponent);
