/**
 * AnimatedEdge — Custom React Flow edge for Modern & Transduction layouts.
 *
 * Features:
 * - Bezier curved path with gradient stroke (source color → target color)
 * - Animated travelling signal dot along the path
 * - Pulsing stroke-width animation
 * - Semi-transparent stroke with glow bleed
 * - Arrow marker at target end
 */
import { memo, useId } from 'react';
import { getBezierPath, type Edge, type EdgeProps } from '@xyflow/react';

interface AnimatedEdgeData extends Record<string, unknown> {
  sourceColor?: string;
  targetColor?: string;
  weight?: number;
  edgeIndex?: number;
}

type AnimatedEdgeModel = Edge<AnimatedEdgeData, 'animatedEdge'>;

function AnimatedEdgeComponent({
  sourceX,
  sourceY,
  targetX,
  targetY,
  sourcePosition,
  targetPosition,
  data,
}: EdgeProps<AnimatedEdgeModel>) {

  const uid = useId().replace(/:/g, '');
  const gradientId = `edge_grad_${uid}`;
  const glowFilterId = `edge_glow_${uid}`;
  const markerId = `edge_arrow_${uid}`;

  const sourceColor = data?.sourceColor || '#06b6d4';
  const targetColor = data?.targetColor || '#14b8a6';
  const weight = data?.weight ?? 0.5;
  const edgeIndex = data?.edgeIndex ?? 0;

  const strokeWidth = Math.max(1.5, weight * 3);
  const dotSize = Math.max(2, weight * 4);
  const animDuration = 2.2 + edgeIndex * 0.3;

  const [edgePath] = getBezierPath({
    sourceX,
    sourceY,
    targetX,
    targetY,
    sourcePosition,
    targetPosition,
    curvature: 0.35,
  });

  return (
    <>
      <defs>
        {/* Gradient along the edge */}
        <linearGradient id={gradientId} x1="0%" y1="0%" x2="100%" y2="0%">
          <stop offset="0%" stopColor={sourceColor} stopOpacity="0.8" />
          <stop offset="100%" stopColor={targetColor} stopOpacity="0.8" />
        </linearGradient>

        {/* Glow filter for the edge */}
        <filter id={glowFilterId} x="-20%" y="-20%" width="140%" height="140%">
          <feGaussianBlur in="SourceGraphic" stdDeviation="2" />
        </filter>

        {/* Arrow marker */}
        <marker
          id={markerId}
          viewBox="0 0 10 10"
          refX="8"
          refY="5"
          markerWidth="6"
          markerHeight="6"
          orient="auto-start-reverse"
        >
          <path d="M 0 0 L 10 5 L 0 10 z" fill={targetColor} opacity="0.7" />
        </marker>
      </defs>

      {/* Glow underlay */}
      <path
        d={edgePath}
        fill="none"
        stroke={sourceColor}
        strokeWidth={strokeWidth + 4}
        strokeOpacity={0.1}
        filter={`url(#${glowFilterId})`}
      />

      {/* Main edge path */}
      <path
        d={edgePath}
        fill="none"
        stroke={`url(#${gradientId})`}
        strokeWidth={strokeWidth}
        strokeLinecap="round"
        markerEnd={`url(#${markerId})`}
        style={{
          animation: `edgePulse${edgeIndex % 3} ${2 + edgeIndex * 0.2}s ease-in-out infinite`,
        }}
      />

      {/* Travelling signal dot */}
      <circle r={dotSize} fill={sourceColor} opacity="0">
        <animateMotion
          dur={`${animDuration}s`}
          repeatCount="indefinite"
          path={edgePath}
          keyPoints="0;1"
          keyTimes="0;1"
          calcMode="spline"
          keySplines="0.4 0 0.2 1"
        />
        <animate
          attributeName="opacity"
          values="0;0.9;0.9;0"
          keyTimes="0;0.1;0.7;1"
          dur={`${animDuration}s`}
          repeatCount="indefinite"
        />
        <animate
          attributeName="r"
          values={`${dotSize * 0.7};${dotSize};${dotSize * 0.7}`}
          dur={`${animDuration * 0.6}s`}
          repeatCount="indefinite"
        />
      </circle>

      {/* Second trailing dot (offset) for richer signal flow */}
      <circle r={dotSize * 0.6} fill={targetColor} opacity="0">
        <animateMotion
          dur={`${animDuration}s`}
          repeatCount="indefinite"
          path={edgePath}
          keyPoints="0;1"
          keyTimes="0;1"
          calcMode="spline"
          keySplines="0.4 0 0.2 1"
          begin={`${animDuration * 0.4}s`}
        />
        <animate
          attributeName="opacity"
          values="0;0.6;0.6;0"
          keyTimes="0;0.1;0.7;1"
          dur={`${animDuration}s`}
          repeatCount="indefinite"
          begin={`${animDuration * 0.4}s`}
        />
      </circle>
    </>
  );
}

export const AnimatedEdge = memo(AnimatedEdgeComponent);
