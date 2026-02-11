/**
 * BrainView — Anatomically accurate mid-sagittal brain visualization for Transduction mode.
 *
 * Renders an anatomically accurate mid-sagittal brain silhouette derived from the
 * cerebroViz scientific brain atlas (Bahl et al., Bioinformatics 2017). The brain
 * is oriented with the frontal lobe on the LEFT and occipital lobe on the RIGHT.
 *
 * Visible anatomical structures:
 *   - Cerebral cortex regions (M1C, S1C, IPC, V1C, MFC, OFC, ITC)
 *   - Cingulate gyrus / corpus callosum region (CNG)
 *   - Thalamus (THA), Hypothalamus (HTH), Amygdala (AMY)
 *   - Hippocampal formation (HIP)
 *   - Striatum / basal ganglia (STR)
 *   - Substantia nigra / VTA (SN)
 *   - Cerebellum (CB)
 *   - Brainstem: Pons (PON) and Medulla (MED)
 *
 * NPE positions are mapped to neuroanatomical locations:
 *   - PrimarySensory (S')  -> Posterior occipital cortex (V1C)
 *   - US                   -> Thalamus (THA)
 *   - AssociativeSensory   -> Parieto-temporal association (IPC)
 *   - Hippocampal (H)      -> Medial temporal lobe (HIP)
 *   - Dopaminergic (D)     -> Ventral tegmental area (SN)
 *   - AssociativeMotor     -> Premotor / supplementary motor area
 *   - PrimaryMotor (M')    -> Precentral gyrus (M1C)
 *
 * This component renders pure SVG -- no React Flow involved.
 */
import { useMemo } from 'react';
import { motion } from 'framer-motion';
import { getDisplayName } from '../../utils/displayNames';
import type { NPE, Connection } from '../../types/ddm';
import {
  brainBackground_PATH,
  CB_PATH,
  MED_PATH,
  PON_PATH,
  CNG_PATH,
  THA_PATH,
  HIP_PATH,
  HTH_PATH,
  AMY_PATH,
  SN_PATH,
  M1C_PATH,
  S1C_PATH,
  IPC_PATH,
  V1C_PATH,
  MFC_PATH,
  OFC_PATH,
  ITC_PATH,
  STR_PATH,
  BRAIN_TRANSFORM,
} from './brainPaths';

/* -- Colour palette per layer -- */
const layerColors: Record<string, string> = {
  US: '#ef4444',
  PrimarySensory: '#3b82f6',
  AssociativeSensory: '#8b5cf6',
  Hippocampal: '#f59e0b',
  AssociativeMotor: '#14b8a6',
  PrimaryMotor: '#10b981',
  Dopaminergic: '#ec4899',
};

/* -- Neuroanatomical positions (in a 600x440 viewBox, frontal left, occipital right) --
 *
 * These positions were computed by transforming the centroids of cerebroViz brain
 * regions through the same flip+scale transform applied to the anatomy paths.
 * Slight manual adjustments ensure good visual spacing between NPE nodes.
 */
const layerPositions: Record<string, { x: number; y: number }> = {
  PrimaryMotor:       { x: 190, y: 80 },    // precentral gyrus (motor strip)
  AssociativeMotor:   { x: 140, y: 130 },   // premotor / SMA
  US:                 { x: 295, y: 215 },   // thalamus (deep center)
  AssociativeSensory: { x: 380, y: 125 },   // parieto-temporal association
  PrimarySensory:     { x: 475, y: 170 },   // occipital / V1
  Hippocampal:        { x: 345, y: 265 },   // medial temporal lobe
  Dopaminergic:       { x: 265, y: 290 },   // VTA midbrain
};

/* -- Cortical region paths for subtle boundary rendering --
 * Rendered as very faint region boundaries to suggest gyri/sulci.
 */
const corticalRegions = [
  { id: 'M1C', path: M1C_PATH },
  { id: 'S1C', path: S1C_PATH },
  { id: 'IPC', path: IPC_PATH },
  { id: 'V1C', path: V1C_PATH },
  { id: 'MFC', path: MFC_PATH },
  { id: 'OFC', path: OFC_PATH },
  { id: 'ITC', path: ITC_PATH },
];

/* -- Component -- */
interface BrainViewProps {
  npes: NPE[];
  connections: Connection[];
  stimulusLabel: string;
  responseLabel: string;
  width?: number;
  height?: number;
}

export function BrainView({
  npes,
  connections,
  stimulusLabel,
  responseLabel,
  width = 600,
  height = 440,
}: BrainViewProps) {

  // Compute NPE positions -- group by layer, spread multiples vertically
  const npePositions = useMemo(() => {
    const byLayer: Record<string, NPE[]> = {};
    npes.forEach(n => {
      if (!byLayer[n.layer]) byLayer[n.layer] = [];
      byLayer[n.layer].push(n);
    });

    const positions: Record<string, { x: number; y: number }> = {};

    Object.entries(byLayer).forEach(([layer, group]) => {
      const base = layerPositions[layer] || { x: 300, y: 220 };
      const n = group.length;
      const SPREAD = 28; // vertical spacing between NPEs in the same layer

      group.forEach((npe, i) => {
        const off = n > 1 ? (i - (n - 1) / 2) * SPREAD : 0;
        positions[npe.name] = { x: base.x, y: base.y + off };
      });
    });

    return positions;
  }, [npes]);

  // Identify sensory/motor NPEs for virtual stimulus/response arrows
  const sensoryNPEs = npes.filter(n => n.layer === 'PrimarySensory');
  const motorNPEs = npes.filter(n => n.layer === 'PrimaryMotor');

  return (
    <div style={{ width, height, position: 'relative', overflow: 'hidden' }}>
      <svg
        viewBox={`0 0 ${width} ${height}`}
        width="100%"
        height="100%"
        style={{ position: 'absolute', inset: 0 }}
      >
        <defs>
          {/* Glow filter for nodes */}
          <filter id="brain_node_glow" x="-100%" y="-100%" width="300%" height="300%">
            <feGaussianBlur in="SourceGraphic" stdDeviation="5" />
          </filter>
          {/* Subtle glow for tracts */}
          <filter id="brain_tract_glow" x="-40%" y="-40%" width="180%" height="180%">
            <feGaussianBlur in="SourceGraphic" stdDeviation="3" />
          </filter>
          {/* Brain interior fill gradient -- warm anatomical tint */}
          <radialGradient id="brainFillGrad" cx="45%" cy="45%" r="55%">
            <stop offset="0%" stopColor="#3b2850" stopOpacity="0.6" />
            <stop offset="50%" stopColor="#1e1530" stopOpacity="0.4" />
            <stop offset="100%" stopColor="#0d0a14" stopOpacity="0.2" />
          </radialGradient>
          {/* Lighter fill for subcortical structures */}
          <radialGradient id="subcorticalFillGrad" cx="50%" cy="50%" r="50%">
            <stop offset="0%" stopColor="#4a3068" stopOpacity="0.4" />
            <stop offset="100%" stopColor="#2d1b45" stopOpacity="0.2" />
          </radialGradient>
          {/* Cerebellum fill -- slightly different texture */}
          <radialGradient id="cerebellumFillGrad" cx="50%" cy="40%" r="60%">
            <stop offset="0%" stopColor="#2d1f48" stopOpacity="0.5" />
            <stop offset="100%" stopColor="#150e24" stopOpacity="0.25" />
          </radialGradient>
          {/* Ambient halo around brain */}
          <radialGradient id="brainAmbient" cx="50%" cy="50%" r="50%">
            <stop offset="0%" stopColor="#a78bfa" stopOpacity="0.03" />
            <stop offset="60%" stopColor="#6d28d9" stopOpacity="0.01" />
            <stop offset="100%" stopColor="#110d18" stopOpacity="0" />
          </radialGradient>
        </defs>

        {/* Background ambient */}
        <rect width={width} height={height} fill="url(#brainAmbient)" />

        {/* ============================================================
            ANATOMICAL BRAIN PATHS
            All anatomy paths are wrapped in a <g> with BRAIN_TRANSFORM
            which flips horizontally (frontal=left) and scales to fit.
            ============================================================ */}
        <g transform={BRAIN_TRANSFORM}>

          {/* -- Brain silhouette (outer glow) -- */}
          <path
            d={brainBackground_PATH}
            fill="none"
            stroke="#a78bfa"
            strokeWidth={5}
            strokeOpacity={0.06}
            filter="url(#brain_tract_glow)"
          />

          {/* -- Brain silhouette (filled) -- */}
          <path
            d={brainBackground_PATH}
            fill="url(#brainFillGrad)"
            stroke="#a78bfa"
            strokeWidth={2}
            strokeOpacity={0.3}
          />

          {/* -- Cortical region boundaries (subtle sulci lines) -- */}
          {corticalRegions.map(region => (
            <path
              key={region.id}
              d={region.path}
              fill="none"
              stroke="#a78bfa"
              strokeWidth={0.8}
              strokeOpacity={0.1}
            />
          ))}

          {/* -- Cingulate gyrus / corpus callosum region -- */}
          <path
            d={CNG_PATH}
            fill="url(#subcorticalFillGrad)"
            stroke="#c4b5fd"
            strokeWidth={1}
            strokeOpacity={0.15}
          />

          {/* -- Thalamus -- */}
          <path
            d={THA_PATH}
            fill="url(#subcorticalFillGrad)"
            stroke="#c4b5fd"
            strokeWidth={0.8}
            strokeOpacity={0.15}
          />

          {/* -- Hippocampus -- */}
          <path
            d={HIP_PATH}
            fill="url(#subcorticalFillGrad)"
            stroke="#c4b5fd"
            strokeWidth={0.8}
            strokeOpacity={0.12}
          />

          {/* -- Hypothalamus -- */}
          <path
            d={HTH_PATH}
            fill="url(#subcorticalFillGrad)"
            stroke="#c4b5fd"
            strokeWidth={0.6}
            strokeOpacity={0.1}
          />

          {/* -- Amygdala -- */}
          <path
            d={AMY_PATH}
            fill="url(#subcorticalFillGrad)"
            stroke="#c4b5fd"
            strokeWidth={0.6}
            strokeOpacity={0.1}
          />

          {/* -- Substantia nigra / VTA -- */}
          <path
            d={SN_PATH}
            fill="url(#subcorticalFillGrad)"
            stroke="#c4b5fd"
            strokeWidth={0.6}
            strokeOpacity={0.1}
          />

          {/* -- Striatum (basal ganglia) -- */}
          <path
            d={STR_PATH}
            fill="none"
            stroke="#c4b5fd"
            strokeWidth={0.6}
            strokeOpacity={0.08}
          />

          {/* -- Cerebellum -- */}
          <path
            d={CB_PATH}
            fill="url(#cerebellumFillGrad)"
            stroke="#a78bfa"
            strokeWidth={1.2}
            strokeOpacity={0.2}
          />

          {/* -- Pons -- */}
          <path
            d={PON_PATH}
            fill="url(#subcorticalFillGrad)"
            stroke="#a78bfa"
            strokeWidth={1}
            strokeOpacity={0.18}
          />

          {/* -- Medulla oblongata -- */}
          <path
            d={MED_PATH}
            fill="url(#subcorticalFillGrad)"
            stroke="#a78bfa"
            strokeWidth={1}
            strokeOpacity={0.18}
          />

        </g>
        {/* End of anatomical brain paths */}

        {/* -- Neural tract connections -- */}
        {connections.map((conn, idx) => {
          const src = npePositions[conn.presynapticNPE];
          const tgt = npePositions[conn.postsynapticNPE];
          if (!src || !tgt) return null;

          const srcLayer = npes.find(n => n.name === conn.presynapticNPE)?.layer || '';
          const tgtLayer = npes.find(n => n.name === conn.postsynapticNPE)?.layer || '';
          const srcColor = layerColors[srcLayer] || '#64748b';
          const tgtColor = layerColors[tgtLayer] || '#64748b';

          // Create a curved path that routes through the interior
          const centerX = 295;
          const centerY = 210;
          const midX = (src.x + tgt.x) / 2;
          const midY = (src.y + tgt.y) / 2;
          // Pull midpoint toward brain center
          const pullStrength = 0.3;
          const cpX = midX + (centerX - midX) * pullStrength;
          const cpY = midY + (centerY - midY) * pullStrength;

          const path = `M ${src.x} ${src.y} Q ${cpX} ${cpY} ${tgt.x} ${tgt.y}`;
          const tractWidth = Math.max(1.2, conn.weight * 3);
          const animDur = 2.5 + idx * 0.35;
          const gradId = `tract_grad_${idx}`;

          return (
            <g key={`${conn.presynapticNPE}-${conn.postsynapticNPE}`}>
              <defs>
                <linearGradient id={gradId} x1="0%" y1="0%" x2="100%" y2="0%">
                  <stop offset="0%" stopColor={srcColor} stopOpacity="0.6" />
                  <stop offset="100%" stopColor={tgtColor} stopOpacity="0.6" />
                </linearGradient>
              </defs>

              {/* Tract glow underlay */}
              <path
                d={path}
                fill="none"
                stroke={srcColor}
                strokeWidth={tractWidth + 5}
                strokeOpacity={0.04}
                filter="url(#brain_tract_glow)"
              />

              {/* Main tract */}
              <path
                d={path}
                fill="none"
                stroke={`url(#${gradId})`}
                strokeWidth={tractWidth}
                strokeLinecap="round"
                strokeOpacity={0.5}
              />

              {/* Travelling neural impulse */}
              <circle r={Math.max(2, conn.weight * 3)} fill={srcColor} opacity="0">
                <animateMotion
                  dur={`${animDur}s`}
                  repeatCount="indefinite"
                  path={path}
                  keyPoints="0;1"
                  keyTimes="0;1"
                  calcMode="spline"
                  keySplines="0.4 0 0.2 1"
                />
                <animate
                  attributeName="opacity"
                  values="0;0.85;0.85;0"
                  keyTimes="0;0.08;0.75;1"
                  dur={`${animDur}s`}
                  repeatCount="indefinite"
                />
              </circle>

              {/* Second pulse (offset) */}
              <circle r={Math.max(1.5, conn.weight * 2)} fill={tgtColor} opacity="0">
                <animateMotion
                  dur={`${animDur}s`}
                  repeatCount="indefinite"
                  path={path}
                  keyPoints="0;1"
                  keyTimes="0;1"
                  calcMode="spline"
                  keySplines="0.4 0 0.2 1"
                  begin={`${animDur * 0.45}s`}
                />
                <animate
                  attributeName="opacity"
                  values="0;0.55;0.55;0"
                  keyTimes="0;0.08;0.75;1"
                  dur={`${animDur}s`}
                  repeatCount="indefinite"
                  begin={`${animDur * 0.45}s`}
                />
              </circle>

              {/* Weight label */}
              <text
                x={cpX}
                y={cpY - 8}
                textAnchor="middle"
                dominantBaseline="auto"
                fill="#94a3b8"
                fontSize={8}
                fontWeight={600}
                opacity={0.4}
              >
                {conn.weight.toFixed(2)}
              </text>
            </g>
          );
        })}

        {/* -- NPE nodes -- */}
        {npes.map((npe, idx) => {
          const pos = npePositions[npe.name];
          if (!pos) return null;
          const color = layerColors[npe.layer] || '#64748b';
          const nodeR = npe.layer === 'Dopaminergic' ? 11
            : npe.layer === 'US' ? 10
            : npe.layer === 'Hippocampal' ? 10
            : 9;

          return (
            <g key={npe.name}>
              {/* Deep glow halo */}
              <motion.circle
                cx={pos.x}
                cy={pos.y}
                r={nodeR + 8}
                fill={color}
                filter="url(#brain_node_glow)"
                initial={{ opacity: 0.12 }}
                animate={{ opacity: [0.08, 0.25, 0.08] }}
                transition={{
                  duration: 3,
                  repeat: Infinity,
                  ease: 'easeInOut',
                  delay: idx * 0.25,
                }}
              />

              {/* Outer ring */}
              <motion.circle
                cx={pos.x}
                cy={pos.y}
                r={nodeR + 2}
                fill="none"
                stroke={color}
                strokeWidth={1}
                initial={{ strokeOpacity: 0.2 }}
                animate={{ strokeOpacity: [0.15, 0.4, 0.15] }}
                transition={{
                  duration: 2.5,
                  repeat: Infinity,
                  ease: 'easeInOut',
                  delay: idx * 0.25,
                }}
              />

              {/* Core node */}
              <circle
                cx={pos.x}
                cy={pos.y}
                r={nodeR}
                fill={color}
                fillOpacity={0.85}
                stroke="rgba(255,255,255,0.15)"
                strokeWidth={1}
              />

              {/* Specular highlight */}
              <circle
                cx={pos.x - nodeR * 0.25}
                cy={pos.y - nodeR * 0.25}
                r={nodeR * 0.3}
                fill="white"
                opacity={0.15}
              />

              {/* Label */}
              <text
                x={pos.x}
                y={pos.y + nodeR + 14}
                textAnchor="middle"
                fill="#e2e8f0"
                fontSize={11}
                fontWeight={700}
                fontStyle="italic"
              >
                {getDisplayName(npe.name)}
              </text>
            </g>
          );
        })}

        {/* -- Virtual stimulus inputs (entering from outside the skull) -- */}
        {sensoryNPEs.map((npe, i) => {
          const pos = npePositions[npe.name];
          if (!pos) return null;
          const color = layerColors['PrimarySensory'] || '#3b82f6';

          // Stimulus enters from outside the brain (to the right of occipital cortex)
          const entryX = pos.x + 75;
          const entryY = pos.y - 10;

          return (
            <g key={`stim_${i}`}>
              {/* Dashed incoming line */}
              <line
                x1={entryX} y1={entryY}
                x2={pos.x + 12} y2={pos.y}
                stroke={color}
                strokeWidth={1.2}
                strokeOpacity={0.3}
                strokeDasharray="5 3"
              />

              {/* Arrow head */}
              <polygon
                points={`${pos.x + 12},${pos.y - 4} ${pos.x + 12},${pos.y + 4} ${pos.x + 4},${pos.y}`}
                fill={color}
                opacity={0.5}
              />

              {/* Incoming pulse animation */}
              <motion.circle
                r={2.5}
                fill={color}
                initial={{ cx: entryX, cy: entryY, opacity: 0.8 }}
                animate={{
                  cx: [entryX, pos.x + 12],
                  cy: [entryY, pos.y],
                  opacity: [0.8, 0],
                }}
                transition={{
                  duration: 1.6,
                  repeat: Infinity,
                  ease: 'easeIn',
                  delay: i * 0.4,
                }}
              />

              {/* Stimulus label */}
              <text
                x={entryX + 5}
                y={entryY - 6}
                textAnchor="start"
                fill={color}
                fontSize={10}
                fontWeight={600}
                fontStyle="italic"
                opacity={0.6}
              >
                {stimulusLabel} {sensoryNPEs.length > 1 ? i + 1 : ''}
              </text>
            </g>
          );
        })}

        {/* -- Virtual response outputs (exiting from motor cortex outward) -- */}
        {motorNPEs.map((npe, i) => {
          const pos = npePositions[npe.name];
          if (!pos) return null;
          const color = layerColors['PrimaryMotor'] || '#10b981';

          // Response exits from motor cortex upward and to the left
          const exitX = pos.x - 60;
          const exitY = pos.y - 50;

          return (
            <g key={`resp_${i}`}>
              {/* Dashed outgoing line */}
              <line
                x1={pos.x - 10} y1={pos.y - 5}
                x2={exitX} y2={exitY}
                stroke={color}
                strokeWidth={1.2}
                strokeOpacity={0.3}
                strokeDasharray="5 3"
              />

              {/* Arrow head at exit */}
              <polygon
                points={`${exitX - 4},${exitY - 2} ${exitX + 2},${exitY - 5} ${exitX},${exitY + 4}`}
                fill={color}
                opacity={0.5}
              />

              {/* Outgoing pulse */}
              <motion.circle
                r={2.5}
                fill={color}
                initial={{ cx: pos.x - 10, cy: pos.y - 5, opacity: 0.8 }}
                animate={{
                  cx: [pos.x - 10, exitX],
                  cy: [pos.y - 5, exitY],
                  opacity: [0.8, 0],
                }}
                transition={{
                  duration: 1.8,
                  repeat: Infinity,
                  ease: 'easeOut',
                  delay: 1.2 + i * 0.4,
                }}
              />

              {/* Response label */}
              <text
                x={exitX - 5}
                y={exitY - 8}
                textAnchor="end"
                fill={color}
                fontSize={10}
                fontWeight={600}
                fontStyle="italic"
                opacity={0.6}
              >
                {responseLabel} {motorNPEs.length > 1 ? i + 1 : ''}
              </text>
            </g>
          );
        })}

        {/* -- Anatomical region labels (subtle) -- */}
        <text x={120} y={55} textAnchor="middle" fill="#a78bfa" fontSize={8} fontWeight={500} opacity={0.2}>
          Frontal Cortex
        </text>
        <text x={490} y={120} textAnchor="middle" fill="#a78bfa" fontSize={8} fontWeight={500} opacity={0.2}>
          Occipital Cortex
        </text>
        <text x={295} y={195} textAnchor="middle" fill="#a78bfa" fontSize={7} fontWeight={500} opacity={0.15}>
          Thalamus
        </text>
        <text x={360} y={290} textAnchor="middle" fill="#a78bfa" fontSize={7} fontWeight={500} opacity={0.15}>
          Hippocampus
        </text>
        <text x={255} y={315} textAnchor="middle" fill="#a78bfa" fontSize={7} fontWeight={500} opacity={0.15}>
          VTA
        </text>
        <text x={465} y={330} textAnchor="middle" fill="#a78bfa" fontSize={7} fontWeight={500} opacity={0.15}>
          Cerebellum
        </text>
        <text x={365} y={385} textAnchor="middle" fill="#a78bfa" fontSize={7} fontWeight={500} opacity={0.15}>
          Brainstem
        </text>
      </svg>
    </div>
  );
}
