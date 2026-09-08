import { useState, useCallback, useMemo, useRef, useEffect } from 'react';
import { motion, AnimatePresence } from 'framer-motion';
import {
  ReactFlow,
  Background,
  Controls,
  useReactFlow,
  ReactFlowProvider,
  applyNodeChanges,
  type Node,
  type Edge,
  type NodeChange,
  type NodeMouseHandler,
  MarkerType,
} from '@xyflow/react';
import '@xyflow/react/dist/style.css';
import { toPng } from 'html-to-image';
import {
  Plus, Trash2, Brain, GitBranch, Settings2, Download, Upload, LayoutGrid, Lock, AlertTriangle,
  Sparkles, ArrowRightLeft,
} from 'lucide-react';
import { Card } from '../components/ui/Card';
import { Badge } from '../components/ui/Badge';
import { useSimStore } from '../stores/useSimStore';
import { useConfirm } from '../components/ui/ConfirmDialog';
import { useToast } from '../components/ui/Toast';
import { getDisplayName } from '../utils/displayNames';
import { downloadArchitectureJSON, parseArchitectureJSON } from '../utils/dataExport';
import { useI18n } from '../i18n';
import type { NPE, NPELayer, NPEType } from '../types/ddm';
import { GlowNode } from '../components/network/GlowNode';
import { AnimatedEdge } from '../components/network/AnimatedEdge';
import { BrainView } from '../components/network/BrainView';
import { PublicationNetwork } from '../components/network/PublicationNetwork';
import { computePublicationLayout, downloadPublicationPng, downloadPublicationSvg } from '../utils/publicationNetwork';

// ── Custom node / edge types for Modern ──
const modernNodeTypes = {
  glowNode: GlowNode,
};

const modernEdgeTypes = {
  animatedEdge: AnimatedEdge,
};

// ── CSS keyframes for edge pulse & dash animations (injected once) ──
const styleId = '__network_anim_styles';
if (typeof document !== 'undefined' && !document.getElementById(styleId)) {
  const style = document.createElement('style');
  style.id = styleId;
  style.textContent = `
    @keyframes edgePulse0 {
      0%, 100% { stroke-opacity: 0.4; }
      50% { stroke-opacity: 0.8; }
    }
    @keyframes edgePulse1 {
      0%, 100% { stroke-opacity: 0.35; }
      50% { stroke-opacity: 0.75; }
    }
    @keyframes edgePulse2 {
      0%, 100% { stroke-opacity: 0.3; }
      50% { stroke-opacity: 0.7; }
    }
    @keyframes dashFlow {
      0% { stroke-dashoffset: 0; }
      100% { stroke-dashoffset: -20; }
    }
    /* Dark canvas overrides for React Flow controls */
    .dark-flow .react-flow__controls {
      background: rgba(15, 23, 42, 0.8) !important;
      border-color: rgba(100, 116, 139, 0.3) !important;
      border-radius: 12px !important;
      backdrop-filter: blur(8px);
    }
    .dark-flow .react-flow__controls button {
      background: transparent !important;
      border-color: rgba(100, 116, 139, 0.2) !important;
      color: #94a3b8 !important;
      fill: #94a3b8 !important;
    }
    .dark-flow .react-flow__controls button:hover {
      background: rgba(100, 116, 139, 0.15) !important;
      color: #e2e8f0 !important;
      fill: #e2e8f0 !important;
    }
    .dark-flow .react-flow__controls svg {
      fill: currentColor !important;
    }
    .dark-flow .react-flow__edge-textbg {
      fill: rgba(15, 23, 42, 0.85) !important;
    }
    .dark-flow .react-flow__edge-text {
      fill: #94a3b8 !important;
    }
  `;
  document.head.appendChild(style);
}

const layerColors: Record<string, string> = {
  US: '#ef4444',
  PrimarySensory: '#3b82f6',
  AssociativeSensory: '#8b5cf6',
  Hippocampal: '#f59e0b',
  AssociativeMotor: '#14b8a6',
  PrimaryMotor: '#10b981',
  Dopaminergic: '#ec4899',
};

// Columns for left-to-right academic layout
const layerColumn: Record<string, number> = {
  US: 0,
  PrimarySensory: 0,
  AssociativeSensory: 1,
  Hippocampal: 1,
  Dopaminergic: 2,
  AssociativeMotor: 2,
  PrimaryMotor: 3,
};

// Row priority within a column (lower = higher on screen)
const layerRow: Record<string, number> = {
  PrimarySensory: 0,
  US: 1,
  AssociativeSensory: 0,
  Hippocampal: 1,
  AssociativeMotor: 0,
  Dopaminergic: 1,
  PrimaryMotor: 0,
};

function getNodeShape(layer: string, type: string): string {
  if (layer === 'US') return '8px';
  if (['PrimarySensory', 'PrimaryMotor'].includes(layer)) return '4px';
  if (type === 'Inhibitory') return '2px';
  if (layer === 'Dopaminergic') return '12px';
  return '50%';
}

/** Compute auto-layout positions based on layer columns/rows */
function computeAutoLayout(npes: NPE[]): Record<string, { x: number; y: number }> {
  const NODE_W = 80;
  const NODE_H = 90;
  const COL_GAP = 200;
  const ROW_GAP = 180;
  const INTRA_GAP = 85;

  // Group NPEs by (column, row)
  const groups: Record<string, NPE[]> = {};
  npes.forEach(npe => {
    const col = layerColumn[npe.layer] ?? 0;
    const row = layerRow[npe.layer] ?? 0;
    const key = `${col}-${row}`;
    if (!groups[key]) groups[key] = [];
    groups[key].push(npe);
  });

  const positions: Record<string, { x: number; y: number }> = {};

  Object.entries(groups).forEach(([key, npesInGroup]) => {
    const [colStr, rowStr] = key.split('-');
    const col = parseInt(colStr);
    const row = parseInt(rowStr);
    const baseX = col * (NODE_W + COL_GAP);
    const baseY = row * (NODE_H + ROW_GAP);

    // Spread multiple NPEs in the same group vertically
    const total = npesInGroup.length;
    npesInGroup.forEach((npe, i) => {
      const yOffset = total > 1 ? (i - (total - 1) / 2) * INTRA_GAP : 0;
      positions[npe.name] = { x: baseX, y: baseY + yOffset };
    });
  });

  return positions;
}

type LayoutMode = 'publication' | 'traditional' | 'modern' | 'transduction';

// ── Node sizes per layer for modern/transduction layouts ──
const modernNodeSize: Record<string, number> = {
  Dopaminergic: 72,
  US: 54,
  Hippocampal: 58,
  AssociativeSensory: 56,
  AssociativeMotor: 56,
  PrimarySensory: 50,
  PrimaryMotor: 50,
};

/**
 * Modern layout — horizontal S-curve "neural pathway".
 *
 * Two parallel streams that weave like a double-helix:
 *   • Top stream (sensory→association→motor): S → S'' → M'' → M'
 *   • Bottom stream (modulatory): US → H → D
 *
 * Each stream follows a sine-wave so the two paths interleave visually,
 * creating a braided neural-pathway aesthetic.
 *
 * Node sizes vary by functional importance (D largest, primary smallest).
 */
function computeModernLayout(npes: NPE[]): Record<string, { x: number; y: number }> {
  const W = 820;
  const H = 420;
  const PAD_X = 60;
  const PAD_Y = 60;
  const INTRA = 70;

  // Horizontal t ∈ [0, 1] for each layer
  const layerT: Record<string, number> = {
    PrimarySensory: 0.00,
    US: 0.05,
    AssociativeSensory: 0.28,
    Hippocampal: 0.40,
    AssociativeMotor: 0.62,
    Dopaminergic: 0.78,
    PrimaryMotor: 1.00,
  };

  // Which stream: top (main pathway) or bottom (modulatory)
  const stream: Record<string, 'top' | 'bottom'> = {
    PrimarySensory: 'top',
    AssociativeSensory: 'top',
    AssociativeMotor: 'top',
    PrimaryMotor: 'top',
    US: 'bottom',
    Hippocampal: 'bottom',
    Dopaminergic: 'bottom',
  };

  // Group by layer
  const byLayer: Record<string, NPE[]> = {};
  npes.forEach(n => {
    if (!byLayer[n.layer]) byLayer[n.layer] = [];
    byLayer[n.layer].push(n);
  });

  const positions: Record<string, { x: number; y: number }> = {};
  const midY = H / 2;

  Object.entries(byLayer).forEach(([layer, group]) => {
    const t = layerT[layer] ?? 0.5;
    const s = stream[layer] ?? 'top';

    const x = PAD_X + t * (W - 2 * PAD_X);

    // Sine wave creates the braiding effect — top stream goes up when bottom goes down
    const wave = Math.sin(t * Math.PI * 1.6) * 55;
    const baseY = s === 'top'
      ? midY - 85 + wave    // upper helix
      : midY + 85 - wave;   // lower helix (inverted)

    // Spread group members vertically
    const total = group.length;
    group.forEach((npe, i) => {
      const off = total > 1 ? (i - (total - 1) / 2) * INTRA : 0;
      positions[npe.name] = {
        x: Math.max(PAD_X, Math.min(W - PAD_X, x)),
        y: Math.max(PAD_Y, Math.min(H - PAD_Y, baseY + off)),
      };
    });
  });

  return positions;
}

// Transduction mode uses BrainView (sagittal brain silhouette SVG) instead of React Flow

const layerOptions: NPELayer[] = [
  'PrimarySensory', 'AssociativeSensory', 'Hippocampal',
  'AssociativeMotor', 'PrimaryMotor', 'Dopaminergic', 'US',
];

const layerDisplayNames: Record<NPELayer, string> = {
  US: 'US (Unconditioned Stimulus)',
  PrimarySensory: "Primary Sensory (S)",
  AssociativeSensory: "Associative Sensory (S'')",
  Hippocampal: 'Hippocampal (H)',
  AssociativeMotor: "Associative Motor (M'')",
  PrimaryMotor: "Primary Motor (M')",
  Dopaminergic: 'Dopaminergic (D)',
};

function NetworkBuilderInner() {
  const {
    npes, connections, addNPE, removeNPE, addConnection, removeConnection,
    appMode, modelKind,
  } = useSimStore();
  const { t, language } = useI18n();
  const { confirm } = useConfirm();
  const toast = useToast();
  const [activeTab, setActiveTab] = useState<'units' | 'connections'>('units');
  const [layoutPreference, setLayoutMode] = useState<LayoutMode>('publication');
  const layoutMode = modelKind !== 'dtd' && layoutPreference === 'publication' ? 'traditional' : layoutPreference;
  const flowRef = useRef<HTMLDivElement>(null);
  const publicationRef = useRef<SVGSVGElement>(null);
  const [isExporting, setIsExporting] = useState(false);
  const { fitView } = useReactFlow();

  // Whether current layout uses the dark canvas + custom components
  const isDarkCanvas = layoutMode === 'modern' || layoutMode === 'transduction';

  // Track custom node positions (set by auto-layout or drag)
  const [customPositions, setCustomPositions] = useState<Record<string, { x: number; y: number }> | null>(null);

  // Lock layout state
  const [isLocked, setIsLocked] = useState(false);

  // Detect disconnected NPEs (no connections at all)
  const disconnectedNPEs = useMemo(() => {
    if (npes.length === 0) return new Set<string>();
    const connected = new Set<string>();
    connections.forEach(c => { connected.add(c.presynapticNPE); connected.add(c.postsynapticNPE); });
    return new Set(npes.filter(n => !connected.has(n.name)).map(n => n.name));
  }, [npes, connections]);

  // Compute positions for the current layout mode (transduction uses OrbitalView, not React Flow)
  const computeLayoutForMode = useCallback((mode: LayoutMode, npeList: NPE[]) => {
    if (mode === 'publication') return Object.fromEntries(Object.entries(computePublicationLayout(npeList, connections))
      .map(([name, point]) => [name, { x: point.x - 30, y: point.y - 30 }]));
    if (mode === 'modern') return computeModernLayout(npeList);
    return computeAutoLayout(npeList);
  }, [connections]);

  // Initialize customPositions on mount / when npes change to prevent first-drag glitch
  useEffect(() => {
    if (npes.length > 0) {
      setCustomPositions(prev => {
        if (prev) {
          // If we already have positions, only add positions for new NPEs
          const autoPositions = computeLayoutForMode(layoutMode, npes);
          const hasAllNpes = npes.every(npe => npe.name in prev);
          if (hasAllNpes) return prev;
          return { ...autoPositions, ...prev };
        }
        return computeLayoutForMode(layoutMode, npes);
      });
    }
  }, [npes, layoutMode, computeLayoutForMode]);

  const handleAutoLayout = useCallback(() => {
    const positions = computeLayoutForMode(layoutMode, npes);
    setCustomPositions(positions);
    // After state update, fit view
    if (layoutMode !== 'publication') setTimeout(() => fitView({ padding: 0.2, duration: 400 }), 50);
  }, [npes, layoutMode, fitView, computeLayoutForMode]);

  const handleLockLayout = useCallback(() => {
    const positions = customPositions || computeLayoutForMode(layoutMode, npes);
    useSimStore.getState().setLockedLayout(positions);
    setIsLocked(true);
    // Reset visual feedback after 2 seconds
    setTimeout(() => setIsLocked(false), 2000);
  }, [customPositions, npes, layoutMode, computeLayoutForMode]);

  const publicationPositions = useMemo(() => {
    const positions = customPositions || computeLayoutForMode('publication', npes);
    return Object.fromEntries(Object.entries(positions).map(([name, point]) => [name, { x: point.x + 30, y: point.y + 30 }]));
  }, [customPositions, computeLayoutForMode, npes]);

  const handlePublicationPositions = useCallback((positions: Record<string, { x: number; y: number }>) => {
    setCustomPositions(Object.fromEntries(Object.entries(positions).map(([name, point]) => [name, { x: point.x - 30, y: point.y - 30 }])));
  }, []);

  const handleDownloadPng = useCallback(async () => {
    if (!flowRef.current || npes.length === 0 || isExporting) return;
    setIsExporting(true);
    try {
      if (layoutMode === 'publication' && publicationRef.current) {
        await downloadPublicationPng(publicationRef.current, 'ddm-network-publication.png', 4);
        return;
      }
      const dataUrl = await toPng(flowRef.current, {
        backgroundColor: layoutMode === 'transduction' ? '#110d18' : layoutMode === 'modern' ? '#0b1120' : '#f8fafc',
        quality: 1.0,
        pixelRatio: 3,
        filter: (node) => !(node instanceof Element && ['react-flow__controls', 'react-flow__minimap', 'react-flow__attribution', 'react-flow__handle'].some(name => node.classList.contains(name))),
      });
      const link = document.createElement('a');
      link.download = 'ddm-network.png';
      link.href = dataUrl;
      link.click();
    } catch (err) {
      console.error('Failed to export PNG:', err);
      toast.error(language === 'es' ? 'No se pudo exportar la red. Inténtalo de nuevo.' : 'Could not export the network. Please try again.');
    } finally {
      setIsExporting(false);
    }
  }, [layoutMode, npes.length, isExporting, toast, language]);

  const handleDownloadSvg = useCallback(async () => {
    if (!publicationRef.current || npes.length === 0) return;
    try {
      await downloadPublicationSvg(publicationRef.current, 'ddm-network-publication.svg');
    } catch {
      toast.error(language === 'es' ? 'No se pudo exportar el SVG.' : 'Could not export the SVG.');
    }
  }, [npes.length, toast, language]);

  const fileInputRef = useRef<HTMLInputElement>(null);

  const handleImportArchitecture = useCallback((event: React.ChangeEvent<HTMLInputElement>) => {
    const file = event.target.files?.[0];
    if (!file) return;

    const reader = new FileReader();
    reader.onload = (e) => {
      const text = e.target?.result as string;
      const data = parseArchitectureJSON(text);
      if (!data) {
        toast.error(t.sim.invalidFile);
        return;
      }
      if (data.validationError) {
        toast.error(data.validationError);
        return;
      }
      const store = useSimStore.getState();
      store.setNPEs(data.npes);
      store.setConnections(data.connections);
      store.setTrials(data.trials);
      store.setContingencies(data.contingencies);
      store.setHasITI(data.hasITI);
      if (data.modelKind) store.setModelKind(data.modelKind);
      // Restore locked layout if present
      if (data.lockedLayout) {
        store.setLockedLayout(data.lockedLayout);
        setCustomPositions(data.lockedLayout);
      } else {
        setCustomPositions(null); // reset layout for new data
      }
      toast.success(t.sim.experimentLoaded);
    };
    reader.readAsText(file);
    event.target.value = '';
  }, [toast, t]);

  // NPE form state
  const [npeName, setNpeName] = useState('');
  const [npeError, setNpeError] = useState('');
  const [npeType, setNpeType] = useState<NPEType>('Excitatory');
  const [npeLayer, setNpeLayer] = useState<NPELayer>('PrimarySensory');
  const [npeActivation, setNpeActivation] = useState(0);
  const [npeTau, setNpeTau] = useState(0.1);
  const [npeKappa, setNpeKappa] = useState(0.1);
  const [npeMu, setNpeMu] = useState(0.2);
  const [npeSigma, setNpeSigma] = useState(0.15);
  const [npeLogis, setNpeLogis] = useState(0.1);

  // Connection form state
  const [connPre, setConnPre] = useState('');
  const [connPost, setConnPost] = useState('');
  const [connWeight, setConnWeight] = useState(0.1);
  const [connAlpha, setConnAlpha] = useState(0.5);
  const [connBeta, setConnBeta] = useState(0.12);
  const [connAlphaP, setConnAlphaP] = useState(0.5);
  const [connBetaP, setConnBetaP] = useState(0.12);

  const handleAddNPE = useCallback(() => {
    if (!npeName.trim()) return;
    if (npes.some(n => n.name === npeName.trim())) {
      setNpeError(t.network.duplicateName);
      return;
    }
    setNpeError('');
    addNPE({
      name: npeName.trim(),
      type: npeType,
      layer: npeLayer,
      activation: npeActivation,
      temporalSummation: npeTau,
      activationDecay: npeKappa,
      mu: npeMu,
      sigma: npeSigma,
      logisSigma: npeLogis,
    });
    setNpeName('');
    setCustomPositions(null); // reset to trigger re-layout
  }, [npeName, npeType, npeLayer, npeActivation, npeTau, npeKappa, npeMu, npeSigma, npeLogis, npes, addNPE, t]);

  const handleAddConnection = useCallback(() => {
    if (!connPre || !connPost || connPre === connPost) return;
    if (connections.some(c => c.presynapticNPE === connPre && c.postsynapticNPE === connPost)) return;
    addConnection({
      presynapticNPE: connPre,
      postsynapticNPE: connPost,
      weight: connWeight,
      alpha: connAlpha,
      beta: connBeta,
      alphaPrime: connAlphaP,
      betaPrime: connBetaP,
    });
  }, [connPre, connPost, connWeight, connAlpha, connBeta, connAlphaP, connBetaP, connections, addConnection]);

  // Build React Flow nodes from store data (only used for traditional + modern; transduction uses OrbitalView)
  const baseNodes = useMemo(() => {
    // In transduction mode, React Flow is not rendered — return empty
    if (layoutMode === 'transduction') return [];

    let positions: Record<string, { x: number; y: number }>;

    if (customPositions) {
      positions = customPositions;
    } else if (layoutMode === 'modern') {
      positions = computeModernLayout(npes);
    } else {
      positions = computeAutoLayout(npes);
    }

    // Layer index for staggered animations (unique per layer)
    const layerIndexMap: Record<string, number> = {
      PrimarySensory: 0, US: 1, AssociativeSensory: 2, Hippocampal: 3,
      AssociativeMotor: 4, Dopaminergic: 5, PrimaryMotor: 6,
    };

    // Build lookup of which NPEs have source/target connections
    const hasSourceSet = new Set(connections.map(c => c.presynapticNPE));
    const hasTargetSet = new Set(connections.map(c => c.postsynapticNPE));

    const realNodes: Node[] = npes.map((npe) => {
      const pos = positions[npe.name] || { x: 0, y: 0 };
      const color = layerColors[npe.layer] || '#64748b';
      const sz = modernNodeSize[npe.layer] || 56;

      if (layoutMode === 'modern') {
        // Modern: SVG glow circles with pulsing animations
        return {
          id: npe.name,
          type: 'glowNode',
          position: pos,
          data: {
            label: getDisplayName(npe.name),
            color,
            size: sz,
            hasSource: hasSourceSet.has(npe.name),
            hasTarget: hasTargetSet.has(npe.name),
            layerIndex: layerIndexMap[npe.layer] ?? 0,
          },
        };
      }

      // Traditional: original grid style with shape-per-layer
      return {
        id: npe.name,
        position: pos,
        data: { label: getDisplayName(npe.name) },
        style: {
          background: color,
          color: '#fff',
          border: `2px solid ${color}`,
          borderRadius: getNodeShape(npe.layer, npe.type),
          width: 60,
          height: 60,
          display: 'flex',
          alignItems: 'center',
          justifyContent: 'center',
          fontSize: '12px',
          fontWeight: 700,
          fontStyle: 'italic' as const,
          boxShadow: `0 0 20px ${color}33`,
        },
      };
    });

    return realNodes;
  }, [npes, connections, customPositions, layoutMode]);

  // Local nodes state for real-time drag updates
  const [flowNodes, setFlowNodes] = useState<Node[]>(baseNodes);

  // Sync flowNodes when baseNodes change (e.g., new NPE added, auto-layout)
  useEffect(() => {
    setFlowNodes(baseNodes);
  }, [baseNodes]);

  // Handle real-time drag via onNodesChange
  const onNodesChange = useCallback((changes: NodeChange[]) => {
    setFlowNodes(prev => applyNodeChanges(changes, prev));
  }, []);

  const edges: Edge[] = useMemo(() => {
    // In transduction mode, React Flow is not rendered
    if (layoutMode === 'transduction') return [];

    // Build a quick NPE-name → layer lookup for edge coloring
    const npeLayerMap: Record<string, string> = {};
    npes.forEach(n => { npeLayerMap[n.name] = n.layer; });

    return connections.map((conn, idx) => {
      const srcColor = layerColors[npeLayerMap[conn.presynapticNPE]] || '#06b6d4';
      const tgtColor = layerColors[npeLayerMap[conn.postsynapticNPE]] || '#14b8a6';

      if (layoutMode === 'modern') {
        // Modern: curved bezier with gradient + double traveling dots
        return {
          id: `${conn.presynapticNPE}-${conn.postsynapticNPE}`,
          source: conn.presynapticNPE,
          target: conn.postsynapticNPE,
          type: 'animatedEdge',
          data: {
            sourceColor: srcColor,
            targetColor: tgtColor,
            weight: conn.weight,
            edgeIndex: idx,
          },
          label: conn.weight.toFixed(2),
          labelStyle: { fill: '#94a3b8', fontSize: 10, fontWeight: 600 },
          labelBgStyle: { fill: 'rgba(15, 23, 42, 0.85)', fillOpacity: 0.9 },
        };
      }

      // Traditional edges
      return {
        id: `${conn.presynapticNPE}-${conn.postsynapticNPE}`,
        source: conn.presynapticNPE,
        target: conn.postsynapticNPE,
        animated: true,
        style: {
          stroke: conn.weight >= 0.999 ? '#ef4444' : '#64748b',
          strokeWidth: Math.max(1.5, Math.abs(conn.weight) * 4),
        },
        markerEnd: { type: MarkerType.ArrowClosed, color: conn.weight >= 0.999 ? '#ef4444' : '#64748b' },
        label: conn.weight.toFixed(2),
        labelStyle: { fill: '#94a3b8', fontSize: 10, fontWeight: 600 },
        labelBgStyle: { fill: '#f8fafc', fillOpacity: 0.9 },
      };
    });
  }, [connections, layoutMode, npes]);

  // Handle node drag to persist positions (skip virtual nodes)
  const onNodeDragStop = useCallback<NodeMouseHandler>((_, node) => {
    if (node.id.startsWith('__virtual_')) return;
    setCustomPositions(prev => ({
      ...(prev || {}),
      [node.id]: node.position,
    }));
  }, []);

  return (
    <div className="max-w-7xl mx-auto space-y-6">
      <div className="flex items-center gap-3">
        <div className="w-10 h-10 rounded-xl bg-cyan-500/10 flex items-center justify-center">
          <Brain size={20} className="text-cyan-600" />
        </div>
        <div>
          <div className="flex items-center gap-2">
            <h1 className="text-2xl font-bold text-slate-800">{t.network.pageTitle}</h1>
            <Badge variant={('muted')}>{modelKind.toUpperCase()}</Badge>
          </div>
          <p className="text-sm text-slate-500">{t.network.pageSubtitle}</p>
        </div>
      </div>

      <div className="grid grid-cols-1 lg:grid-cols-5 gap-6">
        {/* Canvas - 3 cols */}
        <Card className="lg:col-span-3 p-0 overflow-hidden" glow="cyan">
          <div ref={flowRef} className={`h-[550px] transition-colors duration-500 ${isDarkCanvas ? 'dark-flow' : ''}`}
               style={isDarkCanvas ? { background: layoutMode === 'transduction' ? '#110d18' : '#0b1120' } : undefined}>
            {layoutMode === 'publication' ? (
              <PublicationNetwork npes={npes} connections={connections} positions={publicationPositions}
                onPositionsChange={handlePublicationPositions} svgRef={publicationRef} language={language} />
            ) : layoutMode === 'transduction' ? (
              /* ── Transduction: BrainView (sagittal brain silhouette — no React Flow) ── */
              <div style={{ width: '100%', height: '100%', display: 'flex', alignItems: 'center', justifyContent: 'center' }}>
                <BrainView
                  npes={npes}
                  connections={connections}
                  stimulusLabel={t.network.stimulus}
                  responseLabel={t.network.response}
                />
              </div>
            ) : (
              /* ── Traditional / Modern: React Flow ── */
              <ReactFlow
                nodes={flowNodes}
                edges={edges}
                onNodesChange={onNodesChange}
                onNodeDragStop={onNodeDragStop}
                nodeTypes={layoutMode === 'modern' ? modernNodeTypes : undefined}
                edgeTypes={layoutMode === 'modern' ? modernEdgeTypes : undefined}
                fitView
                proOptions={{ hideAttribution: true }}
                style={isDarkCanvas ? { background: 'transparent' } : undefined}
              >
                {layoutMode === 'modern' ? (
                  <>
                    {/* Modern: cold dark navy with cyan/purple ambient glow */}
                    <Background color="rgba(100, 116, 139, 0.08)" gap={24} size={1} />
                    <svg style={{ position: 'absolute', inset: 0, width: '100%', height: '100%', pointerEvents: 'none', zIndex: 0 }}>
                      <defs>
                        <radialGradient id="canvasAmbientModern" cx="50%" cy="50%" r="50%">
                          <stop offset="0%" stopColor="#06b6d4" stopOpacity="0.04" />
                          <stop offset="60%" stopColor="#8b5cf6" stopOpacity="0.02" />
                          <stop offset="100%" stopColor="#0b1120" stopOpacity="0" />
                        </radialGradient>
                      </defs>
                      <rect width="100%" height="100%" fill="url(#canvasAmbientModern)" />
                    </svg>
                  </>
                ) : (
                  <Background color="#e2e8f0" gap={20} />
                )}
                <Controls
                  className={layoutMode === 'modern'
                    ? '!rounded-xl !shadow-md'
                    : '!bg-white !border-slate-200 !rounded-xl !shadow-md'
                  }
                />
              </ReactFlow>
            )}
          </div>
          {/* Layer legend */}
          {layoutMode === 'publication' ? (
            <p className="px-4 py-3 text-xs text-slate-600 border-t border-slate-200 bg-white">
              {language === 'es' ? 'Campos difusos = modulación del aprendizaje, no conexiones. Arrastra las unidades para ajustar la figura. Exportación limpia, sin controles.' : 'Diffuse fields modulate learning; they are not connections. Drag units to arrange your figure. Clean export without controls.'}
            </p>
          ) : <div className={`flex flex-wrap items-center gap-2 px-3 pt-3 pb-1 border-t transition-colors duration-500 ${
            isDarkCanvas
              ? 'border-slate-700/50 bg-slate-900/80'
              : 'border-slate-200 bg-slate-50'
          }`}>
            {Object.entries(layerColors).filter(([layer]) => npes.some(npe => npe.layer === layer)).map(([layer, color]) => (
              <div key={layer} className={`flex items-center gap-1.5 text-xs ${isDarkCanvas ? 'text-slate-400' : 'text-slate-500'}`}>
                <div className="w-3 h-3 rounded-full" style={{ backgroundColor: color, boxShadow: isDarkCanvas ? `0 0 6px ${color}66` : undefined }} />
                {layerDisplayNames[layer as NPELayer] || layer}
              </div>
            ))}
          </div>}
          {/* Action buttons */}
          <div className={`flex flex-wrap items-center gap-2 p-3 transition-colors duration-500 ${isDarkCanvas ? 'bg-slate-900/80' : 'bg-slate-50'}`}>
            <button
              onClick={handleDownloadPng}
              disabled={npes.length === 0 || isExporting}
              className={`flex items-center gap-1.5 px-3 py-1.5 rounded-lg text-xs font-bold border transition-colors ${
                isDarkCanvas
                  ? 'text-violet-300 border-violet-500/30 bg-violet-500/10 hover:bg-violet-500/20'
                  : 'text-violet-600 border-violet-200 bg-violet-50 hover:bg-violet-100'
              }`}
            >
              <Download size={12} /> {t.network.exportPNG}
            </button>
            {layoutMode === 'publication' && <button type="button" onClick={handleDownloadSvg} disabled={npes.length === 0}
              className="px-3 py-1.5 rounded-lg text-xs font-semibold text-slate-700 border border-slate-300 bg-white hover:bg-slate-100 disabled:opacity-40 focus-visible:outline-2 focus-visible:outline-cyan-600">
              SVG
            </button>}

            {/* Layout mode switcher */}
            <div className={`flex gap-0.5 p-0.5 rounded-lg border transition-colors duration-500 ${
              isDarkCanvas
                ? 'bg-slate-800/60 border-slate-700/50'
                : 'bg-slate-200/60 border-slate-200'
            }`}>
              {([
                { mode: 'publication' as LayoutMode, icon: LayoutGrid, label: language === 'es' ? 'Publicación' : 'Publication' },
                { mode: 'traditional' as LayoutMode, icon: LayoutGrid, label: t.network.layoutTraditional },
                { mode: 'modern' as LayoutMode, icon: Sparkles, label: t.network.layoutModern },
                { mode: 'transduction' as LayoutMode, icon: ArrowRightLeft, label: t.network.layoutTransduction },
              ]).map(({ mode, icon: Icon, label }) => (
                <button
                  key={mode}
                  onClick={() => {
                    setLayoutMode(mode);
                    setCustomPositions(null); // trigger re-layout
                    if (mode !== 'publication') setTimeout(() => fitView({ padding: 0.2, duration: 400 }), 50);
                  }}
                  className={`flex items-center gap-1 px-2 py-1 rounded-md text-[11px] font-bold transition-all ${
                    layoutMode === mode
                      ? isDarkCanvas
                        ? 'bg-slate-700 text-cyan-400 shadow-sm shadow-cyan-500/20'
                        : 'bg-white text-cyan-600 shadow-sm'
                      : isDarkCanvas
                        ? 'text-slate-300 hover:text-white'
                        : 'text-slate-400 hover:text-slate-600'
                  }`}
                  title={label}
                >
                  <Icon size={12} />
                  <span className={mode === 'publication' ? '' : 'hidden xl:inline'}>{label}</span>
                </button>
              ))}
            </div>

            <div className="flex-1" />
            <button
              onClick={handleAutoLayout}
              className={`flex items-center gap-1.5 px-3 py-1.5 rounded-lg text-xs font-bold border transition-colors ${
                isDarkCanvas
                  ? 'text-cyan-300 border-cyan-500/30 bg-cyan-500/10 hover:bg-cyan-500/20'
                  : 'text-cyan-600 border-cyan-200 bg-cyan-50 hover:bg-cyan-100'
              }`}
              title={t.network.organizeTitle}
            >
              <LayoutGrid size={12} /> {t.network.organize}
            </button>
            <button
              onClick={handleLockLayout}
              className={`flex items-center gap-1.5 px-3 py-1.5 rounded-lg text-xs font-bold border transition-colors ${
                isLocked
                  ? isDarkCanvas
                    ? 'text-emerald-400 border-emerald-500/30 bg-emerald-500/10'
                    : 'text-emerald-600 border-emerald-300 bg-emerald-50'
                  : isDarkCanvas
                    ? 'text-slate-300 border-slate-600 hover:bg-slate-700/50'
                    : 'text-slate-500 border-slate-200 hover:bg-slate-100'
              }`}
              title={t.network.lockTitle}
            >
              <Lock size={12} /> {isLocked ? t.network.locked : t.network.lock}
            </button>
          </div>
        </Card>

        {/* Editor panel - 2 cols */}
        <div className="lg:col-span-2 space-y-4">
          {/* Tab toggle */}
          <div className="flex gap-1 p-1 rounded-xl bg-slate-100 border border-slate-200">
            <button
              onClick={() => setActiveTab('units')}
              className={`flex-1 flex items-center justify-center gap-2 py-2.5 rounded-lg text-sm font-bold transition-all ${
                activeTab === 'units'
                  ? 'bg-cyan-50 text-cyan-600 shadow-sm'
                  : 'text-slate-500 hover:text-slate-400'
              }`}
            >
              <Brain size={16} /> {t.network.tabUnits}
            </button>
            <button
              onClick={() => setActiveTab('connections')}
              className={`flex-1 flex items-center justify-center gap-2 py-2.5 rounded-lg text-sm font-bold transition-all ${
                activeTab === 'connections'
                  ? 'bg-cyan-50 text-cyan-600 shadow-sm'
                  : 'text-slate-500 hover:text-slate-400'
              }`}
            >
              <GitBranch size={16} /> {t.network.tabConnections}
            </button>
          </div>

          {/* Import/Export */}
          <div className="flex gap-2">
            <button
              onClick={() => {
                const { npes, connections, trials, contingencies, hasITI, lockedLayout } = useSimStore.getState();
                downloadArchitectureJSON(
                  npes, connections, trials, contingencies, hasITI, undefined,
                  { modelKind }, lockedLayout,
                );
              }}
              className="flex-1 flex items-center justify-center gap-1.5 py-2 rounded-lg text-xs font-bold text-slate-500 border border-slate-200 hover:bg-slate-100 transition-colors"
            >
              <Download size={12} /> {t.network.export}
            </button>
            <button
              onClick={() => fileInputRef.current?.click()}
              className="flex-1 flex items-center justify-center gap-1.5 py-2 rounded-lg text-xs font-bold text-slate-500 border border-slate-200 hover:bg-slate-100 transition-colors"
            >
              <Upload size={12} /> {t.network.import}
            </button>
            <input
              ref={fileInputRef}
              type="file"
              accept=".json"
              onChange={handleImportArchitecture}
              className="hidden"
            />
          </div>

          <AnimatePresence mode="wait">
            {activeTab === 'units' ? (
              <motion.div
                key="units"
                initial={{ opacity: 0, x: -10 }}
                animate={{ opacity: 1, x: 0 }}
                exit={{ opacity: 0, x: 10 }}
              >
                <Card className="p-5 space-y-4">
                  <h3 className="text-sm font-bold text-slate-700 flex items-center gap-2">
                    <Plus size={14} className="text-cyan-600" /> {t.network.addNPE}
                  </h3>
                  <div className="space-y-3">
                    <div>
                      <label className="block text-xs font-semibold text-slate-500 mb-1">{t.network.name}</label>
                      <input
                        type="text"
                        value={npeName}
                        onChange={(e) => { setNpeName(e.target.value); setNpeError(''); }}
                        placeholder={t.network.namePlaceholder}
                        className={`w-full px-3 py-2 rounded-lg bg-white border text-slate-800 text-sm focus:outline-none transition-colors ${npeError ? 'border-rose-400 focus:border-rose-500' : 'border-slate-200 focus:border-cyan-500/50'}`}
                      />
                      {npeError && <p className="text-xs text-rose-500 mt-1">{npeError}</p>}
                    </div>
                    <div className="grid grid-cols-2 gap-3">
                      <div>
                        <label className="block text-xs font-semibold text-slate-500 mb-1">{t.network.type}</label>
                        <select
                          value={npeType}
                          onChange={(e) => setNpeType(e.target.value as NPEType)}
                          className="w-full px-3 py-2 rounded-lg bg-white border border-slate-200 text-slate-800 text-sm focus:border-cyan-500/50 focus:outline-none"
                        >
                          <option value="Excitatory">{t.network.excitatory}</option>
                          <option value="Inhibitory">{t.network.inhibitory}</option>
                        </select>
                      </div>
                      <div>
                        <label className="block text-xs font-semibold text-slate-500 mb-1">{t.network.layer}</label>
                        <select
                          value={npeLayer}
                          onChange={(e) => setNpeLayer(e.target.value as NPELayer)}
                          className="w-full px-3 py-2 rounded-lg bg-white border border-slate-200 text-slate-800 text-sm focus:border-cyan-500/50 focus:outline-none"
                        >
                          {layerOptions.map(l => <option key={l} value={l}>{layerDisplayNames[l]}</option>)}
                        </select>
                      </div>
                    </div>

                    {appMode === 'advanced' && (
                    <div className="space-y-3 pt-2 border-t border-slate-200">
                      <div className="flex items-center gap-2 text-xs text-slate-400">
                        <Settings2 size={12} /> {t.network.freeParameters}
                      </div>
                      <div className="grid grid-cols-3 gap-2">
                        {[
                          { label: t.network.activation, value: npeActivation, set: setNpeActivation },
                          { label: t.network.tauSum, value: npeTau, set: setNpeTau },
                          { label: t.network.kappaDecay, value: npeKappa, set: setNpeKappa },
                          { label: t.network.muMean, value: npeMu, set: setNpeMu },
                          { label: t.network.sigmaDev, value: npeSigma, set: setNpeSigma },
                          { label: t.network.logistic, value: npeLogis, set: setNpeLogis },
                        ].map(({ label, value, set }) => (
                          <div key={label}>
                            <label className="block text-[10px] font-semibold text-slate-400 mb-0.5">{label}</label>
                            <input
                              type="number"
                              step="0.01"
                              min={('0')}
                              max={('1')}
                              value={value}
                              onChange={(e) => set(parseFloat(e.target.value) || 0)}
                              className="w-full px-2 py-1.5 rounded-md bg-white border border-slate-200 text-slate-700 text-xs focus:border-cyan-500/50 focus:outline-none"
                            />
                          </div>
                        ))}
                      </div>
                    </div>
                    )}

                    <button
                      onClick={handleAddNPE}
                      disabled={!npeName.trim()}
                      className="w-full py-2.5 rounded-xl bg-gradient-to-r from-cyan-500 to-teal-500 text-white font-bold text-sm disabled:opacity-40 disabled:cursor-not-allowed hover:shadow-lg hover:shadow-cyan-500/20 transition-shadow"
                    >
                      {t.network.addUnit}
                    </button>
                  </div>
                </Card>

                {/* NPE list */}
                <div className="mt-4 space-y-2">
                  {npes.map((npe) => (
                    <motion.div
                      key={npe.name}
                      initial={{ opacity: 0, x: -10 }}
                      animate={{ opacity: 1, x: 0 }}
                      className="flex items-center justify-between p-3 rounded-xl bg-white border border-slate-200"
                    >
                      <div className="flex items-center gap-3">
                        <div
                          className="w-8 h-8 rounded-full flex items-center justify-center text-white text-xs font-bold"
                          style={{ backgroundColor: layerColors[npe.layer] }}
                        >
                          {npe.name.charAt(0)}
                        </div>
                        <div>
                          <div className="flex items-center gap-1.5">
                            <span className="text-sm font-bold text-slate-700 italic">{getDisplayName(npe.name)}</span>
                            {getDisplayName(npe.name) !== npe.name && (
                              <span className="text-xs text-slate-400">({npe.name})</span>
                            )}
                            {disconnectedNPEs.has(npe.name) && (
                              <span title={t.network.disconnectedWarning} className="text-amber-500">
                                <AlertTriangle size={12} />
                              </span>
                            )}
                          </div>
                          <div className="flex gap-1.5 mt-0.5">
                            <Badge variant="info">{npe.layer}</Badge>
                            <Badge variant="muted">{npe.type}</Badge>
                            {disconnectedNPEs.has(npe.name) && (
                              <Badge variant="muted"><span className="text-amber-600">{t.network.disconnected}</span></Badge>
                            )}
                          </div>
                        </div>
                      </div>
                      <button
                        onClick={async () => {
                          const hasConns = connections.some(c => c.presynapticNPE === npe.name || c.postsynapticNPE === npe.name);
                          if (hasConns) {
                            const ok = await confirm({
                              title: t.confirm.deleteNPETitle,
                              message: t.confirm.deleteNPEMessage,
                              confirmLabel: t.common.delete,
                              cancelLabel: t.confirm.cancel,
                              variant: 'danger',
                            });
                            if (!ok) return;
                          }
                          removeNPE(npe.name);
                        }}
                        className="p-1.5 rounded-lg text-slate-400 hover:text-rose-600 hover:bg-rose-50 transition-colors"
                      >
                        <Trash2 size={14} />
                      </button>
                    </motion.div>
                  ))}
                </div>
              </motion.div>
            ) : (
              <motion.div
                key="connections"
                initial={{ opacity: 0, x: 10 }}
                animate={{ opacity: 1, x: 0 }}
                exit={{ opacity: 0, x: -10 }}
              >
                <Card className="p-5 space-y-4">
                  <h3 className="text-sm font-bold text-slate-700 flex items-center gap-2">
                    <Plus size={14} className="text-cyan-600" /> {t.network.addConnection}
                  </h3>
                  <div className="space-y-3">
                    <div className="grid grid-cols-2 gap-3">
                      <div>
                        <label className="block text-xs font-semibold text-slate-500 mb-1">{t.network.source}</label>
                        <select
                          value={connPre}
                          onChange={(e) => setConnPre(e.target.value)}
                          className="w-full px-3 py-2 rounded-lg bg-white border border-slate-200 text-slate-800 text-sm focus:border-cyan-500/50 focus:outline-none"
                        >
                          <option value="">{t.network.select}</option>
                          {npes.map(n => <option key={n.name} value={n.name}>{n.name}</option>)}
                        </select>
                      </div>
                      <div>
                        <label className="block text-xs font-semibold text-slate-500 mb-1">{t.network.target}</label>
                        <select
                          value={connPost}
                          onChange={(e) => setConnPost(e.target.value)}
                          className="w-full px-3 py-2 rounded-lg bg-white border border-slate-200 text-slate-800 text-sm focus:border-cyan-500/50 focus:outline-none"
                        >
                          <option value="">{t.network.select}</option>
                          {npes.map(n => <option key={n.name} value={n.name}>{n.name}</option>)}
                        </select>
                      </div>
                    </div>
                    <div>
                      <label className="block text-xs font-semibold text-slate-500 mb-1">{t.network.weight}: {connWeight.toFixed(2)}</label>
                      <input
                        type="range"
                        min={('0')}
                        max={('1')}
                        step="0.01"
                        value={connWeight}
                        onChange={(e) => setConnWeight(parseFloat(e.target.value))}
                        className="w-full accent-cyan-500"
                      />
                    </div>

                    {appMode === 'advanced' && (
                    <div className="space-y-3 pt-2 border-t border-slate-200">
                      <div className="flex items-center gap-2 text-xs text-slate-400">
                        <Settings2 size={12} /> {t.network.freeParameters}
                      </div>
                      <div className="grid grid-cols-2 gap-2">
                        {[
                          { label: t.network.alphaIncrement, value: connAlpha, set: setConnAlpha, isBeta: false },
                          { label: t.network.betaDecrement, value: connBeta, set: setConnBeta, isBeta: true },
                          { label: t.network.alphaprimeInhibInc, value: connAlphaP, set: setConnAlphaP, isBeta: false },
                          { label: t.network.betaprimeInhibDec, value: connBetaP, set: setConnBetaP, isBeta: true },
                        ].map(({ label, value, set, isBeta }) => (
                          <div key={label}>
                            <label className="block text-[10px] font-semibold text-slate-400 mb-0.5">{label}</label>
                            <input
                              type="number"
                              step="0.01"
                              min="0"
                              max="1"
                              value={isBeta && value === 0.12 ? 0.1 : value}
                              onChange={(e) => {
                                const v = parseFloat(e.target.value);
                                // When user types 0.1 in a beta field, store 0.12 (historical default)
                                if (isBeta && v === 0.1) {
                                  set(0.12);
                                } else {
                                  set(v || 0);
                                }
                              }}
                              className="w-full px-2 py-1.5 rounded-md bg-white border border-slate-200 text-slate-700 text-xs focus:border-cyan-500/50 focus:outline-none"
                            />
                          </div>
                        ))}
                      </div>
                    </div>
                    )}

                    <button
                      onClick={handleAddConnection}
                      disabled={!connPre || !connPost || connPre === connPost}
                      className="w-full py-2.5 rounded-xl bg-gradient-to-r from-cyan-500 to-teal-500 text-white font-bold text-sm disabled:opacity-40 disabled:cursor-not-allowed hover:shadow-lg hover:shadow-cyan-500/20 transition-shadow"
                    >
                      {t.network.addConnection}
                    </button>
                  </div>
                </Card>

                {/* Connection list */}
                <div className="mt-4 space-y-2">
                  {connections.map((conn) => (
                    <motion.div
                      key={`${conn.presynapticNPE}-${conn.postsynapticNPE}`}
                      initial={{ opacity: 0, x: 10 }}
                      animate={{ opacity: 1, x: 0 }}
                      className="flex items-center justify-between p-3 rounded-xl bg-white border border-slate-200"
                    >
                      <div className="flex items-center gap-2">
                        <span className="text-sm font-bold text-slate-700">{conn.presynapticNPE}</span>
                        <span className="text-slate-400">{'\u2192'}</span>
                        <span className="text-sm font-bold text-slate-700">{conn.postsynapticNPE}</span>
                        <Badge variant="info">w={conn.weight.toFixed(2)}</Badge>
                      </div>
                      <button
                        onClick={() => removeConnection(conn.presynapticNPE, conn.postsynapticNPE)}
                        className="p-1.5 rounded-lg text-slate-400 hover:text-rose-600 hover:bg-rose-50 transition-colors"
                      >
                        <Trash2 size={14} />
                      </button>
                    </motion.div>
                  ))}
                </div>
              </motion.div>
            )}
          </AnimatePresence>
        </div>
      </div>
    </div>
  );
}

// Wrap with ReactFlowProvider so useReactFlow() works
export function NetworkBuilder() {
  return (
    <ReactFlowProvider>
      <NetworkBuilderInner />
    </ReactFlowProvider>
  );
}
