import { useId, useMemo, useRef, useState } from 'react';
import type { KeyboardEvent, PointerEvent, RefObject } from 'react';
import type { Connection, NPE } from '../../types/ddm';
import { getDisplayName } from '../../utils/displayNames';
import {
  computePublicationLayout, publicationBounds, publicationNodeRadius, publicationNodeStyle,
} from '../../utils/publicationNetwork';
import type { PublicationBounds, PublicationPoint, PublicationPositions } from '../../utils/publicationNetwork';

export interface PublicationNetworkProps {
  npes: NPE[];
  connections: Connection[];
  /** Centres in SVG document coordinates. These are NOT React Flow node positions. */
  positions?: PublicationPositions;
  onPositionsChange?: (positions: PublicationPositions) => void;
  onSelectUnit?: (name: string) => void;
  selectedUnit?: string | null;
  activations?: Record<string, number>;
  /** Connection keys use the simulation's `${presynapticNPE}-${postsynapticNPE}` convention. */
  weights?: Record<string, number>;
  showValues?: boolean;
  language?: 'es' | 'en';
  svgRef?: RefObject<SVGSVGElement | null>;
  className?: string;
}

interface EdgeGeometry {
  path: string;
  end: PublicationPoint;
  tangent: PublicationPoint;
  label: PublicationPoint;
  bounds: PublicationPoint[];
}

function pointOnUnit(npe: NPE, centre: PublicationPoint, toward: PublicationPoint, margin = 0): PublicationPoint {
  const dx = toward.x - centre.x, dy = toward.y - centre.y;
  const distance = Math.hypot(dx, dy) || 1;
  const ux = dx / distance, uy = dy / distance;
  const radius = publicationNodeRadius(npe);
  const shape = publicationNodeStyle(npe).shape;
  const boundary = shape === 'square' ? radius / Math.max(Math.abs(ux), Math.abs(uy), .01)
    : shape === 'diamond' ? radius * 1.22 / Math.max(Math.abs(ux) + Math.abs(uy), .01)
      : radius;
  return { x: centre.x + ux * (boundary + margin), y: centre.y + uy * (boundary + margin) };
}

function edgeGeometry(connection: Connection, index: number, all: Connection[], npes: NPE[], positions: PublicationPositions): EdgeGeometry | null {
  const from = npes.find(n => n.name === connection.presynapticNPE);
  const to = npes.find(n => n.name === connection.postsynapticNPE);
  if (!from || !to) return null; // A malformed dangling connection is not a real drawable synapse.
  const a = positions[from.name], b = positions[to.name];
  const pathPoint = (p: PublicationPoint) => `${p.x} ${p.y}`;
  const cubic = (c1: PublicationPoint, c2: PublicationPoint) => {
    const start = pointOnUnit(from, a, c1), end = pointOnUnit(to, b, c2, 4);
    return { path: `M ${pathPoint(start)} C ${pathPoint(c1)} ${pathPoint(c2)} ${pathPoint(end)}`, end,
      tangent: { x: end.x - c2.x, y: end.y - c2.y },
      label: { x: (start.x + 3 * c1.x + 3 * c2.x + end.x) / 8, y: (start.y + 3 * c1.y + 3 * c2.y + end.y) / 8 - 11 },
      bounds: [start, c1, c2, end] };
  };
  if (from.name === to.name || Math.hypot(a.x - b.x, a.y - b.y) < 2) {
    const offset = publicationNodeRadius(from) + 90;
    return cubic({ x: a.x - offset, y: a.y - offset * 1.5 }, { x: b.x + offset, y: b.y - offset * 1.5 });
  }
  if (from.layer === 'US' && to.layer === 'PrimaryMotor') {
    const direct = all.slice(0, index).filter(c => npes.some(n => n.name === c.presynapticNPE && n.layer === 'US') && npes.some(n => n.name === c.postsynapticNPE && n.layer === 'PrimaryMotor')).length;
    const right = Math.max(...npes.map(n => positions[n.name].x + publicationNodeRadius(n))) + 62 + direct * 24;
    const bottom = Math.max(...npes.map(n => positions[n.name].y + publicationNodeRadius(n))) + 46 + direct * 21;
    const outputYs = npes.filter(n => n.layer === 'PrimaryMotor').map(n => positions[n.name].y);
    const upper = b.y <= (Math.min(...outputYs) + Math.max(...outputYs)) / 2;
    const portY = b.y + (upper ? -1 : 1) * (publicationNodeRadius(to) + 29);
    const corner = { x: b.x + publicationNodeRadius(to) + 34, y: portY };
    const first = { x: a.x + publicationNodeRadius(from) + 34 + direct * 13, y: bottom };
    const start = pointOnUnit(from, a, first), end = pointOnUnit(to, b, corner, 4);
    const points = [start, first, { x: right, y: bottom }, { x: right, y: portY }, corner, end];
    return { path: points.map((p, i) => `${i === 0 ? 'M' : 'L'} ${pathPoint(p)}`).join(' '), end,
      tangent: { x: end.x - corner.x, y: end.y - corner.y }, label: { x: right - 20, y: (bottom + portY) / 2 }, bounds: points };
  }
  if (to.layer === 'Dopaminergic' && from.layer !== 'US') {
    const incoming = all.filter(c => c.postsynapticNPE === to.name && npes.find(n => n.name === c.presynapticNPE)?.layer !== 'US');
    const rank = incoming.indexOf(connection);
    const side = rank % 2 === 0 ? -1 : 1;
    const offset = 90 + Math.floor(rank / 2) * 36;
    return cubic({ x: a.x + side * offset, y: a.y + (b.y - a.y) * .18 }, { x: b.x + side * offset, y: b.y - (b.y - a.y) * .16 });
  }
  const dx = b.x - a.x, dy = b.y - a.y, distance = Math.hypot(dx, dy);
  const reverse = all.some(c => c.presynapticNPE === to.name && c.postsynapticNPE === from.name);
  const duplicate = all.slice(0, index).filter(c => c.presynapticNPE === from.name && c.postsynapticNPE === to.name).length;
  const obstacle = npes.some(n => {
    if (n.name === from.name || n.name === to.name) return false;
    const p = positions[n.name];
    const t = ((p.x - a.x) * dx + (p.y - a.y) * dy) / (distance * distance);
    return t > .07 && t < .93 && Math.hypot(p.x - (a.x + t * dx), p.y - (a.y + t * dy)) < publicationNodeRadius(n) + 12;
  });
  if (reverse || duplicate || obstacle || Math.abs(dx) < 8) {
    const bend = (obstacle ? 135 : 65) + duplicate * 32;
    const perpendicular = { x: -dy / distance * bend, y: dx / distance * bend };
    return cubic({ x: a.x + dx / 3 + perpendicular.x, y: a.y + dy / 3 + perpendicular.y },
      { x: a.x + dx * 2 / 3 + perpendicular.x, y: a.y + dy * 2 / 3 + perpendicular.y });
  }
  const start = pointOnUnit(from, a, b), end = pointOnUnit(to, b, a, 4);
  return { path: `M ${pathPoint(start)} L ${pathPoint(end)}`, end, tangent: { x: dx, y: dy },
    label: { x: (start.x + end.x) / 2, y: (start.y + end.y) / 2 - 11 }, bounds: [start, end] };
}

function DiffuseField({ nodes, positions, id, blue }: { nodes: NPE[]; positions: PublicationPositions; id: string; blue: boolean }) {
  if (!nodes.length) return null;
  const left = Math.min(...nodes.map(n => positions[n.name].x - publicationNodeRadius(n))) - 31;
  const right = Math.max(...nodes.map(n => positions[n.name].x + publicationNodeRadius(n))) + 42;
  const top = Math.min(...nodes.map(n => positions[n.name].y - publicationNodeRadius(n))) - 40;
  const bottom = Math.max(...nodes.map(n => positions[n.name].y + publicationNodeRadius(n))) + 37;
  const source = nodes.filter(n => n.layer === (blue ? 'Dopaminergic' : 'Hippocampal'));
  const points = (source.length ? source : nodes).map(n => positions[n.name]);
  const cx = points.reduce((sum, p) => sum + p.x, 0) / points.length;
  const cy = points.reduce((sum, p) => sum + p.y, 0) / points.length;
  const color = blue ? '#339CC9' : '#D1A337';
  return <g aria-label={blue ? 'dD: motor and dopaminergic learning modulation' : 'dH: sensory and hippocampal learning modulation'}>
    <defs>
      <radialGradient id={id} cx={`${(cx - left) / (right - left) * 100}%`} cy={`${(cy - top) / (bottom - top) * 100}%`} r="88%">
        <stop offset="0%" stopColor={color} stopOpacity=".28" />
        <stop offset="45%" stopColor={color} stopOpacity=".11" />
        <stop offset="82%" stopColor={color} stopOpacity=".025" />
        <stop offset="100%" stopColor={color} stopOpacity="0" />
      </radialGradient>
      <filter id={`${id}-soft`} x="-35%" y="-35%" width="170%" height="170%" colorInterpolationFilters="sRGB">
        <feGaussianBlur stdDeviation="12" />
      </filter>
    </defs>
    <rect x={left} y={top} width={right - left} height={bottom - top} rx="32" fill={`url(#${id})`} filter={`url(#${id}-soft)`} />
    <text x={left + 13} y={top + 19} fill={blue ? '#368AB9' : '#B88817'} fontSize="18" fontStyle="italic">
      d<tspan baselineShift="sub" fontSize="12">{blue ? 'D,t' : 'H,t'}</tspan>
    </text>
  </g>;
}

/** Native, self-contained SVG: publication artwork never includes editor chrome.
 * The manuscript's palette and diffuse fields are a pinned author requirement.
 * Only explicit connections are drawn; modulating fields never imply D axons.
 */
export function PublicationNetwork({ npes, connections, positions, onPositionsChange, onSelectUnit, selectedUnit, activations, weights, showValues = false, language = 'en', svgRef, className }: PublicationNetworkProps) {
  const localSvgRef = useRef<SVGSVGElement | null>(null);
  const ref = svgRef ?? localSvgRef;
  const uid = useId().replace(/[^a-zA-Z0-9_-]/g, '');
  const [localPositions, setLocalPositions] = useState<PublicationPositions>({});
  const [dragBounds, setDragBounds] = useState<PublicationBounds | null>(null);
  const drag = useRef<{ name: string; pointerId: number; offset: PublicationPoint } | null>(null);
  const resolved = useMemo(() => {
    const defaults = computePublicationLayout(npes, connections);
    npes.forEach(n => {
      const p = positions?.[n.name] ?? localPositions[n.name];
      if (p && Number.isFinite(p.x) && Number.isFinite(p.y)) defaults[n.name] = p;
    });
    return defaults;
  }, [npes, connections, positions, localPositions]);
  const geometry = useMemo(() => connections.map((c, i) => edgeGeometry(c, i, connections, npes, resolved)), [connections, npes, resolved]);
  const contentBounds = useMemo(() => publicationBounds(npes, resolved, geometry.flatMap(g => g ? g.bounds.concat([
    { x: g.label.x - 40, y: g.label.y - 15 }, { x: g.label.x + 40, y: g.label.y + 15 },
  ]) : [])), [npes, resolved, geometry]);
  const box = dragBounds ?? contentBounds;
  const sensory = npes.filter(n => n.layer === 'AssociativeSensory' || n.layer === 'Hippocampal');
  const motor = npes.filter(n => n.layer === 'AssociativeMotor' || n.layer === 'PrimaryMotor' || n.layer === 'Dopaminergic');
  const interactive = Boolean(onPositionsChange || onSelectUnit);
  const title = language === 'es' ? 'Arquitectura DDM — figura para publicación' : 'DDM architecture — publication artwork';
  const description = language === 'es'
    ? 'Los campos difusos indican modulación del aprendizaje, no conexiones sinápticas. Se dibujan únicamente las conexiones de esta red.'
    : 'Diffuse fields indicate learning modulation, not synaptic projections. Only this network’s explicit connections are drawn.';

  const moveUnit = (name: string, p: PublicationPoint) => {
    const next = { ...resolved, [name]: p };
    if (onPositionsChange) onPositionsChange(next);
    else setLocalPositions(next);
  };
  const svgPoint = (event: PointerEvent<SVGGElement>): PublicationPoint | null => {
    const svg = ref.current, matrix = svg?.getScreenCTM();
    if (!svg || !matrix) return null;
    const point = svg.createSVGPoint();
    point.x = event.clientX;
    point.y = event.clientY;
    return point.matrixTransform(matrix.inverse());
  };
  const pointerDown = (event: PointerEvent<SVGGElement>, name: string) => {
    if (event.button !== 0) return;
    onSelectUnit?.(name);
    if (!onPositionsChange) return;
    const point = svgPoint(event);
    if (!point) return;
    event.preventDefault();
    event.currentTarget.focus();
    event.currentTarget.setPointerCapture(event.pointerId);
    drag.current = { name, pointerId: event.pointerId, offset: { x: point.x - resolved[name].x, y: point.y - resolved[name].y } };
    setDragBounds(contentBounds); // Prevent fit-to-content from moving the pointer while dragging.
  };
  const pointerMove = (event: PointerEvent<SVGGElement>) => {
    const active = drag.current;
    if (!active || active.pointerId !== event.pointerId) return;
    const point = svgPoint(event);
    if (point) moveUnit(active.name, { x: point.x - active.offset.x, y: point.y - active.offset.y });
  };
  const endDrag = (event: PointerEvent<SVGGElement>) => {
    if (drag.current?.pointerId !== event.pointerId) return;
    drag.current = null;
    setDragBounds(null);
    if (event.currentTarget.hasPointerCapture(event.pointerId)) event.currentTarget.releasePointerCapture(event.pointerId);
  };
  const keyDown = (event: KeyboardEvent<SVGGElement>, name: string) => {
    if (event.key === 'Enter' || event.key === ' ') { event.preventDefault(); onSelectUnit?.(name); }
    if (!onPositionsChange) return;
    const direction: Record<string, PublicationPoint> = { ArrowLeft: { x: -1, y: 0 }, ArrowRight: { x: 1, y: 0 }, ArrowUp: { x: 0, y: -1 }, ArrowDown: { x: 0, y: 1 } };
    const vector = direction[event.key];
    if (!vector) return;
    event.preventDefault();
    const step = event.shiftKey ? 24 : 8;
    moveUnit(name, { x: resolved[name].x + vector.x * step, y: resolved[name].y + vector.y * step });
  };

  return <svg ref={ref} xmlns="http://www.w3.org/2000/svg" className={className}
    viewBox={`${box.x} ${box.y} ${box.width} ${box.height}`} width="100%" height="100%"
    preserveAspectRatio="xMidYMid meet" role={interactive ? 'group' : 'img'} aria-labelledby={`${uid}-title ${uid}-description`}
    fontFamily="Arial, Helvetica, sans-serif" fill="#202C36" style={{ display: 'block', background: '#FFFFFF', minHeight: 320, touchAction: onPositionsChange ? 'none' : 'auto' }}>
    <title id={`${uid}-title`}>{title}</title>
    <desc id={`${uid}-description`}>{description}</desc>
    <rect x={box.x} y={box.y} width={box.width} height={box.height} fill="#FFFFFF" />
    {!npes.length ? <text x={box.x + box.width / 2} y={box.y + box.height / 2} textAnchor="middle" fill="#576873" fontSize="20">
      {language === 'es' ? 'Añade unidades para construir la figura de tu red.' : 'Add units to build your network figure.'}
    </text> : <>
      <DiffuseField nodes={sensory} positions={resolved} id={`${uid}-sensory-field`} blue={false} />
      <DiffuseField nodes={motor} positions={resolved} id={`${uid}-motor-field`} blue />
      <g fill="none" strokeWidth="2.3" strokeLinecap="round" strokeLinejoin="round">
        {connections.map((c, i) => {
          const g = geometry[i], source = npes.find(n => n.name === c.presynapticNPE);
          if (!g || !source) return null;
          const color = publicationNodeStyle(source).stroke;
          const norm = Math.hypot(g.tangent.x, g.tangent.y) || 1;
          const tx = -g.tangent.y / norm * 8, ty = g.tangent.x / norm * 8;
          return <g key={`${c.presynapticNPE}-${c.postsynapticNPE}-${i}`} data-connection={`${c.presynapticNPE}-${c.postsynapticNPE}`} stroke={color}>
            <title>{`${getDisplayName(c.presynapticNPE)} → ${getDisplayName(c.postsynapticNPE)}`}</title>
            <path d={g.path} />
            {source.type === 'Inhibitory'
              ? <line x1={g.end.x - tx} y1={g.end.y - ty} x2={g.end.x + tx} y2={g.end.y + ty} strokeWidth="3.2" />
              : <circle cx={g.end.x} cy={g.end.y} r="3.7" fill={color} stroke="none" />}
          </g>;
        })}
      </g>
      {showValues && <g fontSize="13" fill="#45535D" textAnchor="middle" pointerEvents="none">
        {connections.map((c, i) => {
          const g = geometry[i], value = weights?.[`${c.presynapticNPE}-${c.postsynapticNPE}`] ?? c.weight;
          return g && Number.isFinite(value) ? <text key={i} x={g.label.x} y={g.label.y} stroke="#FFFFFF" strokeWidth="5" strokeLinejoin="round" paintOrder="stroke fill">{value.toFixed(3)}</text> : null;
        })}
      </g>}
      {npes.map(npe => {
        const p = resolved[npe.name], radius = publicationNodeRadius(npe), style = publicationNodeStyle(npe);
        const selected = selectedUnit === npe.name;
        const activation = activations?.[npe.name] ?? npe.activation;
        return <g key={npe.name} data-unit={npe.name} transform={`translate(${p.x},${p.y})`}
          tabIndex={interactive ? 0 : undefined} role={interactive ? 'button' : undefined}
          aria-label={`${getDisplayName(npe.name)} · ${npe.layer}${showValues && Number.isFinite(activation) ? ` · ${activation.toFixed(3)}` : ''}`}
          onPointerDown={e => pointerDown(e, npe.name)} onPointerMove={pointerMove} onPointerUp={endDrag} onPointerCancel={endDrag} onLostPointerCapture={endDrag}
          onKeyDown={e => keyDown(e, npe.name)}
          style={{ cursor: onPositionsChange ? 'grab' : onSelectUnit ? 'pointer' : 'default', outlineColor: '#0072B2', outlineOffset: 7 }}>
          <title>{`${npe.name} · ${npe.layer} · ${npe.type}`}</title>
          {selected && <circle data-publication-ui="selection" r={radius + 9} fill="none" stroke="#0072B2" strokeWidth="1.6" />}
          <g fill={style.fill} stroke={style.stroke} strokeWidth="2.7">
            {style.shape === 'square' ? <rect x={-radius} y={-radius} width={radius * 2} height={radius * 2} rx="3" />
              : style.shape === 'diamond' ? <polygon points={`0,${-radius * 1.22} ${radius * 1.22},0 0,${radius * 1.22} ${-radius * 1.22},0`} />
                : style.shape === 'hexagon' ? <polygon points={Array.from({ length: 6 }, (_, i) => { const angle = -Math.PI / 2 + i * Math.PI / 3; return `${Math.cos(angle) * radius * 1.12},${Math.sin(angle) * radius * 1.12}`; }).join(' ')} />
                  : <circle r={radius} />}
          </g>
          <text x="0" y="1" dominantBaseline="middle" textAnchor="middle" fontSize="21" fontStyle="italic" fill="#202C36" pointerEvents="none">{getDisplayName(npe.name)}</text>
          {showValues && Number.isFinite(activation) && <text x="0" y={radius + 23} textAnchor="middle" fontSize="13" fill="#45535D" stroke="#FFFFFF" strokeWidth="4" paintOrder="stroke fill" pointerEvents="none">a = {activation.toFixed(3)}</text>}
        </g>;
      })}
      <g fontSize="15" fill="#576873" pointerEvents="none" transform={`translate(${box.x + 25},${box.y + box.height - 30})`}>
        <line x1="0" y1="-5" x2="35" y2="-5" stroke="#576873" strokeWidth="2.3" />
        <circle cx="35" cy="-5" r="3.7" fill="#576873" />
        <text x="46" y="0">{language === 'es' ? 'Excitatoria' : 'Excitatory'}</text>
        <line x1="175" y1="-5" x2="210" y2="-5" stroke="#8E679C" strokeWidth="2.3" />
        <line x1="210" y1="-13" x2="210" y2="3" stroke="#8E679C" strokeWidth="3.2" />
        <text x="222" y="0">{language === 'es' ? 'Inhibitoria' : 'Inhibitory'}</text>
        {box.width > 670 && <text x="390" y="0">{language === 'es' ? 'Campos: modulación del aprendizaje' : 'Fields: learning modulation'}</text>}
      </g>
    </>}
  </svg>;
}

export default PublicationNetwork;
