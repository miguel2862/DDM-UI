import type { Connection, NPE } from '../types/ddm';
import { getDisplayName } from './displayNames';

/** All positions are node CENTRES in SVG document coordinates, not React Flow corners. */
export interface PublicationPoint { x: number; y: number }
export type PublicationPositions = Record<string, PublicationPoint>;
export interface PublicationBounds { x: number; y: number; width: number; height: number }
export interface PublicationNodeStyle {
  fill: string;
  stroke: string;
  shape: 'square' | 'circle' | 'diamond' | 'hexagon';
}

// Exact unit palette of Figure 1 in the author's magnitude/delay manuscript.
export function publicationNodeStyle(npe: NPE): PublicationNodeStyle {
  if (npe.layer === 'US') return { fill: '#FAE2D0', stroke: '#C95C13', shape: 'hexagon' };
  if (npe.type === 'Inhibitory') return { fill: '#EEE2F1', stroke: '#8E679C', shape: 'diamond' };
  switch (npe.layer) {
    case 'PrimarySensory': return { fill: '#DFEDF6', stroke: '#0072B2', shape: 'square' };
    case 'AssociativeSensory': return { fill: '#F0E0ED', stroke: '#A86495', shape: 'circle' };
    case 'Hippocampal': return { fill: '#FBECCB', stroke: '#B88817', shape: 'circle' };
    case 'AssociativeMotor': return { fill: '#D6EEE4', stroke: '#008D79', shape: 'circle' };
    case 'Dopaminergic': return { fill: '#D8EDFB', stroke: '#368AB9', shape: 'circle' };
    default: return { fill: '#E5E8EA', stroke: '#45535D', shape: 'circle' };
  }
}

/** A long user-defined name remains intact; its unit grows instead of clipping it. */
export function publicationNodeRadius(npe: NPE): number {
  const glyphs = Array.from(getDisplayName(npe.name)).length;
  return Math.max(27, 13 + glyphs * 4.3);
}

/** Deterministic DDM layer layout. No connection is inferred, added or removed. */
export function computePublicationLayout(npes: NPE[], connections: Connection[] = []): PublicationPositions {
  const positions: PublicationPositions = {};
  if (!npes.length) return positions;
  const byLayer = (layer: string) => npes.filter(n => n.layer === layer && n.type !== 'Inhibitory');
  const primary = byLayer('PrimarySensory');
  const sensory = byLayer('AssociativeSensory');
  const motor = byLayer('AssociativeMotor');
  const output = byLayer('PrimaryMotor');
  const hippocampal = byLayer('Hippocampal');
  const dopaminergic = byLayer('Dopaminergic');
  const us = byLayer('US');
  const largestRadius = Math.max(...npes.map(publicationNodeRadius));
  const unitGap = Math.max(146, largestRadius * 2 + 55);
  const colGap = Math.max(230, largestRadius * 2 + 110);
  const count = Math.max(primary.length, sensory.length, motor.length, output.length, hippocampal.length + 1, 3);
  const span = (count - 1) * unitGap;
  const top = 112;
  const column = (nodes: NPE[], x: number, y0 = top, y1 = top + span) => {
    nodes.forEach((n, i) => { positions[n.name] = { x, y: nodes.length === 1 ? (y0 + y1) / 2 : y0 + i * (y1 - y0) / (nodes.length - 1) }; });
  };
  column(primary, 100);
  column(sensory, 100 + colGap);
  column(hippocampal, 100 + colGap * 1.48, top + span * .27, top + span * .73);
  // A single H cell sits below a single sensory–motor pathway, as in the
  // manuscript: keeping it on the pathway would force an avoidable crossing.
  if (hippocampal.length === 1 && sensory.length === 1 && motor.length === 1) {
    positions[hippocampal[0].name].y += unitGap * .65;
  }
  column(motor, 100 + colGap * 2.5);
  column(output, 100 + colGap * 3.65);
  us.forEach((n, i) => { positions[n.name] = { x: 100 + i * unitGap, y: top + span + unitGap }; });
  dopaminergic.forEach((n, i) => { positions[n.name] = { x: 100 + colGap * 2.5 + i * unitGap, y: top + span + unitGap * .92 }; });

  // Motor inhibitory interneurons sit between the actual outputs they connect.
  const interneurons = npes.filter(n => n.layer === 'AssociativeMotor' && n.type === 'Inhibitory');
  const occupied: PublicationPoint[] = Object.values(positions);
  interneurons.forEach((n, i) => {
    const neighbours = connections.flatMap(c => c.presynapticNPE === n.name ? [c.postsynapticNPE] : c.postsynapticNPE === n.name ? [c.presynapticNPE] : [])
      .filter(name => output.some(o => o.name === name)).map(name => positions[name]);
    const y = neighbours.length ? neighbours.reduce((sum, p) => sum + p.y, 0) / neighbours.length : top + span / 2;
    const p = { x: 100 + colGap * 3.65 + (i % 2 === 0 ? .34 : -.34) * colGap, y };
    while (occupied.some(q => Math.hypot(q.x - p.x, q.y - p.y) < largestRadius * 2 + 20)) p.y += unitGap * .55;
    positions[n.name] = p;
    occupied.push(p);
  });
  // Nonstandard layers and inhibitory cells in other layers remain visible.
  const otherLayers = [...new Set(npes.filter(n => !positions[n.name]).map(n => n.layer))];
  otherLayers.forEach((layer, i) => column(npes.filter(n => n.layer === layer && !positions[n.name]), 100 + colGap * (5 + i)));
  return positions;
}

export function publicationBounds(npes: NPE[], positions: PublicationPositions, extra: PublicationPoint[] = []): PublicationBounds {
  const bounds = npes.flatMap(n => {
    const p = positions[n.name];
    if (!p || !Number.isFinite(p.x) || !Number.isFinite(p.y)) return [];
    const radius = publicationNodeRadius(n) + 24;
    return [{ x: p.x - radius, y: p.y - radius }, { x: p.x + radius, y: p.y + radius + 20 }];
  }).concat(extra);
  if (!bounds.length) return { x: 0, y: 0, width: 850, height: 460 };
  const left = Math.min(...bounds.map(p => p.x)) - 45;
  const top = Math.min(...bounds.map(p => p.y)) - 64;
  return { x: left, y: top, width: Math.max(380, Math.max(...bounds.map(p => p.x)) - left + 45), height: Math.max(280, Math.max(...bounds.map(p => p.y)) - top + 100) };
}

function exportDimensions(svg: SVGSVGElement): PublicationBounds {
  const box = svg.viewBox.baseVal;
  if (![box.x, box.y, box.width, box.height].every(Number.isFinite) || box.width <= 0 || box.height <= 0) throw new Error('The network has no valid export dimensions.');
  return { x: box.x, y: box.y, width: box.width, height: box.height };
}

/** Serialize artwork only. Editor focus/selection overlays never reach the image. */
export function serializePublicationSvg(svg: SVGSVGElement): string {
  const box = exportDimensions(svg);
  const clone = svg.cloneNode(true) as SVGSVGElement;
  clone.querySelectorAll('[data-publication-ui]').forEach(node => node.remove());
  clone.querySelectorAll('[tabindex]').forEach(node => node.removeAttribute('tabindex'));
  clone.querySelectorAll('[role="button"]').forEach(node => node.removeAttribute('role'));
  clone.removeAttribute('class');
  clone.removeAttribute('style');
  clone.setAttribute('xmlns', 'http://www.w3.org/2000/svg');
  clone.setAttribute('width', String(Math.ceil(box.width)));
  clone.setAttribute('height', String(Math.ceil(box.height)));
  clone.setAttribute('viewBox', `${box.x} ${box.y} ${box.width} ${box.height}`);
  clone.setAttribute('role', 'img');
  if (clone.querySelector('foreignObject, image, script')) throw new Error('Publication export must contain self-contained vector artwork only.');
  return '<?xml version="1.0" encoding="UTF-8"?>\n' + new XMLSerializer().serializeToString(clone);
}

function saveBlob(blob: Blob, filename: string): void {
  const url = URL.createObjectURL(blob);
  const anchor = document.createElement('a');
  anchor.href = url;
  anchor.download = filename;
  document.body.append(anchor);
  anchor.click();
  anchor.remove();
  // WebKit needs the object URL to survive the click's navigation task.
  window.setTimeout(() => URL.revokeObjectURL(url), 1000);
}

export function downloadPublicationSvg(svg: SVGSVGElement, filename = 'DDM-network.svg'): void {
  saveBlob(new Blob([serializePublicationSvg(svg)], { type: 'image/svg+xml;charset=utf-8' }), filename);
}

/** Rasterize the native SVG at a real pixel multiplier, independent of screen size. */
export async function downloadPublicationPng(svg: SVGSVGElement, filename = 'DDM-network.png', scale = 4): Promise<void> {
  if (!Number.isFinite(scale) || scale <= 0) throw new Error('The export scale must be positive.');
  const box = exportDimensions(svg);
  const width = Math.ceil(box.width * scale), height = Math.ceil(box.height * scale);
  if (width > 16384 || height > 16384 || width * height > 64_000_000) throw new Error('This network is too large for PNG at this scale. Use SVG or a smaller export scale.');
  await document.fonts.ready;
  const source = URL.createObjectURL(new Blob([serializePublicationSvg(svg)], { type: 'image/svg+xml;charset=utf-8' }));
  try {
    const image = new Image();
    await new Promise<void>((resolve, reject) => {
      image.onload = () => resolve();
      image.onerror = () => reject(new Error('The network image could not be rendered. Please try SVG export.'));
      image.src = source;
    });
    const canvas = document.createElement('canvas');
    canvas.width = width;
    canvas.height = height;
    const context = canvas.getContext('2d');
    if (!context) throw new Error('PNG export is not supported in this browser.');
    context.fillStyle = '#FFFFFF';
    context.fillRect(0, 0, width, height);
    context.drawImage(image, 0, 0, width, height);
    const blob = await new Promise<Blob>((resolve, reject) => canvas.toBlob(result => result ? resolve(result) : reject(new Error('PNG encoding failed. Please use SVG export.')), 'image/png'));
    saveBlob(blob, filename);
  } finally {
    URL.revokeObjectURL(source);
  }
}
