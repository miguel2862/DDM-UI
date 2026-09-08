import { useState, useMemo, useRef } from 'react';
import {
  LineChart, Line, BarChart, Bar, XAxis, YAxis, CartesianGrid, Tooltip,
  Legend, ResponsiveContainer, ErrorBar, ReferenceLine, ComposedChart,
  Scatter,
} from 'recharts';
import { toPng } from 'html-to-image';
import {
  BarChart3, TrendingUp, Layers, User, Users, Download, Image, CheckSquare, Square,
} from 'lucide-react';
import { Card } from '../components/ui/Card';
import { Badge } from '../components/ui/Badge';
import { PageTransition } from '../components/ui/Skeleton';
import { useSimStore } from '../stores/useSimStore';
import { useToast } from '../components/ui/Toast';
import { downloadResultsCSV, downloadAllNetworksCSV } from '../utils/dataExport';
import { getDisplayName } from '../utils/displayNames';
import { useI18n } from '../i18n';
import type { SimulationResult } from '../types/ddm';

type ChartType = 'activations' | 'weights' | 'aggregate' | 'signals';

interface ChartDataPoint extends Record<string, string | number> {
  trial: number;
  phase: string;
  phaseTrial: string | number;
}

interface TooltipEntry {
  color?: string;
  name?: string | number;
  value?: string | number;
  payload?: ChartDataPoint;
}

interface CustomChartTooltipProps {
  active?: boolean;
  payload?: readonly TooltipEntry[];
}

interface ItalicUnitTickProps {
  x?: number;
  y?: number;
  payload?: { value?: string | number };
}

const EMPTY_STRINGS: string[] = [];
const DTD_SIGNAL_NAMES = ['dVTA', 'dH'];

const COLORS = [
  '#06b6d4', '#14b8a6', '#f59e0b', '#8b5cf6', '#ef4444',
  '#ec4899', '#10b981', '#3b82f6', '#f97316', '#a855f7',
];

const SCATTER_COLORS = [
  '#0e7490', '#0d9488', '#d97706', '#7c3aed', '#dc2626',
  '#db2777', '#059669', '#2563eb', '#ea580c', '#9333ea',
];

// Custom tooltip that shows Phase and Trial info
function CustomChartTooltip({ active, payload }: CustomChartTooltipProps) {
  if (!active || !payload || payload.length === 0) return null;
  const data = payload[0]?.payload;
  return (
    <div style={{
      backgroundColor: '#ffffff',
      border: '1px solid #e2e8f0',
      borderRadius: '12px',
      padding: '10px 14px',
      fontSize: '12px',
      color: '#1e293b',
      boxShadow: '0 4px 12px rgba(0,0,0,0.08)',
    }}>
      <p style={{ fontWeight: 700, marginBottom: 4, color: '#64748b' }}>
        {data?.phase} — Trial {data?.phaseTrial}
      </p>
      {payload.map((entry, i) => (
        <p key={i} style={{ color: entry.color, margin: '2px 0' }}>
          {entry.name}: {typeof entry.value === 'number' ? entry.value.toFixed(4) : entry.value}
        </p>
      ))}
    </div>
  );
}

// Custom XAxis tick that renders unit names in italic with proper display names
function ItalicUnitTick({ x = 0, y = 0, payload }: ItalicUnitTickProps) {
  const displayName = getDisplayName(String(payload?.value ?? ''));
  return (
    <text x={x} y={y} dy={14} textAnchor="middle" fill="#64748b" fontSize={11} fontStyle="italic">
      {displayName}
    </text>
  );
}

// Legend formatter that wraps variable names in italic
function italicLegendFormatter(value: string) {
  return <span style={{ fontStyle: 'italic' }}>{value}</span>;
}

function computeBoundaries(data: ChartDataPoint[]): { trial: number; phase: string }[] {
  const boundaries: { trial: number; phase: string }[] = [];
  let lastPhase = '';
  for (const point of data) {
    if (point.phase !== lastPhase && lastPhase !== '') {
      boundaries.push({ trial: point.trial, phase: point.phase });
    }
    lastPhase = point.phase;
  }
  return boundaries;
}

function getDefaultSelection(units: string[]): string[] {
  const motorUnits = units.filter(unit =>
    /^M\.\d+$/.test(unit) || /^(M'|CR|R)\d*/.test(unit) || /^Motor_/.test(unit)
  );
  return motorUnits.length > 0 ? motorUnits : units.slice(0, 1);
}

function getDefaultSignalSelection(signalNames: string[]): string[] {
  const preferred = (['dVTA', 'dH']);
  const initial = preferred.filter(name => signalNames.includes(name));
  return (initial.length > 0 ? initial : signalNames).slice(0, 5);
}

function resolvePhaseTimestep(
  phase: string,
  overrides: Record<string, number>,
  availableTimesteps: Record<string, number[]>,
): number {
  if (overrides[phase] !== undefined) return overrides[phase];
  const available = availableTimesteps[phase];
  if (available && available.length > 1) return available[available.length - 2];
  if (available && available.length === 1) return available[0];
  return 1;
}

function filterByPhaseTimestep(
  data: SimulationResult[],
  overrides: Record<string, number>,
  availableTimesteps: Record<string, number[]>,
): SimulationResult[] {
  return data.filter(row => {
    const phase = String(row.Phase);
    return Number(row.TimeStep) === resolvePhaseTimestep(phase, overrides, availableTimesteps);
  });
}

export function Results() {
  const { simulationResults, simulationMetadata, selectedNetwork, setSelectedNetwork } = useSimStore();
  const { t } = useI18n();
  const toast = useToast();

  const [tab, setTab] = useState<'individual' | 'general'>('individual');
  const [chartType, setChartType] = useState<ChartType>('activations');
  const [unitSelection, setUnitSelection] = useState<string[] | null>(null);
  const [connectionSelection, setConnectionSelection] = useState<string[] | null>(null);
  const [signalSelection, setSignalSelection] = useState<string[] | null>(null);
  const [phaseSelection, setPhaseSelection] = useState<string | null>(null);
  // Per-phase timestep selection: { "training": 5, "extinction": 5 }
  const [phaseTimestepOverrides, setPhaseTimestepOverrides] = useState<Record<string, number>>({});
  const [aggregateMeasure, setAggregateMeasure] = useState<'mean' | 'median'>('mean');
  const [figureNumber, setFigureNumber] = useState(1);

  const chartRef = useRef<HTMLDivElement>(null);
  const exportWrapperRef = useRef<HTMLDivElement>(null);

  const units = simulationMetadata?.units ?? EMPTY_STRINGS;
  const connNames = simulationMetadata?.connections ?? EMPTY_STRINGS;
  const phases = simulationMetadata?.phases ?? EMPTY_STRINGS;
  const signalNames = simulationMetadata?.signals ?? DTD_SIGNAL_NAMES;
  const effectiveSelectedNetwork = simulationResults && selectedNetwork < simulationResults.length
    ? selectedNetwork
    : 0;

  // Compute available timesteps PER PHASE from actual simulation data
  const phaseAvailableTimesteps = useMemo(() => {
    const result: Record<string, number[]> = {};
    if (!simulationResults || simulationResults.length === 0) return result;
    const firstNet = simulationResults[0];
    if (!firstNet || firstNet.length === 0) return result;

    firstNet.forEach(row => {
      const phase = String(row.Phase);
      const ts = Number(row.TimeStep);
      if (!result[phase]) result[phase] = [];
      if (!result[phase].includes(ts)) result[phase].push(ts);
    });

    // Sort each phase's timesteps
    for (const phase of Object.keys(result)) {
      result[phase].sort((a, b) => a - b);
    }
    return result;
  }, [simulationResults]);

  // Defaults are derived from metadata. A non-null selection, including [], is an
  // explicit user choice and therefore remains stable after "deselect all".
  const defaultSelectedUnits = useMemo(() => getDefaultSelection(units), [units]);
  const defaultSelectedConns = useMemo(() => connNames.slice(0, 3), [connNames]);
  const defaultSelectedSignals = useMemo(
    () => getDefaultSignalSelection(signalNames),
    [signalNames],
  );
  const selectedUnits = unitSelection ?? defaultSelectedUnits;
  const selectedConns = connectionSelection ?? defaultSelectedConns;
  const selectedSignals = signalSelection ?? defaultSelectedSignals;
  const selectedPhase = phaseSelection ?? phases[0] ?? '';

  const setSelectedUnits = (next: string[] | ((previous: string[]) => string[])) => {
    setUnitSelection(previous => typeof next === 'function'
      ? next(previous ?? defaultSelectedUnits)
      : next);
  };
  const setSelectedConns = (next: string[] | ((previous: string[]) => string[])) => {
    setConnectionSelection(previous => typeof next === 'function'
      ? next(previous ?? defaultSelectedConns)
      : next);
  };
  const setSelectedSignals = (next: string[] | ((previous: string[]) => string[])) => {
    setSignalSelection(previous => typeof next === 'function'
      ? next(previous ?? defaultSelectedSignals)
      : next);
  };

  // Helper: get the timestep for a given phase
  const getPhaseTimestep = (phase: string): number => {
    return resolvePhaseTimestep(phase, phaseTimestepOverrides, phaseAvailableTimesteps);
  };

  // Update a single phase's timestep
  const setPhaseTimestep = (phase: string, ts: number) => {
    setPhaseTimestepOverrides(previous => ({ ...previous, [phase]: ts }));
  };

  // Export chart as APA 7 compliant PNG (publication-ready)
  // Strategy: inject APA header/footer into the live exportWrapper (parent of chart),
  // lock wrapper to fixed pixel width so toPng captures everything, then clean up.
  const handleExportChartPNG = async () => {
    if (!exportWrapperRef.current || !chartRef.current) return;
    try {
      // Build localized figure title
      const measureLabel = aggregateMeasure === 'mean' ? t.results.mean : t.results.median;
      const chartTitles: Record<string, string> = {
        activations: tab === 'general'
          ? t.results.figActivationsGeneral.replace('{measure}', measureLabel)
          : t.results.figActivationsIndividual,
        weights: tab === 'general'
          ? t.results.figWeightsGeneral.replace('{measure}', measureLabel)
          : t.results.figWeightsIndividual,
        signals: tab === 'general'
          ? t.results.figSignalsGeneral.replace('{measure}', measureLabel)
          : t.results.figSignalsIndividual,
        aggregate: t.results.figAggregate
          .replace('{measure}', measureLabel)
          .replace('{phase}', selectedPhase),
      };
      const figureTitle = chartTitles[chartType] || 'Chart';
      const currentFigNum = figureNumber;
      const wrapper = exportWrapperRef.current;

      // ── Measure current chart dimensions BEFORE any DOM changes ──
      const wrapperRect = wrapper.getBoundingClientRect();
      const captureWidth = Math.ceil(wrapperRect.width);

      // ── APA 7 header: "Figure X" bold + title italic ──
      const header = document.createElement('div');
      header.setAttribute('data-apa-temp', 'true');
      header.style.cssText = 'padding: 0 0 8px 0;';
      header.innerHTML = `
        <p style="font-size: 14px; font-weight: bold; margin: 0 0 2px 0; font-family: 'Times New Roman', serif; color: #1e293b;">Figure ${currentFigNum}</p>
        <p style="font-size: 13px; font-style: italic; margin: 0; font-family: 'Times New Roman', serif; color: #1e293b;">${figureTitle}</p>
      `;

      // ── APA 7 note (localized) ──
      const timestepInfo = phases.map(p => `${p}: t=${getPhaseTimestep(p)}`).join(', ');
      const networkInfo = tab === 'individual'
        ? t.results.figNetworkIndividual
            .replace('{current}', String(effectiveSelectedNetwork + 1))
            .replace('{total}', String(simulationResults?.length || 1))
        : t.results.figNetworkGeneral
            .replace('{total}', String(simulationResults?.length || 1))
            .replace('{measure}', measureLabel);
      const noteText = t.results.figNote
        .replace('{networkInfo}', networkInfo)
        .replace('{timestepInfo}', timestepInfo);
      const note = document.createElement('div');
      note.setAttribute('data-apa-temp', 'true');
      note.style.cssText = 'padding: 8px 0 0 0; font-size: 11px; color: #475569; font-family: "Times New Roman", serif;';
      note.innerHTML = `<p style="margin: 0;"><span style="font-style: italic;">Note.</span> ${noteText}</p>`;

      // ── Inject header BEFORE chart, note AFTER ──
      wrapper.insertBefore(header, chartRef.current);
      wrapper.appendChild(note);

      // ── Lock wrapper to fixed pixel width + white bg + padding ──
      const origStyles = {
        width: wrapper.style.width,
        padding: wrapper.style.padding,
        background: wrapper.style.background,
        boxSizing: wrapper.style.boxSizing,
      };
      wrapper.style.width = `${captureWidth}px`;
      wrapper.style.padding = '16px 16px 16px 12px';
      wrapper.style.background = '#ffffff';
      wrapper.style.boxSizing = 'content-box';

      // ── Bump SVG text sizes for print readability ──
      // (Italics are handled natively by Recharts via Legend formatter + custom XAxis tick)
      const svgTexts = chartRef.current.querySelectorAll('svg text');
      const origFonts: { el: SVGTextElement; size: string; family: string; fontSize: string }[] = [];
      svgTexts.forEach(textEl => {
        const el = textEl as SVGTextElement;
        origFonts.push({
          el,
          size: el.getAttribute('font-size') || '',
          family: el.style.fontFamily,
          fontSize: el.style.fontSize,
        });
        el.style.fontFamily = 'Arial, Helvetica, sans-serif';
        const currentSize = parseFloat(el.getAttribute('font-size') || el.style.fontSize || '11');
        if (currentSize < 12) {
          el.setAttribute('font-size', '12');
          el.style.fontSize = '12px';
        }
      });

      // Wait for layout to settle
      await new Promise(r => requestAnimationFrame(() => requestAnimationFrame(r)));

      // ── Capture with explicit width to prevent clipping ──
      const totalWidth = wrapper.scrollWidth;
      const totalHeight = wrapper.scrollHeight;
      const dataUrl = await toPng(wrapper, {
        backgroundColor: '#ffffff',
        quality: 1.0,
        pixelRatio: 3,
        width: totalWidth,
        height: totalHeight,
      });

      // ── Cleanup: remove temp elements, restore styles ──
      wrapper.removeChild(header);
      wrapper.removeChild(note);
      wrapper.style.width = origStyles.width;
      wrapper.style.padding = origStyles.padding;
      wrapper.style.background = origStyles.background;
      wrapper.style.boxSizing = origStyles.boxSizing;
      origFonts.forEach(({ el, size, family, fontSize }) => {
        if (size) el.setAttribute('font-size', size);
        else el.removeAttribute('font-size');
        el.style.fontFamily = family;
        el.style.fontSize = fontSize;
      });

      // Auto-increment figure number for next export
      setFigureNumber(prev => prev + 1);

      const link = document.createElement('a');
      link.download = `ddm-figure${currentFigNum}-${chartType}.png`;
      link.href = dataUrl;
      link.click();
    } catch (err) {
      console.error('Failed to export chart PNG:', err);
      toast.error(t.toast.exportPNGError || 'Failed to export chart as PNG. Please try again.');
    }
  };

  // Get the single-network data array to plot (for line charts)
  const plotData = useMemo((): SimulationResult[] | null => {
    if (!simulationResults || simulationResults.length === 0) return null;
    if (tab === 'individual') {
      return simulationResults[effectiveSelectedNetwork] || null;
    }
    // For general: compute mean/median across all networks
    const firstNet = simulationResults[0];
    if (!firstNet || firstNet.length === 0) return null;

    return firstNet.map((refRow, rowIdx): SimulationResult => {
      const result: SimulationResult = {
        Phase: refRow.Phase,
        Trial: refRow.Trial,
        TimeStep: refRow.TimeStep,
      };
      const dataKeys = Object.keys(refRow).filter(k => k !== 'Phase' && k !== 'Trial' && k !== 'TimeStep');
      for (const key of dataKeys) {
        const values = simulationResults
          .map(net => net?.[rowIdx]?.[key])
          .filter((v): v is number => typeof v === 'number');
        if (values.length === 0) {
          result[key] = 0;
        } else if (aggregateMeasure === 'median') {
          const sorted = [...values].sort((a, b) => a - b);
          result[key] = sorted[Math.floor(sorted.length / 2)];
        } else {
          result[key] = values.reduce((a, b) => a + b, 0) / values.length;
        }
      }
      return result;
    });
  }, [simulationResults, tab, effectiveSelectedNetwork, aggregateMeasure]);

  // Activation line chart data -- use a global trial index so phase boundaries are visible
  const activationData = useMemo(() => {
    if (!plotData) return [];
    const filtered = filterByPhaseTimestep(plotData, phaseTimestepOverrides, phaseAvailableTimesteps);
    return filtered.map((row, idx): ChartDataPoint => {
      const point: ChartDataPoint = {
        trial: idx + 1,
        phase: row.Phase,
        phaseTrial: row.Trial,
      };
      selectedUnits.forEach(unit => {
        point[unit] = typeof row[unit] === 'number' ? Number(row[unit].toFixed(4)) : 0;
      });
      return point;
    });
  }, [plotData, phaseTimestepOverrides, phaseAvailableTimesteps, selectedUnits]);

  // Weight line chart data -- use a global trial index
  const weightData = useMemo(() => {
    if (!plotData) return [];
    const filtered = filterByPhaseTimestep(plotData, phaseTimestepOverrides, phaseAvailableTimesteps);
    return filtered.map((row, idx): ChartDataPoint => {
      const point: ChartDataPoint = {
        trial: idx + 1,
        phase: row.Phase,
        phaseTrial: row.Trial,
      };
      selectedConns.forEach(conn => {
        point[conn] = typeof row[conn] === 'number' ? Number(row[conn].toFixed(4)) : 0;
      });
      return point;
    });
  }, [plotData, phaseTimestepOverrides, phaseAvailableTimesteps, selectedConns]);

  // Learning signals data
  const signalData = useMemo(() => {
    if (!plotData) return [];
    const filtered = filterByPhaseTimestep(plotData, phaseTimestepOverrides, phaseAvailableTimesteps);
    return filtered.map((row, idx): ChartDataPoint => {
      const point: ChartDataPoint = {
        trial: idx + 1,
        phase: row.Phase,
        phaseTrial: row.Trial,
      };
      selectedSignals.forEach(signal => {
        point[signal] = typeof row[signal] === 'number' ? Number(row[signal].toFixed(6)) : 0;
      });
      return point;
    });
  }, [plotData, phaseTimestepOverrides, phaseAvailableTimesteps, selectedSignals]);

  // Aggregate bar chart data -- single network (used for individual tab)
  const aggregateDataSingle = useMemo(() => {
    if (!plotData || !selectedPhase) return [];
    const selectedTs = resolvePhaseTimestep(
      selectedPhase,
      phaseTimestepOverrides,
      phaseAvailableTimesteps,
    );
    const filtered = plotData.filter(
      row => row.Phase === selectedPhase && Number(row.TimeStep) === selectedTs
    );
    if (filtered.length === 0) return [];

    return selectedUnits.map(unit => {
      const values = filtered.map(row => typeof row[unit] === 'number' ? row[unit] as number : 0);
      const mean = values.reduce((a: number, b: number) => a + b, 0) / values.length;
      const sorted = [...values].sort((a, b) => a - b);
      const median = sorted[Math.floor(sorted.length / 2)];
      const sd = Math.sqrt(values.reduce((acc: number, v: number) => acc + (v - mean) ** 2, 0) / values.length);
      const se = sd / Math.sqrt(values.length);

      return {
        unit,
        value: aggregateMeasure === 'mean' ? Number(mean.toFixed(4)) : Number(median.toFixed(4)),
        error: Number(se.toFixed(4)),
        sd: Number(sd.toFixed(4)),
      };
    });
  }, [plotData, selectedPhase, phaseTimestepOverrides, phaseAvailableTimesteps, selectedUnits, aggregateMeasure]);

  // Aggregate bar chart data -- general tab: across ALL networks
  const aggregateDataGeneral = useMemo(() => {
    if (!simulationResults || simulationResults.length === 0 || !selectedPhase) return [];
    const selectedTs = resolvePhaseTimestep(
      selectedPhase,
      phaseTimestepOverrides,
      phaseAvailableTimesteps,
    );

    const networkAggregates: { networkIdx: number; values: Record<string, number> }[] = [];

    for (let netIdx = 0; netIdx < simulationResults.length; netIdx++) {
      const netData = simulationResults[netIdx];
      if (!netData) continue;

      const filtered = netData.filter(
        row => row.Phase === selectedPhase && Number(row.TimeStep) === selectedTs
      );
      if (filtered.length === 0) continue;

      const netValues: Record<string, number> = {};
      for (const unit of selectedUnits) {
        const vals = filtered.map(row => typeof row[unit] === 'number' ? row[unit] as number : 0);
        const mean = vals.reduce((a: number, b: number) => a + b, 0) / vals.length;
        const sorted = [...vals].sort((a, b) => a - b);
        const median = sorted[Math.floor(sorted.length / 2)];
        netValues[unit] = aggregateMeasure === 'mean' ? mean : median;
      }

      networkAggregates.push({ networkIdx: netIdx, values: netValues });
    }

    if (networkAggregates.length === 0) return [];

    return selectedUnits.map(unit => {
      const perNetworkValues = networkAggregates.map(na => na.values[unit] || 0);
      const overallMean = perNetworkValues.reduce((a, b) => a + b, 0) / perNetworkValues.length;
      const sorted = [...perNetworkValues].sort((a, b) => a - b);
      const overallMedian = sorted[Math.floor(sorted.length / 2)];
      const sd = Math.sqrt(
        perNetworkValues.reduce((acc, v) => acc + (v - overallMean) ** 2, 0) / perNetworkValues.length
      );
      const se = sd / Math.sqrt(perNetworkValues.length);

      return {
        unit,
        value: aggregateMeasure === 'mean' ? Number(overallMean.toFixed(4)) : Number(overallMedian.toFixed(4)),
        error: Number(se.toFixed(4)),
        sd: Number(sd.toFixed(4)),
        scatterPoints: perNetworkValues.map((v, i) => ({
          value: Number(v.toFixed(4)),
          networkIdx: i,
        })),
        ...Object.fromEntries(perNetworkValues.map((v, i) => [`net${i}`, Number(v.toFixed(4))])),
      };
    });
  }, [simulationResults, selectedPhase, phaseTimestepOverrides, phaseAvailableTimesteps, selectedUnits, aggregateMeasure]);

  // Pick the right aggregate data based on the tab
  const aggregateData = tab === 'general' ? aggregateDataGeneral : aggregateDataSingle;

  // Compute phase boundaries independently for activations and weights
  const activationBoundaries = useMemo(() => computeBoundaries(activationData), [activationData]);
  const weightBoundaries = useMemo(() => computeBoundaries(weightData), [weightData]);
  const signalBoundaries = useMemo(() => computeBoundaries(signalData), [signalData]);

  // Select the appropriate boundaries based on current chart type
  const phaseBoundaries = chartType === 'weights' ? weightBoundaries : chartType === 'signals' ? signalBoundaries : activationBoundaries;

  if (!simulationResults || simulationResults.length === 0) {
    return (
      <PageTransition>
        <div className="max-w-4xl mx-auto flex flex-col items-center justify-center h-[60vh] space-y-4">
          <div className="w-20 h-20 rounded-2xl bg-slate-50 border border-slate-200 flex items-center justify-center">
            <BarChart3 size={36} className="text-slate-400" />
          </div>
          <h2 className="text-xl font-bold text-slate-500">{t.results.noResults}</h2>
          <p className="text-sm text-slate-400 text-center max-w-md">
            {t.results.noResultsDescription}
          </p>
        </div>
      </PageTransition>
    );
  }

  return (
    <PageTransition>
    <div className="max-w-7xl mx-auto space-y-6">
      <div className="flex items-center gap-3">
        <div className="w-10 h-10 rounded-xl bg-violet-50 flex items-center justify-center">
          <BarChart3 size={20} className="text-violet-600" />
        </div>
        <div>
          <h1 className="text-2xl font-bold text-slate-800">{t.results.pageTitle}</h1>
          <p className="text-sm text-slate-500">
            {simulationMetadata?.numNetworks} {t.sim.networksSuffix}, {simulationMetadata?.totalTimesteps} {t.results.timestepsLabel}
          </p>
        </div>
      </div>

      {/* Tab switch */}
      <div className="flex gap-1 p-1 rounded-xl bg-white border border-slate-200 max-w-md">
        <button
          onClick={() => setTab('individual')}
          className={`flex-1 flex items-center justify-center gap-2 py-2.5 rounded-lg text-sm font-bold transition-all ${
            tab === 'individual' ? 'bg-violet-50 text-violet-600' : 'text-slate-500 hover:text-slate-400'
          }`}
        >
          <User size={16} /> {t.results.individual}
        </button>
        <button
          onClick={() => setTab('general')}
          className={`flex-1 flex items-center justify-center gap-2 py-2.5 rounded-lg text-sm font-bold transition-all ${
            tab === 'general' ? 'bg-violet-50 text-violet-600' : 'text-slate-500 hover:text-slate-400'
          }`}
        >
          <Users size={16} /> {t.results.general}
        </button>
      </div>

      <div className="grid grid-cols-1 lg:grid-cols-4 gap-6">
        {/* Filters sidebar */}
        <div className="space-y-4">
          {tab === 'individual' && simulationResults.length > 1 && (
            <Card className="p-4">
              <label className="block text-xs font-bold text-slate-500 mb-2">{t.results.networkLabel}</label>
              <select
                value={effectiveSelectedNetwork}
                onChange={(e) => setSelectedNetwork(parseInt(e.target.value))}
                className="w-full px-3 py-2 rounded-lg bg-white border border-slate-200 text-slate-800 text-sm focus:border-cyan-500/50 focus:outline-none"
              >
                {simulationResults.map((_, i) => (
                  <option key={i} value={i}>{t.results.networkLabel} {i + 1}</option>
                ))}
              </select>
            </Card>
          )}

          <Card className="p-4 space-y-3">
            <label className="block text-xs font-bold text-slate-500">{t.results.chartType}</label>
            {([
              { key: 'activations', label: t.results.activations, icon: TrendingUp },
              { key: 'weights', label: t.results.weights, icon: Layers },
              { key: 'aggregate', label: t.results.aggregate, icon: BarChart3 },
              { key: 'signals', label: t.results.learningSignals, icon: TrendingUp },
            ] as const).map(({ key, label, icon: Icon }) => (
              <button
                key={key}
                onClick={() => setChartType(key)}
                className={`w-full flex items-center gap-2 px-3 py-2 rounded-lg text-sm font-bold transition-all ${
                  chartType === key
                    ? 'bg-cyan-50 text-cyan-600 border border-cyan-200'
                    : 'text-slate-500 border border-transparent hover:bg-slate-100'
                }`}
              >
                <Icon size={14} /> {label}
              </button>
            ))}
          </Card>

          {/* Per-Phase Timestep Selectors */}
          <Card className="p-4 space-y-3">
            <label className="block text-xs font-bold text-slate-500">{t.results.timestepPerPhase}</label>
            {phases.map(phase => {
              const avail = phaseAvailableTimesteps[phase] || [1];
              const currentTs = getPhaseTimestep(phase);
              const maxTs = avail[avail.length - 1] || 1;
              return (
                <div key={phase} className="space-y-1">
                  <div className="flex items-center justify-between">
                    <span className="text-[11px] font-bold text-slate-600">{phase}</span>
                    <Badge variant="info">t = {currentTs}</Badge>
                  </div>
                  {avail.length <= 12 ? (
                    <div className="flex flex-wrap gap-1">
                      {avail.map(ts => (
                        <button
                          key={ts}
                          onClick={() => setPhaseTimestep(phase, ts)}
                          className={`px-1.5 py-0.5 rounded text-[10px] font-bold transition-all ${
                            currentTs === ts
                              ? 'bg-cyan-100 text-cyan-600 border border-cyan-200'
                              : 'text-slate-400 border border-slate-200 hover:bg-slate-100'
                          }`}
                        >
                          {ts}
                        </button>
                      ))}
                    </div>
                  ) : (
                    <input
                      type="range"
                      min={avail[0]}
                      max={maxTs}
                      value={currentTs}
                      onChange={(e) => setPhaseTimestep(phase, parseInt(e.target.value))}
                      className="w-full accent-cyan-500"
                    />
                  )}
                </div>
              );
            })}
          </Card>

          {(chartType === 'activations' || chartType === 'aggregate') && (
            <Card className="p-4 space-y-2">
              <div className="flex items-center justify-between mb-1">
                <label className="block text-xs font-bold text-slate-500">{t.results.units}</label>
                <div className="flex gap-1">
                  <button
                    onClick={() => setSelectedUnits([...units])}
                    className="p-1 rounded text-slate-400 hover:text-cyan-600 transition-colors"
                    title={t.results.selectAll}
                  >
                    <CheckSquare size={13} />
                  </button>
                  <button
                    onClick={() => setSelectedUnits([])}
                    className="p-1 rounded text-slate-400 hover:text-slate-600 transition-colors"
                    title={t.results.deselectAll}
                  >
                    <Square size={13} />
                  </button>
                </div>
              </div>
              <div className="flex flex-wrap gap-1.5">
                {units.map(unit => (
                  <button
                    key={unit}
                    onClick={() => setSelectedUnits(prev =>
                      prev.includes(unit) ? prev.filter(u => u !== unit) : [...prev, unit]
                    )}
                    className={`px-2 py-1 rounded-md text-[11px] font-bold transition-all ${
                      selectedUnits.includes(unit)
                        ? 'bg-cyan-100 text-cyan-600 border border-cyan-200'
                        : 'text-slate-400 border border-slate-200 hover:bg-slate-100'
                    }`}
                  >
                    {unit}
                  </button>
                ))}
              </div>
            </Card>
          )}

          {chartType === 'weights' && (
            <Card className="p-4 space-y-2">
              <div className="flex items-center justify-between mb-1">
                <label className="block text-xs font-bold text-slate-500">{t.results.connectionsLabel}</label>
                <div className="flex gap-1">
                  <button
                    onClick={() => setSelectedConns([...connNames])}
                    className="p-1 rounded text-slate-400 hover:text-teal-600 transition-colors"
                    title={t.results.selectAll}
                  >
                    <CheckSquare size={13} />
                  </button>
                  <button
                    onClick={() => setSelectedConns([])}
                    className="p-1 rounded text-slate-400 hover:text-slate-600 transition-colors"
                    title={t.results.deselectAll}
                  >
                    <Square size={13} />
                  </button>
                </div>
              </div>
              <div className="flex flex-wrap gap-1.5">
                {connNames.map(conn => (
                  <button
                    key={conn}
                    onClick={() => setSelectedConns(prev =>
                      prev.includes(conn) ? prev.filter(c => c !== conn) : [...prev, conn]
                    )}
                    className={`px-2 py-1 rounded-md text-[11px] font-bold transition-all ${
                      selectedConns.includes(conn)
                        ? 'bg-teal-50 text-teal-600 border border-teal-200'
                        : 'text-slate-400 border border-slate-200 hover:bg-slate-100'
                    }`}
                  >
                    {conn}
                  </button>
                ))}
              </div>
            </Card>
          )}

          {chartType === 'signals' && (
            <Card className="p-4 space-y-2">
              <div className="flex items-center justify-between mb-1">
                <label className="block text-xs font-bold text-slate-500">{t.results.learningSignals}</label>
                <button
                  onClick={() => setSelectedSignals(signalNames.slice(0, 6))}
                  className="p-1 rounded text-slate-400 hover:text-violet-600 transition-colors"
                  title={t.results.selectAll}
                >
                  <CheckSquare size={13} />
                </button>
              </div>
              <div className="flex flex-wrap gap-1.5">
                {signalNames.map(signal => (
                  <button
                    key={signal}
                    onClick={() => setSelectedSignals(previous =>
                      previous.includes(signal)
                        ? previous.filter(item => item !== signal)
                        : [...previous, signal].slice(-6)
                    )}
                    className={`px-2 py-1 rounded-md text-[11px] font-bold transition-all ${
                      selectedSignals.includes(signal)
                        ? 'bg-violet-50 text-violet-600 border border-violet-200'
                        : 'text-slate-400 border border-slate-200 hover:bg-slate-100'
                    }`}
                  >
                    {getDisplayName(signal)}
                  </button>
                ))}
              </div>
            </Card>
          )}

          {chartType === 'aggregate' && (
            <Card className="p-4 space-y-3">
              <label className="block text-xs font-bold text-slate-500">{t.results.phase}</label>
              <select
                value={selectedPhase}
                onChange={(e) => setPhaseSelection(e.target.value)}
                className="w-full px-3 py-2 rounded-lg bg-white border border-slate-200 text-slate-800 text-sm focus:border-cyan-500/50 focus:outline-none"
              >
                {phases.map(p => <option key={p} value={p}>{p}</option>)}
              </select>
            </Card>
          )}

          {(chartType === 'aggregate' || tab === 'general') && (
            <Card className="p-4 space-y-3">
              <label className="block text-xs font-bold text-slate-500">{t.results.measure}</label>
              <div className="flex gap-2">
                <button
                  onClick={() => setAggregateMeasure('mean')}
                  className={`flex-1 py-1.5 rounded-lg text-xs font-bold ${
                    aggregateMeasure === 'mean' ? 'bg-cyan-50 text-cyan-600' : 'text-slate-400'
                  }`}
                >
                  {t.results.mean}
                </button>
                <button
                  onClick={() => setAggregateMeasure('median')}
                  className={`flex-1 py-1.5 rounded-lg text-xs font-bold ${
                    aggregateMeasure === 'median' ? 'bg-cyan-50 text-cyan-600' : 'text-slate-400'
                  }`}
                >
                  {t.results.median}
                </button>
              </div>
            </Card>
          )}
        </div>

        {/* Chart area */}
        <div className="lg:col-span-3">
          <Card className="p-6">
            <div className="mb-6 space-y-3">
              <div className="flex items-center justify-between">
                <h3 className="text-lg font-bold text-slate-800">
                  {chartType === 'activations' && (tab === 'general'
                    ? `${t.results.unitActivationsMeasure} (${aggregateMeasure === 'mean' ? t.results.mean : t.results.median}) — ${t.results.allNetworks}`
                    : t.results.unitActivationsOverTrials)}
                  {chartType === 'weights' && (tab === 'general'
                    ? `${t.results.connectionWeightsMeasure} (${aggregateMeasure === 'mean' ? t.results.mean : t.results.median}) — ${t.results.allNetworks}`
                    : t.results.connectionWeightsOverTrials)}
                  {chartType === 'signals' && (tab === 'general'
                    ? `${t.results.learningSignalsMeasure} (${aggregateMeasure === 'mean' ? t.results.mean : t.results.median}) — ${t.results.allNetworks}`
                    : t.results.learningSignalsOverTrials)}
                  {chartType === 'aggregate' && t.results.figAggregate
                    .replace('{measure}', aggregateMeasure === 'mean' ? t.results.mean : t.results.median)
                    .replace('{phase}', selectedPhase)}
                </h3>
                <div className="flex items-center gap-1">
                  {phases.map(p => (
                    <Badge key={p} variant="info">{p}: t={getPhaseTimestep(p)}</Badge>
                  ))}
                </div>
              </div>
              <div className="flex items-center gap-2">
                <div className="flex items-center gap-1">
                  <button
                    onClick={handleExportChartPNG}
                    className="flex items-center gap-1.5 px-3 py-1.5 rounded-lg text-xs font-bold text-violet-600 border border-violet-200 bg-violet-50 hover:bg-violet-100 transition-colors"
                    title={t.results.exportPNG}
                  >
                    <Image size={12} /> PNG
                  </button>
                  <span className="text-[10px] text-slate-400">Fig.</span>
                  <input
                    type="number"
                    min="1"
                    value={figureNumber}
                    onChange={(e) => setFigureNumber(Math.max(1, parseInt(e.target.value) || 1))}
                    className="w-10 px-1 py-0.5 rounded text-[11px] text-center font-bold text-violet-600 bg-violet-50 border border-violet-200 focus:outline-none focus:border-violet-400"
                  />
                </div>
                <div className="w-px h-5 bg-slate-200" />
                <button
                  onClick={() => {
                    if (plotData) {
                      downloadResultsCSV(plotData, `ddm-results-network${effectiveSelectedNetwork + 1}.csv`);
                    }
                  }}
                  className="flex items-center gap-1.5 px-3 py-1.5 rounded-lg text-xs font-bold text-slate-500 border border-slate-200 hover:bg-slate-100 transition-colors"
                  title={t.results.downloadCurrent}
                >
                  <Download size={12} /> {t.results.thisNetwork}
                </button>
                {simulationResults && simulationResults.length > 1 && (
                  <button
                    onClick={() => {
                      if (simulationResults) {
                        downloadAllNetworksCSV(simulationResults, 'ddm-results-all-networks.csv');
                      }
                    }}
                    className="flex items-center gap-1.5 px-3 py-1.5 rounded-lg text-xs font-bold text-cyan-600 border border-cyan-200 bg-cyan-50 hover:bg-cyan-100 transition-colors"
                    title={t.results.downloadAll}
                  >
                    <Download size={12} /> {t.results.allNetworksExport}
                  </button>
                )}
                <button
                  onClick={() => {
                    if (plotData) {
                      const cols = chartType === 'weights' ? selectedConns :
                                   chartType === 'signals' ? selectedSignals : selectedUnits;
                      const filteredData = plotData.map(row => {
                        const filtered: SimulationResult = {
                          Phase: row.Phase,
                          Trial: row.Trial,
                          TimeStep: row.TimeStep,
                        };
                        cols.forEach(col => { filtered[col] = row[col]; });
                        return filtered;
                      });
                      downloadResultsCSV(filteredData, `ddm-selective-export.csv`);
                    }
                  }}
                  className="flex items-center gap-1.5 px-3 py-1.5 rounded-lg text-xs font-bold text-amber-600 border border-amber-200 bg-amber-50 hover:bg-amber-100 transition-colors"
                  title={t.results.downloadSelected}
                >
                  <Download size={12} /> {t.results.selected}
                </button>
              </div>
            </div>

            <div ref={exportWrapperRef}>
            <div ref={chartRef} className="h-[450px]">
              <ResponsiveContainer
                width="100%"
                height="100%"
                minWidth={0}
                minHeight={450}
                initialDimension={{ width: 900, height: 450 }}
              >
                {chartType === 'activations' ? (
                  <LineChart data={activationData}>
                    <CartesianGrid strokeDasharray="3 3" stroke="#e2e8f0" />
                    <XAxis dataKey="trial" stroke="#64748b" fontSize={11} />
                    <YAxis domain={([0, 1])} stroke="#64748b" fontSize={11} />
                    <Tooltip content={<CustomChartTooltip />} />
                    <Legend wrapperStyle={{ fontSize: '12px', color: '#94a3b8' }} formatter={italicLegendFormatter} />
                    {phaseBoundaries.map((b, i) => (
                      <ReferenceLine key={i} x={b.trial} stroke="#94a3b8" strokeDasharray="5 5"
                        label={{ value: b.phase, position: 'top', fill: '#64748b', fontSize: 11 }} />
                    ))}
                    {selectedUnits.map((unit, i) => (
                      <Line
                        key={unit}
                        type="monotone"
                        dataKey={unit}
                        name={getDisplayName(unit)}
                        stroke={COLORS[i % COLORS.length]}
                        strokeWidth={2}
                        dot={false}
                        activeDot={{ r: 4, strokeWidth: 0 }}
                      />
                    ))}
                  </LineChart>
                ) : chartType === 'signals' ? (
                  <LineChart data={signalData}>
                    <CartesianGrid strokeDasharray="3 3" stroke="#e2e8f0" />
                    <XAxis dataKey="trial" stroke="#64748b" fontSize={11} />
                    <YAxis stroke="#64748b" fontSize={11} />
                    <Tooltip content={<CustomChartTooltip />} />
                    <Legend wrapperStyle={{ fontSize: '12px', color: '#94a3b8' }} formatter={italicLegendFormatter} />
                    {phaseBoundaries.map((b, i) => (
                      <ReferenceLine key={i} x={b.trial} stroke="#94a3b8" strokeDasharray="5 5"
                        label={{ value: b.phase, position: 'top', fill: '#64748b', fontSize: 11 }} />
                    ))}
                    {selectedSignals.map((signal, index) => (
                      <Line
                        key={signal}
                        type="monotone"
                        dataKey={signal}
                        stroke={COLORS[index % COLORS.length]}
                        strokeWidth={2}
                        dot={false}
                        name={getDisplayName(signal)}
                        activeDot={{ r: 4, strokeWidth: 0 }}
                      />
                    ))}
                  </LineChart>
                ) : chartType === 'weights' ? (
                  <LineChart data={weightData}>
                    <CartesianGrid strokeDasharray="3 3" stroke="#e2e8f0" />
                    <XAxis dataKey="trial" stroke="#64748b" fontSize={11} />
                    <YAxis domain={([0, 1])} stroke="#64748b" fontSize={11} />
                    <Tooltip content={<CustomChartTooltip />} />
                    <Legend wrapperStyle={{ fontSize: '12px', color: '#94a3b8' }} formatter={italicLegendFormatter} />
                    {phaseBoundaries.map((b, i) => (
                      <ReferenceLine key={i} x={b.trial} stroke="#94a3b8" strokeDasharray="5 5"
                        label={{ value: b.phase, position: 'top', fill: '#64748b', fontSize: 11 }} />
                    ))}
                    {selectedConns.map((conn, i) => {
                      // Connection names are "PreNPE-PostNPE" — display as "Pre'→Post'"
                      const parts = conn.split('-');
                      const displayConn = parts.length === 2
                        ? `${getDisplayName(parts[0])}→${getDisplayName(parts[1])}`
                        : conn;
                      return (
                        <Line
                          key={conn}
                          type="monotone"
                          dataKey={conn}
                          name={displayConn}
                          stroke={COLORS[i % COLORS.length]}
                          strokeWidth={2}
                          dot={false}
                          activeDot={{ r: 4, strokeWidth: 0 }}
                        />
                      );
                    })}
                  </LineChart>
                ) : tab === 'general' && simulationResults && simulationResults.length > 1 ? (
                  <ComposedChart data={aggregateData}>
                    <CartesianGrid strokeDasharray="3 3" stroke="#e2e8f0" />
                    <XAxis dataKey="unit" stroke="#64748b" fontSize={11} tick={<ItalicUnitTick />} />
                    <YAxis domain={([0, 1])} stroke="#64748b" fontSize={11} />
                    <Tooltip content={<CustomChartTooltip />} />
                    <Bar dataKey="value" fill="#06b6d4" radius={[6, 6, 0, 0]} opacity={0.7}>
                      <ErrorBar dataKey="error" width={4} stroke="#64748b" />
                    </Bar>
                    {Array.from({ length: simulationResults.length }, (_, netIdx) => (
                      <Scatter
                        key={`net-${netIdx}`}
                        dataKey={`net${netIdx}`}
                        fill={SCATTER_COLORS[netIdx % SCATTER_COLORS.length]}
                        legendType="none"
                        r={4}
                      />
                    ))}
                  </ComposedChart>
                ) : (
                  <BarChart data={aggregateData}>
                    <CartesianGrid strokeDasharray="3 3" stroke="#e2e8f0" />
                    <XAxis dataKey="unit" stroke="#64748b" fontSize={11} tick={<ItalicUnitTick />} />
                    <YAxis domain={([0, 1])} stroke="#64748b" fontSize={11} />
                    <Tooltip content={<CustomChartTooltip />} />
                    <Bar dataKey="value" fill="#06b6d4" radius={[6, 6, 0, 0]}>
                      <ErrorBar dataKey="error" width={4} stroke="#64748b" />
                    </Bar>
                  </BarChart>
                )}
              </ResponsiveContainer>
            </div>
            </div>{/* closes exportWrapperRef */}

            {/* Legend for scatter points in general aggregate view */}
            {tab === 'general' && chartType === 'aggregate' && simulationResults && simulationResults.length > 1 && (
              <div className="mt-4 flex flex-wrap items-center gap-3 text-xs text-slate-500">
                <span className="font-bold">{t.results.networksLegend}</span>
                {simulationResults.map((_, i) => (
                  <span key={i} className="flex items-center gap-1">
                    <span
                      className="inline-block w-2.5 h-2.5 rounded-full"
                      style={{ backgroundColor: SCATTER_COLORS[i % SCATTER_COLORS.length] }}
                    />
                    {t.results.net} {i + 1}
                  </span>
                ))}
                <span className="ml-2 text-slate-400">| {t.results.barsOverall} {aggregateMeasure === 'mean' ? t.results.mean : t.results.median}, {t.results.dotsPerNetwork} {aggregateMeasure === 'mean' ? t.results.mean : t.results.median}</span>
              </div>
            )}
          </Card>
        </div>
      </div>
    </div>
    </PageTransition>
  );
}
