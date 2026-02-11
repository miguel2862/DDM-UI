export interface VersionEntry {
  version: string;
  date: string;
  changes: string[];
}

export const VERSION_HISTORY: VersionEntry[] = [
  {
    version: '3.0',
    date: 'February 2026',
    changes: [
      'Complete rebuild: React + TypeScript + Tailwind CSS frontend',
      'R Plumber REST API backend (separated from Shiny)',
      'Interactive network visualization with React Flow',
      'Per-phase timestep selectors in Results',
      'Simulation playback with animated network nodes',
      'Phenomenon gallery with 6 verified templates',
      'Recharts-based results explorer with export',
      'Beginner/Expert mode toggle',
      'Update Procedure (async/sync, random/sequential)',
      'Learning signals (dVTA, dH) visualization',
      'EN/ES internationalization support',
      'Model renamed to Diffuse Temporal Discrepancy (DTD), formerly Diffuse Discrepancy Model (DiffDiscM)',
      'Electron desktop packaging (macOS DMG + Windows installer)',
      'Auto-save, error boundary, toast notifications, and keyboard shortcuts',
      'Parameter Sweep with CSV export',
      'APA 7 figure export with localized titles',
      'Cancel simulation with real-time ETA',
      'Navigation blocked during active simulation',
      'Neural animation during simulation progress',
      'Deep import validation for experiment files',
      'Locked layout persistence in exported experiments',
      'Spinner feedback during template loading',
    ],
  },
  {
    version: '2.5',
    date: '2025',
    changes: [
      'R Shiny interface improvements',
      'Added Burgos & Donahoe (2016) blocking/compound/successive templates',
      'Added Burgos (2000) extinction & reacquisition template',
      'Added Alcalá (2017) autoshaped impulsivity template',
      'visNetwork interactive graph visualization',
      'ITI (Inter-Trial Interval) support',
      'Multiple presentation orders (Random, In bulk, Alternated)',
    ],
  },
  {
    version: '2.0',
    date: '2024',
    changes: [
      'Full R Shiny graphical user interface',
      'Dynamic NPE and Connection editors',
      'Trial Designer with timestep grid',
      'Phase/Contingency builder',
      'Results visualization with plotly charts',
      'Network visualization with visNetwork',
      'Multiple network simulations with aggregate results',
      'Threshold type selection (Gaussian/Beta)',
    ],
  },
  {
    version: '1.0',
    date: '2020-2023',
    changes: [
      'Initial R script implementation',
      'Simulate.DBP() core function based on Donahoe, Burgos & Palmer (1993)',
      'Create.Phases() for trial sequencing',
      'Manual parameter entry via R scripts',
      'CSV-based NPE and Connection definitions',
      'RDS trial definitions',
      'Basic R plotting for results',
    ],
  },
];
