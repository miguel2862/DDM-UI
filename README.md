# DDM-UI <img src="images/icon.ico" alt="DDM Simulator Logo" width="120" align="right"/>

> Advancing behavioral science through open simulation

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![Open Science](https://img.shields.io/badge/Open-Science-blue.svg)](https://en.wikipedia.org/wiki/Open_science)
[![Download for macOS](https://img.shields.io/badge/Download-macOS_(arm64)-000000?logo=apple&logoColor=white)](https://github.com/miguel2862/DDM-UI/releases/download/3.0/DDM-UI.dmg)
[![Download for Windows](https://img.shields.io/badge/Download-Windows_(x64)-0078D4?logo=windows&logoColor=white)](https://github.com/miguel2862/DDM-UI/releases/download/3.0/DDM-UI.exe)
[![React](https://img.shields.io/badge/React-18-61DAFB?logo=react&logoColor=white)](https://react.dev)
[![TypeScript](https://img.shields.io/badge/TypeScript-5-3178C6?logo=typescript&logoColor=white)](https://www.typescriptlang.org)
[![R Code (2026)](https://img.shields.io/badge/R_Code_(2026)-Source-blue.svg)](https://github.com/miguel2862/DDM-UI/blob/main/R/DDM_UI%20(2026).R)
[![Online Version](https://img.shields.io/badge/Online_Version-Launch-blue.svg)](https://miguel2862.shinyapps.io/ddm-ui/)

---

## About

The Diffuse Discrepancy Model (DiffDiscM) Simulator is an open-source tool for behavioral research based on Donahoe, Burgos, and Palmer (1993). It provides a connectionist implementation of reinforcement principles for both Pavlovian and operant conditioning.

---

## DDM-UI v3.0 — Standalone Desktop Application

Version 3.0 is a complete ground-up redesign. The original single-file R Shiny application (3,908 lines) has been replaced by a standalone desktop application: a React 18 + TypeScript frontend driving the full R simulation engine through a local Plumber API. Everything is bundled inside the installer — R Portable, all 84+ packages, the frontend, the API. No dependencies. No configuration. Download, install, open.

### Download

| Platform | Installer |
|---|---|
| **macOS** (Apple Silicon) | [**DDM-UI.dmg**](https://github.com/miguel2862/DDM-UI/releases/download/3.0/DDM-UI.dmg) |
| **Windows** (64-bit) | [**DDM-UI.exe**](https://github.com/miguel2862/DDM-UI/releases/download/3.0/DDM-UI.exe) |

> [All releases](https://github.com/miguel2862/DDM-UI/releases/tag/3.0)

Both installers bundle **R Portable** with all required packages pre-installed (plumber, jsonlite, dplyr, igraph, tidygraph, ggraph, visNetwork, httpuv, Rcpp, and every transitive dependency). The user does not need R, RStudio, or any other software.

---

### What's in v3.0

**Frontend stack**: React 18, TypeScript (strict), Tailwind CSS v4, Zustand (state), React Flow (network canvas), Recharts (charts), Framer Motion (animations), Lucide (icons). Dark theme throughout — navy (#0f172a) background, cyan/teal accents.

**Bilingual**: Full English / Spanish interface. Every label, tooltip, template name, status message, and error is translated. Language persists across sessions. Over 200 translated strings.

**Beginner / Advanced toggle**: Beginner mode shows only the essentials (name, type, layer, weight). Advanced mode exposes every free parameter — $\tau$, $\kappa$, $\mu$, $\sigma$, logistic $\sigma$ for units; $\alpha$, $\beta$, $\alpha'$, $\beta'$ for connections; P-update procedure and discrepancy criterion for the simulation.

**Splash screen**: Animated launch screen with gradient network icon, staggered text, loading bar, and a version history modal.

---

### Pages

#### 1. Dashboard

Landing page. Animated SVG hero showing the 7-unit DDM architecture ($S'$, $S''$, $M''$, $M'$, $H$, $D$, $US$) with pulsing nodes, stroking connection lines, and two diffuse discrepancy signal clouds ($S''+H$ hippocampal, $M''+D$ dopaminergic). Four stat cards track the current model (NPEs, connections, trial types, phases). A workflow strip shows progress across the five stages.

**Phenomenon gallery** — 7 pre-built templates, one click each:

| Template | NPEs | Connections | Phases | Description |
|---|---|---|---|---|
| Acquisition | 7 | 6 | 1 | CS+US pairing |
| Extinction | 7 | 6 | 2 | Training → Extinction |
| Spontaneous Recovery | 7 | 6 | 4 | Train → Extinct → Rest → Test |
| Latent Inhibition | 7 | 6 | 2 | Pre-exposure → Training |
| Blocking | 13 | 17 | 3 | $A+$ → $AX+$ → Test |
| Successive | 13 | 17 | 3 | Successive conditioning |
| Autoshaped Impulsivity | 13 | 13 | 1 | SS/LL choice with ITI |

Each template loads the full architecture, trials, and contingencies. Ready to simulate immediately.

#### 2. Network Builder

Visual editor with a React Flow canvas on the left and an editor panel on the right.

**Canvas**: Drag-and-drop NPE nodes colored by layer — US (red), Primary Sensory (blue, $S'$), Associative Sensory (purple, $S''$), Hippocampal (amber, $H$), Associative Motor (teal, $M''$), Primary Motor (green, $M'$), Dopaminergic (pink, $D$). Connection thickness scales with weight. The fixed $US \to D$ connection renders in red. Weight labels display to two decimal places. "Organize" auto-arranges nodes by layer. "Lock" saves positions for playback. "Export PNG" downloads the diagram.

**Editor**: Two tabs — Units and Connections. Add NPEs by name, type (excitatory/inhibitory), and layer. Add connections by source, target, and weight. Advanced mode exposes all parameters. Import/export the full architecture as JSON.

#### 3. Trial Designer

**Trial builder**: Define trial types with 1-10 timesteps. A grid lets you set activation (0–1) for each Primary Sensory and US unit at each timestep, with a learning toggle per row. Bulk buttons: Fill, Clear, Learn On, Learn Off. Separate ITI mode for rest-period trials.

**Contingency builder**: Define phases by selecting trial types, setting counts (dash-separated for mixed, e.g. "100-50"), and choosing order — Random, In Bulk, or Alternated. Each phase can reset activations or insert ITI trials with min/max intervals. Phases reorder with up/down controls.

#### 4. Simulation

**Config**: Networks (1–100), threshold (Gaussian or Beta). Advanced: P-update procedure (async random, async sequential, sync random, sync sequential — each with tooltip) and discrepancy criterion (0.0001–0.1). Save/load the entire experiment as versioned JSON.

**Execution**: The Run button validates the model and shows what's missing. Progress bar with percentage during simulation. Completion banner with elapsed time.

**Trial-by-trial playback**: After simulation, the network visualization replays the entire timeline. Play/Pause, Reset, Skip (+10 steps), speed (1x/2x/5x/10x), draggable slider. Node colors shift with activation — blue (low), yellow (mid), red (high) — with glow proportional to $a_{j,t}$. Each node shows its name and activation to 4 decimals. Connection thickness scales with weight. Phase, trial number, and timestep badges update continuously. Step through acquisition trial by trial, watch extinction unfold moment by moment, observe how discrepancy signals reshape the network.

#### 5. Results

Four chart types:
- **Activations** — line chart of selected units across trials, dashed phase boundaries
- **Weights** — line chart of $w_{i,j}$ across trials
- **Aggregate** — bar chart with SE bars, mean or median per phase. General view overlays each network as a colored scatter dot
- **Learning Signals** — $\bar{d}_D$ (pink) and $\bar{d}_H$ (amber) traces

Individual tab for single-network analysis, General tab for cross-network aggregation. Filter by unit, connection, phase, timestep. Toggle mean/median. Three CSV exports: This Network, All Networks, Selected columns.

#### 6. Parameter Sweep

Sensitivity analysis. Pick a connection or NPE, choose a parameter (weight, $\alpha$, $\beta$, $\mu$, $\sigma$, $\tau$, $\kappa$...), set a range and step count (2–50), run N networks per step (1–20). Output: mean and median activation of a selected motor unit plotted against the swept parameter. Useful for finding thresholds — e.g., the minimum weight where blocking appears, or how $\tau$ affects acquisition speed.

#### 7. Help

In-app documentation.

---

### Quick Start

1. Download the installer for your platform.
2. **macOS**: Open the `.dmg`, drag to Applications. If Gatekeeper blocks it, right-click → Open → confirm.
3. **Windows**: Run the `.exe`, follow prompts. If SmartScreen warns, click "More info" → "Run anyway".
4. Launch DDM-UI. R starts in the background automatically.
5. Click any template on the Dashboard (e.g. Extinction), go to Simulation, click Run.
6. Use playback to step through the network, or go to Results for charts and CSV export.

---

## Other Access Options

### Online Version

Browser-based, no installation: [DDM-UI online](https://miguel2862.shinyapps.io/ddm-ui/). English only. Good for demos and teaching. Limited file workflow, session-based.

### R Version

For researchers who want to script directly in R/RStudio. English only. Tab-based Shiny interface — no visual network builder, no trial-by-trial playback, no parameter sweep, no dark theme. Requires R and manual package installation.

- Latest: [DDM_UI (2026).R](https://github.com/miguel2862/DDM-UI/blob/main/R/DDM_UI%20(2026).R)
- Previous: [DDM_UI (2025).R](R/DDM_UI%20(2025).R)
- Legacy: [DDM_UI.R](R/DDM_UI.R)

```r
shiny::runApp('R/DDM_UI (2026).R')
```

---

## System Requirements

### Desktop (v3.0)

- **macOS**: Apple Silicon (M1/M2/M3/M4), macOS 12+
- **Windows**: Windows 10+ (64-bit)
- 4 GB RAM minimum, 8 GB recommended
- 500 MB disk space
- Nothing else required

### Online

- Modern browser, internet connection

### R Version

- R 4.0+, RStudio

## Comparison

| Feature | v3.0 Desktop | Online | R/RStudio |
|---|---|---|---|
| Installation | Download and run | None | Requires R |
| R required | No (bundled) | No | Yes |
| Language | **EN / ES** | EN | EN |
| Dark theme | Yes | No | No |
| Network editor | Drag-and-drop | Tab-based | Tab-based |
| Trial-by-trial playback | Yes | No | No |
| Parameter sweep | Built-in | No | Manual |
| Templates | 7, one-click | Limited | Limited |
| Beginner / Advanced | Yes | Partial | Partial |
| Save/Load | JSON | Limited | RDS |
| CSV export | Yes | Yes | Yes |
| Offline | Yes | No | Yes |

---

## Reference

- [DDM-UI: A user interface in R for the discrepancy diffuse model in behavioral research](https://link.springer.com/article/10.3758/s13428-025-02648-9)

---

## Core Functions

- `Simulate.DBP()` — simulation engine. Iterates timesteps, updates activations in shuffled order (asynchronous-random), computes discrepancy signals, applies the learning rule.
- `Create.Phases()` — builds the timestep schedule from contingencies/trials, including ITI and presentation order.
- `L(x)` — logistic transform on excitatory and inhibitory inputs.
- `ComputeInputs()` — computes excitatory and inhibitory totals from presynaptic activations and weights.
- `dVTA()` — dopaminergic discrepancy from signed activation changes in $D$ units.
- `dCA1()` — hippocampal discrepancy from absolute activation changes in $H$ units, combined with the dopaminergic signal.
- `Compute.r()` — remaining weight capacity $r_{j,t} = 1 - \sum w_{i,j,t}$.
- `estBetaParams()` — Beta-distributed threshold sampling.

---

## License

MIT License. See [LICENSE](LICENSE).

## Contact

**Miguel Angel Aguayo Mendoza**
miguel.aguayo@academicos.udg.mx | aguayo@iteso.mx
University of Guadalajara — [CEIC Lab](http://www.ceic.cucba.udg.mx/Investigacion/laboratorios?id=13)

## Acknowledgements

- Laboratory for Experimental and Theoretical Research in Learning, Conditioning, and Adaptive Behavior (CEIC), led by Dr. Jose Enrique Burgos Triano.
- Cristiano Valerio Dos Santos, for contributions to model code and logic adaptation from the original framework.

## Troubleshooting

**macOS** — Gatekeeper block: right-click → Open → confirm (once).
**Windows** — SmartScreen warning: "More info" → "Run anyway".
**Slow first launch** — R initializes on first run, wait a few seconds.
**Can't run simulation** — need at least 2 NPEs, 1 connection, 1 trial type, 1 phase.

---

<p align="center">
Advancing behavioral science through open collaboration and simulation
</p>

