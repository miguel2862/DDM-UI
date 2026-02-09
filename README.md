# DDM-UI <img src="images/icon.ico" alt="DDM Simulator Logo" width="120" align="right"/>

> Advancing behavioral science through open simulation

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![Open Science](https://img.shields.io/badge/Open-Science-blue.svg)](https://en.wikipedia.org/wiki/Open_science)
[![Download for macOS](https://img.shields.io/badge/Download-macOS_(arm64)-000000?logo=apple&logoColor=white)](https://github.com/miguel2862/DDM-UI/releases/download/3.0/DDM-UI.dmg)
[![Download for Windows](https://img.shields.io/badge/Download-Windows_(x64)-0078D4?logo=windows&logoColor=white)](https://github.com/miguel2862/DDM-UI/releases/download/3.0/DDM-UI.exe)
[![R Code (2026)](https://img.shields.io/badge/R_Code_(2026)-Download-blue.svg)](https://github.com/miguel2862/DDM-UI/blob/main/R/DDM_UI%20(2026).R)
[![R Code (2025)](https://img.shields.io/badge/R_Code_(2025)-Download-lightgrey.svg)](R/DDM_UI%20(2025).R)
[![Online Version](https://img.shields.io/badge/Online_Version-Launch-blue.svg)](https://miguel2862.shinyapps.io/ddm-ui/)

---

## About

The Diffuse Discrepancy Model (DiffDiscM) Simulator is an open-source tool for behavioral research based on Donahoe, Burgos, and Palmer (1993). It provides a connectionist implementation of reinforcement principles for both Pavlovian and operant conditioning.

---

## DDM-UI v3.0 — Standalone Desktop Application

**The biggest update in DDM-UI history.** Version 3.0 is a complete ground-up redesign: a standalone desktop application that bundles a modern React interface with the full R simulation engine running locally. No dependencies. No installation of R. No configuration. Download, install, and run.

### Download

| Platform | Installer | SHA-256 |
|---|---|---|
| **macOS** (Apple Silicon) | [**DDM-UI.dmg**](https://github.com/miguel2862/DDM-UI/releases/download/3.0/DDM-UI.dmg) | `4322c73a6eb5fdd96184e38d1f2466de72d311e3b72e613baf506ce8d6ad8c42` |
| **Windows** (64-bit) | [**DDM-UI.exe**](https://github.com/miguel2862/DDM-UI/releases/download/3.0/DDM-UI.exe) | `d81186a00636bfb9940b3b7b3243e0e987b4b023d0f9c8ae70539343e9348d35` |

> [All releases](https://github.com/miguel2862/DDM-UI/releases/tag/3.0)

Both installers are **fully self-contained**. They bundle **R Portable** with all 84+ required packages pre-installed (plumber, jsonlite, dplyr, igraph, tidygraph, ggraph, visNetwork, httpuv, Rcpp, and all transitive dependencies). The user does not need R, RStudio, or any other software. Double-click and the simulation engine starts automatically.

---

### What's Inside v3.0

#### Completely New Interface

The entire user interface has been rebuilt from scratch using modern web technologies:

- **React 18 + TypeScript** with strict typing across the entire codebase.
- **Tailwind CSS v4** for a consistent, responsive design system.
- **Framer Motion** for fluid page transitions, staggered list animations, pulsing glows, and interactive hover/tap effects throughout the application.
- **Dark theme** optimized for extended research sessions: deep navy background (#0f172a) with cyan and teal accent colors, carefully tuned contrast ratios across all elements.
- **Lucide icon library** for clean, consistent iconography in every control.

#### Bilingual Interface (English / Spanish)

Full internationalization with one-click language switching. Every label, button, tooltip, description, template name, status message, and error message is translated. The language preference is persisted across sessions via localStorage. Translations cover over 200 UI strings organized by section (dashboard, network, trial, simulation, results, parameter sweep, splash screen).

#### Animated Splash Screen

On launch, an animated splash screen displays for 4.5 seconds with:
- Gradient-filled network icon with animated stroke
- Staggered text entrance animations
- Credits, university, and version information
- An animated loading progress bar
- Version history button opening a modal with the full changelog across all releases

#### Beginner / Advanced Mode

A global toggle in the sidebar switches between:
- **Beginner mode**: exposes only essential parameters (unit name, type, layer, connection source/target/weight). Ideal for students and demonstrations.
- **Advanced mode**: reveals all free parameters for fine-grained control over the simulation (activation, temporal summation, activation decay, mean, standard deviation, logistic sigma for NPEs; alpha, beta, alpha prime, beta prime learning rates for connections; P-update procedure and discrepancy criterion for simulation).

---

### Pages in Detail

DDM-UI v3.0 has **7 integrated pages**, each handling a stage of the modeling workflow:

#### 1. Dashboard

The landing page provides an overview of the current model state and quick access to pre-built templates.

**Hero animation**: A continuously animated SVG visualization of the DDM architecture showing 7 neurocomputational processing elements (S', S'', M'', M', H, D, US) as colored nodes with pulsing opacity and glow effects. Animated connection lines stroke and destroke to represent learning dynamics. Two diffuse discrepancy signal zones (hippocampal S''+H and dopaminergic M''+D) pulse and shift dimensions as soft gradient clouds, visualizing the model's dual-signal learning mechanism.

**Stat cards**: Four cards display the current model configuration (NPEs, connections, trial types, phases) with unique icons and staggered entrance animations.

**Workflow progress**: A visual strip tracks completion across Architecture, Trials, Contingencies, Simulate, and Results. Steps illuminate in emerald as each is completed.

**Phenomenon gallery**: A grid of 7 pre-built conditioning templates, each loadable with a single click:

| Template | NPEs | Connections | Phases | Description |
|---|---|---|---|---|
| **Acquisition** | 7 | 6 | 1 | Simple CS+US pairing |
| **Extinction** | 7 | 6 | 2 | Training followed by extinction |
| **Spontaneous Recovery** | 7 | 6 | 4 | Train, extinguish, rest, test |
| **Latent Inhibition** | 7 | 6 | 2 | Pre-exposure then training |
| **Blocking** | 13 | 17 | 3 | Kamin blocking (A+ then AX+ then test) |
| **Successive** | 13 | 17 | 3 | Successive conditioning |
| **Autoshaped Impulsivity** | 13 | 13 | 1 | Small-sooner vs large-later choice with ITI |

Each template pre-configures the complete network architecture, trial designs, and contingency structure. One click loads everything and the user can immediately run the simulation.

#### 2. Network Builder

A full visual editor for designing the neural network architecture.

**Interactive canvas** (React Flow): Drag-and-drop nodes representing NPEs, connected by weighted edges. Each layer has its own color:
- **US** (red) — unconditioned stimulus
- **Primary Sensory** (blue) — sensory input
- **Associative Sensory** (purple) — sensory associations
- **Hippocampal** (amber) — memory/context
- **Associative Motor** (teal) — motor associations
- **Primary Motor** (green) — motor output
- **Dopaminergic** (pink) — reinforcement signal

Connection lines reflect weight through thickness (1.5x to 4x scaling). The fixed US→D connection (weight = 1.0) is rendered in red. All other connections are gray with animated arrowheads and weight labels.

**Auto-layout**: An "Organize" button arranges all nodes by layer in a structured 4-column academic layout. A "Lock" button saves the current positions for use in Results playback. An "Export PNG" button downloads the network diagram as an image.

**Editor panel**: Two-tab interface for Units and Connections:
- Add/remove NPEs with name, type (excitatory/inhibitory), and layer selection. Advanced mode exposes activation, temporal summation (tau), activation decay (kappa), mean (mu), standard deviation (sigma), and logistic sigma.
- Add/remove connections with source, target, and weight. Advanced mode exposes alpha, beta, alpha prime, and beta prime learning rate parameters. A historical default auto-converts beta values of 0.1 to 0.12.
- **Import/Export**: Save the entire architecture as JSON or load from a previously saved file.

#### 3. Trial Designer

Define trial types and experimental phases.

**Trial builder**: Create trial types with configurable timesteps (1-10). A stimulus activation table lets you set the activation value (0-1) for each Primary Sensory and US unit at each timestep, with a learning checkbox per timestep. Bulk action buttons (Fill, Clear, Learn On, Learn Off) speed up configuration. A separate ITI (inter-trial interval) mode creates rest-period trials.

**Contingency builder**: Define experimental phases by selecting trial types, setting presentation counts, and choosing presentation order:
- **Random**: trials shuffled across all types
- **In bulk**: one trial type completes before the next begins
- **Alternated**: strict interleaving (A, B, A, B...)

Each phase can optionally reset activations between phases or insert ITI trials with configurable min/max intervals. Phases can be reordered with up/down controls.

#### 4. Simulation

Configure and run simulations with live network playback.

**Configuration**: Set number of networks (1-100), threshold type (Gaussian or Beta). Advanced mode reveals the P-update procedure (4 options: async random, async sequential, sync random, sync sequential, each with detailed tooltip explanations) and discrepancy criterion slider.

**Save/Load**: Download the complete experiment configuration as a versioned JSON file (version 3.0), or restore a previously saved experiment.

**Execution**: A single button runs the simulation. The interface validates requirements (minimum 2 NPEs, 1 connection, 1 trial type, 1 phase) and shows specific warnings for missing elements. During execution, an animated progress bar and percentage display track progress. On completion, a success banner shows network count and elapsed time.

**Network playback**: After simulation, an animated React Flow visualization replays the results:
- **Playback controls**: Play/Pause, Reset, Skip Forward (+10 timesteps), speed selection (1x, 2x, 5x, 10x), and a progress slider.
- **Node colors** update in real time based on activation level: blue (< 0.3), yellow (0.3-0.6), red (> 0.6), with glow intensity proportional to activation.
- **Connection thickness** scales with weight (1.5x to 6x).
- **Info badges** display current phase name, trial number, and timestep.

#### 5. Results

Comprehensive data analysis with multiple visualization types.

**Four chart types**:
- **Activations**: Line chart of unit activations over trials with dashed phase-boundary markers
- **Weights**: Line chart of connection weights over trials
- **Aggregate**: Bar chart with standard error bars per phase, showing mean or median activation
- **Learning Signals**: Dopaminergic (dVTA, pink) and hippocampal (dH, amber) signal traces

**Filtering**: Select specific units or connections, choose phases and timesteps (per-phase timestep buttons or sliders for phases with many timesteps), switch between Individual (single network) and General (aggregate across all networks with scatter overlay showing each network as a colored dot) views.

**Statistical measures**: Toggle between mean and median for aggregate and general views.

**Export**: Three export options:
- "This Network": CSV of the current network's full data
- "All Networks": CSV combining all networks
- "Selected": CSV with only the currently selected columns

All charts include interactive tooltips (Phase, Trial, values), legends, and a 10-color palette for multiple series.

#### 6. Parameter Sweep

Systematic sensitivity analysis for exploring how parameter changes affect model output.

Select a target (connection or NPE), choose a parameter to sweep (weight, alpha, beta, etc. for connections; mu, sigma, temporal summation, activation decay, logistic sigma for NPEs), define a min-max range and number of steps (2-50), and set how many networks to run per step (1-20).

The sweep runs the full simulation at each parameter value and plots the results as a line chart with two traces:
- **Mean activation** (solid cyan line) of the selected output unit
- **Median activation** (dashed teal line) of the selected output unit

This reveals parameter sensitivity, optimal ranges, and phase-transition thresholds in the model's behavior.

#### 7. Help

In-app documentation and guidance for using the interface.

---

### Quick Start

1. Download the installer for your platform from the table above.
2. **macOS**: Open the `.dmg` and drag DDM-UI to Applications. On first launch, if Gatekeeper blocks it, right-click the app and select "Open", then confirm.
3. **Windows**: Run the `.exe` installer and follow the prompts. If SmartScreen warns about an unknown publisher, click "More info" then "Run anyway".
4. Launch DDM-UI. The R simulation engine starts automatically in the background.
5. On the Dashboard, click any template (e.g., Extinction) to load a complete experiment, then navigate to Simulation and click Run.

---

## Additional Access Options

### Online Version (Simplified)

Access instantly through your browser without installation:

- [DDM-UI online](https://miguel2862.shinyapps.io/ddm-ui/)
- Responsive design (desktop/tablet/mobile).
- Good for quick demonstrations and teaching.

Limitations: cannot fully use local file workflow; session-based usage is more restricted than local execution.

### R Version

For researchers who prefer working directly in R/RStudio. Note that the R version does not include the bilingual interface, the visual network builder, the live simulation playback, the parameter sweep tool, or the dark theme. It provides the core simulation engine with a tab-based Shiny interface.

- Latest version: [DDM_UI (2026).R](https://github.com/miguel2862/DDM-UI/blob/main/R/DDM_UI%20(2026).R)
- Previous version: [DDM_UI (2025).R](R/DDM_UI%20(2025).R)
- Legacy version: [DDM_UI.R](R/DDM_UI.R)

```r
shiny::runApp('R/DDM_UI (2026).R')
```

---

## System Requirements

### Desktop Application (v3.0)

- **macOS**: Apple Silicon (M1/M2/M3/M4), macOS 12 or later
- **Windows**: Windows 10 or later (64-bit)
- 4 GB RAM minimum (8 GB recommended)
- 500 MB free disk space
- No additional software required

### Online Version

- Modern web browser (Chrome, Edge, Firefox, Safari)
- Internet connection

### R Version

- R 4.0+ and RStudio
- Required R packages (installed automatically on first run)

## Version Comparison

| Feature | v3.0 Desktop | Online | R/RStudio |
|---|---|---|---|
| Installation | Download and run | None | Requires R |
| R required | **No** (bundled inside) | No | Yes |
| Language | **EN / ES** | EN / ES | EN / ES |
| Dark theme | **Yes** | No | No |
| Visual network editor | **Drag-and-drop** | Tab-based | Tab-based |
| Live simulation playback | **Real-time with controls** | No | No |
| Parameter sweep | **Built-in** | No | Manual scripting |
| Pre-built templates | **7 phenomena, one-click** | Limited | Limited |
| Beginner / Advanced mode | **Yes** | Partial | Partial |
| Animated splash & transitions | **Yes** | No | No |
| Save/Load experiments | **JSON** | Limited | RDS files |
| Export results | **CSV** | CSV | Multiple |
| Offline capable | **Yes** | No | Yes |
| Performance | Local processing | Network dependent | Local processing |

---

## Documentation and Theoretical Background

For full conceptual and methodological context:

- [DDM-UI: A user interface in R for the discrepancy diffuse model in behavioral research](https://link.springer.com/article/10.3758/s13428-025-02648-9)

The model follows the Donahoe-Burgos-Palmer architecture and learning logic, with discrepancy-modulated synaptic adaptation.

## Mathematical Formulation

### 1. Logistic transform used in activation updates

$$
L(x,\sigma)=\frac{1}{1+\exp\left(\frac{-x+0.5}{\sigma}\right)}
$$

### 2. Reactivation branch (schematic form)

$$
a_{i,t}=p_{\mathrm{epsp},i,t}+\tau_i\,L(E_{i,t-1},\sigma_i)\,(1-p_{\mathrm{epsp},i,t})-p_{\mathrm{ipsp},i,t}
$$

where $\tau_i$ is temporal summation.

### 3. Dopaminergic discrepancy

$$
d_{D,t}=\frac{1}{N_D}\sum_{k\in D}(a_{k,t}-a_{k,t-1})
$$

### 4. Hippocampal discrepancy

$$
d_{H,t}=\frac{1}{N_H}\sum_{k\in H}|a_{k,t}-a_{k,t-1}|+d_{D,t}(1-d_{H,t-1})
$$

### 5. Weight update rules

For excitatory connections (when discrepancy is below criterion):

$$
w_{ij,t+1}=w_{ij,t}-\beta\,w_{ij,t}\,a_{j,t}\,a_{i,t}
$$

For inhibitory connections:

$$
w_{ij,t+1}=w_{ij,t}-\beta'\,w_{ij,t}\,a_{j,t}\,a_{i,t}
$$

Default values: $\beta=0.1, \quad \beta'=0.1$

## Core Functions (Implementation Map)

- `Simulate.DBP()`: main simulation engine; iterates through timesteps, updates activations, computes discrepancy signals, and updates connection weights.
- `Create.Phases()`: builds the full timestep schedule from contingencies/trials, including optional ITI structure.
- `L(x, sigma)`: logistic transform used in activation and input-related update terms.
- `ComputeInputs()`: computes excitatory and inhibitory input totals from presynaptic activity and current weights.
- `dVTA()`: computes dopaminergic discrepancy from timestep-to-timestep activation change in dopaminergic units.
- `dCA1()`: computes hippocampal discrepancy and combines it with dopaminergic discrepancy as specified by the model logic.
- `Compute.r()`: computes residual capacity terms used in discrepancy-driven weight increment.
- `estBetaParams()`: optional helper for Beta-distributed threshold sampling from mean/deviation parameters.

---

## 2026 Change Log (R Shiny UI)

The 2026 update to the R Shiny version complements the prior release with usability and reproducibility improvements while preserving the core model logic.

- Added a guided top workflow strip and contextual hints by tab.
- Improved path handling with optional directory fields and quick path shortcuts.
- Consolidated help icon behavior (modal guidance now triggers reliably).
- Added beginner/advanced interaction flow for easier onboarding.
- Improved network visualization styling, layout behavior, and export workflow.
- Enforced phase display order in plots according to simulation/contingency order.
- Hardened simulation input handling (including tibble/data.frame compatibility).
- Cleaned plot rendering behavior to stabilize lines, legends, and tooltips.

---

## Open Science and Development

DDM-UI is designed to be extended by the scientific community.

- Modify or extend R code modules.
- Add new graphics or analysis tools.
- Implement alternative file formats/workflows.
- Share improvements via forks and pull requests.

## License

This project is licensed under the MIT License. See [LICENSE](LICENSE).

## Contact

**Miguel Angel Aguayo Mendoza**
Email: miguel.aguayo@academicos.udg.mx or aguayo@iteso.mx
University of Guadalajara

Laboratory website: [CEIC](http://www.ceic.cucba.udg.mx/Investigacion/laboratorios?id=13)

## Acknowledgements

- Laboratory for Experimental and Theoretical Research in Learning, Conditioning, and Adaptive Behavior (CEIC), led by Dr. Jose Enrique Burgos Triano.
- Cristiano Valerio Dos Santos, for contributions to model code and logic adaptation from the original framework.

## Troubleshooting

### Desktop Application (v3.0)

- **macOS**: If Gatekeeper blocks the app, right-click and select "Open", then confirm. This only needs to be done once.
- **Windows**: If SmartScreen shows a warning, click "More info" then "Run anyway". The installer is not code-signed but is safe to use.
- **Slow first launch**: The R engine takes a few seconds to initialize on the first run. Wait for the interface to fully load before interacting.
- **Simulation not starting**: Verify your model has at least 2 NPEs, 1 connection, 1 trial type, and 1 phase. The Run button will display specific warnings about what is missing.

### Online version

- Refresh the browser session if UI controls appear stale.
- Use a current version of Chrome, Edge, or Firefox.
- For full functionality, use the desktop application.

### R version

- Ensure required R packages are installed.
- Run from the repository root so relative paths resolve correctly.
- Use the in-app Help section for field-by-field guidance.

---

<p align="center">
Advancing behavioral science through open collaboration and simulation
</p>
