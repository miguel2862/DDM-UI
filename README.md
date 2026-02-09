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

| Platform | Installer |
|---|---|
| **macOS** (Apple Silicon) | [**DDM-UI.dmg**](https://github.com/miguel2862/DDM-UI/releases/download/3.0/DDM-UI.dmg) |
| **Windows** (64-bit) | [**DDM-UI.exe**](https://github.com/miguel2862/DDM-UI/releases/download/3.0/DDM-UI.exe) |

> [All releases](https://github.com/miguel2862/DDM-UI/releases/tag/3.0)

Both installers are **fully self-contained**. They bundle **R Portable** with all 84+ required packages pre-installed (plumber, jsonlite, dplyr, igraph, tidygraph, ggraph, visNetwork, httpuv, Rcpp, and all transitive dependencies). The user does not need R, RStudio, or any other software. Double-click and the simulation engine starts automatically.

---

### What's Inside v3.0

#### Completely New Interface

The entire user interface has been rebuilt from scratch. The original single-file R Shiny application (3,908 lines of R) has been replaced by a modern **React 18 + TypeScript** frontend communicating with an **R Plumber API** backend:

- **React 18 + TypeScript** with strict typing across the entire codebase.
- **Tailwind CSS v4** for a consistent, responsive design system.
- **Zustand** for centralized state management across all pages.
- **React Flow** (`@xyflow/react`) for interactive, draggable network visualizations with real-time updates.
- **Recharts** for publication-quality data visualization (line, bar, scatter, and composed charts).
- **Framer Motion** for fluid page transitions, staggered list animations, pulsing glows, and interactive hover/tap effects throughout the application.
- **Dark theme** optimized for extended research sessions: deep navy background (#0f172a) with cyan and teal accent colors, carefully tuned contrast ratios across all elements.
- **Lucide icon library** for clean, consistent iconography in every control.

#### Bilingual Interface (English / Spanish)

Full internationalization with one-click language switching. Every label, button, tooltip, description, template name, status message, and error message is translated. The language preference is persisted across sessions via localStorage. Translations cover over 200 UI strings organized by section (dashboard, network, trial, simulation, results, parameter sweep, splash screen).

#### Animated Splash Screen

On launch, an animated splash screen displays with:
- Gradient-filled network icon with animated stroke
- Staggered text entrance animations
- Credits, university, and version information
- An animated loading progress bar
- Version history button opening a modal with the full changelog across all releases

#### Beginner / Advanced Mode

A global toggle in the sidebar switches between:
- **Beginner mode**: exposes only essential parameters (unit name, type, layer, connection source/target/weight). Ideal for students and demonstrations.
- **Advanced mode**: reveals all free parameters for fine-grained control over the simulation. For NPEs: activation, temporal summation ($\tau$), activation decay ($\kappa$), threshold mean ($\mu$), threshold standard deviation ($\sigma$), and logistic sigma. For connections: $\alpha$ (weight gain rate), $\beta$ (weight loss rate), $\alpha'$ and $\beta'$ (inhibitory rates). For simulation: P-update procedure and discrepancy criterion.

---

### Pages in Detail

DDM-UI v3.0 has **7 integrated pages**, each handling a stage of the modeling workflow:

#### 1. Dashboard

The landing page provides an overview of the current model state and quick access to pre-built templates.

**Hero animation**: A continuously animated SVG visualization of the DDM architecture showing 7 neurocomputational processing elements ($S'$, $S''$, $M''$, $M'$, $H$, $D$, $US$) as colored nodes with pulsing opacity and glow effects. Animated connection lines stroke and destroke to represent learning dynamics. Two diffuse discrepancy signal zones (hippocampal $S'' + H$ and dopaminergic $M'' + D$) pulse and shift dimensions as soft gradient clouds, visualizing the model's dual-signal learning mechanism.

**Stat cards**: Four cards display the current model configuration (NPEs, connections, trial types, phases) with unique icons and staggered entrance animations.

**Workflow progress**: A visual strip tracks completion across Architecture, Trials, Contingencies, Simulate, and Results. Steps illuminate in emerald as each is completed.

**Phenomenon gallery**: A grid of 7 pre-built conditioning templates, each loadable with a single click:

| Template | NPEs | Connections | Phases | Description |
|---|---|---|---|---|
| **Acquisition** | 7 | 6 | 1 | Simple CS+US pairing |
| **Extinction** | 7 | 6 | 2 | Training followed by extinction |
| **Spontaneous Recovery** | 7 | 6 | 4 | Train, extinguish, rest, test |
| **Latent Inhibition** | 7 | 6 | 2 | Pre-exposure then training |
| **Blocking** | 13 | 17 | 3 | Kamin blocking ($A+$ then $AX+$ then test) |
| **Successive** | 13 | 17 | 3 | Successive conditioning |
| **Autoshaped Impulsivity** | 13 | 13 | 1 | Small-sooner vs large-later choice with ITI |

Each template pre-configures the complete network architecture, trial designs, and contingency structure. One click loads everything and the user can immediately run the simulation.

#### 2. Network Builder

A full visual editor for designing the neural network architecture.

**Interactive canvas** (React Flow): Drag-and-drop nodes representing NPEs, connected by weighted edges. Each of the 7 layers has its own color and shape:
- **US** (red, hexagon) — unconditioned stimulus
- **Primary Sensory** (blue, rounded square) — sensory input ($S'$)
- **Associative Sensory** (purple, circle) — sensory associations ($S''$)
- **Hippocampal** (amber, circle) — context/novelty detection ($H$)
- **Associative Motor** (teal, circle) — motor associations ($M''$)
- **Primary Motor** (green, rounded square) — motor output ($M'$)
- **Dopaminergic** (pink, circle) — reinforcement signal ($D$)

Connection lines reflect weight through thickness (1.5x to 4x scaling). The fixed $US \to D$ connection (weight = 1.0) is rendered in red. All other connections are gray with animated arrowheads and weight labels displayed to two decimal places.

**Auto-layout**: An "Organize" button arranges all nodes by layer in a structured 4-column academic layout with automatic spacing. A "Lock" button saves the current node positions so they persist into the Results playback visualization. An "Export PNG" button downloads the network diagram as an image file.

**Editor panel**: Two-tab interface for Units and Connections:
- Add/remove NPEs with name, type (excitatory/inhibitory), and layer selection. In advanced mode: activation, $\tau$ (temporal summation), $\kappa$ (activation decay), $\mu$ (threshold mean), $\sigma$ (threshold standard deviation), and logistic $\sigma$.
- Add/remove connections with source, target, and weight. In advanced mode: $\alpha$, $\beta$, $\alpha'$, $\beta'$ learning rate parameters. The historical default of $\beta = 0.12$ is auto-applied.
- **Import/Export**: Save the entire architecture as JSON or load from a previously saved file. The JSON includes NPEs, connections, trials, contingencies, and ITI configuration.

#### 3. Trial Designer

Define trial types and experimental phases with full control over stimulus timing.

**Trial builder**: Create trial types with configurable timesteps (1-10). A stimulus activation table provides a grid where each row is a timestep and each column is a Primary Sensory or US unit. Set activation values (0-1) for each cell, with a learning checkbox per timestep controlling whether weight updates occur. Bulk action buttons (Fill, Clear, Learn On, Learn Off) speed up configuration. A separate ITI (inter-trial interval) mode creates single-timestep rest-period trials with all stimuli at 0.

**Contingency builder**: Define experimental phases by selecting trial types, setting presentation counts (supports different counts per trial type using dash-separated values, e.g., "100-50"), and choosing presentation order:
- **Random**: trials shuffled across all types
- **In bulk**: one trial type completes before the next begins
- **Alternated**: strict interleaving (A, B, A, B...)

Each phase can optionally reset activations between phases or insert ITI trials with configurable min/max intervals and a selectable ITI trial type. Phases can be reordered with up/down controls.

#### 4. Simulation

Configure, run, and visually replay simulations.

**Configuration**: Set number of networks (1-100), threshold type (Gaussian or Beta). In advanced mode: P-update procedure with 4 options (asynchronous random, asynchronous sequential, synchronous random, synchronous sequential — each with detailed tooltip explanations of how activations and weights are updated) and discrepancy criterion slider (0.0001 to 0.1, default 0.001).

**Save/Load**: Download the complete experiment configuration as a versioned JSON file (version 3.0 format) including all NPEs, connections, trials, contingencies, and simulation parameters. Restore a previously saved experiment with a single click.

**Execution**: The Run button validates requirements (minimum 2 NPEs, 1 connection, 1 trial type, 1 phase) and displays specific warnings for missing elements. During execution, an animated progress bar with percentage tracks progress. On completion, a success banner displays network count and elapsed time with options to run again or view results.

**Network playback — trial-by-trial observation**: After simulation, an animated React Flow visualization replays the entire simulation timeline. This is one of the most powerful features in v3.0: you can observe the network evolve across every single timestep, watching activations rise and fall and connections strengthen or weaken in real time.

- **Playback controls**: Play/Pause, Reset (back to timestep 0), Skip Forward (+10 timesteps), speed selection (1x, 2x, 5x, 10x), and a draggable progress slider to jump to any point in the simulation.
- **Node colors** update in real time based on activation level: blue ($a < 0.3$), yellow ($0.3 \leq a \leq 0.6$), red ($a > 0.6$), with glow intensity proportional to activation. Each node displays its name and current activation value to 4 decimal places.
- **Connection thickness** scales dynamically with weight (1.5x to 6x). Fixed connections ($w = 1.0$) remain red.
- **Info badges** display the current phase name, trial number ($T$), and timestep ($t$) at all times.

This allows researchers to step through acquisition, extinction, or any phenomenon moment by moment, observing exactly how the discrepancy signals drive learning across the network.

#### 5. Results

Comprehensive data analysis with multiple visualization types and export options.

**Four chart types**:
- **Activations**: Line chart of unit activations across trials. Dashed vertical reference lines mark phase boundaries with phase labels. Select any combination of units to compare their trajectories.
- **Weights**: Line chart of connection weights across trials. Track how $w_{i,j}$ evolves through training, extinction, rest, and test phases.
- **Aggregate**: Bar chart with standard error bars showing mean or median activation per unit within a selected phase. In General view, each network appears as a colored scatter dot overlaid on the bars.
- **Learning Signals**: Dopaminergic ($\bar{d}_D$, pink) and hippocampal ($\bar{d}_H$, amber) discrepancy signal traces across the entire simulation. Observe how the two signal types diverge during acquisition vs extinction.

**Individual vs General view**: The Individual tab analyzes a single network (selectable by dropdown). The General tab aggregates across all simulated networks with a composed bar + scatter chart where each network is shown as a distinct colored dot, revealing between-network variability.

**Filtering**: Select specific units or connections, choose phases and per-phase timesteps (buttons for phases with few timesteps, sliders for phases with many), and toggle between mean and median measures.

**Export**: Three CSV export options:
- "This Network": full data for the current network
- "All Networks": combined data across all networks
- "Selected": only the currently visible columns (selected units/connections/signals)

All charts include interactive tooltips showing Phase, Trial, and data values, with a 10-color palette for multiple series.

#### 6. Parameter Sweep

Systematic sensitivity analysis for exploring how single parameter changes affect model output.

Select a target element (any connection or NPE), choose a parameter to sweep (weight, $\alpha$, $\beta$, etc. for connections; $\mu$, $\sigma$, $\tau$, $\kappa$, logistic $\sigma$ for NPEs), define a min-max range and number of steps (2-50), and set how many networks to run per step (1-20).

The sweep runs the full simulation at each parameter value and plots the results as a line chart with two traces:
- **Mean activation** (solid cyan line) of the selected output unit
- **Median activation** (dashed teal line) of the selected output unit

This reveals parameter sensitivity, optimal operating ranges, and phase-transition thresholds in the model's behavior — for example, finding the minimum connection weight at which blocking emerges, or how temporal summation ($\tau$) affects acquisition speed.

#### 7. Help

In-app documentation and guidance for using the interface.

---

### Quick Start

1. Download the installer for your platform from the table above.
2. **macOS**: Open the `.dmg` and drag DDM-UI to Applications. On first launch, if Gatekeeper blocks it, right-click the app and select "Open", then confirm.
3. **Windows**: Run the `.exe` installer and follow the prompts. If SmartScreen warns about an unknown publisher, click "More info" then "Run anyway".
4. Launch DDM-UI. The R simulation engine starts automatically in the background.
5. On the Dashboard, click any template (e.g., Extinction) to load a complete experiment, then navigate to Simulation and click Run.
6. After simulation, use the playback controls to step through the network trial by trial, or go to Results for charts and data export.

---

## Additional Access Options

### Online Version (Simplified)

Access instantly through your browser without installation:

- [DDM-UI online](https://miguel2862.shinyapps.io/ddm-ui/)
- Responsive design (desktop/tablet/mobile).
- Good for quick demonstrations and teaching.

Limitations: cannot fully use local file workflow; session-based usage is more restricted than local execution.

### R Version

For researchers who prefer working directly in R/RStudio. The R version provides the core simulation engine with a tab-based Shiny interface but does not include the bilingual (EN/ES) interface, the drag-and-drop visual network builder, the real-time simulation playback, the parameter sweep tool, the dark theme, or the animated transitions. It requires R and all dependencies to be installed manually.

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
| Trial-by-trial playback | **Real-time with controls** | No | No |
| Parameter sweep | **Built-in** | No | Manual scripting |
| Pre-built templates | **7 phenomena, one-click** | Limited | Limited |
| Beginner / Advanced mode | **Yes** | Partial | Partial |
| Animated transitions | **Yes** | No | No |
| Save/Load experiments | **JSON** | Limited | RDS files |
| Export results | **CSV** | CSV | Multiple |
| Offline capable | **Yes** | No | Yes |

---

## Documentation and Theoretical Background

For full conceptual and methodological context:

- [DDM-UI: A user interface in R for the discrepancy diffuse model in behavioral research](https://link.springer.com/article/10.3758/s13428-025-02648-9)

The model follows the Donahoe-Burgos-Palmer architecture and learning logic, with discrepancy-modulated synaptic adaptation.

## Mathematical Formulation

The equations below correspond to Appendices A and B of the published model. They define how units activate and how connections learn.

### Activation Function (Appendix A)

The activation of unit $j$ at moment $t$ is determined by a case-branch equation with two modes: unconditional and conditional.

$$
a_{j,t} =
\begin{cases}
a_{S^*,t}
& \text{if } a_{S^*,t} > 0 \text{ and } j \in D \cup M' \quad \text{(Unconditional)} \\[8pt]

L(\mathit{exc}_{j,t}) + \tau_j \, L(\mathit{exc}_{j,t-1})[1 - L(\mathit{exc}_{j,t})] - L(\mathit{inh}_{j,t})
& \text{if } L(\mathit{exc}) > L(\mathit{inh}) \text{ and } L(\mathit{exc}) \geq \theta_{j,t} \quad \text{(Reactivation)} \\[8pt]

L(\mathit{exc}_{j,t-1}) - \kappa_j \, L(\mathit{exc}_{j,t-1})
& \text{if } L(\mathit{exc}) > L(\mathit{inh}) \text{ and } L(\mathit{exc}) < \theta_{j,t} \quad \text{(Decay)} \\[8pt]

0
& \text{if } L(\mathit{exc}) \leq L(\mathit{inh}) \quad \text{(Deactivation)}
\end{cases}
$$

Where:

- $a_{S^*,t}$ is the activation from a biologically significant stimulus ($S^*$), which directly activates $D$ and $M'$ units
- $\tau_j$ is the temporal summation parameter (default 0.1 for all units)
- $\kappa_j$ is the temporal decay parameter (default 0.1 for all units)
- $\theta_{j,t}$ is a dynamic threshold drawn at each moment from $\mathcal{N}(0.2, 0.15)$

### Excitatory and Inhibitory Input

Each unit receives afferent excitation from $m$ units and inhibition from $n$ units:

$$
\mathit{exc}_{j,t} = \sum_{i=1}^{m} a^+_{i,t} \, w^+_{i,j,t} \qquad \mathit{inh}_{j,t} = \sum_{k=1}^{n} a^-_{k,t} \, w^-_{k,j,t}
$$

### Logistic Transform

All excitatory and inhibitory inputs are passed through a logistic (sigmoid) function before use:

$$
L(x) = \frac{1}{1 + e^{-(x - \mu)/\sigma}}, \quad \mu = 0.5, \; \sigma = 0.1
$$

Note that $L(0) = 0.0006$, meaning there is always a negligible amount of inhibition even in networks with only excitatory units.

### Learning Function (Appendix B)

Connection weights change according to a conditional rule based on the discrepancy magnitude:

$$
\Delta w_{i,j,t} =
\begin{cases}
\alpha_j \, a_{j,t} \, p_{i,t} \, r_{j,t} \, d_t & \text{if } d_t \geq 0.001 \quad \text{(Weight gain)} \\[6pt]
-\beta_j \, a_{i,t} \, a_{j,t} & \text{otherwise} \quad \text{(Weight loss)}
\end{cases}
$$

$$
w_{i,j,t} = w_{i,j,t-1} + \Delta w_{i,j,t}
$$

Where:

- $\alpha_j = 0.5$ (rate of weight gain, for all connections)
- $\beta_j = 0.1$ (rate of weight loss, for all connections)
- $p_{i,t} = \dfrac{a_{i,t} \, w_{i,j,t-1}}{\mathit{exc}_{j,t}}$ (effective synaptic influence from unit $i$)
- $r_{j,t} = 1 - \sum_{i=1}^{n} w_{i,j,t}$ (remaining weight capacity on unit $j$)
- $d_t$ = discrepancy magnitude (defined below)

The term $a_{j,t} \, p_{i,t}$ implements a Hebbian dynamic: co-activation of pre- and post-synaptic units promotes strengthening. The factors $p_{i,t}$ and $r_{j,t}$ implement competitive learning: connections with higher activations and weights gain more, subject to a total weight limit of 1.0 per unit.

All variable weights lie within the open interval $(0.0, 1.0)$, excluding fixed weights ($S^* \to D$ and $S^* \to M'$), which are set to 1.0.

### Discrepancy Signals

The discrepancy $d_t$ depends on the type of post-synaptic unit. It is computed as a temporal difference in activations — not a supervised error signal:

$$
d_t =
\begin{cases}
\bar{d}_{H,t} = \dfrac{1}{n_H} \displaystyle\sum_{k=1}^{n_H} \left| a_{H_k,t} - a_{H_k,t-1} \right| + \bar{d}_{D,t} \left(1 - \bar{d}_{H,t-1}\right) & \text{if } j \in S'' \cup H \\[14pt]
\bar{d}_{D,t} = \dfrac{1}{n_D} \displaystyle\sum_{m=1}^{n_D} \left( a_{D_m,t} - a_{D_m,t-1} \right) & \text{if } j \in M'' \cup D \cup M'
\end{cases}
$$

Key distinctions:

- **Hippocampal discrepancy** ($\bar{d}_H$) uses **absolute** differences: $H$ units detect changes in activation regardless of direction. The term $(1 - \bar{d}_{H,t-1})$ acts as a dynamic saturation mechanism, preventing the signal from exceeding its bound.
- **Dopaminergic discrepancy** ($\bar{d}_D$) uses **signed** differences: $D$ units detect directional mismatches, making motor-pathway connections ($S'' \to M''$, $M'' \to D$, $M'' \to M'$) more susceptible to weight loss under extinction.

These discrepancies are "diffuse" in that the same $d_t$ modulates weight changes across all applicable connections simultaneously.

### The Update Procedure (Appendix C)

All published DTD simulations use an **asynchronous-random** update procedure:

1. At each moment $t$, the list of all updateable units is **shuffled** uniformly at random
2. Activations are updated in that random order, with each new value **replacing the previous one immediately** (asynchronous)
3. After all activations are updated, **discrepancies** are computed
4. Finally, **weights** are updated

This temporal asynchrony is a defining feature of the DTD model. It introduces stochasticity (e.g., if $S''_1$ is updated after $H$, then $H$ sees the previous value of $S''_1$) and is the only update scheme that generates core DTD phenomena such as the interstimulus-interval (ISI) function and its sensitivity to network depth.

---

## Core Functions (Implementation Map)

- `Simulate.DBP()`: main simulation engine; iterates through timesteps, updates activations in shuffled order, computes discrepancy signals, and applies the learning rule.
- `Create.Phases()`: builds the full timestep schedule from contingencies/trials, including optional ITI structure and presentation order (random, bulk, alternated).
- `L(x)`: logistic transform applied to excitatory and inhibitory inputs.
- `ComputeInputs()`: computes $\mathit{exc}_{j,t}$ and $\mathit{inh}_{j,t}$ from presynaptic activations and current weights.
- `dVTA()`: computes dopaminergic discrepancy $\bar{d}_{D,t}$ from signed activation changes in $D$ units.
- `dCA1()`: computes hippocampal discrepancy $\bar{d}_{H,t}$ from absolute activation changes in $H$ units, combined with the dopaminergic signal.
- `Compute.r()`: computes remaining weight capacity $r_{j,t} = 1 - \sum w_{i,j,t}$.
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
