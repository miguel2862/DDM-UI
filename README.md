# DDM-UI <img src="images/icon.ico" alt="DDM Simulator Logo" width="120" align="right"/>

> Advancing behavioral science through open simulation

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![Open Science](https://img.shields.io/badge/Open-Science-blue.svg)](https://en.wikipedia.org/wiki/Open_science)
[![Download for macOS](https://img.shields.io/badge/Download-macOS_(arm64)-000000?logo=apple&logoColor=white)](https://github.com/miguel2862/DDM-UI/releases/download/3.0/DDM-UI-3.0.0-arm64.dmg)
[![Download for Windows](https://img.shields.io/badge/Download-Windows_(x64)-0078D4?logo=windows&logoColor=white)](https://github.com/miguel2862/DDM-UI/releases/download/3.0/DDM-UI.Setup.3.0.0.exe)
[![R Code (2026)](https://img.shields.io/badge/R_Code_(2026)-Download-blue.svg)](https://github.com/miguel2862/DDM-UI/blob/main/R/DDM_UI%20(2026).R)
[![R Code (2025)](https://img.shields.io/badge/R_Code_(2025)-Download-lightgrey.svg)](R/DDM_UI%20(2025).R)
[![Online Version](https://img.shields.io/badge/Online_Version-Launch-blue.svg)](https://miguel2862.shinyapps.io/ddm-ui/)

---

## About

The Diffuse Discrepancy Model (DiffDiscM) Simulator is an open-source tool for behavioral research based on Donahoe, Burgos, and Palmer (1993). It provides a connectionist implementation of reinforcement principles for both Pavlovian and operant conditioning.

---

## DDM-UI v3.0 — Desktop Application

**The biggest update in DDM-UI history.** Version 3.0 is a complete redesign: a standalone desktop application built from the ground up with a modern React interface and the full R simulation engine running locally inside the app. No dependencies, no installation of R, no configuration. Just download, install, and run.

### Download Installers

| Platform | Installer | Size | SHA-256 |
|---|---|---|---|
| **macOS** (Apple Silicon) | [DDM-UI-3.0.0-arm64.dmg](https://github.com/miguel2862/DDM-UI/releases/download/3.0/DDM-UI-3.0.0-arm64.dmg) | ~389 MB | `4322c73a6eb5fdd96184e38d1f2466de72d311e3b72e613baf506ce8d6ad8c42` |
| **Windows** (64-bit) | [DDM-UI Setup 3.0.0.exe](https://github.com/miguel2862/DDM-UI/releases/download/3.0/DDM-UI.Setup.3.0.0.exe) | ~261 MB | `d81186a00636bfb9940b3b7b3243e0e987b4b023d0f9c8ae70539343e9348d35` |

> Both installers are **fully self-contained**: they bundle R Portable with all required packages. The user does not need R, RStudio, or any other software installed. Just double-click and go.

### What's New in v3.0

#### Completely New Interface

- **Modern React frontend** built with TypeScript, Tailwind CSS v4, and Framer Motion animations.
- **Dark theme** designed for extended work sessions (navy-900 background with cyan/teal accents).
- **7 integrated pages**: Dashboard, Network Builder, Trial Designer, Simulation, Results, Parameter Sweep, and Help.
- **Bilingual support**: full English and Spanish interface with one-click language switching (persisted across sessions).

#### Visual Network Builder

- **Drag-and-drop editor** powered by React Flow for designing neural architectures.
- Visually create, connect, and configure NPEs (neurocomputational processing elements) and connections.
- Academic brain-inspired layer layout: US, Primary Sensory, Associative, Motor.
- JSON import/export for saving and sharing network architectures.

#### Trial Designer

- Create trial types with precise timestep-by-timestep stimulus schedules.
- Design experimental phases with multiple trial types and contingencies.
- Visual timeline editor for intuitive trial construction.

#### Live Simulation

- **Real-time execution controls**: play, pause, step-through, and skip.
- Watch NPE activations and connection weights update live during simulation.
- Visual network overlay showing activation levels in real time.

#### Results & Analysis

- **Multi-format charts**: line, bar, and scatter plots with Recharts.
- Filter by phase, trial type, or individual NPEs/connections.
- **CSV export** for further analysis in external tools.
- Interactive tooltips and legends for detailed data inspection.

#### Parameter Sweep

- Systematic parameter exploration across configurable ranges.
- Sweep connection weights, NPE properties, or learning rates.
- Compare results across parameter sets in unified visualizations.

#### Standalone Desktop Application

- **Electron shell** wrapping the React frontend + R Plumber backend.
- **R Portable** bundled inside: all required R packages (plumber, jsonlite, dplyr, igraph, tidygraph, ggraph, visNetwork, and 80+ dependencies) are included.
- No internet connection required after installation.
- Custom application icon on both platforms.

#### Pre-built Templates

- Load classic conditioning phenomena with one click from the Dashboard:
  - Acquisition
  - Extinction
  - Spontaneous Recovery
  - Latent Inhibition
  - Blocking
- Each template pre-configures the full network architecture, trial design, and contingencies.

---

## Additional Access Options

### Online Version (Simplified)

Access instantly through your browser:

- [DDM-UI online](https://miguel2862.shinyapps.io/ddm-ui/)
- No installation required.
- Responsive design (desktop/tablet/mobile).
- Good for quick demonstrations and teaching.

Limitations:

- Cannot fully use local file workflow.
- Session-based usage is more restricted than local execution.

### R Version (Full Scriptable Features)

For researchers who prefer working directly in R/RStudio:

- Latest version: [DDM_UI (2026).R](https://github.com/miguel2862/DDM-UI/blob/main/R/DDM_UI%20(2026).R)
- Previous version: [DDM_UI (2025).R](R/DDM_UI%20(2025).R)
- Legacy version: [DDM_UI.R](R/DDM_UI.R)

Full functionality includes:

- Save/load simulations.
- Import/export NPUs, connections, trials, and contingencies.
- Local persistent storage.
- Reproducible batch runs in R.

### Legacy Windows Installer (R Shiny)

Previous standalone installers for the R Shiny version:

- [English Version (v0.05)](https://drive.google.com/file/d/1_g1aYD9k8oR31n-Mi2L1dPRHYOjSSriN/view?usp=sharing)
- [Spanish Version (v0.04)](https://drive.google.com/file/d/1gy456KA_bwoXmhocAvuYWLrurgJ-OUnx/view?usp=sharing)

---

## System Requirements

### For Desktop Application (v3.0)

- **macOS**: Apple Silicon (M1/M2/M3/M4), macOS 12 or later
- **Windows**: Windows 10 or later (64-bit)
- 4 GB RAM minimum (8 GB recommended)
- 500 MB free disk space
- No additional software required

### For Online Version

- Modern web browser
- Internet connection

### For R Version

- R 4.0+ and RStudio
- Required R packages (installed automatically on first run)

## Version Comparison

| Feature | v3.0 Desktop | Online | R/RStudio |
|---|---|---|---|
| Installation | Download and run | None needed | Requires R |
| R required externally | **No** (bundled) | No | Yes |
| Language | EN / ES | EN / ES | EN / ES |
| Network visual editor | Drag-and-drop | Tab-based | Tab-based |
| Live simulation view | Real-time | No | No |
| Parameter sweep | Built-in | No | Manual scripting |
| Dark theme | Yes | No | No |
| Templates (one-click) | Yes | Limited | Limited |
| File storage | Local | Session-limited | Local |
| CSV export | Yes | Yes | Yes |
| Offline capable | Yes | No | Yes |
| Performance | Local processing | Network dependent | Local processing |

---

## Quick Start

### Desktop Application (v3.0)

1. Download the installer for your platform (see table above).
2. Install:
   - **macOS**: Open the `.dmg`, drag DDM-UI to Applications.
   - **Windows**: Run the `.exe` installer, follow the on-screen prompts.
3. Launch DDM-UI.
4. On the Dashboard, select a template (e.g., Extinction) to get started immediately, or navigate to the Network Builder to design your own architecture.

### R Version

1. Open the repository in RStudio.
2. Run:

```r
shiny::runApp('R/DDM_UI (2026).R')
```

3. Follow tabs in sequence: Home, Network Architecture, Create Trials, Configure Contingencies, Simulate, Individual Results, General Results.

### Example Simulation

To load a pre-configured simulation in the R version:

1. Go to **Simulate**.
2. In **Simulation File Name**, enter: `Extinction_example`
3. In **Simulation Directory Path**, set the root folder that contains the file.
4. Click **Load Simulation**.

Example file: [Simulation example/Extinction_example.rds](Simulation%20example/Extinction_example.rds)

---

## Technical Architecture (v3.0)

The desktop application uses a dual-process architecture:

```
Electron Shell
  ├── React Frontend (TypeScript + Vite + Tailwind CSS v4)
  │     ├── Zustand (state management)
  │     ├── React Flow (network visualization)
  │     ├── Recharts (data visualization)
  │     ├── Framer Motion (animations)
  │     └── i18n (English / Spanish)
  │
  └── R Plumber API (bundled R Portable)
        ├── Simulate.DBP() — simulation engine
        ├── Create.Phases() — phase/trial builder
        └── All R packages (plumber, jsonlite, dplyr, igraph, etc.)
```

The frontend communicates with the R backend via HTTP API calls to `localhost`. The R process starts automatically when the application launches and shuts down when the application closes.

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

## User Interface Overview (R Shiny Version)

1. **Home**: Introduction, quick-start template, and path shortcuts.
2. **Network Architecture**: NPU and connection construction + network graph.
3. **Create Trials**: Trial timing and input schedule definition.
4. **Configure Contingencies**: Phase and trial presentation programming.
5. **Simulate**: Batch execution and save/load of sessions.
6. **Individual Results**: Single-network trajectories.
7. **General Results**: Aggregated patterns across networks.
8. **Help**: In-app guidance and examples.

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

- **macOS**: If the app is blocked by Gatekeeper, right-click the app and select "Open", then confirm.
- **Windows**: If SmartScreen warns about an unknown publisher, click "More info" then "Run anyway". The installer is not code-signed but is safe to use.
- If the simulation does not start, wait a few seconds for the R engine to initialize on first launch.

### Online version

- Refresh the browser session if UI controls look stale.
- Try a current version of Chrome/Edge/Firefox.
- For full data workflow, use the desktop application or R version.

### R version

- Ensure required R packages are installed.
- Run from the repository root so relative paths resolve correctly.
- Use the in-app Help section for field-by-field guidance.

---

<p align="center">
Advancing behavioral science through open collaboration and simulation
</p>
