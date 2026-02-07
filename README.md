# DDM-UI <img src="images/icon.ico" alt="DDM Simulator Logo" width="120" align="right"/>

> Advancing behavioral science through open simulation

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![Open Science](https://img.shields.io/badge/Open-Science-blue.svg)](https://en.wikipedia.org/wiki/Open_science)
[![Latest Installer](https://img.shields.io/badge/Latest_Installer_(EN)-v0.05-green.svg)](https://drive.google.com/file/d/1_g1aYD9k8oR31n-Mi2L1dPRHYOjSSriN/view?usp=sharing)
[![R Code (2026)](https://img.shields.io/badge/R_Code_(2026)-Download-blue.svg)](https://github.com/miguel2862/DDM-UI/blob/main/R/DDM_UI%20(2026).R)
[![R Code (2025)](https://img.shields.io/badge/R_Code_(2025)-Download-lightgrey.svg)](R/DDM_UI%20(2025).R)
[![Online Version](https://img.shields.io/badge/Online_Version-Launch-blue.svg)](https://miguel2862.shinyapps.io/ddm-ui/)

## About

The Diffuse Discrepancy Model (DiffDiscM) Simulator is an open-source tool for behavioral research based on Donahoe, Burgos, and Palmer (1993). It provides a connectionist implementation of reinforcement principles for both Pavlovian and operant conditioning.

## Key Features

### Core model capabilities

- Simulates Pavlovian and operant conditioning.
- Uses neurocomputational units (NPUs) with activation, temporal summation, and decay dynamics.
- Implements discrepancy-driven learning through dopaminergic and hippocampal signals.
- Supports architecture design, trial creation, contingency programming, simulation, and result visualization in one workflow.

### New in 2026

- **Beginner/Advanced modes** to separate essential vs full parameter exposure.
- **One-click example loading workflow** (including root-path selection for local files and quick path shortcuts).
- **Improved UI guidance** with workflow strip + contextual hints across tabs.
- **More reliable help system** (help icons open modal explanations consistently).
- **Refined network visualization** (cleaner styling, better layout handling, export-ready canvas).
- **Plot workflow improvements** (phase order follows contingency sequence, clearer interactive behavior).

## Access Options

### 1. Online Version (Simplified)

Access instantly through your browser:

- [DDM-UI online](https://miguel2862.shinyapps.io/ddm-ui/)
- No installation required.
- Responsive design (desktop/tablet/mobile).
- Good for quick demonstrations and teaching.

Limitations:

- Cannot fully use local file workflow.
- Session-based usage is more restricted than local execution.

### 2. R Version (Full Features)

- Latest version: [DDM_UI (2026).R](https://github.com/miguel2862/DDM-UI/blob/main/R/DDM_UI%20(2026).R)
- Previous version: [DDM_UI (2025).R](R/DDM_UI%20(2025).R)
- Legacy version: [DDM_UI.R](R/DDM_UI.R)

Full functionality includes:

- Save/load simulations.
- Import/export NPUs, connections, trials, and contingencies.
- Local persistent storage.
- Reproducible batch runs in R.

### 3. Windows Installer

1. Download the installer:
   - [English Version (v0.05)](https://drive.google.com/file/d/1_g1aYD9k8oR31n-Mi2L1dPRHYOjSSriN/view?usp=sharing)
   - [Spanish Version (v0.04)](https://drive.google.com/file/d/1gy456KA_bwoXmhocAvuYWLrurgJ-OUnx/view?usp=sharing)
2. Run the installer and follow the on-screen instructions.
3. The installer creates a folder containing:
   - R portable
   - `DDM.R`
   - Required dependencies

## System Requirements

### For Online Version

- Modern web browser
- Internet connection

### For Local Installation

- Windows 10 or later (or R/RStudio on macOS/Linux)
- 4 GB RAM minimum (8 GB recommended)
- 500 MB free disk space

## Version Differences

| Feature | Online | Local R/Installer |
|---------------------|--------------|------------------|
| Installation | None needed | Required |
| File Storage | Session-limited | Local storage |
| Save Configurations | Limited | Yes |
| Load Saved Files | Limited | Yes |
| Results Download | Yes (CSV) | Yes (multiple workflows) |
| Performance | Network dependent | Local processing |
| Accessibility | Any browser | Requires R/installation |

## Quick Start

### Online Version

1. Open [DDM-UI online](https://miguel2862.shinyapps.io/ddm-ui/).
2. Configure your simulation in the browser.
3. Run and inspect results.

### Local R Version

1. Open the repository in RStudio.
2. Run:

```r
shiny::runApp('R/DDM_UI (2026).R')
```

3. Follow tabs in sequence:

- Home
- Network Architecture
- Create Trials
- Configure Contingencies
- Simulate
- Individual Results / General Results

### Example Simulation

To load a pre-configured simulation:

1. Go to **Simulate**.
2. In **Simulation File Name**, enter: `Extinction_example`
3. In **Simulation Directory Path**, set the root folder that contains the file (or use Home shortcuts to auto-fill paths).
4. Click **Load Simulation**.

Example file: [Simulation example/Extinction_example.rds](Simulation%20example/Extinction_example.rds)

## Documentation and Theoretical Background

For full conceptual and methodological context:

- [DDM-UI: A user interface in R for the discrepancy diffuse model in behavioral research](https://link.springer.com/article/10.3758/s13428-025-02648-9)

The model follows the Donahoe-Burgos-Palmer architecture and learning logic, with discrepancy-modulated synaptic adaptation.

## Mathematical Formulation (GitHub-Compatible)

The equations below use GitHub math delimiters (`$` and `$$`).

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

### 5. Weight decrement branch (when discrepancy is below criterion)

For excitatory connections:

$$
w_{ij,t+1}=w_{ij,t}-\beta\,w_{ij,t}\,a_{j,t}\,a_{i,t}
$$

For inhibitory connections:

$$
w_{ij,t+1}=w_{ij,t}-\beta'\,w_{ij,t}\,a_{j,t}\,a_{i,t}
$$

For reporting and examples in this documentation:

$$
\beta=0.1, \quad \beta'=0.1
$$

## Appendix-Based Core Functions (Brief Implementation Map)

This interface follows the appendix logic from the published model. The article contains the full theoretical derivation; here is a concise mapping between equations and implemented functions:

- `Simulate.DBP()`: main simulation engine; iterates through timesteps, updates activations, computes discrepancy signals, and updates connection weights.
- `Create.Phases()`: builds the full timestep schedule from contingencies/trials, including optional ITI structure.
- `L(x, sigma)`: logistic transform used in activation and input-related update terms.
- `ComputeInputs()`: computes excitatory and inhibitory input totals from presynaptic activity and current weights.
- `dVTA()`: computes dopaminergic discrepancy from timestep-to-timestep activation change in dopaminergic units.
- `dCA1()`: computes hippocampal discrepancy and combines it with dopaminergic discrepancy as specified by the model logic.
- `Compute.r()`: computes residual capacity terms used in discrepancy-driven weight increment.
- `estBetaParams()`: optional helper for Beta-distributed threshold sampling from mean/deviation parameters.

In short, the implementation keeps the original flow: build phases → update activations per timestep → compute discrepancy terms → apply learning rule (increment/decrement) → store trajectories for analysis.

## 2026 Change Log (UI and Workflow)

The 2026 update complements the prior release with usability and reproducibility improvements while preserving the core model logic.

- Added a guided top workflow strip and contextual hints by tab.
- Improved path handling with optional directory fields and quick path shortcuts.
- Consolidated help icon behavior (modal guidance now triggers reliably).
- Added beginner/advanced interaction flow for easier onboarding.
- Improved network visualization styling, layout behavior, and export workflow.
- Enforced phase display order in plots according to simulation/contingency order.
- Hardened simulation input handling (including tibble/data.frame compatibility).
- Cleaned plot rendering behavior to stabilize lines, legends, and tooltips.

## User Interface Overview

1. **Home**: Introduction, quick-start template, and path shortcuts.
2. **Network Architecture**: NPU and connection construction + network graph.
3. **Create Trials**: Trial timing and input schedule definition.
4. **Configure Contingencies**: Phase and trial presentation programming.
5. **Simulate**: Batch execution and save/load of sessions.
6. **Individual Results**: Single-network trajectories.
7. **General Results**: Aggregated patterns across networks.
8. **Help**: In-app guidance and examples.

## Open Science and Development

DDM-UI is designed to be extended by the scientific community.

- Modify or extend R code modules.
- Add new graphics or analysis tools.
- Implement alternative file formats/workflows.
- Share improvements via forks and pull requests.

## License

This project is licensed under the MIT License. See [LICENSE](LICENSE).

## Contact

**Miguel Ángel Aguayo Mendoza**  
Email: miguel.aguayo@academicos.udg.mx or aguayo@iteso.mx  
University of Guadalajara

Laboratory website: [CEIC](http://www.ceic.cucba.udg.mx/Investigacion/laboratorios?id=13)

## Acknowledgements

- Laboratory for Experimental and Theoretical Research in Learning, Conditioning, and Adaptive Behavior (CEIC), led by Dr. Jose Enrique Burgos Triano.
- Cristiano Valerio Dos Santos, for contributions to model code and logic adaptation from the original framework.

## Troubleshooting

### Online version

- Refresh the browser session if UI controls look stale.
- Try a current version of Chrome/Edge/Firefox.
- For full data workflow, use the local R version.

### Local version

- Ensure required R packages are installed.
- Run from the repository root so relative paths resolve correctly.
- Use the in-app Help section for field-by-field guidance.

---

<p align="center">
Advancing behavioral science through open collaboration and simulation
</p>


<p align="center">
Advancing behavioral science through open collaboration and simulation
</p>
