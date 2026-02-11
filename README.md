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

DDM-UI is an open-source simulator for the **Diffuse Temporal Discrepancy** model (DTD; Donahoe, Burgos, & Palmer, 1993) — a connectionist model of Pavlovian and operant conditioning grounded in behavioral neuroscience. It lets researchers build neural network architectures, define experimental contingencies, run simulations, and visualize how associative learning unfolds trial by trial.

Whether you study extinction, blocking, latent inhibition, or autoshaped choice, this tool lets you set up the experiment, run the simulation, and see what the model predicts — without writing a single line of code.

---

## The Model

The Diffuse Temporal Discrepancy (DTD) model belongs to a family of biologically inspired neural network models designed to account for behavioral phenomena in conditioning. Unlike purely mathematical models (Rescorla-Wagner, temporal-difference learning), the DTD operates at the level of neural processing elements (NPEs) organized in layers that mirror functional brain systems.

### Architecture

The basic network has **7 NPEs across 6 layers**:

```
                ┌─────────────────────────────────────┐
                │         Network Architecture        │
                └─────────────────────────────────────┘

   Sensory side                              Motor side
  ┌──────────┐      ┌──────────┐      ┌──────────┐
  │   S'     │─────▶│   S''    │─────▶│   M''    │──────▶┌──────────┐
  │ Primary  │      │ Assoc.   │─┐    │ Assoc.   │──┐    │   M'     │
  │ Sensory  │      │ Sensory  │ │    │ Motor    │  │    │ Primary  │
  └──────────┘      └──────────┘ │    └──────────┘  │    │ Motor    │
                                 │                  │    └──────────┘
                                 ▼                  ▼
                            ┌──────────┐      ┌──────────┐
                            │    H     │      │    D     │◀── US
                            │ Hippo-   │      │ Dopami-  │   (weight
                            │ campal   │      │ nergic   │    = 1.0)
                            └──────────┘      └──────────┘
```

- **S'** (Primary Sensory): Receives direct sensory input — the CS.
- **S''** (Associative Sensory): Integrates sensory information and projects forward.
- **H** (Hippocampal): Involved in contextual and configural processing. Modulated by the hippocampal discrepancy signal.
- **M''** (Associative Motor): Bridges sensory associations to motor output. Modulated by the dopaminergic discrepancy signal.
- **M'** (Primary Motor): The behavioral output — conditioned responding.
- **D** (Dopaminergic): Receives a fixed connection from the US (weight = 1.0). Its activation change drives reinforcement.
- **US** (Unconditioned Stimulus): External input representing a biologically significant event.

### How learning works

Learning in the DTD depends on two **discrepancy signals** — diffuse modulatory signals that determine whether synaptic weights increase or decrease:

1. **Dopaminergic discrepancy** — the mean **signed** change in activation across D-layer units. When the US arrives and D units increase their activation, $\bar{d}_D$ is positive — this signals reinforcement. It modulates connections into $M''$, $D$, and $M'$ layers.

$$
\bar{d}\_{D,t} = \frac{1}{n\_D}\sum\_{m=1}^{n\_D}\left(a\_{D\_m,t} - a\_{D\_m,t-1}\right)
$$

2. **Hippocampal discrepancy** — the mean **absolute** change in activation across H-layer units, combined with the dopaminergic signal. It modulates connections into $S''$ and $H$ layers. The factor $(1 - \bar{d}\_{H,t-1})$ acts as a saturation mechanism that prevents the signal from exceeding its upper bound.

$$
\bar{d}\_{H,t} = \frac{1}{n\_H}\sum\_{k=1}^{n\_H}\left|a\_{H\_k,t} - a\_{H\_k,t-1}\right| + \bar{d}\_{D,t}\left(1 - \bar{d}\_{H,t-1}\right)
$$

On each timestep, the learning rule checks whether the relevant discrepancy $d_t$ exceeds a criterion (default: 0.001):

$$
\Delta w\_{i,j,t} = \begin{cases} \alpha \cdot a\_{j,t} \cdot p\_{i,t} \cdot r\_{j,t} \cdot d\_t & \text{if } d\_t \geq 0.001 \quad \text{(weight gain)} \\\ -\beta \cdot a\_{i,t} \cdot a\_{j,t} & \text{otherwise} \quad \text{(weight loss)} \end{cases}
$$

Where $p\_{i,t}$ is the proportional contribution of unit $i$ to the total excitatory input at $j$, and $r\_{j,t} = 1 - \sum w\_{i,j,t}$ is the remaining weight capacity.

This dual mechanism produces the characteristic learning curves seen in conditioning: rapid acquisition when the US is unexpected, slow extinction when it is omitted, spontaneous recovery after a rest interval, and blocking when a redundant predictor adds no new discrepancy.

### What it can simulate

The simulator ships with 7 pre-built templates covering core phenomena:

| Phenomenon | What it shows |
|---|---|
| **Acquisition** | A neutral CS gradually elicits a conditioned response through repeated CS-US pairing |
| **Extinction** | Conditioned responding decreases when the CS is presented without the US |
| **Spontaneous Recovery** | After extinction, responding partially returns following a rest interval |
| **Latent Inhibition** | Pre-exposure to a CS without consequence slows subsequent conditioning |
| **Blocking** | Prior training with A+ prevents learning about X when AX+ is presented |
| **Successive conditioning** | Independent A+ then X+ training — baseline comparison for blocking |
| **Autoshaped Impulsivity** | Smaller-Sooner vs. Larger-Later choice with delay and context stimuli |

Each template loads the full architecture, trials, and contingencies. You can also build your own networks from scratch for any conditioning paradigm.

---

## Demo

https://github.com/miguel2862/DDM-UI/raw/main/video.mp4

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

Systematic sensitivity analysis for any free parameter in the model. This page lets you ask: *"What happens to behavior if I change parameter X from value A to value B?"*

**How it works:**

1. **Choose a target** — either a connection or an NPE.
2. **Choose a parameter** to sweep:
   - For connections: weight, $\alpha$, $\beta$, $\alpha'$, $\beta'$
   - For NPEs: $\mu$, $\sigma$, $\tau$ (temporal summation), $\kappa$ (activation decay), logistic $\sigma$
3. **Set a range** (min and max values) and **number of steps** (2–50). The simulator will divide the range into evenly spaced values.
4. **Set networks per step** (1–20). At each step, it runs N independent networks and averages the results, accounting for stochastic variability from threshold sampling and randomized update order.
5. **Select an output unit** — any non-input NPE (typically a motor unit like M'1) whose mean activation across the last phase will be plotted.
6. **Run**. The sweep launches one full simulation per step per network. A progress bar tracks completion.

**Output**: A line chart plotting the swept parameter (x-axis) against mean and median activation of the output unit (y-axis). This reveals how sensitive the model is to that parameter — for example, the minimum initial weight at which blocking emerges, how temporal summation ($\tau$) affects acquisition speed, or how the discrepancy criterion interacts with decrement rate ($\beta$).

Results can be exported as CSV for further analysis.

#### 7. Help

In-app documentation.

---

## Tutorial: Building an Acquisition Experiment Step by Step

This walkthrough shows how to set up the simplest possible experiment — Pavlovian acquisition — from scratch. The same experiment is available as a one-click template on the Dashboard, but building it manually helps you understand the structure so you can design your own paradigms.

### Step 1: Build the Network (Network Builder)

Create 7 NPEs (neural processing elements):

| NPE | Type | Layer |
|---|---|---|
| US | Excitatory | US |
| D | Excitatory | Dopaminergic |
| S1 | Excitatory | Primary Sensory |
| S''1 | Excitatory | Associative Sensory |
| H1 | Excitatory | Hippocampal |
| M''1 | Excitatory | Associative Motor |
| M'1 | Excitatory | Primary Motor |

Then create 6 connections:

| From → To | Initial Weight | Role |
|---|---|---|
| S1 → S''1 | 0.10 | CS sensory relay |
| S''1 → H1 | 0.10 | Sensory → hippocampal |
| S''1 → M''1 | 0.10 | Sensory → motor association |
| M''1 → D | 0.10 | Motor → dopamine prediction |
| M''1 → M'1 | 0.10 | Association → behavioral output |
| US → D | **1.00** | Fixed — US unconditionally activates D |

All connections use default learning rates: $\alpha = 0.5$, $\beta = 0.1$, $\alpha' = 0.5$, $\beta' = 0.1$.

> The US → D connection must always have weight = 1.0. This is what makes the US biologically significant — it drives the dopaminergic discrepancy signal without requiring learning.

### Step 2: Design the Trial (Trial Designer)

Create a trial type called **"Training"** with **5 timesteps**:

| Timestep | S1 (CS) | US | Learning |
|---|---|---|---|
| 1 | 1.0 | 0.0 | On |
| 2 | 1.0 | 0.0 | On |
| 3 | 1.0 | 0.0 | On |
| 4 | 1.0 | 0.0 | On |
| 5 | 1.0 | **1.0** | On |

The CS (S1) is present throughout the trial. The US appears only on the last timestep — this is the standard delay conditioning arrangement.

### Step 3: Set Up Contingencies (Trial Designer → Contingencies)

Create one phase:

| Phase name | Trial type | Number of trials | Order | ITI |
|---|---|---|---|---|
| Training | Training | 100 | Random | Off |

This gives you 100 CS-US pairings in randomized order.

### Step 4: Run the Simulation (Simulation)

- Networks: **1** (or more if you want to average across stochastic variation)
- Threshold: **Gaussian** (default)
- Click **Run**

The simulation takes a few seconds. Each trial has 5 timesteps, and on each timestep the network computes activations, checks the discrepancy signals, and updates weights.

### Step 5: Interpret the Results (Results)

Go to the **Results** page and select the **Activations** chart. Plot **M'1** (the primary motor unit — the conditioned response).

You should see M'1 activation rise across trials: low at the beginning (the CS doesn't yet predict the US), increasing as weights strengthen, and approaching an asymptote after roughly 40-60 trials. This is the acquisition curve.

Switch to the **Weights** chart to see how S1→S''1, S''1→M''1, and M''1→M'1 grow over training. The US→D weight stays fixed at 1.0.

Switch to **Learning Signals** to see $\bar{d}_D$ — it starts high (the US is unexpected) and decreases as the network learns to predict it.

> **Tip**: This is the same result you get by clicking **"Acquisition"** in the Dashboard gallery. Use the templates for quick results, and the manual setup when you need custom architectures.

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

- Aguayo-Mendoza, M., & Dos Santos, C. V. (2025). [DDM-UI: A user interface in R for the discrepancy diffuse model in behavioral research](https://doi.org/10.3758/s13428-025-02648-9). *Behavior Research Methods*, 57(5), 128.
- Donahoe, J. W., Burgos, J. E., & Palmer, D. C. (1993). A selectionist approach to reinforcement. *Journal of the Experimental Analysis of Behavior*, 60(1), 17–40.

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
