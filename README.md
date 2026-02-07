# DDM-UI <img src="images/icon.ico" alt="DDM Simulator Logo" width="120" align="right"/>

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![Open Science](https://img.shields.io/badge/Open-Science-blue.svg)](https://en.wikipedia.org/wiki/Open_science)
[![Latest Installer (EN)](https://img.shields.io/badge/Latest_Installer_(EN)-v0.05-green.svg)](https://drive.google.com/file/d/1_g1aYD9k8oR31n-Mi2L1dPRHYOjSSriN/view?usp=sharing)
[![R Code](https://img.shields.io/badge/R_Code-Download-blue.svg)](R/DDM_UI%20(2025).R)
[![Online Version](https://img.shields.io/badge/Online_Version-Launch-blue.svg)](https://miguel2862.shinyapps.io/ddm-ui/)

DDM-UI is an R/Shiny interface for the Diffuse Discrepancy Model (DiffDiscM) described by Donahoe, Burgos, and Palmer (1993). It is intended for behavioral research on Pavlovian and operant conditioning under a connectionist framework.

## Reference

Aguayo-Mendoza, M. A. (2025). *DDM-UI: A user interface in R for the discrepancy diffuse model in behavioral research*.
[Behavior Research Methods](https://link.springer.com/article/10.3758/s13428-025-02648-9)

Original SelNet source code (Pascal): [jeburgos-selnet/source-code](https://github.com/jeburgos-selnet/source-code)

## Access Options

### 1. Online version
- [DDM-UI online](https://miguel2862.shinyapps.io/ddm-ui/)
- No local installation
- Suitable for quick exploration
- Limitation: local file import/export is restricted by browser/session context

### 2. Local R version (recommended for full workflow)
- Current script: [R/DDM_UI (2025).R](R/DDM_UI%20(2025).R)
- Legacy script: [R/DDM_UI.R](R/DDM_UI.R)
- Full import/export and persistent local storage
- Compatible with Windows, macOS, and Linux (R/RStudio)

### 3. Windows installer
- [English installer (v0.05)](https://drive.google.com/file/d/1_g1aYD9k8oR31n-Mi2L1dPRHYOjSSriN/view?usp=sharing)
- [Spanish installer (v0.04)](https://drive.google.com/file/d/1gy456KA_bwoXmhocAvuYWLrurgJ-OUnx/view?usp=sharing)

## Quick Start (Local)

1. Open `R/DDM_UI (2025).R` in RStudio.
2. Run:

```r
shiny::runApp('R/DDM_UI (2025).R')
```

3. Use the workflow tabs in order:
- Home
- Network Architecture
- Create Trials
- Configure Contingencies
- Simulate
- Individual Results / General Results

4. Optional: load the example simulation file:
- [Simulation example/Extinction_example.rds](Simulation%20example/Extinction_example.rds)

## Mathematical Core (Implemented Model)

The interface implements the model logic reported in the appendices and the Donahoe-Burgos-Palmer framework, including timestep-wise activation dynamics and discrepancy-driven weight adaptation.

### Activation rule
For each non-input unit \(i\) at time \(t\):

\[
L(x,\sigma)=\frac{1}{1+\exp\left(\frac{-x+0.5}{\sigma}\right)}
\]

- Excitatory and inhibitory inputs are aggregated from presynaptic activity and current weights.
- A stochastic threshold \(\theta_{i,t}\) is sampled (Gaussian or Beta option).
- If excitation dominates inhibition and exceeds threshold, activation is updated with temporal summation.
- Otherwise activation decays or is inhibited to zero, depending on the excitation/inhibition comparison.

A compact form of the reactivation branch is:

\[
a_{i,t}=p_{\mathrm{epsp},i,t}+\tau_i\,L(E_{i,t-1},\sigma_i)\,(1-p_{\mathrm{epsp},i,t})-p_{\mathrm{ipsp},i,t}
\]

where \(\tau_i\) is temporal summation.

### Discrepancy signals
The implementation computes dopaminergic and hippocampal discrepancy terms:

\[
d_{D,t}=\frac{1}{N_D}\sum_{k\in D}(a_{k,t}-a_{k,t-1})
\]

\[
d_{H,t}=\frac{1}{N_H}\sum_{k\in H}|a_{k,t}-a_{k,t-1}|+d_{D,t}(1-d_{H,t-1})
\]

### Learning rule
For connection \(j\rightarrow i\) with weight \(w_{ij,t}\):

If discrepancy is above criterion (\(d_t\geq\delta\)), weights increase proportionally to activity, discrepancy, and remaining capacity term \(r_i\).

If discrepancy is below criterion (\(d_t<\delta\)), weights decrease multiplicatively:

\[
w_{ij,t+1}=w_{ij,t}-\beta\,w_{ij,t}\,a_{j,t}\,a_{i,t}
\]

(and analogously with \(\beta'\) for inhibitory links).

In documentation/examples, decrement parameters are represented as \(\beta=0.1\) and \(\beta'=0.1\).

## 2026 Update Summary

The current update focuses on interface clarity, reproducibility, and robustness while preserving the model structure.

- Full UI text review in English for consistency.
- Help icons now trigger modal explanations reliably.
- Beginner/Advanced mode separation for parameter exposure.
- Path handling improved (optional paths + quick shortcuts).
- Network visualization revised (clean white canvas, export support, stable layout controls).
- Phase order in plots now follows contingency order (e.g., Training before Extinction).
- Simulation pipeline hardened for imported objects (including tibble/data.frame normalization).
- Plot rendering adjusted for stable line visibility and cleaner interactive behavior.

## Repository Structure

- `R/DDM_UI (2025).R`: current application code.
- `R/DDM_UI.R`: previous version for comparison.
- `Simulation example/`: example simulation files.
- `images/`: assets.
- `scripts/`: utility and verification scripts.

## Open Science

This project is released under MIT and intended for extension by the research community.

- Add new simulation templates.
- Extend analysis/visualization modules.
- Adapt file schemas for lab pipelines.
- Contribute improvements through forks and pull requests.

## License

MIT License. See [LICENSE](LICENSE).

## Contact

**Miguel Ángel Aguayo Mendoza**  
Email: miguel.aguayo@academicos.udg.mx / aguayo@iteso.mx  
University of Guadalajara  
Laboratory site: [CEIC](http://www.ceic.cucba.udg.mx/Investigacion/laboratorios?id=13)

## Acknowledgements

- Laboratory for Experimental and Theoretical Research in Learning, Conditioning, and Adaptive Behavior (CEIC).
- Dr. José Enrique Burgos Triano.
- Cristiano Valerio Dos Santos, for adaptation and support in model logic/code integration.
