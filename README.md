# DDM-UI <img src="images/icon.ico" alt="DDM Simulator Logo" width="120" align="right"/>

> Advancing behavioral science through open simulation

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![Open Science](https://img.shields.io/badge/Open-Science-blue.svg)](https://en.wikipedia.org/wiki/Open_science)
[![Latest Installer](https://img.shields.io/badge/Latest_Installer_(EN)-v0.05-green.svg)](https://drive.google.com/file/d/1_g1aYD9k8oR31n-Mi2L1dPRHYOjSSriN/view?usp=sharing)
[![R Code](https://img.shields.io/badge/R_Code-Download-blue.svg)](R/DDM_UI%20(2025).R)
[![Online Version](https://img.shields.io/badge/Online_Version-Launch-blue.svg)](https://miguel2862.shinyapps.io/ddm-ui/)

## 🌟 About

The Diffuse Discrepancy Model (DiffDiscM) Simulator is an innovative, open-source tool designed for researchers in behavioral sciences. Based on the seminal work of Donahoe, Burgos and Palmer (1993), it provides a connectionist interpretation of the unified principle of reinforcement for both operant and Pavlovian conditioning.

## 🚀 Key Features

- 🧠 Simulates Pavlovian and operant conditioning
- 🔬 Grounded in neuroanatomy and neurophysiology principles
- 💡 Models neural processing units (NPUs) using advanced activation and learning rules
- 🔄 Incorporates hippocampal and dopaminergic systems for comprehensive learning simulations

## 💻 Access Options

### 1. Online Version (Simplified)
Access instantly through your browser:
- [DDM-UI online](https://miguel2862.shinyapps.io/ddm-ui/)
- No installation required
- Fully responsive design: works seamlessly on desktop, tablet, and mobile devices
- Access from any modern web browser
- Key limitations:
  - Cannot load files locally
  - Limited to web session storage
  - Best for quick simulations and learning

### 2. R Version (Full Features)
Download and run locally:
- Latest version: [DDM_UI (2025).R](R/DDM_UI%20(2025).R)
- Previous version: [DDM_UI.R](R/DDM_UI.R)
- Full functionality including:
  - Save/load simulations
  - Import/export configurations
  - Persistent storage
- Compatible with Windows, macOS, and Linux
- Run directly through R/RStudio

### 3. Windows Installer

1. Download the installer:
   - [English Version (v0.05)](https://drive.google.com/file/d/1_g1aYD9k8oR31n-Mi2L1dPRHYOjSSriN/view?usp=sharing)
   - [Spanish Version (v0.04)](https://drive.google.com/file/d/1gy456KA_bwoXmhocAvuYWLrurgJ-OUnx/view?usp=sharing)
2. Run the installer and follow the on-screen instructions.
3. The installer will create a folder containing:
   - R portable
   - DDM.R file
   - Other necessary components

## 🖥️ System Requirements

### For Online Version
- Modern web browser
- Internet connection
- No local installation needed

### For Local Installation
- Windows 10 or later
- 4 GB RAM (8 GB recommended)
- 500 MB free disk space

## 🔄 Version Differences

| Feature              | shinyapps.io | Local R/Installer |
|---------------------|--------------|------------------|
| Installation        | None needed  | Required        |
| File Storage        | Web session only | Local storage |
| Save Configurations | No           | Yes             |
| Load Saved Files    | No           | Yes             |
| Results Download    | Yes (CSV)    | Yes (Multiple formats) |
| Performance         | Network dependent | Local processing |
| Accessibility      | Any browser  | Requires R/installation |

## 🏁 Quick Start

### Online Version
1. Enter to [DDM-UI online](https://miguel2862.shinyapps.io/ddm-ui/)
2. Configure your simulation directly in the browser
3. Download results as needed
4. Note: All configurations will be lost when closing the browser

### Local Installation
1. Navigate to the installed directory
2. Run the DDM Simulator shortcut or execute DDM.R using the provided R portable
3. Follow the on-screen instructions to set up your simulation

## 📚 Documentation

The DiffDiscM Simulator comes with comprehensive documentation to help you get started and make the most of its features:

### User Interface Overview

The simulator's interface is divided into several key sections:

1. **Home**: Provides an introduction to the DiffDiscM and its capabilities.
2. **Network Architecture**: Allows you to define and visualize the neural network structure.
3. **Create Trials**: Design experimental trials with specific stimuli and timings.
4. **Configure Contingencies**: Set up the experimental conditions and phases.
5. **Simulate**: Run your designed experiments and observe the results.
6. **Individual Results**: Analyze the outcomes for individual simulations.
7. **General Results**: View aggregated results across multiple simulations.

### Help Section

Within the simulator, you can access the "Help" section, which offers:

- Detailed explanations of each component and parameter
- Step-by-step guides for common tasks
- Troubleshooting tips and FAQs

### Theoretical Background

For an in-depth understanding of the Diffuse Discrepancy Model and its applications, we recommend the following article:

[Autoshaped impulsivity: Some explorations with a neural network model](https://www.sciencedirect.com/science/article/abs/pii/S037663572400055X?via%3Dihub)

This article provides valuable insights into the theoretical foundations of the DiffDiscM. Please note that a more comprehensive article focusing specifically on the model is currently in preparation and will be linked here upon publication.

### Example Simulation

To get started with a pre-configured simulation:

1. Navigate to the 'Simulate' section in the interface.
2. In the 'Simulation File Name' field, enter: `Extinction_example`
3. For the 'Simulation Directory Path', example: `Simulations/first_simulation`
4. Click 'Load Simulation' to begin.

[Download Example Simulation File](Simulation%20example/Extinction_example.rds)

For more detailed instructions and in-depth information, refer to the comprehensive documentation within the simulator.

## 🛠️ Open Science & Development

We embrace the principles of open science. The DiffDiscM Simulator is designed to be modified, extended, and improved by the scientific community:

- Explore and modify the R code: [DDM_UI.R](R/DDM_UI.R)
- Add new functions, graphics, or analysis tools
- Customize the user interface
- Implement new file formats or data structures

We encourage researchers to fork the repository, make improvements, and share their work with the community. Together, we can advance the field of behavioral science through collaborative development and open sharing of knowledge.

The original SelNet code, written in Pascal, is also available at [https://github.com/jeburgos-selnet/source-code](https://github.com/jeburgos-selnet/source-code)

## 📄 License

This project is licensed under the MIT License, promoting open and collaborative science. See the [LICENSE](LICENSE) file for details.

## 📞 Contact

For inquiries or collaboration opportunities:

**Miguel Ángel Aguayo Mendoza**  
📧 Email: miguel.aguayo@academicos.udg.mx or aguayo@iteso.mx  
🏫 University of Guadalajara

Discover more about our research at our [laboratory website](http://www.ceic.cucba.udg.mx/Investigacion/laboratorios?id=13).

## 🙏 Acknowledgements

We extend our heartfelt gratitude to:

- The **Laboratory for Experimental and Theoretical Research in Learning, Conditioning, and Adaptive Behavior**, led by Dr. Jose Enrique Burgos Triano, for their invaluable support and guidance.
- Cristiano Valerio Dos Santos from the Centro de Estudios e Investigaciones en Comportamiento (CEIC) for his significant contribution to the model's code and logic, skillfully adapting the work of Donahoe, Burgos and Palmer (1993) to R.

## 🆘 Troubleshooting

### Online Version Issues
- Clear browser cache if experiencing display problems
- Check internet connection
- Try a different modern browser
- For persistent issues, switch to the local version

### Local Version Issues
If you encounter any issues, please check the Help section within the simulator. For further assistance, contact our support team via email at miguel.aguayo@academicos.udg.mx or aguayo@iteso.mx.

---

<p align="center">
  Advancing behavioral science through open collaboration and simulation
</p>
