# NYU Research: Mathematical Modeling of Multinuclear Scaling in Muscle Cells

This repository contains simulation code developed as part of a research project at **New York University (NYU)** through the **GSTEM 2024 program**, under the mentorship of **Dr. Alex Mogilner**.

The project investigates how signaling molecules released from the **Neuromuscular Junction (NMJ)** influence the spatial organization and positioning of nuclei in multinucleated muscle cells.

## Research Motivation

Multinucleated muscle fibers rely on precise nuclear positioning to maintain efficient gene expression, cellular transport, and muscle function. Disruptions in nuclear spacing have been linked to neuromuscular disorders and muscle degeneration.

Previous biological studies suggest that NMJ-derived signaling molecules play a key role in regulating nuclear dynamics, but the mechanisms by which these signals affect nuclear spacing remain unclear. This project uses **computational modeling** to explore how signal concentration modulates repulsive forces between nuclei, leading to self-organized nuclear distributions.

## Model Overview

The simulation models a **1D muscle cell** containing multiple nuclei and diffusing signaling molecules:

- Nuclei undergo stochastic motion modeled via diffusion  
- Signaling molecules are secreted from a defined NMJ region  
- Signal concentration locally modulates inter-nuclear repulsive forces  
- Higher signal concentration reduces repulsion, allowing nuclei to cluster  
- Boundary conditions enforce confinement within the cell  

The system evolves over time to study steady-state nuclear distributions.

## Key Features of the Simulation

- Random diffusion of nuclei and signaling molecules  
- Signal secretion and probabilistic degradation  
- Signal-dependent, distance-based repulsive forces  
- Boundary constraints to model physical cell limits  
- Time-series tracking of nuclear positions  
- Visualization of nuclear trajectories and spacing  

## Parameters Modeled

- Number of nuclei  
- Cell length  
- Time step size and simulation duration  
- Diffusion coefficients for nuclei and signaling molecules  
- Signal secretion and degradation rates  
- Interaction radius and repulsion strength  

All parameters are explicitly defined and can be adjusted to explore different biological regimes.

## Output and Analysis

- Time-series trajectories of nuclear positions  
- Plots showing nuclear diffusion and clustering over time  
- Analysis of spacing between neighboring nuclei  
- Correlation between signaling molecule density and nuclear proximity  

Simulation outputs were used to generate figures and statistical analyses for the accompanying research paper.

## Technologies Used

- **MATLAB (R2023a)**
  - Stochastic simulation
  - Numerical modeling
  - Data visualization

## Skills Demonstrated

- Computational modeling of biological systems  
- Stochastic processes and diffusion models  
- Simulation-based hypothesis testing  
- Translating biological mechanisms into mathematical rules  
- Scientific computing and data visualization  
- Research-grade code development  

## Academic Context

This work was conducted as part of **NYU GSTEM 2024** and contributed to a formal research paper titled:

**“Mathematical Modeling of Multinuclear Scaling in Muscle Cells”**  
Michelle Lu, Dr. Alex Mogilner  
New York University Center for Data Science

The full paper details the biological background, methodology, results, and implications of the model.

## Notes

- This repository contains research and exploratory code rather than production software  
- The model is intentionally simplified to isolate key biological mechanisms  
- Future extensions could include active transport, variable diffusion, or higher-dimensional models  

## Acknowledgements

Special thanks to **Dr. Alex Mogilner** for mentorship and guidance, and to the **NYU GSTEM program** for the opportunity to conduct this research.
