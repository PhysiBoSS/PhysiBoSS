# PhysiBoSS Gastric Adenocarcinoma (AGS) model 
## Developed by Othmane Hayoun-Mya, Miguel Ponce de León and Arnau Montagud 
### PhysiBoSS version X.XX and PhysiCell version X.XX

This repository contains the PhysiBoSS model built for studying drug synergies in a gastric adenocarcinoma cell line AGS. It integrates the Boolean model built by Flobak et al. (2015) into PhysiCell.

## ODD Protocol description of AGS PhysiBoSS model

For further replicability and standarization of our multiscale ABM, we provide also an ODD description of our model, from the perspective of Pattern-Oriented Modelling (POM).

### Purpose
The model simulates the behavior of gastric adenocarcinoma cells (AGS cell line) in response to various drug treatments (single drugs and their combinations), incorporating intracellular signaling dynamics through Boolean networks into a multiscale agent-based framework, PhysiBoSS.

### Entities, state variables, and scales
- **Entities**: AGS cells (agents) and microenvironment (3D mesh)
- **State variables:**
  - **Cells**: position, volume, cell cycle phase, Boolean network state, pressure
  - **Microenvironment**: Oxygen and drug concentrations
- Spatial scale: 600 × 600 × 80 μm³ domain
- Temporal scale: 4200 simulation minutes

### Process overview and scheduling
1. Microenvironment updates (diffusion of substrates)
2. Cell mechanics updates
3. Cell phenotype updates
4. Boolean network updates
5. Data collection (every 40 minutes)

### Design concepts

#### Basic principles
The model integrates cellular behavior with intracellular signaling using Boolean networks within each agent.

#### Emergence
Cell population dynamics and drug response emerge from individual cell behaviors and intracellular signaling.

#### Adaptation
Cells adapt their behavior based on their microenvironment and internal signaling state.

#### Objectives
Cells aim to survive and proliferate based on their internal signaling state.

#### Sensing
Cells sense drug concentrations in their local microenvironment.

#### Interaction
Cells interact indirectly through substrate consumption and directly through physical contact.

#### Stochasticity
Initial Boolean network states and certain cellular processes incorporate stochasticity.

#### Observation
Cell counts, positions, and phenotypes are recorded every 40 simulation minutes.

### Details

#### Initialization
- Domain: 600 × 600 × 80 μm³, divided into 20 × 20 × 20 μm³ voxels
- Initial cell configuration: 2D disk of 303 μm radius, 2.8 μm cell spacing
- Cell volume: 2495 μm³
- Boolean network initial states: Probability of OFF state is 1 for all nodes except AKT, MEK, PI3K, TAK1, betacatenin, GSK3, and p38alpha, which is at 0.5 probability of ON and 0.5 probability of OFF.

#### Input data
Drug concentrations based on experimental GI₅₀ values.

#### Submodels

1. **Microenvironment**
   - Substrates: oxygen, PI103, PD0325901, AKTi-1,2
   - Oxygen diffusion: 10000 μm²/min
   - Drug diffusion: 600 μm²/min
   - Drug administration: at 1280 minutes, 12-minute pulse

2. **Cell Cycle**
   - Single-phase cycle, duration of 20-24 hours

3. **Cell Death**
   - Apoptosis rate: 5.316 × 10⁻⁵ min⁻¹ (default)
   - Maximum apoptosis rate: 1 × 10⁻⁴ to 1 × 10⁻³ min⁻¹

4. **Boolean network**
   - 75 components representing gastric adenocarcinoma pathways
   - Update frequency: every 10 simulation minutes
   - Steps per update: 50
   - Scaling: 1
   - Inheritance: activated
   - Readouts: FOXO, Caspase8 and Caspase9 nodes for apoptosis; cMYC, TCF and RSK for growth-promoting nodes.

5. **Cell Behavior**
   - Growth and death rates modulated by Boolean network outputs
   - Specific equations using Hill functions (H: Hill index, kₐ half-max value)