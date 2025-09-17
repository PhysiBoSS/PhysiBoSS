# PhysiBoSS 2 - Mutation Tracking and Analysis Branch

**Branch:** `development_mut_oth`  
**Based on:** [PhysiBoSS](https://github.com/PhysiBoSS/PhysiBoSS) v2.2.3  
**PhysiCell base version:** 1.14.2  
**Author:** Othmane Hayou-Mya, Alejandro Madrid (Barcelona Supercomputing Center, BSC-CNS)  
**Purpose:** Extended PhysiBoSS with mutation tracking and lineage analysis

---

## Overview

This repository contains a specialized fork of PhysiBoSS that extends the original framework with comprehensive mutation tracking, lineage analysis, and cellular evolution capabilities. The modifications enable detailed study of cell population dynamics, mutation accumulation, and evolutionary trajectories in agent-based multicellular simulations.

## Key Features Added

### 1. **Mutation Tracking System**
- **Random Node Mutations**: Cells can randomly mutate Boolean network nodes during their lifecycle
- **Mutation Logging**: All mutations are tracked with detailed metadata (time, cell ID, generation, etc.)
- **Lineage Tracing**: Complete parent-child relationships are maintained across cell divisions
- **Export Capabilities**: Comprehensive mutation data export in CSV format for downstream analysis

### 2. **Enhanced Data Output**
- **Mutation CSV Files** (`.mut.csv`): Detailed mutation logs with timestamps and lineage information
- **Boolean Network State Tracking**: Real-time monitoring of intracellular Boolean network states
- **Plotting Scripts**: Custom Python scripts for mutation analysis and visualization
- **Extended MultiCellDS Output**: Enhanced data structure for comprehensive cell state storage

### 3. **Custom Analysis Tools**
- **Plotting Scripts**: Ready-to-use Python scripts for mutation pattern analysis
- **Statistical Analysis**: Tools for studying mutation rates, lineage trees, and population dynamics
- **Visualization**: Automated generation of mutation plots and evolutionary trees

---

## Repository Structure

### Core Modifications

#### **Core PhysiCell Extensions**
```
core/
├── PhysiCell_cell.h          # Extended Cell class with generation/parent tracking
├── PhysiCell_cell.cpp        # Modified cell division and initialization logic
└── PhysiCell_custom.h        # Custom cell behaviors and phenotypes
```

<!-- #### **MaBoSS Integration Enhancements**
```
addons/PhysiBoSS/
├── src/maboss_network.h      # Enhanced MaBoSS network interface
├── src/maboss_network.cpp    # Extended with mutation capabilities
└── MaBoSS/                   # Complete MaBoSS library with compiled binaries
```

#### **Custom Modules**
```
custom_modules/
├── custom.h                  # Mutation function declarations
└── custom.cpp               # Core mutation logic and phase exit functions
```

#### **Configuration Files**
```
config/
├── PhysiCell_settings.xml           # Main simulation configuration
├── cells.csv                        # Cell type definitions
├── model.cfg                        # MaBoSS model configuration
├── model_0.bnd, model_1.bnd        # Boolean network models
├── posiciones.csv                   # Initial cell positions
└── equally_distributed_*.csv       # Parameter distribution files
``` -->

### Analysis and Visualization

#### **Plotting Scripts**
```
plotting_scripts/
├── plot_mutations.py                # Individual mutation analysis
├── plot_mutations_aggregated.py     # Population-level mutation patterns
├── mutations.pdf/.png              # Generated mutation plots
├── mutations_aggregated*.pdf/.png   # Aggregated analysis plots
└── mutations_legend.pdf/.png       # Plot legends and annotations
```

#### **Simulation Output**
```
output/
├── output*.xml                      # Simulation state snapshots
├── output*_mut.csv                  # Mutation tracking data
├── output*_boolean_intracellular.csv # Boolean network states
├── snapshot*.svg                    # Visual cell snapshots
└── *.mat files                     # MATLAB-compatible data export
```

### Sample Projects
Extended the "physiboss_cell_lines" sample project with mutation tracking capabilities as 
a first toy model:


```
sample_projects_intracellular/boolean/physiboss_cell_lines/
├── config/                          # Project-specific configurations
├── custom_modules/                  # Project-specific mutation logic
└── Makefile                        # Build configuration
```

---

## Installation and Usage

### Prerequisites
- GCC compiler (C++11 or later)
- Make build system
- Python 3.x (for analysis scripts)
- Required Python packages: numpy, matplotlib, pandas

### Building the Project
```bash
# Clone this repository
git clone <your-fork-url>
cd mutPhysiBoSS

# Switch to development branch
git checkout development_mut_oth

# Build the project
make clean
make PhysiBoSS-cell-lines; make

# Run a simulation
./project
```

### Configuration
1. **Edit `config/PhysiCell_settings.xml`** for simulation parameters
2. **Modify `config/model.cfg`** for Boolean network settings
3. **Adjust `config/cells.csv`** for cell type definitions
4. **Update `custom_modules/custom.cpp`** for custom mutation logic

---

## Data Analysis

### Mutation Data Format
The mutation tracking system exports data in CSV format with the following structure:

```csv
time,ID,type,parent_ID,generation,mutations
120.0,42,1,17,2,NodeA_1.NodeB_0.NodeC_1
240.0,43,1,42,3,NodeA_1.NodeB_0.NodeC_1.NodeD_0
```

**Fields:**
- `time`: Simulation time when data was recorded
- `ID`: Unique cell identifier
- `type`: Cell type index
- `parent_ID`: ID of parent cell (-1 for initial cells)
- `generation`: Cell generation (0 for initial cells)
- `mutations`: Dot-separated list of mutations (format: `NodeName_NewValue`)

### Visualization Scripts
```bash
# Individual cell mutation analysis
python plotting_scripts/plot_mutations.py

# Population-level mutation patterns
python plotting_scripts/plot_mutations_aggregated.py
```

---

---

## Core Modifications Summary

### Cell Class Extensions
```cpp
class Cell {
    // Added members:
    int generation = 0;        // Cell generation counter
    int parent_ID = -1;        // Parent cell ID for lineage tracking
    std::vector<std::string> mutations;  // Mutation history
    
    // Enhanced division method with lineage tracking
    Cell* divide();
};
```

### Mutation Function
```cpp
void phase_exit_mutation_function(Cell* pCell, Phenotype& phenotype, double dt) {
    // Randomly select and mutate Boolean network nodes
    // Log mutations with metadata
    // Update cell state accordingly
}
```

### Enhanced Data Output
```cpp
// Extended CSV output with lineage and mutation data
void write_mutation_data(std::ostream& os, Cell* pCell, double current_time) {
    os << current_time << "," << pCell->ID << "," << pCell->type 
       << "," << pCell->parent_ID << "," << pCell->generation << ",";
    // Output mutations in dot-separated format
}
```

---

## Documentation

### Core Changes Documentation
See [`core_changes.md`](core_changes.md) for detailed technical documentation of all modifications made to the PhysiCell and PhysiBoSS codebase.

### Original PhysiBoSS Documentation
- [PhysiBoSS User Guide](https://raw.githubusercontent.com/PhysiBoSS/PhysiBoSS/development/documentation/PhysiBoSS_User_Guide.pdf)
- [Tutorial Paper](https://arxiv.org/abs/2406.18371)
- [Reference Paper](https://doi.org/10.1038/s41540-023-00314-4)

---

---

**Note**: This is a development branch with active modifications. For stable production use, consider using tagged releases.
