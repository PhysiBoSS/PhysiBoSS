<img align="right" height="200" src="https://github.com/sysbio-curie/MaBoSS/assets/22156824/b1b9a843-a203-42e6-8f8f-1acd2e2f5285">

# PhysiBoSS 2: a sustainable integration of stochastic Boolean and agent-based modelling frameworks
**Version:** 2.2.3

**Release date:** 13 December 2024

**PhysiCell base version:** 1.14.2

## Overview: 
PhysiBoSS 2.0 is a redesign and reimplementation of PhysiBoSS ([doi:10.1093/bioinformatics/bty766](https://doi.org/10.1093/bioinformatics/bty766)). It has been conceived as an add-on that expands the PhysiCell ([doi:10.1371/journal.pcbi.1005991](https://dx.doi.org/10.1371/journal.pcbi.1005991)) agent-based functionalities with intracellular cell signalling using MaBoSS having a decoupled, maintainable and model-agnostic design. PhysiBoSS 2.0 reproduces simulations reported in the original PhysiBoSS publications and can be used with other Boolean models, for instance to predict drug synergy in a gastric adenocarcinoma cell line.

**Reference paper:** [PhysiBoSS 2.0: a sustainable integration of stochastic Boolean and agent-based modelling frameworks](https://www.nature.com/articles/s41540-023-00314-4)

**Reference paper doi:** [10.1038/s41540-023-00314-4](https://doi.org/10.1038/s41540-023-00314-4)

A tutorial, explaining the new functionalities of PhysiBoSS, including the automated input/output mapping, and more, is available as a preprint. Check it out !

**Tutorial paper**: [https://arxiv.org/abs/2406.18371](https://arxiv.org/abs/2406.18371)

**Tutorial project folder**: [here](https://github.com/PhysiBoSS/PhysiBoSS/tree/master/sample_projects_intracellular/boolean/tutorial)

### How to run a PhysiBoSS sample_project inside PhysiCell:
~~~bash
git clone https://github.com/PhysiBoSS/PhysiBoSS.git
cd PhysiBoSS
make physiboss-tnf-model
make
./spheroid_TNF_model
~~~

### PhysiBoSS documentation:

The documentation of PhysiBoSS is available in the [PhysiBoSS User Guide](https://raw.githubusercontent.com/PhysiBoSS/PhysiBoSS/development/documentation/PhysiBoSS_User_Guide.pdf).


### Key makefile rules, from [PhysiCell repository](https://github.com/MathCancer/PhysiCell):

**`make`**: compiles the current project. If no 
                     project has been defined, it first 
                     populates the cancer heterogeneity 2D 
                     sample project and compiles it 
   
**`make project-name`**: populates the indicated sample project. 
                     Use "make" to compile it. 


   * **PhysiBoSS \[`project-name`\]** choices:
      * physiboss-cell-lines-sample
      * physiboss-tutorial
      * physiboss-tutorial-invasion
      * physiboss-tnf-model

   * **PhysiCell \[`project-name`\]** choices:
      * template 
      * biorobots-sample 
      * cancer-biorobots-sample 
      * cancer-immune-sample
      * celltypes3-sample 
      * heterogeneity-sample 
      * pred-prey-farmer 
      * virus-macrophage-sample 
      * worm-sample
      * ode-energy-sample 
      * cancer-metabolism-sample
      * interaction-sample
      * mechano-sample
      * rules-sample
      * physimess-sample
      * custom-division-sample
      * asymmetric-division-sample
      * immune-function-sample episode-sample

**`make list-projects`** : list all available sample projects 

**`make clean`**         : removes all .o files and the executable, so that the next "make" recompiles the entire project 

**`make data-cleanup`**  : clears out all simulation data 

**`make reset`**         : de-populates the sample project and returns to the original PhysiCell state. Use this when switching to a new PhysiCell sample project. 

**`make save PROJ=name`**: save the current project (including the `Makefile`, `main.cpp`, and everything in `./config` and `./custom_modules/`) in `./user_projects/name`, where `name` is your choice for the project. If the project already exists, overwrite it. 

**`make load PROJ=name`**: load the user project `name` from `./user_projects/name` (including the `Makefile`, `main.cpp`, and everything in `./config` and `./custom_modules/`).  

**`make list-user-projects`**: list all user projects in `./user_projects/`. (Use these names without the trailing `/` in `make load PROJ=name`.)

**`make jpeg`**          : uses ImageMagick to convert the SVG files in the output directory to JPG (with appropriate sizing to make movies). Supply `OUTPUT=foldername` to select a different folder. 

**`make movie`**         : uses ffmpeg to convert the JPG files in the output directory an mp4 movie. Supply `OUTPUT=foldername` to select a different folder, or `FRAMERATE=framerate` to override the frame rate.

## Legacy version

PhysiBoSS 1.0, as described in [PhysiBoSS: a multi-scale agent-based modelling framework integrating physical dimension and cell signalling](https://doi.org/10.1093/bioinformatics/bty766), is accessible at [https://github.com/PhysiBoSS/PhysiBoSSv1](https://github.com/PhysiBoSS/PhysiBoSSv1).

## Acknowledgements

This work has received funding from the Horizon 2020 projects INFORE (ID: 825070) and PerMedCoE (ID: 951773) and from the Horizon Europe project CREXDATA (ID: 101092749). This work was funded in part by the French government under the management of Agence Nationale de la Recherche as part of the “Investissements d’avenir” programme, reference ANR-19-P3IA-0001 (PRAIRIE 3IA Institute).

We thank Anne L'Hévéder for the PhysiBoSS logo.
