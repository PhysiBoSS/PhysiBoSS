/*
 * ags_boolean_model_interface.cpp
 */

#include "../core/PhysiCell.h"
#include "../modules/PhysiCell_standard_modules.h" 
#include "./drug_transport_model.h"

using namespace BioFVM; 
using namespace PhysiCell;

#include "./submodel_data_structures.h" 

void boolean_model_interface_setup();
void update_boolean_model_inputs( Cell* pCell, Phenotype& phenotype, double dt );
void update_cell_from_boolean_model(Cell* pCell, Phenotype& phenotype, double dt);

void ags_bm_interface_main(Cell* pCell, Phenotype& phenotype, double dt); 

void pre_update_intracellular_ags(Cell* pCell, Phenotype& phenotype, double dt);
void post_update_intracellular_ags(Cell* pCell, Phenotype& phenotype, double dt);

std::string get_drug_target(std::string drug_name);

double get_boolean_prosurvival_outputs(Cell* pCell, Phenotype& phenotype);
double get_boolean_antisurvival_outputs(Cell* pCell, Phenotype& phenotype);