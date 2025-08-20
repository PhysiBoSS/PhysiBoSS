#include "../core/PhysiCell.h"
#include "../modules/PhysiCell_standard_modules.h" 
#include "../addons/PhysiBoSS/src/boolean_network.h"
#include "./drug_sensitivity.h"


/**
 *	\boolean network interface 
 *	\brief custom module for prostate example
 * 
 *	\details Modules needed for the prostate example. This custom module can be used to study the inhibition of prostate cell lines.
 *
 *
 *	\date 03/03/2021
 *	\author Annika Meert, BSC-CNS, with code previously developed by Arnau Montagud, Gerard Pradas and Miguel Ponce de Leon, BSC-CNS
 */

using namespace BioFVM; 
using namespace PhysiCell;

void update_custom_variables( Cell* pCell );

void set_boolean_node (Cell* pCell, std::string drug_name, int index, double threshold);
void set_input_nodes(Cell* pCell); 
void from_nodes_to_cell(Cell* pCell, Phenotype& phenotype, double dt);
void boolean_model_interface_main (Cell* pCell, Phenotype& phenotype, double dt);