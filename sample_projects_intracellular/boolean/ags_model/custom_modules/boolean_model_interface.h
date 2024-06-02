/*
 * ags_boolean_model_interface.cpp
 *
 *  Created on: 15 jun. 2020
 *  Author: Miguel Ponce-de-Leon (miguel.ponce@bsc.es)
 *  Contributor: Gerard Pradas
 *  Contributor: Arnau Montagud
 *  Contributor: Thalia Diniaco
 *  Cite as: arXiv:2103.14132 [q-bio.QM]
 *  Description: 
 *      Submodel that work as an interface 
 *      between the Boolean Network (BN) and PhysiCell (PC). 
 *      The submodel run the following steps:
 *      1- updates BN input nodes based on custom cell variables (see receptor model)
 *      2- updates the BN intracellular model by running MaBoSS
 *      3- updates cell state/behaviour based on the state of the BN readout nodes
 *  
 *      The update_monitor_variables funtion is just used to keep track of some relevand
 *      BN nodes' values that are stored as custom variables
 */

#include "../core/PhysiCell.h"
#include "../modules/PhysiCell_standard_modules.h" 

using namespace BioFVM; 
using namespace PhysiCell;

#include "./submodel_data_structures.h" 

void boolean_model_interface_setup();

void update_boolean_model_inputs( Cell* pCell, Phenotype& phenotype, double dt );

void update_cell_from_boolean_model(Cell* pCell, Phenotype& phenotype, double dt);

void ags_bm_interface_main(Cell* pCell, Phenotype& phenotype, double dt);


// Transfer functions
void anti_node_mapping_function( Cell* pCell, std::string drug_name, std::string target_node, double drug_half_max, double drug_Hill_coeff);
void add_noise_to_rates(Cell* pCell);

bool boolean_node_deactivation_prob(double drug_density, double scaling, double GI50 );

double growth_mapping_logistic(double doubling_time, double hill_coeff, double K_half, double S_value);
double apoptosis_mapping_logistic(double basal_apoptosis_rate, double maximum_apoptosis_rate, double hill_coeff, double K_half, double S_value);
double pressure_effect_growth_rate(double pressure, double hill_coeff, double pressure_half);
