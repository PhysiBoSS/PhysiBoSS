/*
 * tnf_boolean_model_interface.cpp
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
#include "../addons/PhysiBoSS/src/maboss_intracellular.h"
#include "../addons/PhysiBoSS/src/maboss_network.h" 

using namespace BioFVM; 
using namespace PhysiCell;

#include "./submodel_data_structures.h" 

void tnf_boolean_model_interface_setup();

void update_boolean_model_inputs( Cell* pCell, Phenotype& phenotype, double dt );

void update_cell_from_boolean_model(Cell* pCell, Phenotype& phenotype, double dt);

void update_phenotype_with_signaling(Cell* pCell, Phenotype& phenotype, double dt);

// helper function to keep updated some cell custom variables
void update_monitor_variables( Cell* pCell );

double pressure_effect_growth_rate(double pressure, double hill_coeff, double pressure_half);

void update_cell_growth_parameters_pressure_based( Cell* pCell, Phenotype& phenotype, double dt );

double Hill_response_function_2( double s, double half_max , double hill_power );

void anti_node_mapping_function( Cell* pCell, std::string drug_name, std::string target_node, double drug_half_max, double drug_Hill_coeff);


// drug transport model

double calculate_diffusion_flux(Cell* pCell, int density_idx, double permeability,std::string drug_name);

// void drug_transport_model_setup();

void drug_transport_model_update( Cell* pCell, Phenotype& phenotype, double dt );

void drug_transport_model_main( double dt );


void pre_update_intracellular(Cell* pCell, Phenotype& phenotype, double dt);
void post_update_intracellular(Cell* pCell, Phenotype& phenotype, double dt);