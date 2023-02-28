/*
 * env_boolean_model_interface.cpp
 *
 *  Created on: 26 feb. 2023
 *      Author: Arnau Montagud
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

// 1, seguir les pressions de les cel.lules
// 2, ixa pressió fa que es tanque els dirichlet (escaló amb Hill molt vertical)
// 3, ixe tancament farà que les cel.lules muiguen

// experiment de ficar sols unes cel.lules en cantó "fruit" de cirurgia

#include "./env_boolean_model_interface.h"

using namespace PhysiCell; 

Submodel_Information env_bm_interface_info;

void tnf_boolean_model_interface_setup()
{
    env_bm_interface_info.name = "Pressure Boolean model interface"; 
    env_bm_interface_info.version = "0.1.0";
	
    env_bm_interface_info.main_function= update_phenotype_with_signaling; 

    // These are just auxiliary variables to keep track of some BN nodes
    env_bm_interface_info.cell_variables.push_back( "Hypoxia_node" );
    env_bm_interface_info.cell_variables.push_back( "Proliferation_node" );
    env_bm_interface_info.cell_variables.push_back( "cell_pressure" );

    env_bm_interface_info.register_model();
}

void update_boolean_model_inputs( Cell* pCell, Phenotype& phenotype, double dt )
{
    // XXX inputs de Fumia: Acidosis, Hypoxia, Nutrients, Carcinogen, GFs, ROS

    // Hypoxia:
   	static int oxygen_substrate_index = pCell->get_microenvironment()->find_density_index( "oxygen" ); 
    // sample the microenvironment to get the pO2 value 
    double pO2 = (pCell->nearest_density_vector())[oxygen_substrate_index]; // PhysiCell_constants::oxygen_index]; 

    // static int nR_EB = pCell->custom_data.find_variable_index( "bound_external_TNFR" ); 
    // static int nTNF_threshold = pCell->custom_data.find_variable_index( "hypoxia_threshold" );
    // This if the step transfer function used to update the state of boolean model inputs
    // using the state of the receptor dynamics model. The continuos value thresholded is
    // the total TNF-recptor complex (doi:10.1016/j.cellsig.2010.08.016)

    if ( pO2 <= pCell->parameters.o2_hypoxic_threshold )
        pCell->phenotype.intracellular->set_boolean_variable_value("Hypoxia", 1);
    else
        pCell->phenotype.intracellular->set_boolean_variable_value("Hypoxia", 0);

    return;
}

void update_cell_from_boolean_model(Cell* pCell, Phenotype& phenotype, double dt)
{	
    // XXX outputs from Fumia: Proliferation, Apoptosis, GSH, GLUT1, VEGF, Lactic_acid, COX412, DNA_Repair

    // static int nTNF_external = microenvironment.find_density_index( "tnf" );      
    // static int nTNF_export_rate = pCell->custom_data.find_variable_index( "TNF_net_production_rate" );
    // static int necrosis_index = phenotype.death.find_death_model_index( PhysiCell_constants::necrosis_death_model ); 
    // static int cycle_start_index = live.find_phase_index( PhysiCell_constants::live ); 
    // static int cycle_end_index = live.find_phase_index( PhysiCell_constants::live );  
    // static float necrosis_rate = pCell->custom_data["necrosis_rate"];
    static int death_decay_idx = pCell->custom_data.find_variable_index( "death_commitment_decay" );
    static int apoptosis_index = phenotype.death.find_death_model_index( PhysiCell_constants::apoptosis_death_model );
    static float apoptosis_rate = pCell->custom_data["apoptosis_rate"];
    static float death_commitment_decay = pCell->custom_data["death_decay_idx"];
    
    // Getting the state of the Boolean model readouts (Readout can be in the XML)
    bool apoptosis = pCell->phenotype.intracellular->get_boolean_variable_value( "Apoptosis" );
    // bool nonACD = pCell->phenotype.intracellular->get_boolean_variable_value( "NonACD" );
    bool survival =pCell->phenotype.intracellular->get_boolean_variable_value( "Proliferation" );
    // bool NFkB = pCell->phenotype.intracellular->get_boolean_variable_value( "NFkB" );
	
    // la lleve perque tinc moooolta apoptosis
    // if ( apoptosis ) {
    //     phenotype.death.rates[apoptosis_index] = apoptosis_rate;
	// 	return;
	// } else {
    //     phenotype.death.rates[apoptosis_index] -= phenotype.death.rates[apoptosis_index] * death_commitment_decay;
    //     if (phenotype.death.rates[apoptosis_index] < 0)
    //         phenotype.death.rates[apoptosis_index] = 0;
    // }

    // if ( nonACD ) {
    //     // pCell->start_death(necrosis_index);
    //     phenotype.death.rates[necrosis_index] = necrosis_rate;
	// 	return;
	// }
    // else {
    //     phenotype.death.rates[necrosis_index] -= phenotype.death.rates[necrosis_index] * death_commitment_decay;
    //     if (phenotype.death.rates[necrosis_index] < 0)
    //         phenotype.death.rates[necrosis_index] = 0;
    // } 
    
    // if ( survival && pCell->phenotype.cycle.current_phase_index() == PhysiCell_constants::Ki67_negative ) 
	if ( survival ) 
    { 
        pCell->phenotype.cycle.advance_cycle(pCell, phenotype, dt); 
    }

    // If NFkB node is active produce some TNF
    // if ( NFkB )	
    // { 
    //     phenotype.secretion.net_export_rates[nTNF_external] = pCell->custom_data[nTNF_export_rate]; 
    // } else 
    // {
    //     phenotype.secretion.net_export_rates[nTNF_external] = 0;
    // }
    // update_cell_and_death_parameters_O2_based(pCell, phenotype, dt);
    return;
}

// custom cell phenotype function to run PhysiBoSS when is needed
void update_phenotype_with_signaling(Cell *pCell, Phenotype &phenotype, double dt)
{
	if( phenotype.death.dead == true )
	{
		pCell->functions.update_phenotype = NULL;
		return;
	}

    if ( pCell->phenotype.intracellular->need_update() )
    {
        // First we update the Boolean Model inputs
        update_boolean_model_inputs(pCell, phenotype, dt );
        
		// Run maboss to update the boolean state of the cell
        pCell->phenotype.intracellular->update();
        
		// update the cell fate based on the boolean outputs
        update_cell_from_boolean_model(pCell, phenotype, dt);

        // Get track of some boolean node values for debugging
        update_monitor_variables(pCell);
    }

    return;
}

void update_monitor_variables( Cell* pCell )
{
    // XXX find pressure of each cell
    
	// static int index_tnf_node  = pCell->custom_data.find_variable_index("tnf_node");
	// static int index_fadd_node = pCell->custom_data.find_variable_index("fadd_node");
	// static int index_nfkb_node = pCell->custom_data.find_variable_index("nfkb_node");
    // pCell->custom_data[index_tnf_node] = pCell->phenotype.intracellular->get_boolean_variable_value("TNF");
	// pCell->custom_data[index_fadd_node] = pCell->phenotype.intracellular->get_boolean_variable_value("FADD");
    // pCell->custom_data[index_nfkb_node] = pCell->phenotype.intracellular->get_boolean_variable_value( "NFkB" ) ;

	static int index_hypoxia_node  = pCell->custom_data.find_variable_index("Hypoxia_node");
	static int index_prolif_node  = pCell->custom_data.find_variable_index("Proliferation_node");
    static int index_cell_pressure = pCell->custom_data.find_variable_index("cell_pressure");
    pCell->custom_data[index_hypoxia_node] = pCell->phenotype.intracellular->get_boolean_variable_value("Hypoxia");
    pCell->custom_data[index_prolif_node] = pCell->phenotype.intracellular->get_boolean_variable_value("Proliferation");
    pCell->custom_data[index_cell_pressure] = pCell->state.simple_pressure; 
    //XXX pCell->phenotype.intracellular->get_boolean_variable_value("TNF");

    return;
}

double pressure_effect_growth_rate(double pressure, double hill_coeff, double pressure_half){

    // Suggestion: Employ the Hill_effect function from PhysiCell instead of this one, just to align with the other mapping functions

    // double pressure_exponential_function = std::pow(6e-03, pressure);
    double pressure_exponential_function =  std::pow(pressure, hill_coeff) / (pressure_half + std::pow(pressure, hill_coeff));
    // if (pressure_exponential_function > 1) pressure_exponential_function = 1.0;
    return pressure_exponential_function;
}

