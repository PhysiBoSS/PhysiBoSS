/*
###############################################################################
# If you use PhysiCell in your project, please cite PhysiCell and the version #
# number, such as below:                                                      #
#                                                                             #
# We implemented and solved the model using PhysiCell (Version x.y.z) [1].    #
#                                                                             #
# [1] A Ghaffarizadeh, R Heiland, SH Friedman, SM Mumenthaler, and P Macklin, #
#     PhysiCell: an Open Source Physics-Based Cell Simulator for Multicellu-  #
#     lar Systems, PLoS Comput. Biol. 14(2): e1005991, 2018                   #
#     DOI: 10.1371/journal.pcbi.1005991                                       #
#                                                                             #
# See VERSION.txt or call get_PhysiCell_version() to get the current version  #
#     x.y.z. Call display_citations() to get detailed information on all cite-#
#     able software used in your PhysiCell application.                       #
#                                                                             #
# Because PhysiCell extensively uses BioFVM, we suggest you also cite BioFVM  #
#     as below:                                                               #
#                                                                             #
# We implemented and solved the model using PhysiCell (Version x.y.z) [1],    #
# with BioFVM [2] to solve the transport equations.                           #
#                                                                             #
# [1] A Ghaffarizadeh, R Heiland, SH Friedman, SM Mumenthaler, and P Macklin, #
#     PhysiCell: an Open Source Physics-Based Cell Simulator for Multicellu-  #
#     lar Systems, PLoS Comput. Biol. 14(2): e1005991, 2018                   #
#     DOI: 10.1371/journal.pcbi.1005991                                       #
#                                                                             #
# [2] A Ghaffarizadeh, SH Friedman, and P Macklin, BioFVM: an efficient para- #
#     llelized diffusive transport solver for 3-D biological simulations,     #
#     Bioinformatics 32(8): 1256-8, 2016. DOI: 10.1093/bioinformatics/btv730  #
#                                                                             #
###############################################################################
#                                                                             #
# BSD 3-Clause License (see https://opensource.org/licenses/BSD-3-Clause)     #
#                                                                             #
# Copyright (c) 2015-2018, Paul Macklin and the PhysiCell Project             #
# All rights reserved.                                                        #
#                                                                             #
# Redistribution and use in source and binary forms, with or without          #
# modification, are permitted provided that the following conditions are met: #
#                                                                             #
# 1. Redistributions of source code must retain the above copyright notice,   #
# this list of conditions and the following disclaimer.                       #
#                                                                             #
# 2. Redistributions in binary form must reproduce the above copyright        #
# notice, this list of conditions and the following disclaimer in the         #
# documentation and/or other materials provided with the distribution.        #
#                                                                             #
# 3. Neither the name of the copyright holder nor the names of its            #
# contributors may be used to endorse or promote products derived from this   #
# software without specific prior written permission.                         #
#                                                                             #
# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" #
# AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE   #
# IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE  #
# ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE   #
# LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR         #
# CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF        #
# SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS    #
# INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN     #
# CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)     #
# ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE  #
# POSSIBILITY OF SUCH DAMAGE.                                                 #
#                                                                             #
###############################################################################
*/

#include "custom.h"
#include "../BioFVM/BioFVM.h"  
#include "../addons/PhysiBoSS/src/maboss_intracellular.h"
using namespace BioFVM;

// declare cell definitions here 

std::vector<bool> nodes;

void create_cell_types( void )
{
	// set the random seed 
	if (parameters.ints.find_index("random_seed") != -1)
	{
		SeedRandom(parameters.ints("random_seed"));
	}
	
	/* 
	   Put any modifications to default cell definition here if you 
	   want to have "inherited" by other cell types. 
	   
	   This is a good place to set default functions. 
	*/ 

	initialize_default_cell_definition(); 
	cell_defaults.phenotype.secretion.sync_to_microenvironment( &microenvironment ); 

	cell_defaults.functions.volume_update_function = standard_volume_update_function;
	cell_defaults.functions.update_velocity = NULL;
	cell_defaults.functions.update_phenotype = NULL; 
	cell_defaults.functions.update_migration_bias = NULL; 
	cell_defaults.functions.pre_update_intracellular = pre_update_intracellular; 
	cell_defaults.functions.post_update_intracellular = post_update_intracellular; 
	cell_defaults.functions.custom_cell_rule = NULL; 
	
	cell_defaults.functions.add_cell_basement_membrane_interactions = NULL; 
	cell_defaults.functions.calculate_distance_to_membrane = NULL; 
	
	cell_defaults.custom_data.add_variable(parameters.strings("node_to_visualize"), "dimensionless", 0.0 ); //for paraview visualization

	/*
	   This parses the cell definitions in the XML config file. 
	*/
	
	initialize_cell_definitions_from_pugixml(); 
	
	/* 
	   Put any modifications to individual cell definitions here. 
	   
	   This is a good place to set custom functions. 
	*/ 

	get_cell_definition("default").functions.update_phenotype = wt_phenotype; 

	// Apply the exit function to specific cell types
    std::vector<std::string> cell_types_with_mutation = {"default"};
    
    for (const auto& cell_type : cell_types_with_mutation)
    {
        Cell_Definition* pCD = find_cell_definition(cell_type);
        if(pCD)
        {
            Cycle_Model* cycle_model = &(pCD->phenotype.cycle.model());
            int live_phase_index = cycle_model->find_phase_index(PhysiCell_constants::live);
            cycle_model->phase_link(live_phase_index, live_phase_index).exit_function = phase_exit_mutation_function;
            
            std::cout << "Applied mutation exit function to " << cell_type << " cells." << std::endl;
        }
        else
        {
            std::cout << "Warning: Cell type " << cell_type << " not found." << std::endl;
        }
    }
	
	/*
	   This builds the map of cell definitions and summarizes the setup. 
	*/

	build_cell_definitions_maps(); 

	/*
	   This intializes cell signal and response dictionaries 
	*/

	setup_signal_behavior_dictionaries();

	/*
	   This summarizes the setup. 
	*/
	
	display_cell_definitions( std::cout ); 


	return; 
}

void setup_microenvironment( void )
{
	// set domain parameters 
	
	// put any custom code to set non-homogeneous initial conditions or 
	// extra Dirichlet nodes here. 
	
	// initialize BioFVM 
	
	initialize_microenvironment(); 	
	
	return; 
}

void setup_tissue( void )
{
	// load cells from your CSV file
	load_cells_from_pugixml(); 	
}

void pre_update_intracellular( Cell* pCell, Phenotype& phenotype, double dt )
{
	if (PhysiCell::PhysiCell_globals.current_time >= 100.0 
		&& pCell->phenotype.intracellular->get_parameter_value("$time_scale") == 0.0
	){
		pCell->phenotype.intracellular->set_parameter_value("$time_scale", 0.1);
	}

}



void post_update_intracellular( Cell* pCell, Phenotype& phenotype, double dt )
{
	color_node(pCell);

	// std::cout << "Cell " << pCell->ID << " generation " << pCell->generation << " parent " << pCell->parent_ID << std::endl;

	// ADD NODE C RESISTANCE
	update_cell_from_boolean_model(pCell, phenotype, dt);
	
}

void update_cell_from_boolean_model(Cell* pCell, Phenotype& phenotype, double dt)
{
	static int death_decay_idx = pCell->custom_data.find_variable_index( "death_commitment_decay" );
	static int apoptosis_index = phenotype.death.find_death_model_index( PhysiCell_constants::apoptosis_death_model );
	static float apoptosis_rate = pCell->custom_data["apoptosis_rate"];
	static float death_commitment_decay = pCell->custom_data["death_decay_idx"];
	bool C_node = pCell->phenotype.intracellular->get_boolean_variable_value( "C" );

	static int density_idx = microenvironment.find_density_index("drug");
	double drug_density_ext = pCell->nearest_density_vector()[density_idx]; // A density (mM)
	double density_int = pCell->phenotype.molecular.internalized_total_substrates[density_idx];
	// density_int /= cell_volume; // divide int tot substrate to get concentration


	float basal_apoptosis = 1e-06;
	float final_apoptosis_rate = basal_apoptosis + apoptosis_rate*drug_density_ext;



	if ( C_node )
	{
		pCell-> phenotype.death.rates[apoptosis_index] = basal_apoptosis;

	} else {

		pCell-> phenotype.death.rates[apoptosis_index] = final_apoptosis_rate;
	}


	return;
}

std::vector<std::string> my_coloring_function( Cell* pCell )
{
	std::vector< std::string > output( 4 , "rgb(0,0,0)" );
	
	if ( !pCell->phenotype.intracellular->get_boolean_variable_value( parameters.strings("node_to_visualize") ) )
	{
		output[0] = "rgb(255,0,0)";
		output[2] = "rgb(125,0,0)";
		
	}
	else{
		output[0] = "rgb(0, 255,0)";
		output[2] = "rgb(0, 125,0)";
	}
	
	return output;
}

void color_node(Cell* pCell){
	std::string node_name = parameters.strings("node_to_visualize");
	pCell->custom_data[node_name] = pCell->phenotype.intracellular->get_boolean_variable_value(node_name);
}

void wt_phenotype( Cell* pCell, Phenotype& phenotype, double dt )
{ 
	// Default phenotype, doesn't do anything special

	if( phenotype.death.dead == true )
	{
		pCell->functions.update_phenotype = NULL; 
		return; 
	}

	// Add info on parent cell ID and generation
	// change_custom_data_var(pCell, "parent_ID", pCell->parent_ID);
	// change_custom_data_var(pCell, "generation", pCell->generation);

	// First check O2 availability
	// update_cell_and_death_parameters_O2_based(pCell, phenotype, dt);

	return; 
}

void add_custom_cycle_function()
{

	for( int i=0; i < all_cells->size() ; i++ )
	{
		// to call each cell, the pointer is *all_cells)[i] instead of pCell
		(*all_cells)[i]->phenotype.cycle.model().phases[0].entry_function = phase_exit_mutation_function;
	}
	
	return;
}


void phase_exit_mutation_function( Cell* pCell, Phenotype& phenotype, double dt )
{

	// Choose a random node
	static std::default_random_engine generator;
    
	
	// Only proceed if the cell has a MaBoSS model
    if (pCell->phenotype.intracellular &&
        pCell->phenotype.intracellular->intracellular_type == "maboss")
    {
        MaBoSSIntracellular* maboss_model = static_cast<MaBoSSIntracellular*>(pCell->phenotype.intracellular);
        std::vector<std::string> node_names = maboss_model->maboss.get_all_node_names();

        if (!node_names.empty())
        {			
			if (uniform_random() < parameters.doubles("mutation_rate_threshold")){

				std::uniform_int_distribution<size_t> node_dist(0, node_names.size() - 1);
				size_t node_idx = node_dist(generator);
				std::string node_name = node_names[node_idx];

				if (node_name == "A" || node_name == "B"){

					// Flip its value
					bool current_value = maboss_model->maboss.get_node_value(node_name);
					bool new_value = !current_value;
					maboss_model->maboss.set_node_value(node_name, new_value);

					// Record the mutation in the cell's mutations vector
					std::string mutation_record = node_name + "_" + std::to_string(new_value);
					pCell->custom_data.mutations.push_back(mutation_record);

					// std::cout << "Cell " << pCell->ID << ": flipped node " << node_name
					//           << " from " << current_value << " to " << new_value << std::endl;

				}
			}

        }
        else
        {
            std::cout << "Cell " << pCell->ID << " has no nodes in its MaBoSS model!" << std::endl;
        }
    }
    else
    {
        std::cout << "Cell " << pCell->ID << " does not have a MaBoSS model!" << std::endl;
    }
}


double get_custom_data_variable(Cell* pCell, std::string variable_name){
	int tmp_variable_idx = pCell->custom_data.find_variable_index(variable_name);
	double tmp_variable_value = pCell->custom_data[tmp_variable_idx];
	return tmp_variable_value;
}

void change_custom_data_var(Cell* pCell, std::string variable_name, double variable_new_value){
	int tmp_variable_idx = pCell->custom_data.find_variable_index(variable_name);
	pCell->custom_data[tmp_variable_idx] = variable_new_value;
	return;
}