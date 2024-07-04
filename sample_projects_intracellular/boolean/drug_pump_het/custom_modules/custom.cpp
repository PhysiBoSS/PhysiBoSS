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
# Copyright (c) 2015-2021, Paul Macklin and the PhysiCell Project             #
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

#include "./custom.h"
#include "./drug_transport_model.h"

void create_cell_types( void )
{
	// set the random seed 
	SeedRandom( parameters.ints("random_seed") );  
	
	/* 
	   Put any modifications to default cell definition here if you 
	   want to have "inherited" by other cell types. 
	   
	   This is a good place to set default functions. 
	*/ 
	
	initialize_default_cell_definition(); 
	cell_defaults.phenotype.secretion.sync_to_microenvironment( &microenvironment ); 
	
	cell_defaults.functions.volume_update_function = standard_volume_update_function;
	cell_defaults.functions.update_velocity = standard_update_cell_velocity;

	cell_defaults.functions.update_migration_bias = NULL; 
	cell_defaults.functions.update_phenotype = NULL; // update_cell_and_death_parameters_O2_based; 
	cell_defaults.functions.custom_cell_rule = NULL; 
	cell_defaults.functions.contact_function = NULL; 
	
	cell_defaults.functions.add_cell_basement_membrane_interactions = NULL; 
	cell_defaults.functions.calculate_distance_to_membrane = NULL;

	// @oth: TO ADD here

	cell_defaults.functions.pre_update_intracellular =  NULL;
	cell_defaults.functions.post_update_intracellular = NULL;
	
	/*
	   This parses the cell definitions in the XML config file. 
	*/
	
	initialize_cell_definitions_from_pugixml(); 

	/*
	   This builds the map of cell definitions and summarizes the setup. 
	*/
		
	build_cell_definitions_maps(); 

	/*
	   This intializes cell signal and response dictionaries 
	*/

	setup_signal_behavior_dictionaries(); 	

	/*
       Cell rule definitions 
	*/

	setup_cell_rules(); 

	/* 
	   Put any modifications to individual cell definitions here. 
	   
	   This is a good place to set custom functions. 
	*/ 
	
	cell_defaults.functions.update_phenotype = phenotype_function; 
	cell_defaults.functions.custom_cell_rule = custom_function; 
	cell_defaults.functions.contact_function = contact_function; 
	
	/*
	   This builds the map of cell definitions and summarizes the setup. 
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
	
	// load cells from your CSV file (if enabled)
	load_cells_from_pugixml();
	setup_cells_disk(); // add cells based on user params argument

	add_lognorm_distro( "pump_activated" );

	// Add distro of specific variables
	add_param_from_distro( "mutation_rate" );




	return; 
}

void setup_cells_disk( void )
{
	double Xmin = microenvironment.mesh.bounding_box[0]; 
	double Ymin = microenvironment.mesh.bounding_box[1]; 
	double Zmin = microenvironment.mesh.bounding_box[2]; 

	double Xmax = microenvironment.mesh.bounding_box[3]; 
	double Ymax = microenvironment.mesh.bounding_box[4]; 
	double Zmax = microenvironment.mesh.bounding_box[5]; 
	
	if( default_microenvironment_options.simulate_2D == true )
	{
		Zmin = 0.0; 
		Zmax = 0.0; 
	}
	
	double Xrange = Xmax - Xmin; 
	double Yrange = Ymax - Ymin; 
	double Zrange = Zmax - Zmin; 
	
	// create some of each type of cell 
	
	Cell* pC;
	
	for( int k=0; k < cell_definitions_by_index.size() ; k++ )
	{
		Cell_Definition* pCD = cell_definitions_by_index[k]; 
		std::cout << "Placing cells of type " << pCD->name << " ... " << std::endl; 
		for( int n = 0 ; n < parameters.ints("number_of_cells") ; n++ )
		{
			std::vector<double> position = {0,0,0}; 
			position[0] = Xmin + UniformRandom()*Xrange; 
			position[1] = Ymin + UniformRandom()*Yrange; 
			position[2] = Zmin + UniformRandom()*Zrange; 
			
			pC = create_cell( *pCD ); 
			pC->assign_position( position );
		}
	}
	std::cout << std::endl; 
}

std::vector<std::string> my_coloring_function( Cell* pCell )
{ return paint_by_number_cell_coloring(pCell); }

void phenotype_function( Cell* pCell, Phenotype& phenotype, double dt )
{ 
	if( phenotype.death.dead == true )
	{ return; }

	// First check O2 availability
	update_cell_and_death_parameters_O2_based(pCell, phenotype, dt);

	// Custom cycle exit in the phenotype
	PhysiCell::Cycle_Model cycle_model = cell_defaults.phenotype.cycle.model();

	// #TODO @oth: why are we using the same idx for the cell cycle? it should be automated to detect the start and end
	static int start_phase_idx = phenotype.cycle.model().find_phase_index( PhysiCell_constants::live );
	static int end_phase_idx = phenotype.cycle.model().find_phase_index( PhysiCell_constants::live );
	cycle_model.phase_link(start_phase_idx, end_phase_idx).exit_function = phase_exit_mutation_function;

	// Pressure affects negatively growth rate
	double growth_rate = phenotype.cycle.data.transition_rate(start_phase_idx, end_phase_idx);
	double pressure_scaling = get_pressure_scaling( pCell );
	growth_rate *= (1 - pressure_scaling);
	if (growth_rate < 0)
		growth_rate = 0;
	phenotype.cycle.data.transition_rate(start_phase_idx, end_phase_idx) = growth_rate;

	// Pump presence
	bool pump_state = static_cast<bool>(get_custom_data_variable(pCell, "pump_activated"));
	double pump_expression = get_custom_data_variable(pCell, "pump_expression_rate") * dt;
	if( uniform_random() < pump_expression && pump_state == true)
	{ 
		// Activate pump
		change_custom_data_var(pCell, "pump_activated", 1.0);
	}


	return; 

}

double get_pressure_scaling(Cell* pCell)
{
	// PRESSURE-BASED CONTACT-INHIBITION
	// @oth: #TODO Encapsulate this in a function? Or add it through PhysiCell rules?
	// Check relative pressure to eith er number of neighbor cells or set max logistic function to number of neighbor cells
	// pressure threshold set to 1, above this value there is no growth

	double cell_pressure = pCell->state.simple_pressure; 
	double hill_coeff_pressure = get_custom_data_variable(pCell, "hill_coeff_pressure");
	double pressure_half = get_custom_data_variable(pCell, "pressure_half");
	double scaling = Hill_response_function(cell_pressure, pressure_half, hill_coeff_pressure);
	// std::cout << "I am " << pCell->ID << " and my growth rate is " << growth_rate << std::endl;

	return scaling;
}


// @oth: TODO - for this function


void phase_exit_mutation_function( Cell* pCell, Phenotype& phenotype, double dt )
{
	// std::cout << "applying mutation at exit phase on cell " << pCell->ID << std::endl;


	bool pump_state = get_custom_data_variable(pCell, "pump_activated");
	double expression = get_custom_data_variable(pCell, "pump_expression_rate") * dt; // Gillespie algorithm
	double mutation = get_custom_data_variable(pCell, "mutation_rate") * dt;

	if(uniform_random() < mutation){ // Monte Carlo, flipping a coin to see if you mutate or not 

		double mutation_type = uniform_random();
		int apoptosis_model_index = phenotype.death.find_death_model_index( "Apoptosis" );


		if (mutation_type > 0.75){
			pCell-> phenotype.death.rates[apoptosis_model_index] *= 1.05;
			change_custom_data_var(pCell, "pump_expression_rate", expression*0.8);
			return;
			} // Deleterous
		
		if (mutation_type > 0.5 && mutation_type < 0.75){ 
			// Increase the rate of expression
			change_custom_data_var(pCell, "pump_expression_rate", expression*2.0);
			return;
			} // favourable
		
		return; // If neutral  

	}

	return;
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



void custom_function( Cell* pCell, Phenotype& phenotype , double dt )
{ return; } 

void contact_function( Cell* pMe, Phenotype& phenoMe , Cell* pOther, Phenotype& phenoOther , double dt )
{ return; } 

void treatment_function () 
{
	if (PhysiCell::parameters.bools.find_index("treatment") != -1) 
	{
		int treatment_substrate_index = BioFVM::microenvironment.find_density_index(PhysiCell::parameters.strings("treatment_substrate"));

		if (PhysiCell::parameters.bools("treatment")){
		
			if (
				(((int)PhysiCell::PhysiCell_globals.current_time) % PhysiCell::parameters.ints("treatment_period")) == 0 
				&& !BioFVM::microenvironment.get_substrate_dirichlet_activation(treatment_substrate_index)
			)
			{
				std::cout << PhysiCell::parameters.strings("treatment_substrate") << " activation at t=" << PhysiCell::PhysiCell_globals.current_time << std::endl;
				BioFVM::microenvironment.set_substrate_dirichlet_activation(treatment_substrate_index, true);	
			}

			if (
				(((int)PhysiCell::PhysiCell_globals.current_time) % PhysiCell::parameters.ints("treatment_period")) == PhysiCell::parameters.ints("treatment_duration") 
				&& BioFVM::microenvironment.get_substrate_dirichlet_activation(treatment_substrate_index)
			)
			{
				std::cout << PhysiCell::parameters.strings("treatment_substrate") << " inactivation at t=" << PhysiCell::PhysiCell_globals.current_time << std::endl;
				BioFVM::microenvironment.set_substrate_dirichlet_activation(treatment_substrate_index, false);	
			}
			
		} else if ( BioFVM::microenvironment.get_substrate_dirichlet_activation(treatment_substrate_index) ){
			std::cout << PhysiCell::parameters.strings("treatment_substrate") << " inactivation (NO TREATMENT) at t=" << PhysiCell::PhysiCell_globals.current_time << std::endl;
			BioFVM::microenvironment.set_substrate_dirichlet_activation(treatment_substrate_index, false);	
		}
	}
}


// AUXILIAR FUNCTIONS
std::string get_drug_target(std::string drug_name){
    std::string param_name = drug_name + "_target";
    std::string drug_target = parameters.strings(param_name);
    return drug_target;
}

void add_param_from_distro( std::string variable_name )
{
	std::string var_mean_name =  variable_name + "_mean";
	std::string var_sd_name =  variable_name + "_sd";

	double param_mean = parameters.doubles( var_mean_name );
	double param_sd = parameters.doubles( var_sd_name );

	// generate distro
	double param_distro = NormalRandom(param_mean, param_sd);

	// Iterate over all cells to add value from the distro to their custom data

	for( int i=0; i < all_cells->size() ; i++ )
	{
		// to call each cell, the pointer is *all_cells)[i] instead of pCell
		change_custom_data_var( (*all_cells)[i], variable_name, param_distro);
	}

}

void add_lognorm_distro( std::string param_name )
{
	// @oth: FOR NOW ONLY EMPLOYED FOR DISTRO OF THE PUMP in setup_tissue(), but can be extended to other params

	// Add log-normal distro for pump presence
	std::string param_name_sd = param_name + "_sd";
	double param_sd = parameters.doubles(param_name_sd);

	// Generate a random number from a distribution
	std::random_device rd;
    std::default_random_engine generator(rd());
	double rand_from_gaussian = NormalRandom(0.0, 3.0);
	double logsample = std::exp(rand_from_gaussian);
	bool pumpsate = static_cast<bool>(logsample); // need to booleanize based on distro
	// std::lognormal_distribution<double> distribution(0.0, param_name_sd); // Always mean = 0
	// double number = distribution(generator);

	
	for (int i = 0; i < (*all_cells).size(); i++)
	{
		// to call each cell, the pointer is *all_cells)[i] instead of pCell
		change_custom_data_var( (*all_cells)[i], param_name, pumpsate);
	}

	return;

}





// INTRACELLULAR FUNCTIONS -- @oth Not yet added, first start with just the PhysiCell model