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

#include "./custom.h"

// declare cell definitions here
void create_cell_types(void)
{
	// use the same random seed so that future experiments have the
	// same initial histogram of oncoprotein, even if threading means
	// that future division and other events are still not identical
	// for all runs

	SeedRandom(parameters.ints("random_seed")); // or specify a seed here

	// housekeeping
	std::cout << cell_defaults.name << std::endl;

	cell_defaults.functions.volume_update_function = standard_volume_update_function;
	cell_defaults.functions.update_velocity = standard_update_cell_velocity;
	cell_defaults.functions.update_migration_bias = NULL;
	cell_defaults.functions.custom_cell_rule = NULL;

	cell_defaults.custom_data.add_variable("time", "min", PhysiCell_globals.current_time );
	cell_defaults.custom_data.add_variable("total_live_cells", "dimensionless", 0.0 );

	cell_defaults.custom_data.add_variable("reactivation_rate", "1/min", 0.0 );
	cell_defaults.custom_data.add_variable("mutation_rate", "1/generation", 0.0 );
	
	cell_defaults.functions.add_cell_basement_membrane_interactions = NULL;
	cell_defaults.functions.calculate_distance_to_membrane = NULL;
	cell_defaults.functions.set_orientation = NULL;

	// cell_defaults.functions.update_phenotype = NULL;
	cell_defaults.functions.update_phenotype = tumor_cell_phenotype_with_signaling;

	initialize_cell_definitions_from_pugixml();

	drug_transport_model_setup();
	boolean_model_interface_setup();

	submodel_registry.display(std::cout);

	build_cell_definitions_maps();
	display_cell_definitions(std::cout);

	return;
}

void setup_microenvironment(void)
{
	// make sure to override and go back to 2D
	// if (default_microenvironment_options.simulate_2D == true)
	// {
	// 	std::cout << "Warning: overriding XML config option and setting to 3D!" << std::endl;
	// 	default_microenvironment_options.simulate_2D = false;
	// }

	// microenvironment.add_density( "drug_X", "mM" );
	// microenvironment.diffusion_coefficients[1] = 0.0; 
	// microenvironment.decay_rates[1] = 0.0;
	

	// initialize BioFVM
	initialize_microenvironment();



	return;
}


void setup_tissue( void )
{
	// place a cluster of tumor cells at the center 
	
	double cell_radius = cell_defaults.phenotype.geometry.radius; 
	double cell_spacing = parameters.doubles("cell_spacing") * cell_radius; 
	
	double tumor_radius = parameters.doubles( "tumor_radius" ); 
	// Parameter<double> temp; 
	// int i = parameters.doubles.find_index( "tumor_radius" ); 
	
	Cell* pCell = NULL; 

	std::cout << "about to call basic_2D_disk_setup" << std::endl;

	basic_2D_disk_setup(pCell, tumor_radius, cell_spacing);

	// PATHWAY-BASED RESISTANCE HETEROGENEITY
	// double het_mean = parameters.doubles("heterogeneity_mean");
	// double het_sd = parameters.doubles("heterogeneity_sd");
	// add_heterogeneity(pCell, het_mean, het_sd, tumor_radius, cell_spacing);

	return; 
}

void update_cell_gowth_parameters_pressure_based( Cell* pCell, Phenotype& phenotype, double dt ) 
{
	// supported cycle models:
		// advanced_Ki67_cycle_model= 0;
		// basic_Ki67_cycle_model=1
		// live_cells_cycle_model = 5; 
	
	if( phenotype.death.dead == true )
	{ return; }
	
	// set up shortcuts to find the Q and K(1) phases (assuming Ki67 basic or advanced model)
	static bool indices_initiated = false; 
	static int start_phase_index; // Q_phase_index; 
	static int end_phase_index; // K_phase_index;
	static int necrosis_index; 
	
	static int oxygen_substrate_index = pCell->get_microenvironment()->find_density_index( "oxygen" ); 
	
	if( indices_initiated == false )
	{
		// Ki67 models
		
		if( phenotype.cycle.model().code == PhysiCell_constants::advanced_Ki67_cycle_model || 
			phenotype.cycle.model().code == PhysiCell_constants::basic_Ki67_cycle_model )
		{
			start_phase_index = phenotype.cycle.model().find_phase_index( PhysiCell_constants::Ki67_negative );
			necrosis_index = phenotype.death.find_death_model_index( PhysiCell_constants::necrosis_death_model ); 
			
			if( phenotype.cycle.model().code == PhysiCell_constants::basic_Ki67_cycle_model )
			{
				end_phase_index = 
					phenotype.cycle.model().find_phase_index( PhysiCell_constants::Ki67_positive );
				indices_initiated = true; 
			}
			if( phenotype.cycle.model().code == PhysiCell_constants::advanced_Ki67_cycle_model )
			{
				end_phase_index = 
					phenotype.cycle.model().find_phase_index( PhysiCell_constants::Ki67_positive_premitotic );
				indices_initiated = true; 
			}
		}
		
		// live model 
			
		if( phenotype.cycle.model().code == PhysiCell_constants::live_cells_cycle_model )
		{
			start_phase_index = phenotype.cycle.model().find_phase_index( PhysiCell_constants::live );
			necrosis_index = phenotype.death.find_death_model_index( PhysiCell_constants::necrosis_death_model ); 
			end_phase_index = phenotype.cycle.model().find_phase_index( PhysiCell_constants::live );
			indices_initiated = true; 
		}
		
		// cytometry models 
		
		if( phenotype.cycle.model().code == PhysiCell_constants::flow_cytometry_cycle_model || 
			phenotype.cycle.model().code == PhysiCell_constants::flow_cytometry_separated_cycle_model )
		{
			start_phase_index = phenotype.cycle.model().find_phase_index( PhysiCell_constants::G0G1_phase );
			necrosis_index = phenotype.death.find_death_model_index( PhysiCell_constants::necrosis_death_model ); 
			end_phase_index = phenotype.cycle.model().find_phase_index( PhysiCell_constants::S_phase );
			indices_initiated = true; 
		}	

		if( phenotype.cycle.model().code == PhysiCell_constants::cycling_quiescent_model )
		{
			start_phase_index = phenotype.cycle.model().find_phase_index( PhysiCell_constants::quiescent );
			necrosis_index = phenotype.death.find_death_model_index( PhysiCell_constants::necrosis_death_model ); 
			end_phase_index = phenotype.cycle.model().find_phase_index( PhysiCell_constants::cycling );
			indices_initiated = true; 
		}
		
	}
	
	// don't continue if we never "figured out" the current cycle model. 
	if( indices_initiated == false )
	{
		return; 
	}

	// this multiplier is for linear interpolation of the oxygen value 
	double multiplier = 1.0;
	
	// now, update the appropriate cycle transition rate 


	// PRESSURE-BASED CONTACT-INHIBITION
	// Check relative pressure to either number of neighbor cells or set max logistic function to number of neighbor cells
	// pressure threshold set to 1, above this value there is no growth

	double p = pCell->state.simple_pressure; 
    double hill_coeff_pressure = parameters.doubles("hill_coeff_pressure");
    double pressure_half = parameters.doubles("pressure_half");
    double scaling = pressure_effect_growth_rate(p, hill_coeff_pressure, pressure_half );
	// std::cout << "scaling is: " << scaling << std::endl;

	double rate = phenotype.cycle.data.transition_rate(0, 0);
	rate *= (1 - scaling);
	if (rate < 0)
		rate = 0;
	
	phenotype.cycle.data.transition_rate(start_phase_index, end_phase_index) = rate;

	// adding noise after division

	std::string drug_X_target_node = parameters.strings("drug_X_target");
	// if (drug_X_target_node != "none"){ add_noise_after_division(pCell, phenotype, drug_X_target_node);}
	
	std::string drug_Y_target_node = parameters.strings("drug_Y_target");
	// if (drug_Y_target_node != "none"){ add_noise_after_division(pCell, phenotype, drug_Y_target_node);}

	// Adding a custom phase exit function -- the smooth option for adding a mutation
	phenotype.cycle.model().phase_link(start_phase_index, end_phase_index).exit_function = my_mutation_function;

}

void my_mutation_function( Cell* pCell, Phenotype& phenotype, double dt )
{
	
	// std::cout << "applying mutation at exit phase on cell " << pCell->ID << std::endl;

	static int reactivation_value_idx = pCell->custom_data.find_variable_index("reactivation_value");
	static int mutation_value_idx = pCell->custom_data.find_variable_index("mutation_value");
	static int apoptosis_model_index = phenotype.death.find_death_model_index( "Apoptosis" );

	double mutation_value = pCell->custom_data.variables[mutation_value_idx].value * dt; // Gillespie
	double reactivation_value = pCell->custom_data.variables[reactivation_value_idx].value * dt; // Gillespie

	// double mutation_prob_mean = parameters.doubles("mutation_prob_mean");
	// double mutation_prob_sd = parameters.doubles("mutation_prob_sd");
	// double mutation_value = NormalRandom(mutation_prob_mean, mutation_prob_sd );

	pCell->custom_data.variables[mutation_value_idx].value = mutation_value;

	if(uniform_random() < mutation_value){ // Monte Carlo 

		double mutation_type = uniform_random();

		if (mutation_value > 0.75){  
			pCell-> phenotype.death.rates[apoptosis_model_index] *= 1.05;
			return;
			} // Deleterous
		
		if (mutation_value > 0.5){ 
			pCell->custom_data.variables[reactivation_value_idx].value *= 2.0;
			return;
			} // favourable
		
		return; // If neutral 

	}

	return;
}


// void add_noise_after_division(Cell *pCell, Phenotype& phenotype, std::string drug_name){

// 	// ADDING NOISE TO SPECIFIC CUSTOM DATA AFTER EACH DIVISION
// 	// This is done through flagging agents with elapsed time in phase 0

// 	double time_in_phase = pCell->phenotype.cycle.data.elapsed_time_in_phase;
// 	static int reactivation_value_idx = pCell->custom_data.find_variable_index("reactivation_value");


// 	double mutation_prob_mean = parameters.doubles("mutation_prob_mean");
// 	double mutation_prob_sd = parameters.doubles("mutation_prob_sd");
// 	double mutation_rate = NormalRandom(mutation_prob_mean, mutation_prob_sd ) * dt; 

// 	// std::cout << time_in_phase << std::endl;

	
// 	// Mutate after division (Implementation method B)
// 	// if (( uniform_random() < pCell->custom_data.variables[reactivation_prob_idx].value ) & (drug_X_target_node != "none") )

// 	if (( time_in_phase == 0) & (uniform_random() < mutation_value)){ 

// 		// std::cout << "this cell " << pCell->ID << " has time elapsed in phase " << phenotype.cycle.data.elapsed_time_in_phase << std::endl;
// 		std::cout << "this cell " << pCell->ID << " is mutating " << std::endl;

// 		std::string drug_target_node = parameters.strings(drug_name);

// 		double mutation_type = uniform_random();

// 		if (mutation_value < 0.5){ return; } // Neutral mutation
// 		if ( (mutation_value > 0.5) & (mutation_value < 0.5) ){ pCell->phenotype.intracellular->set_boolean_variable_value(drug_target_node, 1);} // 
// 		if (mutation_value = 0){ return;} 

// 		// this allows for maintaining some sort of reactivation value
// 		// pCell->custom_data.variables[reactivation_value_idx].value += NormalRandom(reactivation_noise_mean, reactivation_noise_sd) * 0.1;
// 		// if (pCell->custom_data.variables[reactivation_value_idx].value < 0){ pCell->custom_data.variables[reactivation_value_idx].value = 0;}
// 		// if (pCell->custom_data.variables[reactivation_value_idx].value > 1){ pCell->custom_data.variables[reactivation_value_idx].value = 1;}


// 	}

// 	return;
// }


// custom cell phenotype function to run PhysiBoSS when is needed
void tumor_cell_phenotype_with_signaling(Cell *pCell, Phenotype &phenotype, double dt)
{
	
	// std::cout << "Choosing penotype with signalling " << std::endl;

	if (phenotype.death.dead == true)
	{
		pCell->functions.update_phenotype = NULL;
		return;
	}
	
	update_cell_and_death_parameters_O2_based(pCell, phenotype, dt);
	// drug_transport_model_main( dt );
	ags_bm_interface_main(pCell, phenotype, dt);
	update_cell_gowth_parameters_pressure_based(pCell, phenotype, dt);

	// Adding time to custom data
	static int time_index = pCell->custom_data.find_variable_index("time");
	pCell->custom_data.variables[time_index].value = PhysiCell_globals.current_time;

	// Adding total number of cells
	static int total_live_cells_index = pCell->custom_data.find_variable_index("total_live_cells");
	float cells = 0.0;

	#pragma omp parallel for 
	for( int i=0; i < (*all_cells).size() ; i++ )
	{
		Cell* pCell = (*all_cells)[i]; 
		if( pCell->phenotype.death.dead == true)
		{ continue; }

		cells += 1; 
	}

	pCell->custom_data.variables[total_live_cells_index].value = cells;
}


// cell coloring function for ploting the svg files
std::vector<std::string> my_coloring_function(Cell *pCell)
{
	// start with live coloring
	// std::vector<std::string> output = false_cell_coloring_live_dead(pCell);
	std::vector<std::string> output = simple_cell_coloring(pCell);
	static int reactivation_prob_idx = pCell->custom_data.find_variable_index("reactivation_rate");
	static int mutation_prob_idx = pCell->custom_data.find_variable_index("mutation_rate");

	if(pCell->custom_data.variables[reactivation_prob_idx].value > 1.0){
		pCell->custom_data.variables[reactivation_prob_idx].value == 1.0;
	}

	if(pCell->custom_data.variables[mutation_prob_idx].value > 1.0){
		pCell->custom_data.variables[mutation_prob_idx].value == 1.0;
	}

	int reactivation_prob = (int) round(255 * pCell->custom_data.variables[reactivation_prob_idx].value);
	int mutation_prob = (int) round(255 * pCell->custom_data.variables[mutation_prob_idx].value);

	char szTempString[128];
	sprintf(szTempString, "rgb(%u,%u,%u)", 0, reactivation_prob, 0); // green shade depicts reactivation probability
	output[0].assign(szTempString);

	sprintf(szTempString, "rgb(%u,%u,%u)", 0, 0, mutation_prob); // green shade depicts reactivation probability
	output[1].assign(szTempString);

	// sprintf( szTempString , "rgb(%u,%u,%u)", (int)round(output[0][0]/2.0) , (int)round(output[0][1]/2.0) , (int)round(output[0][2]/2.0) );
	output[2].assign(szTempString);
	output[3].assign(szTempString);



	// if not, dead colors 
	if (pCell->phenotype.cycle.current_phase().code == PhysiCell_constants::apoptotic )  // Apoptotic - Red
	{
		output[0] = "rgb(255,0,0)";
		output[2] = "rgb(125,0,0)";
	}
	
	// Necrotic - Brown
	if( pCell->phenotype.cycle.current_phase().code == PhysiCell_constants::necrotic_swelling || 
		pCell->phenotype.cycle.current_phase().code == PhysiCell_constants::necrotic_lysed || 
		pCell->phenotype.cycle.current_phase().code == PhysiCell_constants::necrotic )
	{
		output[0] = "rgb(250,138,38)";
		output[2] = "rgb(139,69,19)";
	}	

	return output;
}


// Funtion to read init files created with PhysiBoSSv2
std::vector<init_record> read_init_file(std::string filename, char delimiter, bool header)
{
	// File pointer
	std::fstream fin;
	std::vector<init_record> result;

	// Open an existing file
	fin.open(filename, std::ios::in);

	// Read the Data from the file
	// as String Vector
	std::vector<std::string> row;
	std::string line, word;

	if (header)
		getline(fin, line);

	do
	{
		row.clear();

		// read an entire row and
		// store it in a string variable 'line'
		getline(fin, line);

		// used for breaking words
		std::stringstream s(line);

		// read every column data of a row and
		// store it in a string variable, 'word'
		while (getline(s, word, delimiter))
		{

			// add all the column data
			// of a row to a vector
			row.push_back(word);
		}

		init_record record;
		record.x = std::stof(row[2]);
		record.y = std::stof(row[3]);
		record.radius = std::stof(row[5]);
		record.phase = std::stoi(row[13]);
		record.z = std::stof(row[4]);
		record.elapsed_time = std::stod(row[14]);

		result.push_back(record);
	} while (!fin.eof());

	return result;
}


void set_input_nodes(Cell* pCell) {}

void from_nodes_to_cell(Cell* pCell, Phenotype& phenotype, double dt) {}



void color_node(Cell* pCell){
	std::string node_name = parameters.strings("node_to_visualize");
	pCell->custom_data[node_name] = pCell->phenotype.intracellular->get_boolean_variable_value(node_name);
}

void inject_density_sphere(int density_index, double concentration, double membrane_lenght) // it's lengTH
{
	// Inject given concentration on the extremities only
	#pragma omp parallel for
	for (int n = 0; n < microenvironment.number_of_voxels(); n++)
	{
		auto current_voxel = microenvironment.voxels(n);
		std::vector<double> cent = {current_voxel.center[0], current_voxel.center[1], current_voxel.center[2]};
		
		microenvironment.density_vector(n)[density_index] = concentration;

		// if (current_voxel.center[2] >= 2)

		// std::cout << norm(cent) << std::endl;

		// if ((membrane_lenght - norm(cent)) <= 0)
		// 	microenvironment.density_vector(n)[density_index] = concentration;

	}
}

void remove_density(int density_index)
{
	for (int n = 0; n < microenvironment.number_of_voxels(); n++)
		microenvironment.density_vector(n)[density_index] = 0;
}



double total_live_cell_count()
{
	double out = 0.0;

	for (int i = 0; i < (*all_cells).size(); i++)
	{
		if ((*all_cells)[i]->phenotype.death.dead == false && (*all_cells)[i]->type == 0)
		{
			out += 1.0;
		}
	}

	return out;
}

double total_dead_cell_count()
{
	double out = 0.0;

	for (int i = 0; i < (*all_cells).size(); i++)
	{
		if ((*all_cells)[i]->phenotype.death.dead == true && (*all_cells)[i]->phenotype.death.current_death_model_index == 0)
		{
			out += 1.0;
		}
	}

	return out;
}

double total_necrosis_cell_count()
{
	double out = 0.0;

	for (int i = 0; i < (*all_cells).size(); i++)
	{
		if ((*all_cells)[i]->phenotype.death.dead == true && (*all_cells)[i]->phenotype.death.current_death_model_index == 1)
		{
			out += 1.0;
		}
	}

	return out;
}

void add_reactivation_prob(Cell* pCell, double reactivation_mean, double reactivation_sd){

	static int reactivation_prob_idx = pCell->custom_data.find_variable_index("reactivation_rate");
	double reactivation_rate_min = 0.0;
	// double resistance_max = 1.0;

	pCell->custom_data.variables[reactivation_prob_idx].value = NormalRandom( reactivation_mean, reactivation_sd );
	if( pCell->custom_data.variables[reactivation_prob_idx].value < reactivation_rate_min ){ pCell->custom_data.variables[reactivation_prob_idx].value = reactivation_rate_min; }

	std::cout << "This cell has " << pCell->custom_data.variables[reactivation_prob_idx].value << " probab. of reactivating different pathways" << std::endl;
}

void add_mutation_rate(Cell* pCell, double mutation_rate_mean, double mutation_rate_sd){

	static int mutation_prob_idx = pCell->custom_data.find_variable_index("mutation_rate");
	double mutation_rate_min = 0.0;
	double mutation_rate_max = 1.0;

	pCell->custom_data.variables[mutation_prob_idx].value = NormalRandom( mutation_rate_mean, mutation_rate_sd );
	if( pCell->custom_data.variables[mutation_prob_idx].value < mutation_rate_min ){ pCell->custom_data.variables[mutation_prob_idx].value = mutation_rate_min; }

	std::cout << "This cell has " << pCell->custom_data.variables[mutation_prob_idx].value << " mutation rate " << std::endl;
}


// void add_pump_heterogeneity(Cell* pCell, std::string variable_name, double variable_mean, double variable_sd, double variable_min, double variable_max){

// 	// This function adds, to each agent, a value for the Km, the k2 and the total amount of transporter (Enzyme),
// 	// drawn from a given distribution in the config XML.

// 	std::string variable_name_value = variable_name + "_value";

// 	static int variable_idx = pCell->custom_data.find_variable_index(variable_name_value);
	
// 	pCell->custom_data.variables[variable_idx].value = NormalRandom( variable_mean, variable_sd );
// 	if( pCell->custom_data.variables[variable_idx].value < variable_min){ pCell->custom_data.variables[variable_idx].value = variable_min; }
// 	if( pCell->custom_data.variables[variable_idx].value > variable_max ){ pCell->custom_data.variables[variable_idx].value = variable_max; }

// 	std::cout << "The variable " << variable_name << " has a prob of " << pCell->custom_data.variables[variable_idx].value << " in cell " << pCell->ID << std::endl;
// }

void basic_2D_disk_setup( Cell* pCell, double tumor_radius, double cell_spacing){

	// double pump_total_Enzyme_mean = parameters.doubles("pump_total_Enzyme_mean");
	// double pump_total_Enzyme_sd = parameters.doubles("pump_total_Enzyme_sd");
	// double pump_k_2_mean = parameters.doubles("pump_k_2_mean");
	// double pump_k_2_sd = parameters.doubles("pump_k_2_sd");
	// double pump_kM_mean = parameters.doubles("pump_kM_mean");
	// double pump_kM_sd = parameters.doubles("pump_kM_sd");

	// This basic setup also accounts for heterogeneity in the population.

	double reactivation_rate_mean = parameters.doubles("reactivation_rate_mean");
	double reactivation_rate_sd = parameters.doubles("reactivation_rate_sd");

	double mutation_rate_mean = parameters.doubles("mutation_rate_mean");
	double mutation_rate_sd = parameters.doubles("mutation_rate_sd");
	
	double x = 0.0; 
	double x_outer = tumor_radius; 
	double y = 0.0; 

	int n = 0; 
	while( y < tumor_radius )
	{
		x = 0.0; 
		if( n % 2 == 1 )
		{ x = 0.5*cell_spacing; }
		x_outer = sqrt( tumor_radius*tumor_radius - y*y ); 
		
		while( x < x_outer )
		{
			pCell = create_cell();
			pCell->assign_position( x , y , 0.0 );
			add_reactivation_prob(pCell, reactivation_rate_mean, reactivation_rate_sd);
			add_mutation_rate(pCell, mutation_rate_mean, mutation_rate_sd);

	
			if( fabs( y ) > 0.01 )
			{
				pCell = create_cell(); 
				pCell->assign_position( x , -y , 0.0 );
				add_reactivation_prob(pCell, reactivation_rate_mean, reactivation_rate_sd);
				add_mutation_rate(pCell, mutation_rate_mean, mutation_rate_sd);

				
			}
			
			if( fabs( x ) > 0.01 )
			{ 
				pCell = create_cell(); 
				pCell->assign_position( -x , y , 0.0 );
				add_reactivation_prob(pCell, reactivation_rate_mean, reactivation_rate_sd);
				add_mutation_rate(pCell, mutation_rate_mean, mutation_rate_sd);


				if( fabs( y ) > 0.01 )
				{
					pCell = create_cell();  
					pCell->assign_position( -x , -y , 0.0 );
					add_reactivation_prob(pCell, reactivation_rate_mean, reactivation_rate_sd);
					add_mutation_rate(pCell, mutation_rate_mean, mutation_rate_sd);

				}
			}
			x += cell_spacing; 
			
		}
		
		y += cell_spacing * sqrt(3.0)/2.0; 
		n++; 
	}


	get_heterogeneity_summary(pCell);

}

void get_heterogeneity_summary(Cell* pCell){

	static int reactivation_prob_idx = pCell->custom_data.find_variable_index("reactivation_value");

	double sum = 0.0; 
	double min = 9e9; 
	double max = -9e9; 
	for( int i=0; i < all_cells->size() ; i++ )
	{
		double r = (*all_cells)[i]->custom_data[reactivation_prob_idx]; 
		sum += r;
		if( r < min )
		{ min = r; } 
		if( r > max )
		{ max = r; }
	}
	double mean = sum / ( all_cells->size() + 1e-15 ); 
	// compute standard deviation 
	sum = 0.0; 
	for( int i=0; i < all_cells->size(); i++ )
	{
		sum +=  ( (*all_cells)[i]->custom_data[reactivation_prob_idx] - mean )*( (*all_cells)[i]->custom_data[reactivation_prob_idx] - mean ); 
	}
	double standard_deviation = sqrt( sum / ( all_cells->size() - 1.0 + 1e-15 ) ); 
	
	std::cout << std::endl << "Heterogeneity summary: " << std::endl
			  << "===================" << std::endl; 
	std::cout << "mean: " << mean << std::endl; 
	std::cout << "standard deviation: " << standard_deviation << std::endl; 
	std::cout << "[min max]: [" << min << " " << max << "]" << std::endl << std::endl; 

}




// OLD CODE


	// output[3].assign("black");

	// std::cout << "this cell " << pCell->ID <<" has a prob of " << reactivation_prob << std::endl;

	// int oncoprotein = (int) round( 0.5 * pCell->custom_data[oncoprotein_i] * 255.0 ); 
	// 	char szTempString [128];
	// 	sprintf( szTempString , "rgb(%u,%u,%u)", oncoprotein, oncoprotein, 255-oncoprotein );
	// 	output[0].assign( szTempString );
	// 	output[1].assign( szTempString );

	// 	sprintf( szTempString , "rgb(%u,%u,%u)", (int)round(output[0][0]/2.0) , (int)round(output[0][1]/2.0) , (int)round(output[0][2]/2.0) );
	// 	output[2].assign( szTempString );
		

	// // dead cells
	// if (pCell->phenotype.death.dead == false)
	// {
	// 	// static double V_cell = pCell->phenotype.volume.total; 
	// 	// static int drug_X_index = microenvironment.find_density_index("drug_X");
	// 	// static int drug_Y_index = microenvironment.find_density_index("drug_Y");

	// 	// double I_drug_X = pCell->phenotype.molecular.internalized_total_substrates[drug_X_index] / V_cell;
	// 	// double I_drug_Y = pCell->phenotype.molecular.internalized_total_substrates[drug_Y_index] / V_cell;

	// 	// float activation_threshold = pCell->custom_data.find_variable_index("activation threshold");

	// }