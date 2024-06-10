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

#define _USE_MATH_DEFINES
#include <cmath>
#include <sstream>
#include "./custom.h"
#include "../BioFVM/BioFVM.h"  

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
	cell_defaults.functions.custom_cell_rule = NULL; // @oth: add here custom cell cycle function?
	cell_defaults.functions.contact_function = NULL; 
	
	cell_defaults.functions.add_cell_basement_membrane_interactions = NULL; 
	cell_defaults.functions.calculate_distance_to_membrane = NULL; 

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

	// This sets the pre and post simulation functions
	cell_defaults.functions.pre_update_intracellular =  pre_update_intracellular_ags;
	cell_defaults.functions.post_update_intracellular = post_update_intracellular_ags;
	cell_defaults.functions.update_phenotype = update_cell_gowth_parameters_pressure_based;

	drug_transport_model_setup();
	boolean_model_interface_setup();

	cell_defaults.custom_data.add_variable("time", "min", PhysiCell_globals.current_time );
	cell_defaults.custom_data.add_variable("total_live_cells", "dimensionless", 0.0 );

	cell_defaults.custom_data.add_variable("reactivation_rate", "1/min", 0.0 );
	cell_defaults.custom_data.add_variable("mutation_rate", "1/generation", 0.0 );

	submodel_registry.display(std::cout);
	/*
	   This builds the map of cell definitions and summarizes the setup. 
	*/

	display_cell_definitions(std::cout);

	// @oth: moved cell cycle function outside of here to the phenotype function
	return;
}

void setup_microenvironment(void)
{
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

// @othmane [WIP] -- Copied from cancer_invasion model, left empty so far -- should this go here?

// void pre_update_intracellular(Cell* pCell, Phenotype& phenotype, double dt){
// 	return;
// }



void update_cell_gowth_parameters_pressure_based( Cell* pCell, Phenotype& phenotype, double dt ) 
{
	
	if( phenotype.death.dead == true )
	{ return; }

	// Custom cycle exit in the phenotype
	PhysiCell::Cycle_Model cycle_model = cell_defaults.phenotype.cycle.model();
	static int oxygen_substrate_index = pCell->get_microenvironment()->find_density_index( "oxygen" ); 
	static int start_phase_idx = phenotype.cycle.model().find_phase_index( PhysiCell_constants::live );
	static int end_phase_idx = phenotype.cycle.model().find_phase_index( PhysiCell_constants::live );
	static int necrosis_idx = phenotype.death.find_death_model_index( PhysiCell_constants::necrosis_death_model );

	cycle_model.phase_link(start_phase_idx, end_phase_idx).exit_function = my_mutation_function;

	
	// this multiplier is for linear interpolation of the oxygen value 	
	// PRESSURE-BASED CONTACT-INHIBITION
	// Check relative pressure to eith er number of neighbor cells or set max logistic function to number of neighbor cells
	// pressure threshold set to 1, above this value there is no growth

	double cell_pressure = pCell->state.simple_pressure; 
    double hill_coeff_pressure = parameters.doubles("hill_coeff_pressure");
    double pressure_half = parameters.doubles("pressure_half");

	// @oth: Discarded -- use PC Hill function instead
	// double scaling = pressure_effect_growth_rate(cell_pressure, hill_coeff_pressure, pressure_half); // @oth: Discarded -- use PC Hill function instead

	double scaling = Hill_response_function(cell_pressure, pressure_half, hill_coeff_pressure);
	double growth_rate = phenotype.cycle.data.transition_rate(start_phase_idx, end_phase_idx);

	// Pressure affects negatively
	
	growth_rate *= (1 - scaling);
	if (growth_rate < 0)
		growth_rate = 0;
	
	phenotype.cycle.data.transition_rate(start_phase_idx, end_phase_idx) = growth_rate;


	return;
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

void inject_density_sphere(int density_index, double concentration, double membrane_length) 
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


// FUNCTIONS TO PLOT CELLS

std::vector<std::string> my_coloring_function_for_stroma( double concentration, double max_conc, double min_conc )
{
	 return paint_by_density_percentage( concentration,  max_conc,  min_conc); 

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