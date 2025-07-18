

#include "./custom.h"

/**
 *	\main prostate custom
 *	\brief Custom module file for prostate example
 * 
 *	\details Modules needed for the prostate example. 
 *
 *
 *	\date 14/12/2020
 *	\author Annika Meert, BSC-CNS, with code previously developed by Arnau Montagud, Gerard Pradas and Miguel Ponce de Leon, BSC-CNS
 */

// declare cell definitions here 

void create_cell_types( void )
{
	// use the same random seed so that future experiments have the 
	// same initial histogram of oncoprotein, even if threading means 
	// that future division and other events are still not identical 
	// for all runs 
	
	SeedRandom( parameters.ints("random_seed") ); // or specify a seed here 
	
	/* 
	   Put any modifications to default cell definition here if you 
	   want to have "inherited" by other cell types. 
	   
	   This is a good place to set default functions. 
	*/ 

	cell_defaults.functions.volume_update_function = standard_volume_update_function;
	cell_defaults.functions.update_velocity = standard_update_cell_velocity;
	
	cell_defaults.functions.pre_update_intracellular = pre_update_intracellular;
	cell_defaults.functions.post_update_intracellular = post_update_intracellular;
	
	cell_defaults.functions.update_migration_bias = NULL; 
	cell_defaults.functions.update_phenotype = tumor_cell_phenotype_with_signaling;
	cell_defaults.functions.custom_cell_rule = NULL; 
	
	cell_defaults.functions.add_cell_basement_membrane_interactions = NULL; 
	cell_defaults.functions.calculate_distance_to_membrane = NULL; 
	
	cell_defaults.functions.set_orientation = NULL;

	/*
	   This parses the cell definitions in the XML config file. 
	*/
	initialize_cell_definitions_from_pugixml();

	cell_defaults.phenotype.secretion.sync_to_microenvironment( &microenvironment );

	// set molecular properties for each drug defined in PhysiCell_settings.xml
	// starts with the second density because the first density is oxygen
	for (int i = 0; i < microenvironment.number_of_densities(); i++) 
	{
		if (microenvironment.density_names[i] != "oxygen") {
			cell_defaults.phenotype.molecular.fraction_released_at_death[i] = 0.0;
		}
	}

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

double get_decay_rate(double half_life){
	// natural logarithm of 2 / half-life = k (reaction rate coefficient)
	double decay_rate = log(2)/(half_life);
	return decay_rate;
}

void setup_microenvironment( void )
{
	// make sure to override and go back to 2D 
	if( default_microenvironment_options.simulate_2D == true )
	{
		std::cout << "Warning: overriding XML config option and setting to 3D!" << std::endl; 
		default_microenvironment_options.simulate_2D = false; 
	}	

	// set intial conditions and dirichlet boundary conditions for the drugs; vector already contains the condition for oxygen
	double oxygen_condition = 160.0;
	vector<double> condition_vector = {oxygen_condition};
	vector<bool> activation_vector {1};
	vector<double> decay_vector {0.1};
	for (int i = 0; i < microenvironment.number_of_densities(); i++)
	{
		std::string drug_name = microenvironment.density_names[i];
		if (drug_name != "oxygen") {
			// int current_drug_level= parameters.ints("current_concentration_level_" + drug_name);
			// int total_drug_levels = parameters.ints("total_concentration_levels");
			string cell_line = parameters.strings("cell_line");
			int simulation_mode = parameters.ints("simulation_mode");

			//drug_conc can either be an IC value or an actual drug concentration
			string drug_conc = parameters.strings("drug_concentration_" + drug_name);
			double drug_concentration;
			// check if drug_conc contains the string "IC"
			if(drug_conc.find("IC") != string::npos) {
				drug_concentration = get_drug_concentration_from_IC(cell_line, drug_name, drug_conc, simulation_mode);
			
			}
			else {
				drug_concentration = stod(drug_conc);
				}
			
			// double drug_concentration = get_drug_concentration_from_level(cell_line, drug_name, current_drug_level, total_drug_levels, simulation_mode);
			condition_vector.push_back(drug_concentration);
			activation_vector.push_back(1);

			// double half_life = get_value(half_lives, drug_name);
			// double decay_rate = get_decay_rate(half_life);
			// decay_vector.push_back(decay_rate);

		}
	}
	default_microenvironment_options.Dirichlet_activation_vector = activation_vector;
	default_microenvironment_options.Dirichlet_condition_vector = condition_vector;
	default_microenvironment_options.Dirichlet_zmax_values = condition_vector;
	default_microenvironment_options.Dirichlet_zmin_values = condition_vector;
	default_microenvironment_options.Dirichlet_ymax_values = condition_vector;
	default_microenvironment_options.Dirichlet_ymin_values = condition_vector;
	default_microenvironment_options.Dirichlet_xmax_values = condition_vector;
	default_microenvironment_options.Dirichlet_xmin_values = condition_vector;


	// initialize BioFVM 
	initialize_microenvironment(); 	

	// microenvironment.decay_rates = decay_vector;
	
	return; 
}

void setup_tissue( void )
{
	Cell* pC;

	std::vector<init_record> cells = read_init_file(parameters.strings("init_cells_filename"), ';', true);
	
	for (int i = 0; i < cells.size(); i++)
	{
		float x = cells[i].x;
		float y = cells[i].y;
		float z = cells[i].z;
		float radius = cells[i].radius;
		int phase = cells[i].phase;
		double elapsed_time = cells[i].elapsed_time;

		double random_num_1 = (double) rand()/RAND_MAX;
		double random_num_2 = (double) rand()/RAND_MAX;

		if (PhysiCell::parameters.ints("simulation_mode") == 0)
		{
			// single inhibition - just one drug is present 
			if (random_num_1 < PhysiCell::parameters.doubles("prop_drug_resistant_" + microenvironment.density_names[1]))
			{
				// cell is sensitive to the drug
				pC = create_cell(get_cell_definition(microenvironment.density_names[1] + "_resistant"));
			}
			else 
			{
				// cell is not sensitive to the drug
				pC = create_cell(get_cell_definition(microenvironment.density_names[1] + "_sensitive"));
			}
		}
		else if (PhysiCell::parameters.ints("simulation_mode") == 1)
		{
			// double inhibition - two drugs are present - we have 4 cell strains 
			if (random_num_1 < PhysiCell::parameters.doubles("prop_drug_resistant_" + microenvironment.density_names[1]))
			{
				if (random_num_2 < PhysiCell::parameters.doubles("prop_drug_resistant_" + microenvironment.density_names[2]))
				{
					// cell is resistant to both drugs
					pC = create_cell(get_cell_definition(microenvironment.density_names[1] + "_resistant"));
				}
				else 
				{
					// cell is only resistant to the first drug
					pC = create_cell(get_cell_definition(microenvironment.density_names[2] + "_sensitive"));
				}
			}
			else
			{
				if (random_num_2 < PhysiCell::parameters.doubles("prop_drug_resistant_" + microenvironment.density_names[2]))
				{
					// cell is only resistant to the second drug
					pC = create_cell(get_cell_definition(microenvironment.density_names[2] + "_resistant"));
				}
				else
				{
					// cell is resistant to no drug
					pC = create_cell(get_cell_definition(microenvironment.density_names[1] + "_sensitive"));
				}
				
			}
			
		}
		else
		{
			pC = create_cell(get_cell_definition("default"));
		}
 
		pC->assign_position( x, y, z );
		// pC->set_total_volume(sphere_volume_from_radius(radius));
		
		// pC->phenotype.cycle.data.current_phase_index = phase;
		pC->phenotype.cycle.data.elapsed_time_in_phase = elapsed_time;	
		
		update_custom_variables(pC);
	}

	return; 
}

// custom cell phenotype function to run PhysiBoSS when is needed
void tumor_cell_phenotype_with_signaling( Cell* pCell, Phenotype& phenotype, double dt)
{
	update_cell_and_death_parameters_O2_based(pCell, phenotype, dt);

	// update motility state variable
	static int index_motility_state = pCell->custom_data.find_variable_index("motility_state");
	pCell->custom_data.variables.at(index_motility_state).value = int(pCell->phenotype.motility.is_motile);
	
	boolean_model_interface_main (pCell, phenotype, dt);
}

std::vector<std::string> prolif_apoptosis_coloring( Cell* pCell )
{
	std::vector<std::string> output;
	if (pCell->phenotype.cycle.current_phase().code == PhysiCell_constants::apoptosis_death_model)
	{
		//apoptotic cells are colored red
		output = {"crimson", "black","darkred", "darkred"};
	}

	else if (pCell->phenotype.cycle.current_phase().code == PhysiCell_constants::necrosis_death_model)
	{
		//necrotic cells are colored brown
		output = {"peru", "black","saddlebrown", "saddlebrown"};
	}

	else if (PhysiCell::parameters.ints("simulation_mode") == 0) 
	{
		std::string drug_name = microenvironment.density_names[1];
		if (pCell->type_name == drug_name + "_sensitive")
		{
			//drug sensitive living cells are colored blue
			output = {"deepskyblue", "black", "darkblue", "darkblue"};
		} 
		else 
		{
			//drug resistant living cells are colored green
			output = {"limegreen", "black", "darkgreen", "darkgreen"};
		}
	}
	else if (PhysiCell::parameters.ints("simulation_mode") == 1) 
	{
		// // color living cells just in one color 
		// output = {"limegreen", "black", "darkgreen", "darkgreen"};

		// In case we want to color all 4 strains differently:
		// double inhibitions --> 4 cell strains 
		std::string drug1_name = microenvironment.density_names[1];
		std::string drug2_name = microenvironment.density_names[2]; 
		if (pCell->type_name == drug1_name + "_sensitive")
		{
			//cells that are sensitive to both drugs are colored blue
			output = {"deepskyblue", "black", "darkblue", "darkblue"};
		}
		else if (pCell->type_name == drug2_name + "_resistant")
		{
			// cells that are just sensitive to the first drug 
			output = {"limegreen", "black", "darkgreen", "darkgreen"};

		}
		else if (pCell->type_name == drug2_name + "_sensitive")
		{
			// cells are just sensitive to the second drug
			output = {"gold", "black", "orange", "orange"};
		}
		else 
		{
			// cells aren't sensitive to any drug
			output = {"mediumorchid", "black", "mediumpurple", "mediumpurple"};
		}
	}
	else 
	{
		// no drug simulation - living cells are colored green
		output = {"limegreen", "black", "darkgreen", "darkgreen"};
	}
	return output;

}

// ***********************************************************
// * NOTE: Funtion to read init files created with PhysiBoSS *
// ***********************************************************
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

	if(header)
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
		while (getline(s, word, delimiter)) { 

			// add all the column data 
			// of a row to a vector 
			row.push_back(word); 
		}

		init_record record;
		record.x = std::stof(row[2]);
		record.y = std::stof(row[3]);
		record.z = std::stof(row[4]);
		record.radius = std::stof(row[5]);
		record.phase = std::stoi(row[13]);
		record.elapsed_time = std::stod(row[14]);

		result.push_back(record);
	} while (!fin.eof());
	
	return result;
}