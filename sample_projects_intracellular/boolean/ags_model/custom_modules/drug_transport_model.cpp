/*
	
	ags_receptor_dynamics.cpp

	each drug that enters does so with its own mechanism
*/

#include "./drug_transport_model.h" 
#include "./boolean_model_interface.h" 

using namespace PhysiCell; 

Submodel_Information ags_receptor_info;


void drug_transport_model_setup()
{
    ags_receptor_info.name = "AGS model"; 
	ags_receptor_info.version = "0.0.1";
	
    ags_receptor_info.main_function = drug_transport_model_update; 

	// what custom data do I need?
	ags_receptor_info.cell_variables.push_back( "activation_threshold" );

	// drug_X variables
	ags_receptor_info.cell_variables.push_back( "drug_X_permeability" );
	ags_receptor_info.cell_variables.push_back( "drug_X_external_density" );
	ags_receptor_info.cell_variables.push_back( "drug_X_internal_density" );

	// drug_Y variables
	ags_receptor_info.cell_variables.push_back( "drug_Y_permeability" );
	ags_receptor_info.cell_variables.push_back( "drug_Y_external_density" );
	ags_receptor_info.cell_variables.push_back( "drug_Y_internal_density" );

	// oxygen
	ags_receptor_info.cell_variables.push_back( "oxygen_external_density" );
	ags_receptor_info.cell_variables.push_back( "oxygen_internal_density" );

	ags_receptor_info.register_model();

	return;
}


double calculate_diffusion_flux(Cell* pCell, int density_idx, double permeability, std::string drug_name){
	
	// Cell and voxel geometries
	static double voxel_volume = microenvironment.mesh.dV;

	double cell_volume = pCell->phenotype.volume.total;
	double cell_radius = pCell->phenotype.geometry.radius;
	double cell_surface = 4 * M_PI * std::pow(cell_radius, 2);

	double density_ext = pCell->nearest_density_vector()[density_idx];	 // A density (mM)
	if(density_ext < 0.0){ density_ext = 0.0;}
	// if(density_ext > parameters.doubles("drug_X_pulse_concentration")){density_ext =  parameters.doubles("drug_X_pulse_concentration");}

	double density_int = pCell->phenotype.molecular.internalized_total_substrates[density_idx];
	density_int /= cell_volume; // divide int tot substrate to get concentration
	if(density_int < 0.0){ density_int = 0.0;}
	// if(density_int > parameters.doubles("drug_X_pulse_concentration")){density_int =  parameters.doubles("drug_X_pulse_concentration");}
	
	// Simple Diffusion entry of drug
	double flux = 0.0;
	flux = permeability * (density_int - density_ext) * cell_surface; // amol/min
	// then map to custom data
	std::string drug_external_density = drug_name + "_external_density";
	std::string drug_internal_density = drug_name + "_internal_density";

	

	// Adding an active bomb -- simplified model from the MSc Thesis
	// Also includes a trade-off with cell growth 

	double k_2 = parameters.doubles("pump_k_2_value");
	double total_Enzyme = parameters.doubles("pump_total_Enzyme_value");
	double Km = parameters.doubles("pump_Km_value");

	double pumping_rate = ((k_2 * total_Enzyme) * density_int) / ( Km + density_int); 
	pumping_rate *= cell_volume; // Convert mM/min to amol/min

	flux += pumping_rate;

	// Add effect on growth of pumping	
	double rate = pCell->phenotype.cycle.data.transition_rate(0, 0);
	rate *= 1/(pumping_rate + 1);
	if (rate < 0)
		rate = 0;
	
	pCell->phenotype.cycle.data.transition_rate(0, 0) = rate;

	// // Drug X
	static int drug_external_density_custom_data_idx = pCell->custom_data.find_variable_index(drug_external_density);
	static int drug_internal_density_custom_data_idx = pCell->custom_data.find_variable_index(drug_internal_density);
	pCell->custom_data.variables[drug_external_density_custom_data_idx].value = density_ext;
	pCell->custom_data.variables[drug_internal_density_custom_data_idx].value = density_int;


	return flux;
}



void drug_transport_model_update( Cell* pCell, Phenotype& phenotype, double dt )
{	
	if( phenotype.death.dead == true )
	{ return; }

	// Fetch microenvironment variables
	static int drug_X_idx = microenvironment.find_density_index("drug_X");
	static int drug_Y_idx = microenvironment.find_density_index("drug_Y");
	
	double drug_X_permeability = parameters.doubles("drug_X_permeability");
	double drug_Y_permeability = parameters.doubles("drug_Y_permeability");

	// pCell->phenotype.secretion.uptake_rates[drug_X_idx] = drug_X_permeability;
	// pCell->phenotype.secretion.uptake_rates[drug_Y_idx] = drug_X_permeability;
	// pCell->phenotype.secretion.secretion_rates[drug_X_idx] = drug_Y_permeability;
	// pCell->phenotype.secretion.secretion_rates[drug_Y_idx] = drug_Y_permeability;

	double drug_X_flux = calculate_diffusion_flux(pCell, drug_X_idx, drug_X_permeability, "drug_X");
	std::cout << "drug_X_flux " << drug_X_flux << std::endl;
	pCell->phenotype.secretion.net_export_rates[drug_X_idx] = drug_X_flux;

	double drug_Y_flux = calculate_diffusion_flux(pCell, drug_Y_idx, drug_Y_permeability, "drug_Y");
	pCell->phenotype.secretion.net_export_rates[drug_Y_idx] = drug_Y_flux;

	return;
}

void drug_transport_model_main( double dt )
{

	#pragma omp parallel for 
	for( int i=0; i < (*all_cells).size() ; i++ )
	{
		Cell* pCell = (*all_cells)[i]; 


		if( pCell->phenotype.death.dead == true || pCell->is_out_of_domain == true)
		{ continue; }


		int cell_microenv_index = pCell->get_current_voxel_index();
		bool DC_node = microenvironment.is_dirichlet_node(cell_microenv_index);
		bool agent_out_domain = pCell->is_out_of_domain;

		if (cell_microenv_index >= 0 && DC_node == false) // Avoid segfault of NER
		{
			// pCell->phenotype.secretion.advance( pCell, pCell->phenotype, dt );
			drug_transport_model_update( pCell, pCell->phenotype , dt );
			update_diffusion_monitor_vars(pCell);
		}
		
		
		pCell->set_internal_uptake_constants(dt);

		// pCell->phenotype.secretion.advance( pCell, pCell->phenotype, dt );

		// Adding time to custom data
		static int time_index = pCell->custom_data.find_variable_index("time");
		pCell->custom_data.variables[time_index].value = PhysiCell_globals.current_time;

	}
	

	return;
}



/*
	There must be a better way of automatising this...

		Previously, the mapping of the internal and external density values for each drug were automatically
		mapped to their respective custom data (see lines 106-111 of this script). For some reason, PhysiCell does not like this, so there was no 
		way of plotting both substrate densities at the same time, for debugging purposes. 

*/

void update_diffusion_monitor_vars( Cell* pCell ){

	// std::cout << "UPDATING diffusion monitor vars at " << PhysiCell_globals.current_time << std::endl;


	// Fetch microenvironment variables
	static int drug_X_idx = microenvironment.find_density_index("drug_X");
	static int drug_Y_idx = microenvironment.find_density_index("drug_Y");
	static double cell_volume = pCell->phenotype.volume.total; 

	// // Drug X
	static int drug_X_external_density_custom_data_idx = pCell->custom_data.find_variable_index("drug_X_external_density");
	static int drug_X_internal_density_custom_data_idx = pCell->custom_data.find_variable_index("drug_X_internal_density");

	float drug_X_external_density = pCell->nearest_density_vector()[drug_X_idx];
	float drug_X_internal_density = pCell->phenotype.molecular.internalized_total_substrates[drug_X_idx];	
	drug_X_internal_density /= cell_volume;

	// std::cout << "MAPPING density ext " << drug_X_external_density << std::endl;
	// std::cout << "MAPPING density int " << drug_X_internal_density << std::endl;

	pCell->custom_data.variables[drug_X_external_density_custom_data_idx].value = drug_X_external_density;
	pCell->custom_data.variables[drug_X_internal_density_custom_data_idx].value = drug_X_internal_density;

	// Drug Y 
	static int drug_Y_external_density_custom_data_idx = pCell->custom_data.find_variable_index("drug_Y_external_density");
	static int drug_Y_internal_density_custom_data_idx = pCell->custom_data.find_variable_index("drug_Y_internal_density");
	float drug_Y_external_density = pCell->nearest_density_vector()[drug_Y_idx];
	float drug_Y_internal_density = pCell->phenotype.molecular.internalized_total_substrates[drug_Y_idx];
	drug_Y_internal_density /= cell_volume;

	pCell->custom_data.variables[drug_Y_external_density_custom_data_idx].value = drug_Y_external_density;
	pCell->custom_data.variables[drug_Y_internal_density_custom_data_idx].value = drug_Y_internal_density;
		
	// map oxygen to custom data -- TODO: Encapsulate this in a function
	static int oxy_idx = microenvironment.find_density_index("oxygen");
	float oxy_int = pCell->phenotype.molecular.internalized_total_substrates[oxy_idx] / cell_volume;
	float oxy_ext = pCell->nearest_density_vector()[oxy_idx];
	static int oxy_int_index = pCell->custom_data.find_variable_index("oxygen_internal_density");
	static int oxy_ext_index = pCell->custom_data.find_variable_index("oxygen_external_density");
	pCell->custom_data.variables[oxy_int_index].value = oxy_int;
	pCell->custom_data.variables[oxy_ext_index].value = oxy_ext;


	return;

}

