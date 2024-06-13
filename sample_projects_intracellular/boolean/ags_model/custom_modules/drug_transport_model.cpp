/*
	
	Here are the drug transport dynamics
	
	Currently, both drugs enter through simple diffusion.
	In the future, it could be changed to include more complicated mechanisms.

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

double get_custom_data_variable(Cell* pCell, std::string variable_name){
	static int tmp_variable_idx = pCell->custom_data.find_variable_index(variable_name);
	double tmp_variable_value = pCell->custom_data[tmp_variable_idx];
	return tmp_variable_value;
}

void update_diffusion_monitor_vars( Cell* pCell ){

	// std::cout << "UPDATING diffusion monitor vars at " << PhysiCell_globals.current_time << std::endl;

	int n_densities = 2;
    std::string densities[n_densities] = { "drug_X", "drug_Y", "oxygen"};
  
    for (int i = 0; i < n_densities; i++){
        std::string density_name = densities[i];

		static int density_idx = microenvironment.find_density_index(density_name); 
		static int density_ext_custom_var_idx = pCell->custom_data.find_variable_index(density_name + "_external_density");
		static int density_int_custom_var_idx = pCell->custom_data.find_variable_index(density_name + "_internal_density");

		double cell_volume = pCell->phenotype.volume.total;
		double density_ext = pCell->nearest_density_vector()[density_idx];
		double total_density_int = pCell->phenotype.molecular.internalized_total_substrates[density_idx];

		if (density_ext > 0.0){ std::cout << "for drug " << density_name << " ext is " << density_ext << std::endl; }

		pCell->custom_data.variables[density_ext_custom_var_idx].value = density_ext;
		pCell->custom_data.variables[density_int_custom_var_idx].value = total_density_int / cell_volume;
	}

	return;
}


// @oth: #TODO encapsulate this in an external function
// double pump_model()

	// Adding an active bomb -- simplified model from the MSc Thesis
	// Also includes a trade-off with cell growth 

	// double k_2 = parameters.doubles("pump_k_2_value");
	// double total_Enzyme = parameters.doubles("pump_total_Enzyme_value");
	// double Km = parameters.doubles("pump_Km_value");

	// double pumping_rate = ((k_2 * total_Enzyme) * density_int) / ( Km + density_int); 
	// pumping_rate *= cell_volume; // Convert mM/min to amol/min

	// flux += pumping_rate;

	// // Add effect on growth of pumping	
	// double rate = pCell->phenotype.cycle.data.transition_rate(0, 0);
	// rate *= 1/(pumping_rate + 1);
	// if (rate < 0)
	// 	rate = 0;
	
	// pCell->phenotype.cycle.data.transition_rate(0, 0) = rate;


double calculate_diffusion_flux(Cell* pCell, int density_idx, double permeability){
	
	// Cell and voxel geometries
	static double voxel_volume = microenvironment.mesh.dV;

	double cell_volume = pCell->phenotype.volume.total;
	double cell_radius = pCell->phenotype.geometry.radius;
	double cell_surface = 4 * M_PI * std::pow(cell_radius, 2);

	double density_ext = pCell->nearest_density_vector()[density_idx];	 // A density (mM)
	double total_substrate = pCell->phenotype.molecular.internalized_total_substrates[density_idx];
	double density_int = total_substrate / cell_volume; // divide int tot substrate to get concentration
	if(density_int < 0.0)
	{ density_int = 0.0;}

	// Simple Diffusion entry of drug
	static double diffusion_constant = permeability * cell_surface;
	double flux = diffusion_constant * (density_int - density_ext) ; // amol/min

	return flux;
}



void drug_transport_model_update( Cell* pCell, Phenotype& phenotype, double dt )
{	
	if( phenotype.death.dead == true )
	{ return; }
	int n_drugs = 2;
    std::string drugs[n_drugs] = { "drug_X", "drug_Y" };
  
    for (int i = 0; i < n_drugs; i++){
        std::string drug_name = drugs[i];

		static int drug_idx = microenvironment.find_density_index(drug_name);
		static std::string param_permea = drug_name + "_permeability";
		double drug_permeability = get_custom_data_variable(pCell, param_permea);
		double drug_flux = calculate_diffusion_flux(pCell, drug_idx, drug_permeability);

		pCell->phenotype.secretion.net_export_rates[drug_idx] = drug_flux;
		pCell->set_internal_uptake_constants(dt);
	}
	
	return;
}

void drug_transport_model_main( double dt )
{
	#pragma omp parallel for 
	for( int i=0; i < (*all_cells).size() ; i++ )
	{
		Cell* pCell = (*all_cells)[i]; 

		if( pCell->phenotype.death.dead == true )
		{ continue; }
		if( pCell->is_out_of_domain == true )
		{ continue; }
		if( microenvironment.is_dirichlet_node(pCell->get_current_voxel_index()) )
		{ continue; }
		
		drug_transport_model_update( pCell, pCell->phenotype , dt );

		update_diffusion_monitor_vars(pCell);
	}
	
	return;
}

