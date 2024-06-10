

	// Following lines of code do not work properly, do not know why
	// std::string substrate_external_idx = drug_name + "_external_density";
	// std::string substrate_internal_idx = drug_name + "_internal_density";
	// static int external_density_custom_data_idx = pCell->custom_data.find_variable_index(substrate_external_idx);
	// static int internal_density_custom_data_idx = pCell->custom_data.find_variable_index(substrate_internal_idx);
	// pCell->custom_data.variables[external_density_custom_data_idx].value = density_ext;
	// pCell->custom_data.variables[internal_density_custom_data_idx].value = density_int;




	// -- Option A 
	// First, check for the equilibirum concentration

	// double flux_net_rate = std::abs(flux*diffusion_dt); // in amol
	// std::string drug_initial_external_conc = drug_name + "_pulse_concentration";

	// float initial_external_density = parameters.doubles(drug_initial_external_conc); // in mM
	// float initial_external_net_amount = initial_external_density * voxel_volume; // in mmol
	// float actual_external_net_amount = density_ext * voxel_volume;

	// float initial_internal_density = 0.0; // mM
	// float initial_internal_net_amount = initial_internal_density * cell_volume; // amol
	// float actual_internal_net_amount = density_int * cell_volume;

	// float total_net_amount = initial_external_net_amount + initial_internal_net_amount;
	// float actual_total_net_amount = actual_external_net_amount + actual_internal_net_amount;

	
	// // At equilibirum (net amols)
	// float eq_external_density = (total_net_amount / 2) / voxel_volume;
	// float eq_internal_density = (total_net_amount / 2) / cell_volume;

	// float dynamic_eq_external_density = (actual_total_net_amount / 2) / voxel_volume;
	// float dynamic_eq_internal_density = (actual_total_net_amount / 2) / cell_volume;

	// // Check if net export rate step will pass equilibirum
	// float internal_density_at_t1 = density_int + (flux_net_rate / cell_volume);
	// float external_density_at_t1 = density_ext + (flux_net_rate / voxel_volume);

	// if (internal_density_at_t1 < 0.0){ internal_density_at_t1 = 0.0; }
	// if (external_density_at_t1 < 0.0){ external_density_at_t1 = 0.0; }

	// std::cout << "Unadjusted flux is: " << flux << std::endl;
	// std::cout << "Int. dens at t1 is: " << internal_density_at_t1 << std::endl;
	// std::cout << "Eq. int density: " << dynamic_eq_internal_density << std::endl;

	// if ( flux < 0.0 && internal_density_at_t1 > dynamic_eq_internal_density )
	// {
	// 	// Introduce the difference between the internal density and 
	// 	flux = - ( density_int - dynamic_eq_internal_density ); // remember: negative is uptake
	// 	// std::cout << "Adjusting entry flux" << std::endl;
	// 	flux = 0.0;
	// }

	// if ( flux > 0.0 && external_density_at_t1 > dynamic_eq_external_density ){
	// 	// Introduce the difference between the internal density and 
	// 	flux = (density_ext - dynamic_eq_external_density);
	// 	// std::cout << "Adjusting exit flux" << std::endl;
	// 	flux = 0.0;
	// }
	

	// -- Option B: Using only the delta of net amount, not concentration

	// if( flux < 0.0 && std::abs(flux) >= density_ext * voxel_volume * 0.00001 ){
	// 	std::cout << "Adjusting exit flux" << std::endl;
	// 	flux = - (density_ext * voxel_volume) - (density_int * cell_volume); // amol
	// } else if ( flux > 0.0 && flux >= density_int * cell_volume * 0.00001){
	// 	std::cout << "Adjusting entry flux" << std::endl;
	// 	flux = (density_int * cell_volume) - (density_ext * voxel_volume); // amol
	// }

	// for( int i=0; i < (*all_cells).size() ; i++ ){}
	
	// pCell->phenotype.secretion.net_export_rates[ density_idx ] = flux;

	// OPTION C: Using same adjusting as for the 

	// if ( flux < 0.0 ) { // Uptake
	// 	flux = - std::min( std::abs(flux), std::abs(density_ext * voxel_volume - density_int * cell_volume) );
	// }
	// else if ( flux > 0.0 ) { // Secretion
	// 	flux = std::min(flux, std::abs(density_ext * voxel_volume - density_int * cell_volume) );
	// }

	// std::cout << "flux is: " << flux << std::endl;