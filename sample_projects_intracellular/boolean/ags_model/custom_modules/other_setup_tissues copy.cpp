// void setup_tissue( void ){
	
// 	double Xmin = microenvironment.mesh.bounding_box[0]; 
// 	double Ymin = microenvironment.mesh.bounding_box[1]; 
// 	double Zmin = microenvironment.mesh.bounding_box[2]; 

// 	double Xmax = microenvironment.mesh.bounding_box[3]; 
// 	double Ymax = microenvironment.mesh.bounding_box[4]; 
// 	double Zmax = microenvironment.mesh.bounding_box[5]; 
	
// 	double max_radius = parameters.doubles("tumor_radius");
// 	if( Xmax > max_radius )
// 	{ Xmax = max_radius; }
// 	if( Xmin < -max_radius )
// 	{ Xmin = -max_radius; }
	
// 	if( Ymax > max_radius )
// 	{ Ymax = max_radius; }
// 	if( Ymin < -max_radius )
// 	{ Ymin = -max_radius; }

// 	if( Zmax > max_radius )
// 	{ Zmax = max_radius; }
// 	if( Zmin < -max_radius )
// 	{ Zmin = -max_radius; }
	
// 	if( default_microenvironment_options.simulate_2D == true )
// 	{
// 		Zmin = 0.0; 
// 		Zmax = 0.0; 
// 	}

// 	double Xrange = Xmax - Xmin; 
// 	double Yrange = Ymax - Ymin; 
// 	double Zrange = Zmax - Zmin; 

// 	Cell* pC;
	
// 	for( int n = 0 ; n < parameters.ints("number_of_A") ; n++ )
// 	{
// 		std::vector<double> position = {0,0,0}; 
		
// 		double r = max_radius + 1; 
// 		while( r > max_radius )
// 		{
// 			position[0] = Xmin + UniformRandom()*Xrange; 
// 			position[1] = Ymin + UniformRandom()*Yrange; 
// 			position[2] = Zmin + UniformRandom()*Zrange; 
			
// 			r = norm( position ); 
// 		}
		
// 		pC = create_cell( get_cell_definition("default") ); 
// 		pC->assign_position( position );
// 		for( int k=0 ; k < pC->phenotype.death.rates.size() ; k++ )
// 		{ pC->phenotype.death.rates[k] = 0.0; }
// 	}
// }

// void setup_tissue( void )
// {
// 	// place a cluster of tumor cells at the center 
	
// 	double cell_radius = cell_defaults.phenotype.geometry.radius * 2.0 * 0.8; 
// 	double cell_spacing = 1.8 * 2.0 * cell_radius; 
	
// 	double tumor_radius = 
// 		parameters.doubles("tumor_radius"); // 250.0; 
	
// 	Cell* pCell = NULL; 
	
// 	// std::vector<std::vector<double>> positions = create_cell_sphere_positions(cell_radius, tumor_radius); 
// 	std::vector<std::vector<double>> positions = create_cell_circle_positions(cell_radius, tumor_radius); 
// 	std::cout << "creating " << positions.size() << " closely-packed tumor cells ... " << std::endl; 

// 	// #pragma omp parallel for 
// 	for( int i=0; i < positions.size(); i++ )
// 	{
		
// 		pCell = create_cell(); 
// 		pCell->assign_position( positions[i] );
// 		std::cout << "Creating cell in position " << positions[i] << std::endl;

// 		pCell->phenotype.intracellular->start();

// 		const int max_iter = 1;

// 		#pragma omp parallel for 
// 		for (int j=0; j < max_iter; j++){
// 			pCell->phenotype.intracellular->update(); }
		
// 		// pCell->phenotype.intracellular->print_current_nodes();
		
	
// 		update_monitor_variables(pCell);
	
// 	}

// 	return; 
// }