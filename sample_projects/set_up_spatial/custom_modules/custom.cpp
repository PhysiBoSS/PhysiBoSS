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
	// cell_defaults.functions.update_velocity = standard_update_cell_velocity;
	cell_defaults.functions.update_velocity = custom_update_velocity;

	cell_defaults.functions.update_migration_bias = NULL; 
	cell_defaults.functions.update_phenotype = NULL; // update_cell_and_death_parameters_O2_based; 
	cell_defaults.functions.custom_cell_rule = NULL; 
	cell_defaults.functions.contact_function = NULL; 
	
	cell_defaults.functions.add_cell_basement_membrane_interactions = NULL; 
	cell_defaults.functions.calculate_distance_to_membrane = NULL; 
	
	/*
	   This parses the cell definitions in the XML config file. 
	*/
	
	initialize_cell_definitions_from_pugixml(); 
	
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
		
	build_cell_definitions_maps(); 
	display_cell_definitions( std::cout ); 
	
	return; 
}

void setup_microenvironment( void )
{
	// make sure ot override and go back to 2D 
	if( default_microenvironment_options.simulate_2D == false )
	{
		std::cout << "Warning: overriding XML config option and setting to 2D!" << std::endl; 
		default_microenvironment_options.simulate_2D = true; 
	}
			
	initialize_microenvironment();
	double o2_conc;
	if ( parameters.bools("read_init") )
	{
		o2_conc = parameters.doubles("o2_conc");
	}
	else
	{
		double o2_conc = 38.06;
	}
	std::vector<double> dirichlet_o2( 1 , o2_conc );
	// std::ifstream file( "/home/alejandro/Escritorio/MasterBioinformatica/TFM/TFM2022AlejandroMadrid/pruebas/Stutilityvsscanpy/setup_physi/setup_endo_norm.csv", std::ios::in );
	// if( !file )
	// { 
	// 	std::cout << "Error: " << "set_upf_endo_norm.csv" << " not found during cell loading. Quitting." << std::endl; 
	// 	exit(-1);
	// }
	// std::string line;
	// while (std::getline(positions, line))
	// {
	// 	std::vector<double> data;
	// 	csv_to_vector( line.c_str() , data ); 
	// 	if( data.size() != 3 )
	// 	{
	// 		std::cout << "Error! Importing cells from a CSV file expects each row to be x,y,z." << std::endl;
	// 		exit(-1);
	// 	};
	// 	/*std::vector<double> position = { data[0] , data[1] , data[2] }; */	
	// 	microenvironment.add_dirichlet_node( microenvironment.voxel_index(data[0],data[1],data[2]) , dirichlet_o2 );
	// }

	std::vector<std::vector<double>> positions;
	if ( parameters.bools("read_init") )
	{
		std::string csv_fname = parameters.strings("blood_source_file");
		positions = read_cells_positions(csv_fname, '\t', true);
	}

	for (int i = 0; i < positions.size(); i++)
	{
		int x = (positions[i][0]);
		int y = ( positions[i][1]);
		int z = ( positions[i][2]);
		microenvironment.add_dirichlet_node( microenvironment.voxel_index(x,y,z) , dirichlet_o2 );
		// microenvironment.add_dirichlet_node( microenvironment.voxel_index(positions[0],positions[1],positions[2]) , dirichlet_o2 );
	}
	// file.close(); 	
	microenvironment.set_substrate_dirichlet_activation( 0, true ); 

	/* here we try to set up the ecm */

	return; 
}	

std::vector<std::vector<double>> read_cells_positions(std::string filename, char delimiter, bool header)
{
	// File pointer
	std::fstream fin;
	std::vector<std::vector<double>> positions;

	// Open an existing file
	fin.open(filename, std::ios::in);

	// Read the Data from the file
	// as String Vector
	std::vector<std::string> row;
	std::string line, word;

	if (header)
	{ getline(fin, line); }

	do
	{
		row.clear();
		// read an entire row and
		// store it in a string variable 'line'
		getline(fin, line);
		// used for breaking words
		std::stringstream s(line);

		while (getline(s, word, delimiter))
		{ 
			row.push_back(word); 
		}

		std::vector<double> tempPoint(3,0.0);
		tempPoint[0]= std::stof(row[0]);
		tempPoint[1]= std::stof(row[1]);
		tempPoint[2]= std::stof(row[2]);

		positions.push_back(tempPoint);
	} while (!fin.eof());

	return positions;
}

void setup_tissue( void )
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
	
	// load cells from your CSV file (if enabled)

	load_cells_from_pugixml();



	
	return; 
}


std::vector<std::string> regular_colors( Cell* pCell )
{
	static int A_type = get_cell_definition( "1" ).type; 
	static int B_type = get_cell_definition( "2" ).type; 
	static int C_type = get_cell_definition( "3" ).type; 
	
	// start with flow cytometry coloring 
	
	std::vector<std::string> output = {"black" , "black" , "black" , "black"} ;

	// color live C 

	if( pCell->type == A_type )
	{
		 output[0] = "blue";  
		 output[2] = "blue";  
	}
	
	// color live B

	if( pCell->type == B_type )
	{
		 output[0] = "orange";  
		 output[2] = "orange";  
	}
	
	// color live C

	if( pCell->type == C_type )
	{
		 output[0] = "green";  
		 output[2] = "green";  
	}
	

	return output; 
}

void set_substrate_density( void )
{

	// double minim = 20;
	double maxim = 30;

	// std::ifstream file( "/home/alejandro/Escritorio/MasterBioinformatica/TFM/TFM2022AlejandroMadrid/pruebas/Stutilityvsscanpy/setup_physi/setup_ecm.csv", std::ios::in );
	// if( !file )
	// { 
	// 	std::cout << "Error: " << "set_ecm.csv" << " not found during cell loading. Quitting." << std::endl; 
	// 	exit(-1);
	// }

	std::vector<std::vector<double>> positions;
	// if ( parameters.bools("read_init") )	{
		std::string csv_fname = parameters.strings("ecm_density_file");
		positions = read_cells_positions(csv_fname, '\t', true);
	// }

	for (int i = 0; i < positions.size(); i++)
	{
		int x = (positions[i][0]);
		int y = ( positions[i][1]);
		int z = ( positions[i][2]);
		microenvironment.density_vector(microenvironment.voxel_index(x,y,z))[microenvironment.find_density_index("ecm")] = maxim;
	}

	// std::string line;
	// while (std::getline(positions, line))
	// {
	// 	std::vector<double> data2;
	// 	csv_to_vector( line.c_str() , data2 ); 		
	// 	if( data2.size() != 3 )
	// 	{
	// 		std::cout << "Error! Importing cells from a CSV file expects each row to be x,y,z." << std::endl;
	// 		exit(-1);
	// 	};
		/*std::vector<double> position = { data[0] , data[1] , data[2] }; */
		//std::cout << microenvironment.density_vector(microenvironment.voxel_index(data2[0],data2[1],data2[2]))[microenvironment.find_density_index("ecm")] << "kweeeeeeeeeeee" << std::endl;
		// microenvironment.density_vector(microenvironment.voxel_index(data2[0],data2[1],data2[2]))[microenvironment.find_density_index("ecm")] = maxim;
		//std::cout << microenvironment.density_vector(microenvironment.voxel_index(data2[0],data2[1],data2[2]))[microenvironment.find_density_index("ecm")] << "kweeeeeeeeeeee" << std::endl;
		//std::cout << microenvironment.find_density_index("ecm") << std::endl;
	// }

	// file.close(); 	
	return;
	// std::cout << microenvironment.number_of_voxels() << "\n";
	// #pragma omp parallel for
	// for (int n = 0; n < microenvironment.number_of_voxels(); n++)
	// {
	// 	auto current_voxel = microenvironment.voxels(n);
	// 	double t_norm = norm(current_voxel.center);

	// 	if ((radius - t_norm) <= 0)
	// 		microenvironment.density_vector(n)[density_index] = current_value(min, max, uniform_random());
	// }
}


void phenotype_function( Cell* pCell, Phenotype& phenotype, double dt )
{ return; }

void custom_function( Cell* pCell, Phenotype& phenotype , double dt )
{ return; } 

void contact_function( Cell* pMe, Phenotype& phenoMe , Cell* pOther, Phenotype& phenoOther , double dt )
{ return; } 

double current_value( double min, double max, double percent )
{ return (min + (max-min) * percent); };

double add_ecm_interaction(Cell* pC, int index_ecm, int index_voxel )
{
	// Check if there is ECM material in given voxel
	// double dens2 = pC->get_microenvironment()->density_vector(index_voxel)[index_ecm];
	double dens = pC->get_microenvironment()->nearest_density_vector(index_voxel)[index_ecm];
	double ecmrad = sqrt(3.0) * pC->get_microenvironment()->mesh.dx / 2.0;
	// if voxel is "full", density is 1
	dens = std::min( dens, 1.0 ); 
	if ( dens > EPSILON )
	{
		// Distance between agent center and ECM voxel center
		pC->displacement = pC->position - pC->get_container()->underlying_mesh.voxels[index_voxel].center;
		double distance = norm(pC->displacement);
		// Make sure that the distance is not zero
		distance = std::max(distance, EPSILON);
		
		double dd = pC->phenotype.geometry.radius + ecmrad;  
		double dnuc = pC->phenotype.geometry.nuclear_radius + ecmrad;  

		double tmp_r = 0;
		// Cell overlap with ECM node, add a repulsion term
		if ( distance < dd )
		{
			// repulsion stronger if nucleii overlap, see Macklin et al. 2012, 2.3.1
			if ( distance < dnuc )
			{
				double M = 1.0;
				double c = 1.0 - dnuc/dd;
				c *= c;
				c -= M;
				tmp_r = c*distance/dnuc + M;
				pC->custom_data["nucleus_deform"] += (dnuc-distance);
			}
			else
			{
				tmp_r = ( 1 - distance / dd );
				tmp_r *= tmp_r;
			}

			tmp_r *= dens * PhysiCell::parameters.doubles("cell_ecm_repulsion");
		}

		// // Cell adherence to ECM through integrins
		// double max_interactive_distance = (PhysiCell::parameters.doubles("max_interaction_factor")*pC->phenotype.geometry.radius) + ecmrad;
		// if ( distance < max_interactive_distance ) 
		// {	
		// 	double temp_a = 1 - distance/max_interactive_distance; 
		// 	temp_a *= temp_a; 
		// 	/* \todo change dens with a maximal density ratio ? */
		// 	pC->custom_data["ecm_contact"] += dens * (max_interactive_distance-distance);
		// 	// temp_a *= dens * ( static_cast<Cell*>(this) )->integrinStrength();

		// 	double temp_integrins = get_integrin_strength( pC->custom_data["integrin"] );

		// 	temp_a *= dens * temp_integrins;
			
		// 	tmp_r -= temp_a;
		// }
		
		/////////////////////////////////////////////////////////////////
		// if (tmp_r==0) 
		// 	return;

		tmp_r/=distance;

	
		double zeroo = 0;

		if (tmp_r != 0){
			
			pC->velocity = zeroo * pC->displacement;
			pC -> phenotype.motility.motility_vector.assign( 0, 0);

		}


		
		// std::cout << "JJJJJJJJJJ" << pC->velocity << std::endl;

		return tmp_r;
	}
}

void custom_update_velocity( Cell* pCell, Phenotype& phenotype, double dt)
{
	pCell->custom_data["ecm_contact"] = 0;
	pCell->custom_data["nucleus_deform"] = 0;
	
	if( pCell->functions.add_cell_basement_membrane_interactions )
	{
		pCell->functions.add_cell_basement_membrane_interactions(pCell, phenotype,dt);
	}
	
	pCell->state.simple_pressure = 0.0; 
	
	//First check the neighbors in my current voxel
	for( auto neighbor : pCell->get_container()->agent_grid[pCell->get_current_mechanics_voxel_index()] )
	{
		pCell->add_potentials( neighbor );
	}

	pCell->state.simple_pressure = 0.0; 
	pCell->state.neighbors.clear(); // new 1.8.0

	std::vector<Cell*>::iterator neighbor;
	std::vector<Cell*>::iterator end = pCell->get_container()->agent_grid[pCell->get_current_mechanics_voxel_index()].end();
	for(neighbor = pCell->get_container()->agent_grid[pCell->get_current_mechanics_voxel_index()].begin(); neighbor != end; ++neighbor)
	{
		pCell->add_potentials(*neighbor);
	}
	std::vector<int>::iterator neighbor_voxel_index;
	std::vector<int>::iterator neighbor_voxel_index_end = 
		pCell->get_container()->underlying_mesh.moore_connected_voxel_indices[pCell->get_current_mechanics_voxel_index()].end();

	for( neighbor_voxel_index = 
		pCell->get_container()->underlying_mesh.moore_connected_voxel_indices[pCell->get_current_mechanics_voxel_index()].begin();
		neighbor_voxel_index != neighbor_voxel_index_end; 
		++neighbor_voxel_index )
	{
		if(!is_neighbor_voxel(pCell, pCell->get_container()->underlying_mesh.voxels[pCell->get_current_mechanics_voxel_index()].center, pCell->get_container()->underlying_mesh.voxels[*neighbor_voxel_index].center, *neighbor_voxel_index))
			continue;
		end = pCell->get_container()->agent_grid[*neighbor_voxel_index].end();
		for(neighbor = pCell->get_container()->agent_grid[*neighbor_voxel_index].begin();neighbor != end; ++neighbor)
		{
			pCell->add_potentials(*neighbor);
		}
	}

	int ecm_index = BioFVM::microenvironment.find_density_index("ecm");
	if ( ecm_index >= 0 ){
		double tmp_r2 = add_ecm_interaction( pCell, ecm_index, pCell->get_current_mechanics_voxel_index() );
		double owo = 0;

		if ( tmp_r2 != 0 ){
			// pCell -> phenotype.motility.motility_vector.assign( 0, 0);
			// std::cout << pCell -> velocity << "uwuwuwuwuw" << std::endl;
			pCell->velocity = owo * pCell->displacement;
			
		} else {

			
			pCell->update_motility_vector(dt);
	//std::cout << pCell->state.simple_pressure << " \n ";
			pCell->velocity += phenotype.motility.motility_vector;
			
		}
		
		return;
	}
/*
	if (pCell->custom_data["freezed"] > 2){
		return ;
	}
*/
	


	
	return; 
}









void SVG_plot_ecm( std::string filename , Microenvironment& M, double z_slice , double time, std::vector<std::string> (*cell_coloring_function)(Cell*), std::string sub )
{




	double X_lower = M.mesh.bounding_box[0];
	double X_upper = M.mesh.bounding_box[3];
 
	double Y_lower = M.mesh.bounding_box[1]; 
	double Y_upper = M.mesh.bounding_box[4]; 

	double plot_width = X_upper - X_lower; 
	double plot_height = Y_upper - Y_lower; 

	double font_size = 0.025 * plot_height; // PhysiCell_SVG_options.font_size; 
	double top_margin = font_size*(.2+1+.2+.9+.5 ); 

	// open the file, write a basic "header"
	std::ofstream os( filename , std::ios::out );
	if( os.fail() )
	{ 
		std::cout << std::endl << "Error: Failed to open " << filename << " for SVG writing." << std::endl << std::endl; 

		std::cout << std::endl << "Error: We're not writing data like we expect. " << std::endl
		<< "Check to make sure your save directory exists. " << std::endl << std::endl
		<< "I'm going to exit with a crash code of -1 now until " << std::endl 
		<< "you fix your directory. Sorry!" << std::endl << std::endl; 
		exit(-1); 
	} 
	
	Write_SVG_start( os, plot_width , plot_height + top_margin );

	// draw the background 
	Write_SVG_rect( os , 0 , 0 , plot_width, plot_height + top_margin , 0.002 * plot_height , "white", "white" );

	// write the simulation time to the top of the plot
 
	char* szString; 
	szString = new char [1024]; 
 
	int total_cell_count = all_cells->size(); 
 
	double temp_time = time; 

	std::string time_label = formatted_minutes_to_DDHHMM( temp_time ); 
 
	sprintf( szString , "Current time: %s, z = %3.2f %s", time_label.c_str(), 
		z_slice , PhysiCell_SVG_options.simulation_space_units.c_str() ); 
	Write_SVG_text( os, szString, font_size*0.5,  font_size*(.2+1), 
		font_size, PhysiCell_SVG_options.font_color.c_str() , PhysiCell_SVG_options.font.c_str() );
	sprintf( szString , "%u agents" , total_cell_count ); 
	Write_SVG_text( os, szString, font_size*0.5,  font_size*(.2+1+.2+.9), 
		0.95*font_size, PhysiCell_SVG_options.font_color.c_str() , PhysiCell_SVG_options.font.c_str() );
	
	delete [] szString; 


	// add an outer "g" for coordinate transforms 
	
	os << " <g id=\"tissue\" " << std::endl 
	   << "    transform=\"translate(0," << plot_height+top_margin << ") scale(1,-1)\">" << std::endl; 
	   
	// prepare to do mesh-based plot (later)
	
	double dx_stroma = M.mesh.dx; 
	double dy_stroma = M.mesh.dy; 
	double dz_stroma = M.mesh.dz;
	
	os << "  <g id=\"ECM\">" << std::endl; 
  
	int ratio = 1; 
	double voxel_size = dx_stroma / (double) ratio ; 
  
	double half_voxel_size = voxel_size / 2.0; 
	double normalizer = 78.539816339744831 / (voxel_size*voxel_size*voxel_size); 
	
	// color the dark background 

	int sub_index = M.find_density_index(sub);

	// to add for loop to look for mac conc.

	double max_conc = 0.0;
	//look for the mac conc among all the substrates
	#pragma omp parallel for
	for (int n = 0; n < M.number_of_voxels(); n++)
	{
		double concentration = M.density_vector(n)[sub_index];
		if (concentration > max_conc)
			max_conc = concentration;
	}

	//double max_conc = default_microenvironment_options.initial_condition_vector[sub_index];

	if(max_conc == 0){

		std::cout << "it is not possible to correctly print the substrate, make sure to indicate the max value of your substrate in 'initial_condition' in the microenv section of your xml" << std::endl;

		max_conc = 1.0;

	};

	for (int n = 0; n < M.number_of_voxels(); n++)
	{
		auto current_voxel = M.voxels(n);
		int z_center = current_voxel.center[2];
		double z_displ = z_center -  dz_stroma/2;
		
		double z_compare = z_displ;

		if (default_microenvironment_options.simulate_2D == true){
		z_compare = z_center;
		};

		if (z_slice == z_compare){
			int x_center = current_voxel.center[0];
			int y_center = current_voxel.center[1];
			
			double x_displ = x_center -  dx_stroma/2;
			double y_displ = (y_center - dy_stroma) +  dy_stroma/2;
			//std::cout <<  x_displ - X_lower << "  __  " << y_displ - Y_lower << "\n" ;

			std::vector< std::string > output( 4 , "black" );

			double concentration = M.density_vector(n)[sub_index];

			int color = (int) round( (concentration / max_conc) * 255.0 );
			if(color > 255){
				color = 255;
			}
			char szTempString [128];
			sprintf( szTempString , "rgb(%u,234,197)", 255 - color);
			output[0].assign( szTempString );

			// std::vector< std::string > output( 4 , "black" );
			// int color = (int) round( concentration * 255.0 );
			// if(color > 255){
			// 	color = 255;
			// }
			//char szTempString [128];

			int green = 255 - color;
			int blue = 255 - color;
			int red = 255 - color;


			//sprintf( szTempString , "rgb(%u, %u, %u)", red, green, blue );
			//output[0].assign( szTempString );

			Write_SVG_rect( os , x_displ - X_lower , y_displ - Y_lower, dx_stroma, dy_stroma , 0 , "none", output[0] );
		}

	}

	// Write_SVG_rect( os , 0 , 0 , plot_width, plot_height , 0 , "none", "black" );

 
 // color in the background ECM
/* 
 if( ECM.TellRows() > 0 )
 {
  // find the k corresponding to z_slice
  
  
  
  Vector position; 
  *position(2) = z_slice; 
  

  // 25*pi* 5 microns^2 * length (in source) / voxelsize^3
  
  for( int j=0; j < ratio*ECM.TellCols() ; j++ )
  {
   // *position(1) = *Y_environment(j); 
   *position(1) = *Y_environment(0) - dy_stroma/2.0 + j*voxel_size + half_voxel_size; 
   
   for( int i=0; i < ratio*ECM.TellRows() ; i++ )
   {
    // *position(0) = *X_environment(i); 
    *position(0) = *X_environment(0) - dx_stroma/2.0 + i*voxel_size + half_voxel_size; 
	
    double E = evaluate_Matrix3( ECM, X_environment , Y_environment, Z_environment , position );	
	double BV = normalizer * evaluate_Matrix3( OxygenSourceHD, X_environment , Y_environment, Z_environment , position );
	if( isnan( BV ) )
	{ BV = 0.0; }

	vector<string> Colors;
	Colors = hematoxylin_and_eosin_stroma_coloring( E , BV );
	Write_SVG_rect( os , *position(0)-half_voxel_size-X_lower , *position(1)-half_voxel_size+top_margin-Y_lower, 
	voxel_size , voxel_size , 1 , Colors[0], Colors[0] );
   
   }
  }
 
 }
*/
	os << "  </g>" << std::endl; 

 
	// plot intersecting cells 
	os << "  <g id=\"cells\">" << std::endl; 


	double lowest_pressure = 1.0;
	double max_pressure = 1.0;
	for( int i=0 ; i < total_cell_count ; i++ )
	{
		Cell* pC = (*all_cells)[i]; // global_cell_list[i]; 

		double pressure = pC->custom_data["padhesion"];

		if (pressure < lowest_pressure)
			lowest_pressure = pressure;
		if (pressure > max_pressure)
			max_pressure = pressure;

		static std::vector<std::string> Colors; 
		if( fabs( (pC->position)[2] - z_slice ) < pC->phenotype.geometry.radius )
		{
			double r = pC->phenotype.geometry.radius ; 
			double rn = pC->phenotype.geometry.nuclear_radius ; 
			double z = fabs( (pC->position)[2] - z_slice) ; 
   
			Colors = cell_coloring_function( pC ); 

			os << "   <g id=\"cell" << pC->ID << "\">" << std::endl; 
  
			// figure out how much of the cell intersects with z = 0 
   
			double plot_radius = sqrt( r*r - z*z ); 

			Write_SVG_circle( os, (pC->position)[0]-X_lower, (pC->position)[1]-Y_lower, 
				plot_radius , 0.5, Colors[1], Colors[0] ); 

			// plot the nucleus if it, too intersects z = 0;
			if( fabs(z) < rn && PhysiCell_SVG_options.plot_nuclei == true )
			{   
				plot_radius = sqrt( rn*rn - z*z ); 
			 	Write_SVG_circle( os, (pC->position)[0]-X_lower, (pC->position)[1]-Y_lower, 
					plot_radius, 0.5, Colors[3],Colors[2]); 
			}					  
			os << "   </g>" << std::endl;
		}
	}

	//std::cout << max_pressure << "  --  " << lowest_pressure;

	os << "  </g>" << std::endl; 
	
	// end of the <g ID="tissue">
	os << " </g>" << std::endl; 
 
	// draw a scale bar
 
	double bar_margin = 0.025 * plot_height; 
	double bar_height = 0.01 * plot_height; 
	double bar_width = PhysiCell_SVG_options.length_bar; 
	double bar_stroke_width = 0.001 * plot_height; 
	
	std::string bar_units = PhysiCell_SVG_options.simulation_space_units; 
	// convert from micron to mm
	double temp = bar_width;  

	if( temp > 999 && std::strstr( bar_units.c_str() , PhysiCell_SVG_options.mu.c_str() )   )
	{
		temp /= 1000;
		bar_units = "mm";
	}
	// convert from mm to cm 
	if( temp > 9 && std::strcmp( bar_units.c_str() , "mm" ) == 0 )
	{
		temp /= 10; 
		bar_units = "cm";
	}
	
	szString = new char [1024];
	sprintf( szString , "%u %s" , (int) round( temp ) , bar_units.c_str() );
 
	Write_SVG_rect( os , plot_width - bar_margin - bar_width  , plot_height + top_margin - bar_margin - bar_height , 
		bar_width , bar_height , 0.002 * plot_height , "rgb(255,255,255)", "rgb(0,0,0)" );
	Write_SVG_text( os, szString , plot_width - bar_margin - bar_width + 0.25*font_size , 
		plot_height + top_margin - bar_margin - bar_height - 0.25*font_size , 
		font_size , PhysiCell_SVG_options.font_color.c_str() , PhysiCell_SVG_options.font.c_str() ); 
	
	delete [] szString; 

	// plot runtime 
	szString = new char [1024]; 
	RUNTIME_TOC(); 
	std::string formatted_stopwatch_value = format_stopwatch_value( runtime_stopwatch_value() );
	Write_SVG_text( os, formatted_stopwatch_value.c_str() , bar_margin , top_margin + plot_height - bar_margin , 0.75 * font_size , 
		PhysiCell_SVG_options.font_color.c_str() , PhysiCell_SVG_options.font.c_str() );
	delete [] szString; 

	// draw a box around the plot window
	Write_SVG_rect( os , 0 , top_margin, plot_width, plot_height , 0.002 * plot_height , "rgb(0,0,0)", "none" );
	
	// close the svg tag, close the file
	Write_SVG_end( os ); 
	os.close();
 
	return; 
}