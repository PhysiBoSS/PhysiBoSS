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
	SeedRandom(parameters.ints("random_seed"));

	initialize_default_cell_definition();

	/*  This parses the cell definitions in the XML config file.  */
	

	//  This sets the custom phenotype function
	cell_defaults.functions.update_phenotype = update_phenotype_with_signaling;

	cell_defaults.functions.update_velocity = custom_update_velocity;

	// esto nuevo

	cell_defaults.functions.volume_update_function = standard_volume_update_function;

	//  This initializes the the TNF receptor model

	tnf_boolean_model_interface_setup();
	submodel_registry.display(std::cout);

	// // Needs to initialize one of the receptor state to the total receptor value
	// cell_defaults.custom_data["unbound_external_TNFR"] = cell_defaults.custom_data["TNFR_receptors_per_cell"];
	// cell_defaults.custom_data["bound_external_TNFR"] = 0;
	// cell_defaults.custom_data["bound_internal_TNFR"] = 0;
	
	
	initialize_cell_definitions_from_pugixml();
	
	build_cell_definitions_maps();
	display_cell_definitions(std::cout);

	return;
}


void setup_microenvironment(void)
{
	initialize_microenvironment();
	
	double o2_conc;
	double drug_conc;
	double ecm_conc;
	
	std::vector<double> dirichlet_o2( 3 ); // number of sustrates, creates the vector of 3 positions

	o2_conc = parameters.doubles("o2_conc1");
	drug_conc = parameters.doubles("drug_conc");
	ecm_conc = parameters.doubles("ecm_conc");
		
	dirichlet_o2[0] = o2_conc;
	dirichlet_o2[1] = drug_conc;
	dirichlet_o2[2] = ecm_conc;
	
	std::vector<std::vector<double>> positions;

	std::string csv_fname = parameters.strings("blood_source_file");
	positions = read_cells_positions(csv_fname, '\t', true);

	for (int i = 0; i < positions.size(); i++)
	{
		int x = (positions[i][0]);
		int y = ( positions[i][1]);
		int z = ( positions[i][2]);
		microenvironment.add_dirichlet_node( microenvironment.voxel_index(x,y,z) , dirichlet_o2 );
		microenvironment.set_substrate_dirichlet_activation( 1,  microenvironment.voxel_index(x,y,z),false );
		// microenvironment.add_dirichlet_node( microenvironment.voxel_index(data[0],data[1],data[2]) , dirichlet_drug );

		//microenvironment.set_substrate_dirichlet_activation( 2,  microenvironment.voxel_index(x,y,z),true );
		// microenvironment.add_dirichlet_node( microenvironment.voxel_index(positions[0],positions[1],positions[2]) , dirichlet_o2 );
	}
	
	return;
}

void change_dirichlet_nodes ( void )
{


	std::vector<std::vector<double>> positions;

	std::string csv_fname = parameters.strings("blood_source_file");
	positions = read_cells_positions(csv_fname, '\t', true);

	for (int i = 0; i < positions.size(); i++)
	{
		int x = (positions[i][0]);
		int y = ( positions[i][1]);
		int z = ( positions[i][2]);
		
		microenvironment.set_substrate_dirichlet_activation( 1,  microenvironment.voxel_index(x,y,z),true );
		// microenvironment.add_dirichlet_node( microenvironment.voxel_index(data[0],data[1],data[2]) , dirichlet_drug );

		//microenvironment.set_substrate_dirichlet_activation( 2,  microenvironment.voxel_index(x,y,z),true );
		// microenvironment.add_dirichlet_node( microenvironment.voxel_index(positions[0],positions[1],positions[2]) , dirichlet_o2 );
	}
	
	std::cout << "voy" << std::endl;
	return;
}	
void setup_tissue(void)
{


	
	// load cells from your CSV file (if enabled)

	load_cells_from_pugixml();


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

std::vector<std::vector<double>> create_cell_sphere_positions(double cell_radius, double sphere_radius)
{
	std::vector<std::vector<double>> cells;
	int xc=0,yc=0,zc=0;
	double x_spacing= cell_radius*sqrt(3);
	double y_spacing= cell_radius*2;
	double z_spacing= cell_radius*sqrt(3);
	
	std::vector<double> tempPoint(3,0.0);
	// std::vector<double> cylinder_center(3,0.0);
	
	for(double z=-sphere_radius;z<sphere_radius;z+=z_spacing, zc++)
	{
		for(double x=-sphere_radius;x<sphere_radius;x+=x_spacing, xc++)
		{
			for(double y=-sphere_radius;y<sphere_radius;y+=y_spacing, yc++)
			{
				tempPoint[0]=x + (zc%2) * 0.5 * cell_radius;
				tempPoint[1]=y + (xc%2) * cell_radius;
				tempPoint[2]=z;
				
				if(sqrt(norm_squared(tempPoint))< sphere_radius)
				{ cells.push_back(tempPoint); }
			}
			
		}
	}
	return cells;
	
}

std::vector<std::vector<double>> create_cell_disc_positions(double cell_radius, double disc_radius)
{	 
	double cell_spacing = 0.95 * 2.0 * cell_radius; 
	
	double x = 0.0; 
	double y = 0.0; 
	double x_outer = 0.0;

	std::vector<std::vector<double>> positions;
	std::vector<double> tempPoint(3,0.0);
	
	int n = 0; 
	while( y < disc_radius )
	{
		x = 0.0; 
		if( n % 2 == 1 )
		{ x = 0.5 * cell_spacing; }
		x_outer = sqrt( disc_radius*disc_radius - y*y ); 
		
		while( x < x_outer )
		{
			tempPoint[0]= x; tempPoint[1]= y;	tempPoint[2]= 0.0;
			positions.push_back(tempPoint);			
			if( fabs( y ) > 0.01 )
			{
				tempPoint[0]= x; tempPoint[1]= -y;	tempPoint[2]= 0.0;
				positions.push_back(tempPoint);
			}
			if( fabs( x ) > 0.01 )
			{ 
				tempPoint[0]= -x; tempPoint[1]= y;	tempPoint[2]= 0.0;
				positions.push_back(tempPoint);
				if( fabs( y ) > 0.01 )
				{
					tempPoint[0]= -x; tempPoint[1]= -y;	tempPoint[2]= 0.0;
					positions.push_back(tempPoint);
				}
			}
			x += cell_spacing; 
		}		
		y += cell_spacing * sqrt(3.0)/2.0; 
		n++; 
	}
	return positions;
}

void inject_density_sphere(int density_index, double concentration, double membrane_lenght)
{
	// Inject given concentration on the extremities only
	#pragma omp parallel for
	for (int n = 0; n < microenvironment.number_of_voxels(); n++)
	{
		auto current_voxel = microenvironment.voxels(n);
		std::vector<double> cent = {current_voxel.center[0], current_voxel.center[1], current_voxel.center[2]};

		if ((membrane_lenght - norm(cent)) <= 0)
			microenvironment.density_vector(n)[density_index] = concentration;
	}
}

void remove_density(int density_index)
{
	for (int n = 0; n < microenvironment.number_of_voxels(); n++)
		microenvironment.density_vector(n)[density_index] = 0;
}

std::vector<std::string> my_coloring_function(Cell *pCell)
{
	// start with live coloring
	std::vector<std::string> output = false_cell_coloring_live_dead(pCell);

	return output;
}

std::vector<std::string> regular_colors( Cell* pCell )
{
	static int A_type = get_cell_definition( "green" ).type; 
	
	static int B_type = get_cell_definition( "blue" ).type;

	static int C_type = get_cell_definition( "orange" ).type; 
	
	// start with flow cytometry coloring 
	
	std::vector<std::string> output = {"black" , "black" , "black" , "black"} ;

	// color live C 

	if( pCell->type == A_type )
	{
		 output[0] = "green";  
		 output[2] = "green";  
	}
	
	if( pCell->type == B_type )
	{
		 output[0] = "blue";  
		 output[2] = "blue";  
	}

	if( pCell->type == C_type )
	{
		 output[0] = "orange";  
		 output[2] = "orange";  
	}
	
	return output; 
}

double total_live_cell_count()
{
        double out = 0.0;

        for( int i=0; i < (*all_cells).size() ; i++ )
        {
                if( (*all_cells)[i]->phenotype.death.dead == false && (*all_cells)[i]->type == 0 )
                { out += 1.0; }
        }

        return out;
}

double total_dead_cell_count()
{
        double out = 0.0;

        for( int i=0; i < (*all_cells).size() ; i++ )
        {
                if( (*all_cells)[i]->phenotype.death.dead == true && (*all_cells)[i]->phenotype.death.current_death_model_index == 0 )
                { out += 1.0; }
        }

        return out;
}

double total_necrosis_cell_count()
{
        double out = 0.0;

        for( int i=0; i < (*all_cells).size() ; i++ )
        {
                if( (*all_cells)[i]->phenotype.death.dead == true && (*all_cells)[i]->phenotype.death.current_death_model_index == 1 )
                { out += 1.0; }
        }

        return out;
}

void set_substrate_density( void )
{

	std::vector<std::vector<double>> positions;
	// if ( parameters.bools("read_init") )	{
		std::string csv_fname = parameters.strings("ecm_density_file");
		positions = read_cells_positions(csv_fname, '\t', true);
	// }
	double o2_conc;
	double drug_conc;
	double ecm_conc;
	if ( parameters.bools("read_init") )
	{
		o2_conc = parameters.doubles("o2_conc1");
		drug_conc = parameters.doubles("drug_conc");
		ecm_conc = parameters.doubles("ecm_conc");
	}
	else
	{
		double o2_conc = 30; //deuria de ser 38.06
	}
	std::vector<double> vector_substrates( 3 ); // number of sustrates, creates the vector of 3 positions

	vector_substrates[0] = o2_conc;
	vector_substrates[1] = drug_conc;
	vector_substrates[2] = ecm_conc;

	for (int i = 0; i < positions.size(); i++)
	{
		int x = (positions[i][0]);
		int y = ( positions[i][1]);
		int z = ( positions[i][2]);
		
		microenvironment.density_vector(microenvironment.voxel_index(x,y,z))[microenvironment.find_density_index("ecm")] = vector_substrates[2];
	}
	
	//microenvironment.density_vector(microenvironment.voxel_index(15,10,0))[microenvironment.find_density_index("drug_A")] = vector_substrates[2];
	return;
}

double current_value( double min, double max, double percent )
{ return (min + (max-min) * percent); };

// Calculate repulsion/adhesion between agent and ecm according to its local density, from Ruscone
void add_ecm_interaction(Cell* pC, int index_ecm, int index_voxel )
{
	// Check if there is ECM material in given voxel
	//double dens2 = get_microenvironment()->density_vector(index_voxel)[index_ecm];
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

		if(tmp_r==0)
			return;
		tmp_r/=distance;

		pC->velocity += tmp_r * pC->displacement;
	}
}

//from Ruscone
void custom_update_velocity( Cell* pCell, Phenotype& phenotype, double dt)
{
	// pCell->custom_data["ecm_contact"] = 0;
	// pCell->custom_data["nucleus_deform"] = 0;
	// pCell->custom_data["TGFbeta_contact"] = 0;
	// pCell->custom_data["cell_contact"] = 0;
	
	if( pCell->functions.add_cell_basement_membrane_interactions )
	{
		pCell->functions.add_cell_basement_membrane_interactions(pCell, phenotype,dt);
	}
	
	pCell->state.simple_pressure = 0.0;  //XXX per que?
	
	//First check the neighbors in my current voxel
	for( auto neighbor : pCell->get_container()->agent_grid[pCell->get_current_mechanics_voxel_index()] )
	{
		pCell->add_potentials( neighbor );
	}

	int ecm_index = BioFVM::microenvironment.find_density_index("ecm");
	if ( ecm_index >= 0 )
		add_ecm_interaction( pCell, ecm_index, pCell->get_current_mechanics_voxel_index() );
		// add_TGFbeta_interaction(pCell, pCell->get_current_mechanics_voxel_index());

	for (auto neighbor_voxel_index : pCell->get_container()->underlying_mesh.moore_connected_voxel_indices[pCell->get_current_mechanics_voxel_index()])
	{
		if(!is_neighbor_voxel(pCell, pCell->get_container()->underlying_mesh.voxels[pCell->get_current_mechanics_voxel_index()].center, pCell->get_container()->underlying_mesh.voxels[neighbor_voxel_index].center, neighbor_voxel_index))
			continue;

		if ( ecm_index >= 0 )
			add_ecm_interaction( pCell, ecm_index, neighbor_voxel_index );
			// add_TGFbeta_interaction(pCell, pCell->get_current_mechanics_voxel_index());
		for( auto other_neighbor : pCell->get_container()->agent_grid[neighbor_voxel_index] )
		{
			pCell->add_potentials(other_neighbor);
		}
	}

	pCell->update_motility_vector(dt);
	// std::cout << phenotype.motility.motility_vector << "  ";
	//std::cout << pCell->state.simple_pressure << " \n ";
	pCell->velocity += phenotype.motility.motility_vector;
	
	return; 
}
