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
	if (parameters.ints.find_index("random_seed") != -1)
	{
		SeedRandom(parameters.ints("random_seed"));
	}
	
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


	cell_defaults.functions.custom_cell_rule = custom_function; 
	cell_defaults.functions.contact_function = contact_function; 

	Cell_Definition* pCD = find_cell_definition( "green");


	pCD->functions.custom_cell_rule = custom_function; 
	pCD->functions.contact_function = contact_function; 
	pCD->functions.update_velocity = standard_update_cell_velocity; 
	pCD->functions.pre_update_intracellular =  pre_update_intracellular;
	pCD->functions.post_update_intracellular = post_update_intracellular;

	pCD = find_cell_definition( "blue");

	pCD->functions.custom_cell_rule = custom_function; 
	pCD->functions.contact_function = contact_function;
	pCD->functions.update_velocity = standard_update_cell_velocity; 
	pCD->functions.pre_update_intracellular =  pre_update_intracellular;
	pCD->functions.post_update_intracellular = post_update_intracellular;


	pCD = find_cell_definition( "orange");

	pCD->functions.custom_cell_rule = custom_function; 
	pCD->functions.contact_function = contact_function;
	pCD->functions.update_velocity = standard_update_cell_velocity; 
	pCD->functions.pre_update_intracellular =  pre_update_intracellular;
	pCD->functions.post_update_intracellular = post_update_intracellular;
	
	/*
	   This builds the map of cell definitions and summarizes the setup. 
	*/
		
	display_cell_definitions( std::cout ); 
	
	return; 
}

void setup_microenvironment( void )
{
	initialize_microenvironment();
	
	double o2_conc;
	double drug_conc;
	
	std::vector<double> dirichlet_o2( 2 ); // number of sustrates, creates the vector of 3 positions

	o2_conc = parameters.doubles("o2_conc1");
	drug_conc = parameters.doubles("drug_conc");
		
	dirichlet_o2[0] = o2_conc;
	dirichlet_o2[1] = drug_conc;
	
	std::vector<std::vector<double>> positions;

	std::string csv_fname = parameters.strings("blood_source_file");
	positions = read_cells_positions(csv_fname, '\t', true);

	for (int i = 0; i < positions.size(); i++)
	{
		int x = (positions[i][0]);
		int y = ( positions[i][1]);
		int z = ( positions[i][2]);
		microenvironment.add_dirichlet_node( microenvironment.voxel_index(x,y,z) , dirichlet_o2 );
		//microenvironment.set_substrate_dirichlet_activation( 1,  microenvironment.voxel_index(x,y,z),true );
		// microenvironment.add_dirichlet_node( microenvironment.voxel_index(data[0],data[1],data[2]) , dirichlet_drug );

		//microenvironment.set_substrate_dirichlet_activation( 2,  microenvironment.voxel_index(x,y,z),true );
		// microenvironment.add_dirichlet_node( microenvironment.voxel_index(positions[0],positions[1],positions[2]) , dirichlet_o2 );
	}
	
	return;
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

void change_dirichlet_nodes ( void )
{

	double o2_conc;
	double drug_conc;
	
	std::vector<double> dirichlet_o2( 2 ); // number of sustrates, creates the vector of 3 positions

	o2_conc = parameters.doubles("o2_conc1");
	drug_conc = parameters.doubles("drug_conc");
		
	dirichlet_o2[0] = o2_conc;
	dirichlet_o2[1] = drug_conc;
	
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
	
	std::cout << "Adding drug." << std::endl;
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


void custom_function( Cell* pCell, Phenotype& phenotype , double dt )
{ 	

	pCell->custom_data["cell_contact"] = 0.0;
	for( int j=0; j < pCell->nearby_interacting_cells().size(); j++ )
    {
        Cell* pTest = pCell->nearby_interacting_cells()[j]; 
		contact_function(pCell, phenotype, pTest, pTest->phenotype, dt);
    }
	
	for( int j=0; j < pCell->state.spring_attachments.size(); j++ )
    {
        Cell* pTest = pCell->state.spring_attachments[j]; 
        contact_function(pCell, phenotype, pTest, pTest->phenotype, dt);
    }

	// ADDING ECM PHYSICAL INTERACTION AND ADHESION

	pCell->custom_data["ecm_contact"] = 0.0;
	pCell->custom_data["nucleus_deform"] = 0.0;
	//std::cout << pCell->custom_data["nucleus_deform"] << std::endl;

	int ecm_index = BioFVM::microenvironment.find_density_index("ecm");
	if ( ecm_index >= 0 ){
		add_ecm_interaction( pCell, ecm_index, pCell->get_current_voxel_index() );
		//add_TGFbeta_interaction(pCell, pCell->get_current_mechanics_voxel_index());
		std::vector<int>::iterator neighbor_voxel_index;
		std::vector<int>::iterator neighbor_voxel_index_end = 
		microenvironment.mesh.moore_connected_voxel_indices[pCell->get_current_voxel_index()].end();

		for( neighbor_voxel_index = 
			microenvironment.mesh.moore_connected_voxel_indices[pCell->get_current_voxel_index()].begin();
			neighbor_voxel_index != neighbor_voxel_index_end; 
			++neighbor_voxel_index )
		{
			add_ecm_interaction( pCell, ecm_index, *neighbor_voxel_index );
			
		}

		/* pCell->update_motility_vector(dt); 
		pCell->velocity += phenotype.motility.motility_vector; */
	}
	
	return; 	
} 

// This function is never called because I am not using "dynamic attachment" but I am using the "dynamic SPRING adhesion"
void contact_function( Cell* pMe, Phenotype& phenoMe , Cell* pOther, Phenotype& phenoOther , double dt )
{ 

	std::vector<double> displacement = pOther->position;
	displacement -= pMe->position;
	double distance = norm( displacement ); 
			
	double max_distance = pMe->phenotype.geometry.radius + 
				pOther->phenotype.geometry.radius; 
	max_distance *=  pMe->phenotype.mechanics.relative_maximum_adhesion_distance;  //parameters.doubles("max_interaction_factor"); 

			//std::cout << max_distance << " - " << distance << "\n";

	double interaction_distance = max_distance - distance;

	if (interaction_distance > 0){

		double perc_distance = distance / pMe->phenotype.geometry.radius ;
		pMe->custom_data["cell_contact"] += perc_distance;
			}
	else {
		detach_cells_as_spring(pMe, pOther);
	}

	return; 
} 




void add_ecm_interaction(Cell* pC, int index_ecm, int index_voxel )
{
	// Check if there is ECM material in given voxel
	//double dens2 = get_microenvironment()->density_vector(index_voxel)[index_ecm];
	double dens = pC->get_microenvironment()->nearest_density_vector(index_voxel)[index_ecm];
	double ecmrad = sqrt(3.0) * pC->get_microenvironment()->mesh.dx * 0.5;
	// if voxel is "full", density is 1
	dens = std::min( dens, 1.0 ); 
	if ( dens > EPSILON )
	{
		// Distance between agent center and ECM voxel center
		pC->displacement = pC->position - microenvironment.mesh.voxels[index_voxel].center;
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

		
		
		/////////////////////////////////////////////////////////////////
		if(tmp_r==0)
			return;
		tmp_r/=distance;

		axpy( &pC->velocity , tmp_r , pC->displacement ); 
	}

}






std::vector<std::string> my_coloring_function( Cell* pCell )
{ return paint_by_number_cell_coloring(pCell); }

std::vector<std::string> regular_colors( Cell* pCell )
{
	static int A_type = get_cell_definition( "green" ).type; 
	static int B_type = get_cell_definition( "blue" ).type; 
	static int C_type = get_cell_definition( "orange" ).type; 
	
	// start with flow cytometry coloring 
	
	std::vector<std::string> output = {"black" , "black" , "black" , "black"} ;

	// color live a
		
	if( pCell->type == A_type )
	{
		 output[0] = "rgb(34,139,34)";  
		 output[2] = "rgb(34,139,34)";  
	}
	
	// color live B

	if( pCell->type == B_type )
	{
		 output[0] = "rgb(0,0,255)";  
		 output[2] = "rgb(0,0,255)";  
	}
	
	// color live C

	if( pCell->type == C_type )
	{
		 output[0] = "rgb(255,140,0)";  
		 output[2] = "rgb(255,140,0)";  
	}
	
	if( pCell->phenotype.death.dead == true )
	{
		// Necrotic - Brown
		if( pCell->phenotype.cycle.current_phase().code == PhysiCell_constants::necrotic_swelling || 
			pCell->phenotype.cycle.current_phase().code == PhysiCell_constants::necrotic_lysed || 
			pCell->phenotype.cycle.current_phase().code == PhysiCell_constants::necrotic )
		{
			output[2] = "chocolate";
		}
		else
		{
			output[2] = "black"; 
		}
	}

	return output; 
}

void phenotype_function( Cell* pCell, Phenotype& phenotype, double dt )
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
