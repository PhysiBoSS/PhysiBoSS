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
# Copyright (c) 2015-2022, Paul Macklin and the PhysiCell Project             #
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
# POSSIBILITY OF SUCH DAMAGE.                  F                               #
#                                                                             #
###############################################################################
*/

#include <cstdio>
#include <cstdlib>
#include <iostream>
#include <ctime>
#include <cmath>
#include <omp.h>
#include <fstream>

#include "./core/PhysiCell.h"
#include "./modules/PhysiCell_standard_modules.h" 

// put custom code modules here! 

#include "./custom_modules/custom.h" 
#include "./addons/PhysiBoSS/src/maboss_intracellular.h"

using namespace BioFVM;
using namespace PhysiCell;

int main( int argc, char* argv[] )
{

	std::cout << "Entering main..." << std::endl;

	// load and parse settings file(s)
	
	bool XML_status = false; 
	char copy_command [1024]; 
	char copy_command_2 [1024];
	if( argc > 1 )
	{
		XML_status = load_PhysiCell_config_file( argv[1] ); 
		sprintf( copy_command , "cp %s %s/PhysiCell_settings.xml" , argv[1] , PhysiCell_settings.folder.c_str() ); 
		sprintf( copy_command_2 , "cp %s %s" , argv[1] , PhysiCell_settings.folder.c_str() ); 
	}
	else
	{
		XML_status = load_PhysiCell_config_file( "./config/PhysiCell_settings.xml" );
		sprintf( copy_command , "cp ./config/PhysiCell_settings.xml %s" , PhysiCell_settings.folder.c_str() ); 
	}
	if( !XML_status )
	{ exit(-1); }
	
	// copy config file to output directry 
	system( copy_command ); 
	
	if ( argc > 1 )
	{
		system( copy_command_2 );
	}

	// OpenMP setup
	omp_set_num_threads(PhysiCell_settings.omp_num_threads);
	
	// PNRG setup 
	SeedRandom( parameters.ints("random_seed") ); // or specify a seed here
	
	// time setup 
	std::string time_units = "min"; 

	/* Microenvironment setup */ 
	
	setup_microenvironment(); // modify this in the custom code 


	
	/* PhysiCell setup */ 

	// set mechanics voxel size, and match the data structure to BioFVM
	double mechanics_voxel_size = 30; 
	Cell_Container* cell_container = create_cell_container_for_microenvironment( microenvironment, mechanics_voxel_size );

	
	/* Users typically start modifying here. START USERMODS */ 

	create_cell_types();
	// std::cout << "create cells OK" << std::endl;



	// update monitoring variables too in setup_tissue();
	setup_tissue();
	// std::cout << "setup tissue OK" << std::endl;



	// @oth: this is deprecated, delete it
	// Run the MaBoSS model prior to the simulation  
	// needs to be done before setup_tissue(), because that's where you call the update variables functions
	// int max_iter = 1;
	// #pragma omp parallel for 
	// for( int i=0; i < (*all_cells).size(); i++ )
	// {
	// 	// std::cout << "Running MaBoSS on agent:" << i << std::endl;
	// 	#pragma omp parallel for 
	// 	for (int j=0; j < max_iter; j++){
	// 		(*all_cells)[i]->phenotype.intracellular->update();
	// 		(*all_cells)[i]->phenotype.intracellular->print_current_nodes();
	// 		// std::cout << "Iteration " << j << std::endl;
	// 	}
	// }

	std::cout << "beginning maboss first save..." << std::endl;
	char filename2[1024];
	sprintf( filename2 , "%s/BM_initial_readouts.txt" , PhysiCell_settings.folder.c_str() ); 
	MaBoSSIntracellular::save( filename2, *PhysiCell::all_cells );
	std::cout << "maboss first save OK" << std::endl;


	/* Users typically stop modifying here. END USERMODS */ 
	
	// set MultiCellDS save options 

	set_save_biofvm_mesh_as_matlab( true ); 
	set_save_biofvm_data_as_matlab( true ); 
	set_save_biofvm_cell_data( true ); 
	set_save_biofvm_cell_data_as_custom_matlab( true );
	
	// save a simulation snapshot 
	
	char filename[1024];
	sprintf( filename , "%s/initial" , PhysiCell_settings.folder.c_str() ); 
	save_PhysiCell_to_MultiCellDS_v2( filename , microenvironment , PhysiCell_globals.current_time ); 
	
	sprintf( filename , "%s/states_initial.csv", PhysiCell_settings.folder.c_str());
	MaBoSSIntracellular::save( filename, *PhysiCell::all_cells);
	
	// save a quick SVG cross section through z = 0, after setting its 
	// length bar to 200 microns 

	PhysiCell_SVG_options.length_bar = 200; 

	// for simplicity, set a pathology coloring function 
	
	std::vector<std::string> (*cell_coloring_function)(Cell*) = paint_by_number_cell_coloring; 
	std::vector<std::string> (*ECM_coloring_function)(double, double, double) = my_coloring_function_for_stroma; 
	
	sprintf( filename , "%s/initial.svg" , PhysiCell_settings.folder.c_str() ); 
	SVG_plot( filename , microenvironment, 0.0 , PhysiCell_globals.current_time, cell_coloring_function );
	
	sprintf( filename , "%s/legend.svg" , PhysiCell_settings.folder.c_str() ); 
	create_plot_legend( filename , cell_coloring_function ); 

	add_software_citation( "PhysiBoSS" , PhysiBoSS_Version , PhysiBoSS_DOI, PhysiBoSS_URL); 
	
	display_citations(); 
	
	// set the performance timers 
	BioFVM::RUNTIME_TIC();
	BioFVM::TIC();
	
	std::ofstream report_file;
	if( PhysiCell_settings.enable_legacy_saves == true )
	{	
		sprintf( filename , "%s/simulation_report.txt" , PhysiCell_settings.folder.c_str() ); 
		
		report_file.open(filename); 	// create the data log file 
		report_file<<"simulated time\tnum cells\tnum division\tnum death\twall time"<<std::endl;
	}

	// main loop 

	try 
	{		
		while( PhysiCell_globals.current_time < PhysiCell_settings.max_time + 0.1*diffusion_dt )
		{
			// save data if it's time. 
			if( fabs( PhysiCell_globals.current_time - PhysiCell_globals.next_full_save_time ) < 0.01 * diffusion_dt )
			{
				std::cout << "Time to save" << std::endl;
				display_simulation_status( std::cout ); 
				if( PhysiCell_settings.enable_legacy_saves == true )
				{	
					log_output( PhysiCell_globals.current_time , PhysiCell_globals.full_output_index, microenvironment, report_file);
				}
				
				if( PhysiCell_settings.enable_full_saves == true )
				{	
					sprintf( filename , "%s/output%08u" , PhysiCell_settings.folder.c_str(),  PhysiCell_globals.full_output_index ); 
					
					save_PhysiCell_to_MultiCellDS_v2( filename , microenvironment , PhysiCell_globals.current_time ); 

					sprintf( filename , "%s/states_%08u.csv", PhysiCell_settings.folder.c_str(), PhysiCell_globals.full_output_index);
					
					MaBoSSIntracellular::save( filename, *PhysiCell::all_cells );
				}
				
				PhysiCell_globals.full_output_index++; 
				PhysiCell_globals.next_full_save_time += PhysiCell_settings.full_save_interval;
			}
			
			// save SVG plot if it's time
			if( fabs( PhysiCell_globals.current_time - PhysiCell_globals.next_SVG_save_time  ) < 0.01 * diffusion_dt )
			{
				if( PhysiCell_settings.enable_SVG_saves == true )
				{	
					sprintf( filename , "%s/snapshot%08u.svg" , PhysiCell_settings.folder.c_str() , PhysiCell_globals.SVG_output_index ); 
					SVG_plot( filename , microenvironment, 0.0 , PhysiCell_globals.current_time, cell_coloring_function, ECM_coloring_function);
					
					PhysiCell_globals.SVG_output_index++; 
					PhysiCell_globals.next_SVG_save_time  += PhysiCell_settings.SVG_save_interval;
				}
			}

			/*
			  Custom add-ons could potentially go here. 
			*/	

			// @oth: #TODO: complete function from template_BM
			treatment_function();

			// std::cout << "added treatment function" << std::endl;
			
			// update the microenvironment
			microenvironment.simulate_diffusion_decay( diffusion_dt );

			// std::cout << "simulated diffusion" << std::endl;


			// update te AGS receptor model of each cell 
			drug_transport_model_main( diffusion_dt );


			// std::cout << "about to run pc" << std::endl;
			
			// run PhysiCell 
			((Cell_Container *)microenvironment.agent_container)->update_all_cells( PhysiCell_globals.current_time );
			
			/*
			  Custom add-ons could potentially go here. 
			*/
			
			PhysiCell_globals.current_time += diffusion_dt;
		}

		if( PhysiCell_settings.enable_legacy_saves == true )
		{			
			log_output(PhysiCell_globals.current_time, PhysiCell_globals.full_output_index, microenvironment, report_file);
			report_file.close();
		}
	}
	catch( const std::exception& e )
	{ // reference to the base of a polymorphic object
		std::cout << e.what(); // information from length_error printed
	}
	
	// save a final simulation snapshot 
	
	sprintf( filename , "%s/final" , PhysiCell_settings.folder.c_str() ); 
	save_PhysiCell_to_MultiCellDS_v2( filename , microenvironment , PhysiCell_globals.current_time ); 
	
	sprintf( filename , "%s/states_final.csv", PhysiCell_settings.folder.c_str());
	MaBoSSIntracellular::save( filename, *PhysiCell::all_cells );
	
	sprintf( filename , "%s/final.svg" , PhysiCell_settings.folder.c_str() ); 
	SVG_plot( filename , microenvironment, 0.0 , PhysiCell_globals.current_time, cell_coloring_function, ECM_coloring_function);
	
	// timer 
	std::cout << std::endl << "Total simulation runtime: " << std::endl; 
	BioFVM::display_stopwatch_value( std::cout , BioFVM::runtime_stopwatch_value() ); 

	return 0; 
}