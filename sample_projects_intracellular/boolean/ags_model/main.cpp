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

#include <cstdio>
#include <cstdlib>
#include <iostream>
#include <ctime>
#include <cmath>
#include <omp.h>
#include <fstream>

#include "./core/PhysiCell.h"
#include "./modules/PhysiCell_standard_modules.h" 
#include "./addons/PhysiBoSS/src/maboss_intracellular.h"
#include "./addons/PhysiBoSS/src/maboss_network.h"

// put custom code modules here! 
#include "./custom_modules/custom.h" 


using namespace BioFVM;
using namespace PhysiCell;

int main( int argc, char* argv[] )
{

	std::cout << "Entering main..." << std::endl;

	// load and parse settings file(s)
	
	bool XML_status = false; 
	if( argc > 1 )
	{ XML_status = load_PhysiCell_config_file( argv[1] ); }
	else
	{ XML_status = load_PhysiCell_config_file( "./config/PhysiCell_settings.xml" ); }
	if( !XML_status )
	{ exit(-1); }

	// OpenMP setup
	omp_set_num_threads(PhysiCell_settings.omp_num_threads);
	
	// PNRG setup 
	SeedRandom( parameters.ints("random_seed") ); // or specify a seed here
	
	
	// time setup 
	std::string time_units = "min"; 

	/* Microenvironment setup */ 
	
	setup_microenvironment(); // modify this in the custom code 
	std::cout << "setup OK" << std::endl;

	// User parameters 
	double drug_X_pulse_period = parameters.doubles("drug_X_pulse_period");
	double drug_X_pulse_duration = parameters.doubles("drug_X_pulse_duration");
	double drug_X_pulse_concentration = parameters.doubles("drug_X_pulse_concentration");
	double membrane_lenght = parameters.doubles("membrane_length"); // radious around which the drug_X pulse is injected

	double drug_X_pulse_timer = drug_X_pulse_period;
	double drug_X_pulse_injection_timer = -1;


	double drug_Y_pulse_period = parameters.doubles("drug_Y_pulse_period");
	double drug_Y_pulse_duration = parameters.doubles("drug_Y_pulse_duration");
	double drug_Y_pulse_concentration = parameters.doubles("drug_Y_pulse_concentration");

	double drug_Y_pulse_timer = drug_Y_pulse_period;
	double drug_Y_pulse_injection_timer = -1;

	// Drugs density index
	static int drug_X_ix = microenvironment.find_density_index("drug_X");	
	static int drug_Y_ix = microenvironment.find_density_index("drug_Y");

	// this is to emulate PhysiBoSSv1 TNF experiment
	// bool seed_tnf = true; // Was set to false in TNF model (?)

	// do small diffusion steps alone to initialize densities
	// if ( seed_tnf )
	// {
	// 	inject_density_sphere(drug_X_ix, drug_X_pulse_concentration, membrane_lenght);
	// 	inject_density_sphere(drug_Y_ix, drug_Y_pulse_concentration, membrane_lenght);

	// 	for ( int i = 0; i < 25; i ++ )
	// 		microenvironment.simulate_diffusion_decay( diffusion_dt );
	// }
	
	/* PhysiCell setup */ 
 	
	// set mechanics voxel size, and match the data structure to BioFVM
	double mechanics_voxel_size = 30; 
	Cell_Container* cell_container = create_cell_container_for_microenvironment( microenvironment, mechanics_voxel_size );
	
	/* Users typically ssetup_microenvironmenttart modifying here. START USERMODS */ 


	
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


	create_cell_types();
	std::cout << "create cells OK" << std::endl;



	// update monitoring variables too in setup_tissue();
	setup_tissue();
	std::cout << "setup tissue OK" << std::endl;




	// main loop
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
	save_PhysiCell_to_MultiCellDS_xml_pugi( filename , microenvironment , PhysiCell_globals.current_time ); 
	
	// save a quick SVG cross section through z = 0, after setting its 
	// length bar to 200 microns 

	PhysiCell_SVG_options.length_bar = 200; 

	// for simplicity, set a pathology coloring function 
	
	std::vector<std::string> (*cell_coloring_function)(Cell*) = my_coloring_function;
	// std::vector<std::string> (*cell_coloring_function)(Cell*) = heterogeneity_coloring_function;
	
	sprintf( filename , "%s/initial.svg" , PhysiCell_settings.folder.c_str() ); 
	SVG_plot( filename , microenvironment, 0.0 , PhysiCell_globals.current_time, cell_coloring_function );
	
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
					
					//Count Necrotic Apoptotic Alive cells
					// Producer
					std::string message;
					std::string topic_name = "cells";
					double timepoint = PhysiCell_globals.current_time;
					int alive_no,necrotic_no,apoptotic_no=0;
					alive_no = total_live_cell_count();
					necrotic_no = total_necrosis_cell_count();
					apoptotic_no = total_dead_cell_count();
					pid_t pid_var = getpid();
					// MessageBuilder builder(topic_name);
					message = std::to_string(pid_var) + ';' + std::to_string(timepoint) + ';' + std::to_string(alive_no) + ';' + std::to_string(apoptotic_no) + ';' + std::to_string(necrotic_no) + ';';
					
					
					// message = '{' + 'process_id' + std::to_string(pid_var) + ',' + 'timepoint' + std::to_string(timepoint) + ',' 
					// + 'alive' + std::to_string(alive_no) + ',' + 'apoptotic' + std::to_string(apoptotic_no) + ',' 
					// + 'necrotic' + std::to_string(necrotic_no) + '}';
					// Define the configuration structure
					// Configuration config = { { "metadata.broker.list", "localhost:9092" } };
				    // Create the producer
				    // BufferedProducer<std::string> producer(config);
					//Produce a message
					// The message that will be sent
					// std::cout << "Message to Kafka: " << message << std::endl;
					// std::string str(message);
					// builder.partition(0).payload(str);
				    // producer.add_message(builder);
				    // producer.flush();
				}
				
				if( PhysiCell_settings.enable_full_saves == true )
				{	
					sprintf( filename , "%s/output%08u" , PhysiCell_settings.folder.c_str(),  PhysiCell_globals.full_output_index ); 
					
					save_PhysiCell_to_MultiCellDS_xml_pugi( filename , microenvironment , PhysiCell_globals.current_time ); 
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
					SVG_plot( filename , microenvironment, 0.0 , PhysiCell_globals.current_time, cell_coloring_function );
					
					PhysiCell_globals.SVG_output_index++; 
					PhysiCell_globals.next_SVG_save_time  += PhysiCell_settings.SVG_save_interval;
				}
			}

			/*
			  Custom add-ons could potentially go here. 
			*/	


			/*
				Slight modification on the use of the inject_density and remove_density() functions in TNF
					- Just a single addition of substrate
					- remove_desntiy used to avoid any substrate prior to time of addition
					- Controlled by two parameters: time of addition and duration of pulse

			*/

			if ( PhysiCell_globals.current_time >= drug_X_pulse_timer && PhysiCell_globals.current_time <= drug_X_pulse_timer + drug_X_pulse_duration  )
				inject_density_sphere(drug_X_ix, drug_X_pulse_concentration, membrane_lenght);
				
			if (PhysiCell_globals.current_time < drug_X_pulse_timer)
				remove_density(drug_X_ix);
			
			if ( PhysiCell_globals.current_time >= drug_Y_pulse_timer && PhysiCell_globals.current_time <= drug_Y_pulse_timer + drug_Y_pulse_duration  )
				inject_density_sphere(drug_Y_ix, drug_Y_pulse_concentration, membrane_lenght);

			if (PhysiCell_globals.current_time < drug_Y_pulse_timer)
				remove_density(drug_Y_ix);

			
			// update the microenvironment
			microenvironment.simulate_diffusion_decay( diffusion_dt );

			
			// update te AGS receptor model of each cell 
			drug_transport_model_main( diffusion_dt );

			
			// run PhysiCell 
			((Cell_Container *)microenvironment.agent_container)->update_all_cells( PhysiCell_globals.current_time );
			
			
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
	save_PhysiCell_to_MultiCellDS_xml_pugi( filename , microenvironment , PhysiCell_globals.current_time ); 
	
	sprintf( filename , "%s/final.svg" , PhysiCell_settings.folder.c_str() ); 
	SVG_plot( filename , microenvironment, 0.0 , PhysiCell_globals.current_time, cell_coloring_function );

	
	// timer 
	
	std::cout << std::endl << "Total simulation runtime: " << std::endl; 
	BioFVM::display_stopwatch_value( std::cout , BioFVM::runtime_stopwatch_value() ); 

	return 0; 
}