
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


/**
 *	\main main-template_BM
 *	\brief Main file of the template examples with Boolena models 
 * 
 *	\details Compiling this main file runs the template-BM sample_project. This file and its custom modules can be used as templates to generate other PhysiBoSS examples.
 *
 *	\date 19/10/2020
 *	\author Arnau Montagud, BSC-CNS, with code previously developed by Gerard Pradas and Miguel Ponce de Leon, BSC-CNS
 */

// put custom code modules here! 

#include "./custom_modules/custom.h" 
#include "./custom_modules/custom_main.h"


using namespace BioFVM;
using namespace PhysiCell;

int main( int argc, char* argv[] )
{
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

	// User parameters
	double time_add_tnf = parameters.ints("time_add_tnf");
	double time_put_tnf = 0;
	double duration_add_tnf = parameters.ints("duration_add_tnf");
	double time_tnf_next = 0;
	double time_remove_tnf = parameters.ints("time_remove_tnf");
	double concentration_tnf = parameters.doubles("concentration_tnf") * microenvironment.voxels(0).volume * 0.000001;
	double time_add_ainhib = parameters.ints("time_add_ainhib");
	double time_put_ainhib = 0;
	double duration_add_ainhib = parameters.ints("duration_add_ainhib");
	double time_ainhib_next = 0;
	double time_remove_ainhib = parameters.ints("time_remove_ainhib");
	double concentration_ainhib = parameters.doubles("concentration_ainhib") * microenvironment.voxels(0).volume * 0.000001;
	double membrane_lenght = parameters.ints("membrane_length");
	double time_add_binhib = parameters.ints("time_add_binhib");
	double time_put_binhib = 0;
	double duration_add_binhib = parameters.ints("duration_add_binhib");
	double time_binhib_next = 0;
	double time_remove_binhib = parameters.ints("time_remove_binhib");
	double concentration_binhib = parameters.doubles("concentration_binhib") * microenvironment.voxels(0).volume * 0.000001;

	// do small diffusion steps alone to initialize densities
	int k = microenvironment.find_density_index("tnf");
	if ( k >= 0 ) 
		inject_density_sphere(k, concentration_tnf, membrane_lenght);
	for ( int i = 0; i < 25; i ++ )
		microenvironment.simulate_diffusion_decay( 5*diffusion_dt );
	
	/* PhysiCell setup */ 
 	
	// set mechanics voxel size, and match the data structure to BioFVM
	double mechanics_voxel_size = 30; 
	Cell_Container* cell_container = create_cell_container_for_microenvironment( microenvironment, mechanics_voxel_size );
	
	/* Users typically start modifying here. START USERMODS */ 
	
	create_cell_types();
	
	setup_tissue();

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
	
	sprintf( filename , "%s/initial.svg" , PhysiCell_settings.folder.c_str() ); 
	SVG_plot( filename , microenvironment, 0.0 , PhysiCell_globals.current_time, cell_coloring_function );
	
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
				display_simulation_status( std::cout ); 
				if( PhysiCell_settings.enable_legacy_saves == true )
				{	
					log_output( PhysiCell_globals.current_time , PhysiCell_globals.full_output_index, microenvironment, report_file);
				}
				
				if( PhysiCell_settings.enable_full_saves == true )
				{	
					sprintf( filename , "%s/output%08u" , PhysiCell_settings.folder.c_str(),  PhysiCell_globals.full_output_index ); 
					
					save_PhysiCell_to_MultiCellDS_xml_pugi( filename , microenvironment , PhysiCell_globals.current_time ); 
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
				
			//Ainhib
			if ( PhysiCell_globals.current_time >= time_put_ainhib )
			{
				time_ainhib_next = PhysiCell_globals.current_time + duration_add_ainhib;
				time_put_ainhib += time_add_ainhib;
			}

			if ( PhysiCell_globals.current_time >= time_remove_ainhib )
			{
				int k = microenvironment.find_density_index("Ainhib");
				if ( k >= 0 )
					remove_density(k);
				time_remove_ainhib += PhysiCell_settings.max_time;
			}

			if ( PhysiCell_globals.current_time <= time_ainhib_next )
			{
				int k = microenvironment.find_density_index("Ainhib");
				if ( k >= 0 )
					inject_density_sphere(k, concentration_ainhib, membrane_lenght);
			}
			
			//Binhib
			if ( PhysiCell_globals.current_time >= time_put_binhib )
			{
				time_binhib_next = PhysiCell_globals.current_time + duration_add_binhib;
				time_put_binhib += time_put_binhib;
			}

			if ( PhysiCell_globals.current_time >= time_remove_binhib )
			{
				int k = microenvironment.find_density_index("Binhib");
				if ( k >= 0 )
					remove_density(k);
				time_remove_binhib += PhysiCell_settings.max_time;
			}

			if ( PhysiCell_globals.current_time <= time_binhib_next )
			{
				int k = microenvironment.find_density_index("Binhib");
				if ( k >= 0 )
					inject_density_sphere(k, concentration_binhib, membrane_lenght);
			}
			
			// update the microenvironment
			microenvironment.simulate_diffusion_decay( diffusion_dt );
			
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
