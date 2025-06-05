/*
 * tnf_boolean_model_interface.cpp
 *
 *  Created on: 15 jun. 2020
 *  Author: Miguel Ponce-de-Leon (miguel.ponce@bsc.es)
 *  Contributor: Gerard Pradas
 *  Contributor: Arnau Montagud
 *  Contributor: Thalia Diniaco
 *  Cite as: arXiv:2103.14132 [q-bio.QM]
 *  Description: 
 *      Submodel that work as an interface 
 *      between the Boolean Network (BN) and PhysiCell (PC). 
 *      The submodel run the following steps:
 *      1- updates BN input nodes based on custom cell variables (see receptor model)
 *      2- updates the BN intracellular model by running MaBoSS
 *      3- updates cell state/behaviour based on the state of the BN readout nodes
 *  
 *      The update_monitor_variables funtion is just used to keep track of some relevand
 *      BN nodes' values that are stored as custom variables
 */


#include "./tnf_boolean_model_interface.h"

using namespace PhysiCell; 


double Hill_response_function_2( double s, double half_max , double hill_power )
{ 
    // newer. only one expensive a^b operation. 45% less computationl expense. 

	// give an early exit possibility to cut cost on "empty" rules
	if(s < 1e-16 ) // maybe also try a dynamic threshold: 0.0001 * half_max 
	{ return 0.0; } 

	// operations to reduce a^b operations and minimize hidden memory allocation / deallocation / copy operations. 
	// Hill = (s/half_max)^hill_power / ( 1 + (s/half_max)^hill_power  )
	double temp = s; // s 
	temp /= half_max; // s/half_max 
	double temp1 = pow(temp,hill_power); // (s/half_max)^h 
	temp = temp1;  // (s/half_max)^h 
	temp +=1 ;  // (1+(s/half_max)^h ); 
	temp1 /= temp; // (s/half_max)^h / ( 1 + s/half_max)^h) 
    if(temp1 < 1e-16)
        temp1 = 0.0;
	return temp1; 
}



double calculate_diffusion_flux(Cell* pCell, int density_idx, double permeability,std::string drug_name){
	
	// Cell and voxel geometries
	static double voxel_volume = microenvironment.mesh.dV;

	double cell_volume = pCell->phenotype.volume.total;
	double cell_radius = pCell->phenotype.geometry.radius;
	double cell_surface = 4 * M_PI * std::pow(cell_radius, 2);

	double density_ext = pCell->nearest_density_vector()[density_idx];	 // A density (mM)
	
	if (density_ext < 0){
	
		pCell->nearest_density_vector()[density_idx] = 0;

	}
	
	double density_int = pCell->phenotype.molecular.internalized_total_substrates[density_idx];
    //std::cout << density_int << "antes" << std::endl;
	density_int /= cell_volume; // divide int tot substrate to get concentration
    //std::cout << density_int << "despues" << std::endl;

	double flux = 0.0; // açò igual no cal, no ?? 
	flux = permeability * (density_int - density_ext) * cell_surface;
    // std::cout << flux << "flux" << std::endl;


	if( flux < 0.0 && std::abs(flux) >= density_ext * voxel_volume ){
		flux = (density_int * cell_volume) - (density_ext * voxel_volume); // amol
	} else if ( flux > 0.0 && flux >= density_int * cell_volume){
		flux = (density_int * cell_volume) - (density_ext * voxel_volume); // amol
	}

	// Then map to the custom data for post-simulation analysis
	std::string substrate_external_idx = drug_name + "_external_density";
	std::string substrate_internal_idx = drug_name + "_internal_density";

	int external_density_custom_data_idx = pCell->custom_data.find_variable_index(substrate_external_idx);
	int internal_density_custom_data_idx = pCell->custom_data.find_variable_index(substrate_internal_idx);

	pCell->custom_data.variables[external_density_custom_data_idx].value = density_ext;
	pCell->custom_data.variables[internal_density_custom_data_idx].value = density_int;

	//std::cout << drug_name << substrate_external_idx << " : " << density_ext << std::endl; 
	//std::cout << drug_name << substrate_internal_idx << " : " << density_int << std::endl; 

	return flux;
}

void drug_transport_model_update( Cell* pCell, Phenotype& phenotype, double dt )
{
	if( phenotype.death.dead == true )
	{ return; }

	// Fetch microenvironment variables
	static int drug_X_idx = microenvironment.find_density_index("drug_A");
	static int drug_Y_idx = microenvironment.find_density_index("oxygen");

	double drug_X_permeability = parameters.doubles("drug_X_permeability");
	//double drug_Y_permeability = parameters.doubles("drug_Y_permeability");

	double drug_X_flux = calculate_diffusion_flux(pCell, drug_X_idx, drug_X_permeability,"drug_A");
	double drug_Y_flux = calculate_diffusion_flux(pCell, drug_Y_idx, drug_X_permeability, "oxygen");

	pCell->phenotype.secretion.net_export_rates[drug_X_idx] = drug_X_flux;
	pCell->phenotype.secretion.net_export_rates[drug_Y_idx] = drug_Y_flux;
	

	return;
}


void update_boolean_model_inputs( Cell* pCell, Phenotype& phenotype, double dt )
{

    if( pCell->phenotype.death.dead == true )
	    { return; }
	
	
	//std::cout << "eweeeeeeeeeeeeeeeeeeeeeeeeeeeeeee" << std::endl;
    static int oxygen_substrate_index = pCell->get_microenvironment()->find_density_index( "oxygen" );
    static int drug_substrate_index = pCell->get_microenvironment()->find_density_index( "drug_A" ); 


    double pO2 = (pCell->nearest_density_vector())[oxygen_substrate_index]; // PhysiCell_constants::oxygen_index]; 
    double conc_drug = (pCell->nearest_density_vector())[drug_substrate_index];


    if ( pO2 <= pCell->parameters.o2_hypoxic_threshold ) {
 
       pCell->phenotype.intracellular-> set_boolean_variable_value("Hypoxia", 1);
       //std::cout << "poco o2" << std::endl;        
    }
    else {
    
        pCell->phenotype.intracellular->set_boolean_variable_value("Hypoxia", 0);
        //std::cout << "mucho o2" << std::endl;            
    }
	
	anti_node_mapping_function(pCell, "drug_A", parameters.strings("drug_X_target"), parameters.doubles("drug_X_half_max"), parameters.doubles("drug_X_Hill_coeff"));

    return;
}

void anti_node_mapping_function( Cell* pCell, std::string drug_name, std::string target_node, double drug_half_max, double drug_Hill_coeff){

    double cell_volume = pCell->phenotype.volume.total;
	static int drug_idx = microenvironment.find_density_index( drug_name );

    double drug_int = pCell->phenotype.molecular.internalized_total_substrates[drug_idx];
    drug_int /= cell_volume; // Convert to density (mM)

    double target_inactivate_p = Hill_response_function_2(drug_int, drug_half_max, drug_Hill_coeff );

    //std::cout << "drug_int: " << drug_int << std::endl;
    //std::cout << "aver  num:   " << target_inactivate_p  << std::endl;
    //std::cout << "For drug " << drug_name << " uniform_random value " << uniform_random() << " target inactivation prob " << target_inactivate_p << std::endl;

    if (target_node != "none"){
        if ( uniform_random() <  target_inactivate_p ){ // Added normalization by maximum rand value
            pCell->phenotype.intracellular->set_boolean_variable_value(target_node, 1);
            bool nodo = pCell->phenotype.intracellular->get_boolean_variable_value("TNFalpha");
            //std::cout << nodo << "jj" <<std::endl;
            //std::cout << "1" << std::endl;

        } else { 
            pCell->phenotype.intracellular->set_boolean_variable_value(target_node, 0);
            //std::cout << "iwi" << std::endl;
        }
    }

    return;
}

void update_cell_from_boolean_model(Cell* pCell, Phenotype& phenotype, double dt)
{	

    static int death_decay_idx = pCell->custom_data.find_variable_index( "death_commitment_decay" );
    static int apoptosis_index = phenotype.death.find_death_model_index( PhysiCell_constants::apoptosis_death_model );
    static float apoptosis_rate = pCell->custom_data["apoptosis_rate"];
    static float death_commitment_decay = pCell->custom_data["death_decay_idx"];
    
    // // Getting the state of the Boolean model readouts (Readout can be in the XML)
    bool apoptosis = pCell->phenotype.intracellular->get_boolean_variable_value( "Apoptosis" );
    // // bool nonACD = pCell->phenotype.intracellular->get_boolean_variable_value( "NonACD" );
    bool survival = pCell->phenotype.intracellular->get_boolean_variable_value( "Proliferation" );

    if ( survival ) 
    {   
        pCell->phenotype.cycle.advance_cycle(pCell, phenotype, dt); 
        //std::cout << "survival" << std::endl; 
    }

    if ( apoptosis ) 
    { 
        pCell-> phenotype.death.rates[apoptosis_index] =  apoptosis_rate * 50;
    }


    return;
}


void update_monitor_variables(Cell* pCell )
{

	static int index_hypoxia_node  = pCell->custom_data.find_variable_index("Hypoxia_node");
	static int index_prolif_node  = pCell->custom_data.find_variable_index("Proliferation_node");
    static int index_cell_pressure = pCell->custom_data.find_variable_index("cell_pressure");
    static int index_raf_node = pCell->custom_data.find_variable_index("TNFalpha_node");
   
    pCell->custom_data[index_hypoxia_node] = pCell->phenotype.intracellular->get_boolean_variable_value("Hypoxia");
    pCell->custom_data[index_prolif_node] = pCell->phenotype.intracellular->get_boolean_variable_value("Proliferation");
    pCell->custom_data[index_raf_node] = pCell->phenotype.intracellular->get_boolean_variable_value("TNFalpha");

    pCell->custom_data[index_cell_pressure] = pCell->state.simple_pressure; 

    return;
}

void update_cell_growth_parameters_pressure_based( Cell* pCell, Phenotype& phenotype, double dt ) 
{
	// supported cycle models:
		// advanced_Ki67_cycle_model= 0;
		// basic_Ki67_cycle_model=1
		// live_cells_cycle_model = 5; 
	
	if( phenotype.death.dead == true )
	{ return; }
	
	// set up shortcuts to find the Q and K(1) phases (assuming Ki67 basic or advanced model)
	static bool indices_initiated = false; 
	static int start_phase_index; // Q_phase_index; 
	static int end_phase_index; // K_phase_index;
	static int necrosis_index; 
	
	static int oxygen_substrate_index = pCell->get_microenvironment()->find_density_index( "oxygen" ); 
	
	if( indices_initiated == false )
	{
		// Ki67 models
		
		if( phenotype.cycle.model().code == PhysiCell_constants::advanced_Ki67_cycle_model || 
			phenotype.cycle.model().code == PhysiCell_constants::basic_Ki67_cycle_model )
		{
			start_phase_index = phenotype.cycle.model().find_phase_index( PhysiCell_constants::Ki67_negative );
			necrosis_index = phenotype.death.find_death_model_index( PhysiCell_constants::necrosis_death_model ); 
			
			if( phenotype.cycle.model().code == PhysiCell_constants::basic_Ki67_cycle_model )
			{
				end_phase_index = 
					phenotype.cycle.model().find_phase_index( PhysiCell_constants::Ki67_positive );
				indices_initiated = true; 
			}
			if( phenotype.cycle.model().code == PhysiCell_constants::advanced_Ki67_cycle_model )
			{
				end_phase_index = 
					phenotype.cycle.model().find_phase_index( PhysiCell_constants::Ki67_positive_premitotic );
				indices_initiated = true; 
			}
		}
		
		// live model 
			
		if( phenotype.cycle.model().code == PhysiCell_constants::live_cells_cycle_model )
		{
			start_phase_index = phenotype.cycle.model().find_phase_index( PhysiCell_constants::live );
			necrosis_index = phenotype.death.find_death_model_index( PhysiCell_constants::necrosis_death_model ); 
			end_phase_index = phenotype.cycle.model().find_phase_index( PhysiCell_constants::live );
			indices_initiated = true; 
		}
		
		// cytometry models 
		
		if( phenotype.cycle.model().code == PhysiCell_constants::flow_cytometry_cycle_model || 
			phenotype.cycle.model().code == PhysiCell_constants::flow_cytometry_separated_cycle_model )
		{
			start_phase_index = phenotype.cycle.model().find_phase_index( PhysiCell_constants::G0G1_phase );
			necrosis_index = phenotype.death.find_death_model_index( PhysiCell_constants::necrosis_death_model ); 
			end_phase_index = phenotype.cycle.model().find_phase_index( PhysiCell_constants::S_phase );
			indices_initiated = true; 
		}	

		if( phenotype.cycle.model().code == PhysiCell_constants::cycling_quiescent_model )
		{
			start_phase_index = phenotype.cycle.model().find_phase_index( PhysiCell_constants::quiescent );
			necrosis_index = phenotype.death.find_death_model_index( PhysiCell_constants::necrosis_death_model ); 
			end_phase_index = phenotype.cycle.model().find_phase_index( PhysiCell_constants::cycling );
			indices_initiated = true; 
		}
		
	}
	
	// don't continue if we never "figured out" the current cycle model. 
	if( indices_initiated == false )
	{
		return; 
	}

	// this multiplier is for linear interpolation of the oxygen value 
	double multiplier = 1.0;
	
	// now, update the appropriate cycle transition rate 

	// Check relative pressure to either number of neighbor cells or set max logistic function to number of neighbor cells
	// pressure threshold set to 1, above this value there is no growth

	double p = pCell->state.simple_pressure; 
    double hill_coeff_pressure = parameters.doubles("hill_coeff_pressure");
    double pressure_half = parameters.doubles("pressure_half");
    double scaling = pressure_effect_growth_rate(p, hill_coeff_pressure, pressure_half );
	// std::cout << "scaling is: " << scaling << std::endl;

	double rate = phenotype.cycle.data.transition_rate(0, 0);
	rate *= (1 - scaling);
	if (rate < 0)
		rate = 0;
	
	phenotype.cycle.data.transition_rate(start_phase_index, end_phase_index) = rate;
}

double pressure_effect_growth_rate(double pressure, double hill_coeff, double pressure_half){

    // Suggestion: Employ the Hill_effect function from PhysiCell instead of this one, just to align with the other mapping functions

    // double pressure_exponential_function = std::pow(6e-03, pressure);
    double pressure_exponential_function =  std::pow(pressure, hill_coeff) / (pressure_half + std::pow(pressure, hill_coeff));
    // if (pressure_exponential_function > 1) pressure_exponential_function = 1.0;
    return pressure_exponential_function;
}

void pre_update_intracellular(Cell* pCell, Phenotype& phenotype, double dt)
{
    if( phenotype.death.dead == true )
	{
		pCell->functions.update_phenotype = NULL;
		return;
	}
    // Update MaBoSS input nodes based on the environment and cell state
    update_boolean_model_inputs(pCell, phenotype, dt);
    
    return;
}


void post_update_intracellular(Cell* pCell, Phenotype& phenotype, double dt)
{
    if( phenotype.death.dead == true )
	{
		pCell->functions.update_phenotype = NULL;
		return;
	}
    
    // update the cell fate based on the boolean outputs
    update_cell_from_boolean_model(pCell, phenotype, dt);
    
    update_cell_and_death_parameters_O2_based(pCell, phenotype, dt);
    update_cell_growth_parameters_pressure_based(pCell, phenotype, dt);


    // Get track of some boolean node values for debugging
    // @oth: Probably not needed anymore with pcdl
    update_monitor_variables(pCell);

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

		int cell_microenv_index = pCell->get_current_voxel_index();
		bool DC_node = microenvironment.is_dirichlet_node(cell_microenv_index);
		bool agent_out_domain = pCell->is_out_of_domain;

		if (cell_microenv_index >= 0 && DC_node == false && agent_out_domain == false) // Avoid segfault of NER
			drug_transport_model_update( pCell, pCell->phenotype , dt );


		// Adding time to custom data
		static int time_index = pCell->custom_data.find_variable_index("time");
		pCell->custom_data.variables[time_index].value = PhysiCell_globals.current_time;

	}

	return;
}
