#include "./boolean_model_interface.h"
#include <math.h>

using namespace PhysiCell; 

Submodel_Information bm_interface_info;

// Auxiliar functions

std::string get_drug_target(std::string drug_name){
    std::string param_name = drug_name + "_target";
    std::string drug_target = parameters.strings(param_name);
    return drug_target;
}

void boolean_model_interface_setup()
{
    bm_interface_info.name = "AGS Boolean model interface"; 
	bm_interface_info.version = "0.0.2";
	
    bm_interface_info.main_function = ags_bm_interface_main; 

	// These are just auxiliary variables to keep track of some BN nodes

    bm_interface_info.cell_variables.push_back( "mek_node" );
    bm_interface_info.cell_variables.push_back( "pi3k_node" );
    bm_interface_info.cell_variables.push_back( "tak1_node" );
    bm_interface_info.cell_variables.push_back( "akt_node" );

    bm_interface_info.cell_variables.push_back( "anti_mek_node" );
    bm_interface_info.cell_variables.push_back( "anti_pi3k_node" );
    bm_interface_info.cell_variables.push_back( "anti_tak1_node" );
    bm_interface_info.cell_variables.push_back( "anti_akt_node" );


    bm_interface_info.cell_variables.push_back( "prosurvival_b1_node" );
    bm_interface_info.cell_variables.push_back( "prosurvival_b2_node" );
    bm_interface_info.cell_variables.push_back( "prosurvival_b3_node" );

    bm_interface_info.cell_variables.push_back( "antisurvival_b1_node" );
    bm_interface_info.cell_variables.push_back( "antisurvival_b2_node" );
    bm_interface_info.cell_variables.push_back( "antisurvival_b3_node" );


    // Could add here output of transfer functions
	bm_interface_info.register_model();
}


// @othmane: New PhysiBoSS, this might not be even needed, can be tracked through the BM states CSV in the output folder

void update_monitor_variables(Cell* pCell ) 
{
	static int mek_node_ix = pCell->custom_data.find_variable_index("mek_node");
	static int akt_node_ix = pCell->custom_data.find_variable_index("akt_node");
	static int pi3k_node_ix = pCell->custom_data.find_variable_index("pi3k_node");
	static int tak1_node_ix = pCell->custom_data.find_variable_index("tak1_node");

    static int anti_mek_node_ix = pCell->custom_data.find_variable_index("anti_mek_node");
    static int anti_akt_node_ix = pCell->custom_data.find_variable_index("anti_akt_node");
    static int anti_pi3k_node_ix = pCell->custom_data.find_variable_index("anti_pi3k_node");
    static int anti_tak1_node_ix = pCell->custom_data.find_variable_index("anti_tak1_node");

    static int antisurvival_b1_ix = pCell->custom_data.find_variable_index("antisurvival_b1_node");
    static int antisurvival_b2_ix = pCell->custom_data.find_variable_index("antisurvival_b2_node");
    static int antisurvival_b3_ix = pCell->custom_data.find_variable_index("antisurvival_b3_node");

    static int prosurvival_b1_ix = pCell->custom_data.find_variable_index("prosurvival_b1_node");
    static int prosurvival_b2_ix = pCell->custom_data.find_variable_index("prosurvival_b2_node");
    static int prosurvival_b3_ix = pCell->custom_data.find_variable_index("prosurvival_b3_node");

	pCell->custom_data[mek_node_ix] = pCell->phenotype.intracellular->get_boolean_variable_value( "MEK" );
    pCell->custom_data[akt_node_ix] = pCell->phenotype.intracellular->get_boolean_variable_value( "AKT" );
    pCell->custom_data[pi3k_node_ix] = pCell->phenotype.intracellular->get_boolean_variable_value( "PI3K" );
    pCell->custom_data[tak1_node_ix] = pCell->phenotype.intracellular->get_boolean_variable_value( "TAK1" );

    pCell->custom_data[anti_mek_node_ix] = pCell->phenotype.intracellular->get_boolean_variable_value( "anti_MEK" );
    pCell->custom_data[anti_akt_node_ix] = pCell->phenotype.intracellular->get_boolean_variable_value( "anti_AKT" );
    pCell->custom_data[anti_pi3k_node_ix] = pCell->phenotype.intracellular->get_boolean_variable_value( "anti_PI3K" );
    pCell->custom_data[anti_tak1_node_ix] = pCell->phenotype.intracellular->get_boolean_variable_value( "anti_TAK1" );

    pCell->custom_data[antisurvival_b1_ix] = pCell->phenotype.intracellular->get_boolean_variable_value( "Antisurvival_b1" );
    pCell->custom_data[antisurvival_b2_ix] = pCell->phenotype.intracellular->get_boolean_variable_value( "Antisurvival_b2" );
    pCell->custom_data[antisurvival_b3_ix] = pCell->phenotype.intracellular->get_boolean_variable_value( "Antisurvival_b3" );

    pCell->custom_data[prosurvival_b1_ix] = pCell->phenotype.intracellular->get_boolean_variable_value( "Prosurvival_b1" );
    pCell->custom_data[prosurvival_b2_ix] = pCell->phenotype.intracellular->get_boolean_variable_value( "Prosurvival_b2" );
    pCell->custom_data[prosurvival_b3_ix] = pCell->phenotype.intracellular->get_boolean_variable_value( "Prosurvival_b3" );

    return;
}

// Functions used to update the Boolean model just before running MaBoSS

double calculate_drug_effect(Cell* pCell, std::string drug_name){
    
	std::string p_half_max_name   = drug_name + "_half_max";
    std::string p_hill_coeff_name = drug_name + "_Hill_coeff";
    
    static int drug_idx         = microenvironment.find_density_index( drug_name );
    static int p_half_max_idx   = pCell->custom_data.find_variable_index(p_half_max_name);
    static int p_hill_coeff_idx = pCell->custom_data.find_variable_index(p_hill_coeff_name);
	
    double cell_volume   = pCell->phenotype.volume.total;
    double ic_drug_total = pCell->phenotype.molecular.internalized_total_substrates[drug_idx];
    double ic_drug_conc  = ic_drug_total / cell_volume; // Convert to concentration

    std::cout << pCell->custom_data[p_half_max_idx] << std::endl;

    double p_half_max    = pCell->custom_data[p_half_max_idx];
    double p_hill_coeff  = pCell->custom_data[p_hill_coeff_idx];

    return Hill_response_function(ic_drug_conc, p_half_max, p_hill_coeff);
}

void update_boolean_model_inputs( Cell* pCell, Phenotype& phenotype, double dt )
{
    if( pCell->phenotype.death.dead == true )
	{ return; }

    int n_drugs = 2;
    std::string drugs[n_drugs] = { "drug_X", "drug_Y" };
  
    for (int i = 0; i < n_drugs; i++){
        std::string drug_name = drugs[i];
        std::string target_node = get_drug_target(drug_name);
        double drug_effect = calculate_drug_effect(pCell, drug_name);    
        if (drug_effect > 0){ // why 0 here?
            if ( uniform_random() < drug_effect )
                pCell->phenotype.intracellular->set_boolean_variable_value(target_node, 1);
            else
                pCell->phenotype.intracellular->set_boolean_variable_value(target_node, 0);
        }
    }
    return;
}

void pathway_reactivation( Cell* pCell, Phenotype& phenotype, double dt ){

    static int reactivation_prob_idx = pCell->custom_data.find_variable_index("reactivation_value");
    double p_reactivation = pCell->custom_data.variables[reactivation_prob_idx].value;
    
    int n_drugs = 2;
    std::string drugs[n_drugs] = { "drug_X", "drug_Y" };
    for (int i = 0; i < n_drugs; i++){
        std::string drug_name = drugs[i];
        std::string target_node = get_drug_target(drug_name);
        if (target_node == "none" )
            continue;
        if ( uniform_random() < p_reactivation )
            pCell->phenotype.intracellular->set_boolean_variable_value(target_node, 1);
    }
    return;
}

void pre_update_intracellular_ags(Cell* pCell, Phenotype& phenotype, double dt)
{
    if( phenotype.death.dead == true )
	{
		pCell->functions.update_phenotype = NULL;
		return;
	}
    // Update MaBoSS input nodes based on the environment and cell state
    update_boolean_model_inputs(pCell, phenotype, dt);

    // This function can be use for the reactivation mechanisms
    // pathway_reactivation( Cell* pCell, Phenotype& phenotype, double dt )
    
    return;
}

// Functions used to update the cell agent just after MaBoSS was updated

double growth_mapping_logistic(double doubling_time, double hill_coeff, double K_half, double S_value){

    // double growth_logistic_function = doubling_time / (1 + std::exp(- log10(readout_value) * scaling));
    double growth_logistic_function;
    growth_logistic_function = (doubling_time * std::pow(S_value, hill_coeff ) ) / (K_half + std::pow(S_value, hill_coeff) ) ;

    return growth_logistic_function;
}

double apoptosis_mapping_logistic(double basal_apoptosis_rate, double maximum_apoptosis_rate, double hill_coeff, double K_half, double S_value){
    double apoptosis_mapping_function;
    apoptosis_mapping_function = (maximum_apoptosis_rate * std::pow(S_value, hill_coeff) ) / (K_half + std::pow(S_value, hill_coeff) ) ;

    return apoptosis_mapping_function +  basal_apoptosis_rate;
}

void update_cell_from_boolean_model(Cell* pCell, Phenotype& phenotype, double dt){

    if( pCell->phenotype.death.dead == true )
	{ return; }

    static int apoptosis_model_index = phenotype.death.find_death_model_index( "Apoptosis" );
    static int necrosis_model_index = phenotype.death.find_death_model_index( "Necrosis" );

    bool casp37_b1 = pCell->phenotype.intracellular->get_boolean_variable_value( "Caspase37_b1" );
    bool casp37_b2 = pCell->phenotype.intracellular->get_boolean_variable_value( "Caspase37_b2" );
    bool FOXO = pCell->phenotype.intracellular->get_boolean_variable_value( "FOXO" );

    bool antisurvival_b1 = pCell->phenotype.intracellular->get_boolean_variable_value( "Antisurvival_b1" );
    bool antisurvival_b2 = pCell->phenotype.intracellular->get_boolean_variable_value( "Antisurvival_b2" );
    bool antisurvival_b3 = pCell->phenotype.intracellular->get_boolean_variable_value( "Antisurvival_b3" );

    double anti_w1 = parameters.doubles("w1_apoptosis");
    double anti_w2 = parameters.doubles("w2_apoptosis") + anti_w1;
    double anti_w3 = 1 - (anti_w1 + anti_w2);
    double S_anti = (anti_w1*antisurvival_b1) + (anti_w2 * antisurvival_b2) + (anti_w3 * antisurvival_b3);
    double S_anti_real = (anti_w1*casp37_b1) + (anti_w2 * casp37_b2) + (anti_w3 * FOXO);

    bool prosurvival_b1 = pCell->phenotype.intracellular->get_boolean_variable_value( "Prosurvival_b1" );
    bool prosurvival_b2 = pCell->phenotype.intracellular->get_boolean_variable_value( "Prosurvival_b2" );
    bool prosurvival_b3 = pCell->phenotype.intracellular->get_boolean_variable_value( "Prosurvival_b3" );

    bool cMYC = pCell->phenotype.intracellular->get_boolean_variable_value( "cMYC" );
    bool CCND_b1 = pCell->phenotype.intracellular->get_boolean_variable_value( "CCND1_b1" );
    bool CCND_b2 = pCell->phenotype.intracellular->get_boolean_variable_value( "CCND1_b2" );

    double pro_w1 = parameters.doubles("w1_growth");
    double pro_w2 = parameters.doubles("w2_growth") + pro_w1;
    double pro_w3 = 1 - (pro_w1 + pro_w2);
    double S_pro = (pro_w1*prosurvival_b1) + (pro_w2 * prosurvival_b2) + (pro_w3 * prosurvival_b3);
    double S_pro_real = (pro_w1*cMYC) + (pro_w2 * CCND_b1) + (pro_w3 * CCND_b2);

    // Connect output from model to actual cell variables

    double apoptosis_rate_basal = parameters.doubles("apoptosis_rate_basal");
    double maximum_apoptosis_rate = parameters.doubles("max_apoptosis_rate");
    double hill_coeff_apoptosis = parameters.doubles("hill_coeff_apoptosis");
    double K_half_apoptosis = parameters.doubles("K_half_apoptosis");
    
    double apoptosis_value =  apoptosis_mapping_logistic(apoptosis_rate_basal, maximum_apoptosis_rate, 
                                                            hill_coeff_apoptosis, K_half_apoptosis, S_anti);
    
    double apoptosis_value_Hill = Hill_response_function(S_anti, K_half_apoptosis, 
                                                            hill_coeff_apoptosis);


    double basal_growth_rate = parameters.doubles("basal_growth_rate");
    double hill_coeff_growth = parameters.doubles("hill_coeff_growth");
    double K_half_growth = parameters.doubles("K_half_growth");
    double growth_value = growth_mapping_logistic(basal_growth_rate, hill_coeff_growth, K_half_growth, S_pro);
    double growth_value_Hill = Hill_response_function(S_pro, K_half_growth, hill_coeff_growth); // Max value is 1 

    // Another option is to use the growth_value_Hill as a probability

    if ( uniform_random() < growth_value_Hill ){ 
        phenotype.cycle.data.transition_rate(0, 0) = basal_growth_rate;
    } else { 
        phenotype.cycle.data.transition_rate(0, 0) = 0;  // Best to actually add an arrest function
    }

    pCell-> phenotype.death.rates[apoptosis_model_index] = (apoptosis_value_Hill * maximum_apoptosis_rate) + apoptosis_rate_basal;
    
    return;
}

void post_update_intracellular_ags(Cell* pCell, Phenotype& phenotype, double dt)


{
    if( phenotype.death.dead == true )
	{
		pCell->functions.update_phenotype = NULL;
		return;
	}
    
    // update the cell fate based on the boolean outputs
    update_cell_from_boolean_model(pCell, phenotype, dt);
    
    // Get track of some boolean node values for debugging
    // Probably not needed anymore
    update_monitor_variables(pCell);

    return;
}

//  @oth: added the main interface function

void ags_bm_interface_main (Cell* pCell, Phenotype& phenotype, double dt){
    
	if( phenotype.death.dead == true )
	{
		pCell->functions.update_phenotype = NULL;
		return;
	}
}