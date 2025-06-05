#include "./boolean_model_interface.h"
#include <math.h>

using namespace PhysiCell; 

Submodel_Information bm_interface_info;



void pre_update_intracellular_ags(Cell* pCell, Phenotype& phenotype, double dt)
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
    // @oth: Probably not needed anymore with pcdl
    update_monitor_variables(pCell);

    return;
}





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

    // bm_interface_info.cell_variables.push_back( "prosurvival_b1_node" );
    // bm_interface_info.cell_variables.push_back( "prosurvival_b2_node" );
    // bm_interface_info.cell_variables.push_back( "prosurvival_b3_node" );

    // bm_interface_info.cell_variables.push_back( "antisurvival_b1_node" );
    // bm_interface_info.cell_variables.push_back( "antisurvival_b2_node" );
    // bm_interface_info.cell_variables.push_back( "antisurvival_b3_node" );


    // Could add here output of transfer functions
	bm_interface_info.register_model();
}



// @oth: This is not really needed, only used for the setup function above
void ags_bm_interface_main (Cell* pCell, Phenotype& phenotype, double dt){
    
	if( phenotype.death.dead == true )
	{
		pCell->functions.update_phenotype = NULL;
		return;
	}
}


// @oth: New PhysiBoSS, this might not be even needed, can be tracked through the BM states CSV in the output folder
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

    // static int antisurvival_b1_ix = pCell->custom_data.find_variable_index("antisurvival_b1_node");
    // static int antisurvival_b2_ix = pCell->custom_data.find_variable_index("antisurvival_b2_node");
    // static int antisurvival_b3_ix = pCell->custom_data.find_variable_index("antisurvival_b3_node");

    // static int prosurvival_b1_ix = pCell->custom_data.find_variable_index("prosurvival_b1_node");
    // static int prosurvival_b2_ix = pCell->custom_data.find_variable_index("prosurvival_b2_node");
    // static int prosurvival_b3_ix = pCell->custom_data.find_variable_index("prosurvival_b3_node");

	pCell->custom_data[mek_node_ix] = pCell->phenotype.intracellular->get_boolean_variable_value( "MEK" );
    pCell->custom_data[akt_node_ix] = pCell->phenotype.intracellular->get_boolean_variable_value( "AKT" );
    pCell->custom_data[pi3k_node_ix] = pCell->phenotype.intracellular->get_boolean_variable_value( "PI3K" );
    pCell->custom_data[tak1_node_ix] = pCell->phenotype.intracellular->get_boolean_variable_value( "TAK1" );

    pCell->custom_data[anti_mek_node_ix] = pCell->phenotype.intracellular->get_boolean_variable_value( "anti_MEK" );
    pCell->custom_data[anti_akt_node_ix] = pCell->phenotype.intracellular->get_boolean_variable_value( "anti_AKT" );
    pCell->custom_data[anti_pi3k_node_ix] = pCell->phenotype.intracellular->get_boolean_variable_value( "anti_PI3K" );
    pCell->custom_data[anti_tak1_node_ix] = pCell->phenotype.intracellular->get_boolean_variable_value( "anti_TAK1" );

    // pCell->custom_data[antisurvival_b1_ix] = pCell->phenotype.intracellular->get_boolean_variable_value( "Antisurvival_b1" );
    // pCell->custom_data[antisurvival_b2_ix] = pCell->phenotype.intracellular->get_boolean_variable_value( "Antisurvival_b2" );
    // pCell->custom_data[antisurvival_b3_ix] = pCell->phenotype.intracellular->get_boolean_variable_value( "Antisurvival_b3" );

    // pCell->custom_data[prosurvival_b1_ix] = pCell->phenotype.intracellular->get_boolean_variable_value( "Prosurvival_b1" );
    // pCell->custom_data[prosurvival_b2_ix] = pCell->phenotype.intracellular->get_boolean_variable_value( "Prosurvival_b2" );
    // pCell->custom_data[prosurvival_b3_ix] = pCell->phenotype.intracellular->get_boolean_variable_value( "Prosurvival_b3" );

    return;
}

// Functions used to update the Boolean model just before running MaBoSS

double calculate_drug_effect(Cell* pCell, std::string drug_name){

    /*

    This function calculates the effect of a drug on a specific node in the Boolean model.

    The input is the drug internalized concentration, the output is the Probability of turning OFF a specific node
    This Hill function is shaped XML parameters for the Hill index and the half-max

    */
    
	std::string p_half_max_name   = drug_name + "_half_max";
    std::string p_hill_coeff_name = drug_name + "_Hill_coeff";
    
    static int drug_idx         = microenvironment.find_density_index( drug_name );
    static int p_half_max_idx   = pCell->custom_data.find_variable_index(p_half_max_name);
    static int p_hill_coeff_idx = pCell->custom_data.find_variable_index(p_hill_coeff_name);
	
    double cell_volume   = pCell->phenotype.volume.total;
    double ic_drug_total = pCell->phenotype.molecular.internalized_total_substrates[drug_idx];
    double ic_drug_conc  = ic_drug_total / cell_volume; // Convert to concentration

    // std::cout << pCell->custom_data[p_half_max_idx] << std::endl;

    double p_half_max    = pCell->custom_data[p_half_max_idx];
    double p_hill_coeff  = pCell->custom_data[p_hill_coeff_idx];

    return Hill_response_function(ic_drug_conc, p_half_max, p_hill_coeff);
}

// @oth: Instead of this, we could have the input rules
void update_boolean_model_inputs( Cell* pCell, Phenotype& phenotype, double dt )
{
    /*
    This function updates the Boolean model inputs based on the drug concentrations.a
    After the drug effect is calculated, a Gillespie process is used to update the Boolean model.
    */

    if( pCell->phenotype.death.dead == true )
	{ return; }

    int n_drugs = 2;
    std::string drugs[n_drugs] = { "drug_X", "drug_Y" };
  
    for (int i = 0; i < n_drugs; i++){
        std::string drug_name = drugs[i];
        std::string target_node = get_drug_target(drug_name);

        // Single drug cases, drug_Y is not used so it's set to "null" in the config file
        if (target_node == "null") {
            return; 
        }

        double drug_effect = calculate_drug_effect(pCell, drug_name);

        if (drug_effect > 0){ // why 0 here?
            // Apply Gillespie only for the target_node obtained
            if (uniform_random() < drug_effect) // Gillespie
                pCell->phenotype.intracellular->set_boolean_variable_value(target_node, 1);
            else
                pCell->phenotype.intracellular->set_boolean_variable_value(target_node, 0);
        } 
    }
    return;
}



// @oth: added these functions to compute the readout nodes from the BM
double get_boolean_antisurvival_outputs(Cell* pCell, Phenotype& phenotype){

    /*
    This function computes the output of the apoptosis pathway in the Boolean model.
    It does a convex combination of the three output nodes.
    */

    // Antisurvival node outputs
    bool FOXO = pCell->phenotype.intracellular->get_boolean_variable_value( "FOXO" );
    bool casp8 = pCell->phenotype.intracellular->get_boolean_variable_value( "Caspase8" );
    bool casp9 = pCell->phenotype.intracellular->get_boolean_variable_value( "Caspase9" );
    double anti_w1 = get_custom_data_variable(pCell, "w1_apoptosis");
    double anti_w2 = get_custom_data_variable(pCell, "w2_apoptosis");
    double anti_w3 = get_custom_data_variable(pCell, "w3_apoptosis");
    double total_anti = anti_w1 + anti_w2 + anti_w3;
    double anti_w1_scaled = anti_w1 / total_anti;
    double anti_w2_scaled = anti_w2 / total_anti;
    double anti_w3_scaled = anti_w3 / total_anti;
    double S_anti_real = (anti_w1_scaled*casp8) + (anti_w2_scaled * casp9) + (anti_w3_scaled * FOXO);

    // For the Control case 
    if (total_anti == 0.0){
        return 0.0;
    } else {
        return S_anti_real;
    }
}

double get_boolean_prosurvival_outputs(Cell* pCell, Phenotype& phenotype){

    /*
    This function computes the output of the prosurvival pathway in the Boolean model.
    It does a convex combination of the three output nodes.
    */

    // bool CCND1 = pCell->phenotype.intracellular->get_boolean_variable_value( "CCND1" );
    bool cMYC = pCell->phenotype.intracellular->get_boolean_variable_value( "cMYC" );
    bool TCF = pCell->phenotype.intracellular->get_boolean_variable_value( "TCF" );
    bool RSK = pCell->phenotype.intracellular->get_boolean_variable_value( "RSK" );

    double pro_w1 = get_custom_data_variable(pCell, "w1_growth");
    double pro_w2 = get_custom_data_variable(pCell, "w2_growth");
    double pro_w3 = get_custom_data_variable(pCell, "w3_growth");
    double total_pro = pro_w1 + pro_w2 + pro_w3;
    double pro_w1_scaled = pro_w1 / total_pro;
    double pro_w2_scaled = pro_w2 / total_pro;
    double pro_w3_scaled = pro_w3 / total_pro;

    double S_pro_real = (pro_w1_scaled * cMYC) + (pro_w2_scaled * TCF) + (pro_w3_scaled * RSK);

    // For the Control curve, when all weights are 0, the output should be 1.0 (maximum growth rate)
    // This avoids the -nan value as an input for the Hill function
    if (total_pro == 0.0){
        return 1.0;
    } else {
        return S_pro_real;
    }
}


void update_cell_from_boolean_model(Cell* pCell, Phenotype& phenotype, double dt){

    /*
    This function updates the cell variables based on the Boolean model outputs.
    Maps the value of the convex combination of the output nodes to the cell variables.
    Prosurvival and Antisurvival are mapped to the growth rate and apoptosis rate respectively.

    */

    if( pCell->phenotype.death.dead == true )
	{ return; }

    // Effect on the apoptosis rate
    static int apoptosis_model_index = phenotype.death.find_death_model_index( "Apoptosis" );
    static int necrosis_model_index = phenotype.death.find_death_model_index( "Necrosis" );
    // Connect output from model to actual cell variables
    double apoptosis_rate_basal = get_custom_data_variable(pCell, "apoptosis_rate_basal");
    double maximum_apoptosis_rate =  get_custom_data_variable(pCell, "max_apoptosis_rate");
    double hill_coeff_apoptosis = get_custom_data_variable(pCell, "hill_coeff_apoptosis");
    double K_half_apoptosis = get_custom_data_variable(pCell, "K_half_apoptosis");
    double S_anti_real = get_boolean_antisurvival_outputs(pCell, phenotype);
    // sigmoidal mapping
    double apoptosis_value_Hill = maximum_apoptosis_rate * (Hill_response_function(S_anti_real, K_half_apoptosis, hill_coeff_apoptosis));
    apoptosis_value_Hill += apoptosis_rate_basal;
    pCell-> phenotype.death.rates[apoptosis_model_index] = apoptosis_value_Hill;

    // Effect on the growth rate
    double basal_growth_rate = get_custom_data_variable(pCell, "basal_growth_rate");
    double hill_coeff_growth = get_custom_data_variable(pCell, "hill_coeff_growth");
    double K_half_growth = get_custom_data_variable(pCell, "K_half_growth");
    double S_pro_real = get_boolean_prosurvival_outputs(pCell, phenotype);
    // sigmoidal mapping
    double growth_value_Hill =  (Hill_response_function(S_pro_real, K_half_growth, hill_coeff_growth));
    growth_value_Hill *= basal_growth_rate; // Max value is the basal growth rate
    double growth_value_Hill_arrested = Hill_response_function(S_pro_real, K_half_growth, hill_coeff_growth); // Max is 1, min is 0
    // Old approach: Just map directly the value of the Hill function to the growth rate
    pCell->phenotype.cycle.data.transition_rate(0, 0) = growth_value_Hill;

    return;
}


