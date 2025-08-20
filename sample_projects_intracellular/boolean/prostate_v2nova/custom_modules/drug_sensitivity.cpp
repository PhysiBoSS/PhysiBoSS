
#include "./drug_sensitivity.h"

/**
 *	\drug sensitivity
 *	\brief drug sensitivity module for prostate example
 * 
 *	\details Module needed for the prostate example. 
 *
 *
 *	\date 15/02/2020
 *	\author Annika Meert, BSC-CNS
 */

using namespace std;


const vector<pair<string, string>> drug_targets = {
    { "Ipatasertib", "AKT"},
    { "Afuresertib", "AKT"},
    { "Afatinib", "EGFR"},
    { "Erlotinib", "EGFR"},
    { "Ulixertinib", "ERK"},
    { "Luminespib", "HSPs"},
    { "Trametinib", "MEK1_2"},
    { "Selumetinib", "MEK1_2"},
    { "Pictilisib", "PI3K"},
    { "Alpelisib", "PI3K"},
    { "BIBR1532", "TERT"}
};

const vector<pair<string, int>> half_lives = {
    { "Ipatasertib", 2748},
    { "Afuresertib", 2448},
    { "Afatinib", 2220},
    { "Erlotinib", 2172},
    { "Ulixertinib", 105},
    { "Luminespib", 7200},
    { "Trametinib", 5760},
    { "Selumetinib", 822},
    { "Pictilisib", 1062},
    { "Alpelisib", 822},
    
};

string get_value (const vector<pair<string, string>> dict, string key) {
    vector< pair<string, string>>::const_iterator dict_iterator = find_if( dict.begin(), dict.end(),[&key](const pair < string, string>& element){ return element.first  == key;} );
    return (*dict_iterator).second;
}

int get_value (const vector<pair<string, int>> dict, string key) {
    vector< pair<string, int>>::const_iterator dict_iterator = find_if( dict.begin(), dict.end(),[&key](const pair < string, int>& element){ return element.first  == key;} );
    return (*dict_iterator).second;
}


vector<double> get_drug_sensitivity_values (string drug_name) {
    static double max_conc = parameters.doubles(drug_name+"_maxc");
    static double xmid = parameters.doubles(drug_name+"_xmid");
    static double scale = parameters.doubles(drug_name+"_scal");

    return {max_conc, xmid, scale};
}


double get_drug_concentration_from_IC (string cell_line, string drug_name, string IC_value, int simulation_mode) {
    std::stringstream ss;
    ss << IC_value[3] << "." << IC_value[2];
    double inhibition_value = stod(ss.str());
    std::cout << "Inhibition value for "<< drug_name << " is: " << inhibition_value << '\n';
    double final_viability = 1.0 - inhibition_value;
    
   // call functions to retrieve data from datastructure cell line and drug
    vector<double> drug_sens_vals = get_drug_sensitivity_values(drug_name);
    // call linear_mixed_model_function
    double max_conc = drug_sens_vals[0];
    double xmid = drug_sens_vals[1];
    double scale = drug_sens_vals[2];

    // get the drug concentration for the cell viability
    double x = get_x_for_cell_viability(xmid, scale, final_viability);
    double drug_conc = get_conc_from_x(x, max_conc);

    return drug_conc;
}

// returns x: the concentration scaled for 9 different concentrations
double get_x_from_conc(double x_conc, double max_conc) {
    double x = (log (x_conc / max_conc) / log (2) ) + 9;
    return(x);
}

// returns conc: concentration in micromolar 
double get_conc_from_x (double x, double max_conc) {
    double x_conc = max_conc * pow(2 , (x - 9));
    return(x_conc);
}

// returns the natural logarithm for x
double get_lx_from_x (double x, double max_conc) {
    double lx = log (get_conc_from_x(x, max_conc));
    return(lx);
} 

// two parameter logistic function used to model GDSC dose response data 
// this model is described further in published in Vis, D.J. et al Pharmacogenomics 2016, 17(7):691-700) 
// double logist3 (double x, double xmid, double scale) {
//     //https://github.com/CancerRxGene/gdscIC50/blob/master/R/nlme_fit.R

//     return;
// }


// this function then has to be called whenever we translate concentration into maboss
// check in the datastructure for cell line and drug for the important values and enter concentration
// then we receive a y_hat that we can use to calculate the maboss activity
double get_cell_viability_for_x(double x, double max_conc, double xmid, double scale) 
{
    // double y_hat = 1 - logist3(x, xmid, scale);
    double y_hat = 1/(1 + pow( exp(1), (x- xmid) /scale));
    return(y_hat);
}


double get_x_for_cell_viability (double xmid, double scale, double cell_viability) {
    double x = (log((1/cell_viability) - 1) * scale) + xmid;
    return x;
}

double get_cell_viability_for_drug_conc (double drug_conc, string cell_line, string drug_name) 
{
    //std::cout << "Concentration of " << drug_name << ": " << drug_conc << std::endl;
    // call functions to retrieve data from datastructure cell line and drug
    vector<double> drug_sens_vals = get_drug_sensitivity_values(drug_name);
    // call linear_mixed_model_function
    double max_conc = drug_sens_vals[0];
    double xmid = drug_sens_vals[1];
    double scale = drug_sens_vals[2];
    // transfor drug concentration into lx
    double x = get_x_from_conc(drug_conc, max_conc);
    double y_hat = get_cell_viability_for_x(x, max_conc, xmid, scale);
    // return cell viability value y_hat
    //std::cout << "Cell viability for " << drug_name << ": " << y_hat << std::endl;
    return y_hat;
}



