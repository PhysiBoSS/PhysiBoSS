#include "./drug_sensitivity.h"
#include <cmath>
#include <limits>

namespace PhysiCell {

const std::unordered_map<std::string, std::string> drug_targets = {
    {"Ipatasertib", "AKT"},
    {"Afuresertib", "AKT"},
    {"Afatinib", "EGFR"}
    // ... other mappings
};

const std::unordered_map<std::string, int> half_lives = {
    {"Ipatasertib", 2748},
    {"Afuresertib", 2448}
    // ... other mappings
};

DrugSensitivityParams get_drug_sensitivity_values(const std::string& drug_name) {
    try {
        return {
            parameters.doubles(drug_name + "_maxc"),
            parameters.doubles(drug_name + "_xmid"),
            parameters.doubles(drug_name + "_scal")
        };
    } catch (const std::exception& e) {
        throw std::runtime_error("Failed to get sensitivity values for drug: " + drug_name);
    }
}

void ConcentrationCalculator::validate_concentration(double conc) {
    if (conc < 0 || std::isnan(conc) || std::isinf(conc)) {
        throw std::invalid_argument("Invalid concentration value");
    }
}

void ConcentrationCalculator::validate_max_concentration(double max_conc) {
    if (max_conc <= 0 || std::isnan(max_conc) || std::isinf(max_conc)) {
        throw std::invalid_argument("Invalid maximum concentration value");
    }
}

double ConcentrationCalculator::get_x_from_conc(double x_conc, double max_conc) {
    validate_concentration(x_conc);
    validate_max_concentration(max_conc);
    
    if (x_conc > max_conc) {
        throw std::invalid_argument("Concentration cannot exceed maximum concentration");
    }
    
    return (std::log2(x_conc / max_conc)) + 9.0;
}

double get_cell_viability_for_drug_conc(
    double drug_conc, 
    const std::string& cell_line, 
    const std::string& drug_name) 
{
    ConcentrationCalculator::validate_concentration(drug_conc);
    
    auto params = get_drug_sensitivity_values(drug_name);
    double x = ConcentrationCalculator::get_x_from_conc(drug_conc, params.max_conc);
    
    // Logistic function with bounds checking
    double viability = 1.0 / (1.0 + std::exp((x - params.xmid) / params.scale));
    return std::max(0.0, std::min(1.0, viability));
}

} // namespace PhysiCell