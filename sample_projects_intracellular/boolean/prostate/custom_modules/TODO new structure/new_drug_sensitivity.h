#ifndef __DRUG_SENSITIVITY_H__
#define __DRUG_SENSITIVITY_H__

#include <string>
#include <unordered_map>
#include <stdexcept>
#include "../core/PhysiCell.h"
#include "../modules/PhysiCell_standard_modules.h" 

namespace PhysiCell {

// Use unordered_map for O(1) lookup instead of vector<pair>
extern const std::unordered_map<std::string, std::string> drug_targets;
extern const std::unordered_map<std::string, int> half_lives;

// Add error handling and const correctness
inline std::string get_drug_target(const std::string& drug_name) {
    try {
        return drug_targets.at(drug_name);
    } catch (const std::out_of_range& e) {
        throw std::runtime_error("Drug target not found for: " + drug_name);
    }
}

inline int get_half_life(const std::string& drug_name) {
    try {
        return half_lives.at(drug_name);
    } catch (const std::out_of_range& e) {
        throw std::runtime_error("Half life not found for: " + drug_name);
    }
}

// Improve drug sensitivity value retrieval with error handling
struct DrugSensitivityParams {
    double max_conc;
    double xmid;
    double scale;
};

DrugSensitivityParams get_drug_sensitivity_values(const std::string& drug_name);

// Improve cell viability calculation with bounds checking
double get_cell_viability_for_drug_conc(
    double drug_conc, 
    const std::string& cell_line, 
    const std::string& drug_name
);

// Add validation for concentration calculations
class ConcentrationCalculator {
public:
    static double get_x_from_conc(double x_conc, double max_conc);
    static double get_conc_from_x(double x, double max_conc);
    static double get_lx_from_x(double x, double max_conc);
    
private:
    static void validate_concentration(double conc);
    static void validate_max_concentration(double max_conc);
};

} // namespace PhysiCell

#endif // __DRUG_SENSITIVITY_H__