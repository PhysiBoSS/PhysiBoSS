#ifndef __BOOLEAN_MODEL_INTERFACE_H__
#define __BOOLEAN_MODEL_INTERFACE_H__

#include "../core/PhysiCell.h"
#include "./drug_sensitivity.h"
#include <memory>

namespace PhysiCell {

class BooleanModelInterface {
public:
    static void update_custom_variables(Cell* pCell);
    static void set_boolean_node(Cell* pCell, const std::string& drug_name, int drug_index, double threshold);
    static void set_input_nodes(Cell* pCell);
    
private:
    static void handle_single_inhibition(Cell* pCell);
    static void handle_double_inhibition(Cell* pCell);
};

} // namespace PhysiCell

#endif