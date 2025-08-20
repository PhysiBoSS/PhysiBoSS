
#include "../core/PhysiCell.h"
#include "../modules/PhysiCell_standard_modules.h" 
// #include "../addons/PhysiBoSS/src/maboss_network.h"
#include "drug_sensitivity.h"
#include "boolean_model_interface.h"

/**
 *	\main drug custom
 *	\brief Custom module file for  example
 * 
 *	\details Modules needed for the drug example. This custom module can be used to study the inhibition of cell lines with inhibitors.
 *
 *
 *	\date 20/08/2025
 *	\author Arnau Montagud, Annika Meert, Gerard Pradas and Miguel Ponce de Leon, BSC-CNS
 */

using namespace BioFVM; 
using namespace PhysiCell;

void create_cell_types( void );
double get_decay_rate(double half_life);

// set up the BioFVM microenvironment 
void setup_microenvironment( void ); 

void setup_tissue( void ); 
void setup_tissue_resistant( void ); 

void phenotype_function( Cell* pCell, Phenotype& phenotype, double dt );
void pre_update_intracellular(Cell* pCell, Phenotype& phenotype, double dt);
void post_update_intracellular(Cell* pCell, Phenotype& phenotype, double dt);

// custom cell phenotype functions could go here 
void tumor_cell_phenotype_with_signaling( Cell* pCell, Phenotype& phenotype, double dt );

//  needed for setup_tissue_resistant
struct init_record
{
	float x;
	float y;
	float z;
	float radius;
	int phase;
	double elapsed_time;
};

// custom pathology coloring function 
// std::vector<std::string> prolif_apoptosis_coloring( Cell* );

std::vector<init_record> read_init_file(std::string filename, char delimiter, bool header);

inline float sphere_volume_from_radius(float radius) {return 4/3 * PhysiCell_constants::pi * std::pow(radius, 3);}
