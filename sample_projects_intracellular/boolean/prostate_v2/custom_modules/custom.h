
#include "../core/PhysiCell.h"
#include "../modules/PhysiCell_standard_modules.h" 
#include "drug_sensitivity.h"
#include "boolean_model_interface.h"

/**
 *	\main drug_AGS custom
 *	\brief Custom module file for drug_AGS example
 * 
 *	\details Modules needed for the drug_AGS example. This custom module can be used to study the inhibition of AGS cell lines with AKT, beta-catenin and TAK inhibitors.
 *
 *
 *	\date 19/10/2020
 *	\author Arnau Montagud, BSC-CNS, with code previously developed by Gerard Pradas and Miguel Ponce de Leon, BSC-CNS
 */

using namespace BioFVM; 
using namespace PhysiCell;

struct init_record
{
	float x;
	float y;
	float z;
	float radius;
	int phase;
	double elapsed_time;
};

// setup functions to help us along 
void create_cell_types( void );
void setup_tissue( void ); 

// set up the BioFVM microenvironment 
double get_decay_rate(double half_life);
void setup_microenvironment( void ); 

// custom pathology coloring function 
std::vector<std::string> prolif_apoptosis_coloring( Cell* );

// custom cell phenotype functions could go here 
void tumor_cell_phenotype_with_signaling( Cell* pCell, Phenotype& phenotype, double dt );

std::vector<init_record> read_init_file(std::string filename, char delimiter, bool header);

inline float sphere_volume_from_radius(float radius) {return 4/3 * PhysiCell_constants::pi * std::pow(radius, 3);}
