#!/usr/bin/env python3
# coding: utf-8

# Tool: physiboss_drugsim.py
# Date: January 2021
# Author: Annika Meert
# E-mail: annikameert@protonmail.com

#%% Imports 
import os, sys
import argparse
import re
import shutil
from itertools import combinations, product
from lxml import etree
import subprocess
import pandas as pd

current_wd = os.getcwd()
arg = sys.argv
parser = argparse.ArgumentParser()

parser.add_argument("-p", "--project", required=True, help="Name of the tissue collection project that the drug simulations are based on (ex. 'prostate' or 'drug_AGS_template').")
parser.add_argument("--cell_line", default="LNCaP", choices=["22Rv1", "BHP1", "DU145", "PC3", "LNCaP", "VCaP", "AGS"], help = "Cell line to be simulated.")
parser.add_argument("-d", "--drugs", required=True, help="Names of simulated drugs, comma separated (ex. 'Ipatasertib, Afatinib').")
parser.add_argument("-r", "--drug_rest", default="0", help="Proportion of drug resistances in the cells, between 0 and 1, comma separated (ex: '0, 0.2, 0.4, 0.6, 0.8, 1.0').")
parser.add_argument("-m", "--mode", default="both", choices=['single', 'double', 'both'], help="Mode of simulation for drug inhibition: single, double or both.")
parser.add_argument("--drug-concs", nargs='+', default=["IC10", "IC30", "IC50", "IC70", "IC90"], help="Drug concentrations to be simulated either in IC values e.g.: IC10 IC50 or in uM e.g: 0.1 2.3 .")
# parser.add_argument("--levels", default=3, type=int, help="[Deprecated: Use --drug-concs instead!] Number of levels for the drug simulation.")
# parser.add_argument("-i", "--input_cond", default='00', nargs='?', choices=['00', 'AR', 'AR_EGF', 'EGF'], help="Initial condition for drug simulation.")
parser.add_argument("-cl", "--cluster", default=False ,type=bool, help="Use of cluster or not. This flag creates the bash-script to run with SLURM workload manager.") 
parser.add_argument("--replicates", default=1, type=int, help="Number of replicates for each simulation.")

# example usage (to obtain the setup for the showcasing project in my master thesis):
# python3 addons/drugsims/physiboss_drugsim.py -p prostate --cell_line LNCaP -d Ipatasertib,Afatinib,Ulixertinib,Luminespib,Selumetinib,Pictilisib --mode both --drug-concs IC10 IC30 IC50 IC70 IC90 --replicates 10 -cl True
args = parser.parse_args()
#%% Process Arguments

# Dictionary for drug-node pairs
drug_node_pairs_prostate = {
    "Ipatasertib": "AKT",
    "Afuresertib": "AKT",
    "Afatinib": "EGFR",
    "Erlotinib": "EGFR",
    "Ulixertinib": "ERK",
    "Luminespib": "HSPs",
    "Trametinib": "MEK1_2",
    "Selumetinib": "MEK1_2",
    "Pictilisib":"PI3K",
    "Alpelisib": "PI3K",
    "BIBR1532": "TERT"
}

drug_node_pairs_AGS = {
    "GSK3" : "GSK3",
    "p38alpha"  : "p38alpha",
    "betacatenin" : "betacatenin",
    "TAK1" : "TAK1",
    "PI3K" : "PI3K",
    "MEK" : "MEK",
    "AKT" : "AKT"
}

cell_line_ids = {
    "LNCaP":907788,
    "BPH1": 924105,
    "DU145": 905935,
    "22Rv1": 924100,
    "VCaP": 1299075,
    "PC3": 905934
}

drug_ids = {
    "Ipatasertib": 1924,
    "Afuresertib": 1912,
    "Afatinib": 1032,
    "Erlotinib": 1168,
    "Ulixertinib": 2047,
    "Luminespib":1559,
    "Trametinib": 1372,
    "Selumetinib": 1736,
    "Pictilisib":1058,
    "Alpelisib": 1560,
    "BIBR1532": 2043
}

####################################################################
# Process arguments
####################################################################

print("Arguments:\n")
print(args)
print("\n")

project=args.project
cluster = args.cluster
cell_line = args.cell_line
# levels = args.levels
drug_concs = args.drug_concs
replicates = args.replicates

# check if the project is in sample projects 
if not os.path.isdir("sample_projects/" + project):
    print(project + " folder is not found")
    sys.exit(1)

drugs=args.drugs
drug_list = str(drugs.replace (" ", "").replace(",", ", ")).split(", ")
print("Drugs to be administered: "+str(drug_list).replace("[", "").replace("]","").replace("'",""))
# convert druglist into node-list
if (project == "prostate") :
    node_list = [drug_node_pairs_prostate.get(item,item)  for item in drug_list]
else :
    node_list = [drug_node_pairs_AGS.get(item,item)  for item in drug_list]


# specify boolean model path
#input_cond = args.input_cond
bool_model_path_dir = "{}/{}/{}/{}".format("sample_projects", project, "config", "boolean_network")
if (project == "prostate") :
    bool_model_filename = cell_line + "_mut_RNA_00"
else:
    bool_model_filename = cell_line
bool_model = "{}/{}".format(bool_model_path_dir, bool_model_filename)


node_list_1 = [x.split(", ") for x in node_list]

# set drug sensitivity path
drug_csv_path =  "{}/{}/{}/{}_{}".format("sample_projects", project, "config", project,"drug_sensitivity.csv")
drug_dataframe = pd.read_csv(drug_csv_path)

# create two list containing the value and the name of the chosen drug resistance levels 
drug_rest = args.drug_rest
drug_rest_value_list = [float(i) for i in drug_rest.replace(" ","").replace(",",", ").split(", ")]
drug_rest_name_list = [str(i).replace(".", "_") for i in drug_rest_value_list]
print("Applied drug resistance proportions: "+str(drug_rest_value_list).replace("[","").replace("]",""))
print("Applied drug concentrations: " + str(drug_concs).replace("[","").replace("]",""))
# set base paths for output and project folders 
sample_project_path = "sample_projects"
prostate_path = "{}/{}".format(sample_project_path, project)
project_name = "{}_{}_{}".format("physiboss_drugsim", project, cell_line)
project_path = "{}/{}".format(sample_project_path, project_name) 

####################################################################
# Function definitions
####################################################################

# adds a drug to a network by modifying and adding corresponding nodes in .bnd and .cfg files 
def add_nodes_to_network(bool_model, nodelist):
    # creating new .bnd file containing modified nodes and anti-nodes for all drugs
    bnd_file = "{}.{}".format(bool_model,"bnd")
    bnd_file_content = open(bnd_file).read()

    bnd_file_content = bnd_file_content.replace(" {", "{")
    bnd_file_content = bnd_file_content.replace("\n{", "{")
    bnd_file_content = bnd_file_content.replace("{", "\n{")

    bnd_nodes = bnd_file_content.split("Node")[1:]
    for item_index, item in enumerate(bnd_nodes):
        bnd_nodes[item_index] = "Node" + item 

    # filter nodes that are selected in the node_list 
    bnd_nodes_list = [re.findall(r'Node.*',line) for line in open(str(bnd_file))]
    bnd_nodes_list = list(filter(None, bnd_nodes_list))
    bnd_nodes_list = [[w.replace("Node ", "").replace(" {", "") for w in line] for line in bnd_nodes_list]

    # check if all selected nodes are present in the .bnd file 
    nodes_missing_bnd = [item for item in node_list_1 if item not in bnd_nodes_list]
    if len(nodes_missing_bnd) > 0:
        print("Nodes NOT in bnd file: " + str(nodes_missing_bnd).replace("[", "").replace("]", "").replace("'",""))

    new_nodes = str()
    for drug in nodelist:
        print("Modifying bnd file: Adding inhibitor node: " + str(drug))
        string = "Node {}\n".format(drug)
        item_count = 0
        for item_index, item in enumerate(bnd_nodes):
            if string in (bnd_nodes[item_index]):
                print("{} found".format(drug))
                item_count += 1
                for line_index, line in enumerate(bnd_nodes[item_index].split("\n")):
                    if " logic" in line:
                        new_nodes += line.split(";")[0].split(" = ")[0] + " = (" + line.split(";")[0].split(" = ")[1] + ") AND NOT anti_{};\n".format(drug)
                    else:
                        new_nodes += line + "\n"
                bnd_nodes.pop(item_index)
        if item_count == 0:
            print("{} not found".format(drug))
    
    for last_item in bnd_nodes:
        new_nodes += last_item 
    
    for inhibitor_to_add in nodelist:
        new_nodes += "Node anti_{}".format(inhibitor_to_add) + " \n{\n" + "\tlogic = (anti_{});\n".format(inhibitor_to_add) + "\trate_up = @logic ? 0 : 0;\n\trate_down = @logic ? 0 : 0;\n}\n\n\n"

    new_bnd_name = "{}_{}.{}".format(bool_model, "all_nodes", "bnd")
    fw1 = open(new_bnd_name,"w")
    fw1.write(new_nodes)
    fw1.close()

    # create new .cfg file setting anti-nodes on 0 and drug-nodes on random 
    cfg_file = "{}.{}".format(bool_model,"cfg")
    cfg_file_content = open(cfg_file).read()

    cfg_file_content = cfg_file_content.replace("[1]", " [1]").replace("[0]", " [0]").replace("] ,", "],").replace("=", " = ").replace("  ", " ")
    new_cfg = str()

    for line_index, line in enumerate(cfg_file_content.split("\n")):
        if ".istate = TRUE" in line:
            line = "[" + line.replace(".istate = TRUE", "].istate = 1 [1], 0 [0]")
        elif ".istate = FALSE" in line:
            line = "[" + line.replace(".istate = FALSE", "].istate = 0 [1], 1 [0]")
        elif "].istate = 0 [0], 1 [1]" in line:
            line = line.replace("].istate = 0 [0], 1 [1]", "].istate = 1 [1], 0 [0]")
        elif "].istate = 1 [0], 0 [1]" in line:
            line = line.replace("].istate = 1 [0], 0 [1]", "].istate = 0 [1], 1 [0]")
        else:
            line = line
        new_cfg += line + "\n"  

    b1 = [re.findall(r'.*istate.*',line) for line in open(str(bool_model+".cfg"))]
    b2 = list(filter(None, b1))
    cfg_nodes = [[re.sub(r"].istate.*","",w) for w in line] for line in b2]
    cfg_nodes = [[re.sub(r"\[","",w) for w in line] for line in cfg_nodes]
    missing_nodes_cfg = [item for item in bnd_nodes_list if item not in cfg_nodes]
    if (len(missing_nodes_cfg) > 0):
        missing_line = str([re.sub(r"$",".istate = 0.5 [0] , 0.5 [1];",str(line)) for line in missing_nodes_cfg]).replace("'","").replace('"[',"").replace('"]',"").replace(";",";\n").replace('", ','[')
        new_cfg += missing_line
    drugs_in_cfg = [item for item in nodelist if [item] in cfg_nodes]

            
    name_cfg = "{}_{}.cfg".format(bool_model,"all_nodes")
    fw2 = open(name_cfg, "w")
    fw2.write(new_cfg)
    fw2.close()

    for node in nodelist:
        f = open(name_cfg, "r")
        lines = f.readlines()
        f.close()
        nf = open(name_cfg, "w")
        for line in lines:
            if "[{}].istate = ".format(node) in line and not line.startswith("//"):
                new_line = line
                new_line = new_line = "[{}].istate = {} [1], {} [0];\n[anti_{}].istate = {} [1], {} [0];\n".format(node, "0.5", "0.5", node, "0", "1")
                nf.write(str(new_line))
            else:
                nf.write(str(line))
        nf.close()

def get_param_from_drug_df(drug_name, drug_df, param):
    drug_id = drug_ids.get(drug_name)
    cell_line_id = cell_line_ids.get(cell_line)
    filtered_df = drug_df[(drug_df["DRUG_ID_lib"] == drug_id) & (drug_df["CL"] == cell_line_id)]
    param = filtered_df.iloc[0][param]
    return param

# adds a drug to a physicell xml file
def add_drug_to_xml(drug, drug_concs, rest, path_to_xml, xml_output_path, model_name, mode, output_path):
    parser = etree.XMLParser(remove_blank_text=True)
    root = etree.parse(path_to_xml, parser).getroot()

    # set the output directory for the current run
    save = root.find("save")
    output_folder = save.find("folder")
    output_folder.text = output_path

    # insert drug as a density in the microenvironment 
    # 1. childs 
    microenv_setup = root.find('microenvironment_setup')

    # 2. childs
    variables = microenv_setup.findall('variable')
    variable_count = 0
    for variable in variables:
        variable_count += 1
    new_variable = etree.SubElement(microenv_setup, "variable")
    new_variable.set('name', drug)
    new_variable.set('units', 'mmol')
    new_variable.set("ID", str(variable_count))
    physical_parameter_set = etree.SubElement(new_variable, "physical_parameter_set")
    # initial_cond = etree.SubElement(new_variable, "initial_condition")
    # initial_cond.set('units', 'mmHg')
    # initial_cond.text = str(0.0)
    # dich = etree.SubElement(new_variable, "Dirichlet_boundary_condition")
    # dich.set('units', 'mmHg')
    # dich.set('enabled', 'true')
    # dich.text = str(0.0)

    # 3. childs 
    diffusion_coeff = etree.SubElement(physical_parameter_set, "diffusion_coefficient")
    diffusion_coeff.set('units', 'micron^2/min')
    diffusion_coeff.text = str(1200.0)
    decay_rate = etree.SubElement(physical_parameter_set, "decay_rate")
    decay_rate.set('units', '1/min')
    decay_rate.text = str(0.0275)
    
    cell_definitions = root.find('cell_definitions')
    # insert drug as a secreting substance for the default strain 
    cell_definition = cell_definitions.find('cell_definition')
    phenotype = cell_definition.find('phenotype')
    if (phenotype != None):
        secretion = phenotype.find('secretion')
    if (secretion != None):
        substrate = etree.SubElement(secretion, "substrate")
        substrate.set("name", drug)
        secretion_rate = etree.SubElement(substrate, "secretion_rate")
        secretion_rate.set("units", "1/min")
        secretion_rate.text = str(0.0)
        secretion_target = etree.SubElement(substrate, "secretion_target")
        secretion_target.set("units", "substrate density")
        secretion_target.text = str(1.0)
        uptake_rate = etree.SubElement(substrate, "uptake_rate")
        uptake_rate.set("units", "1/min")
        uptake_rate.text = str(0.0025)
        net_export_rate = etree.SubElement(substrate, "net_export_rate")
        net_export_rate.set("units", "total substrate/min")
        net_export_rate.text = str(0.0)

    # insert drug in custom data 
    custom_data = cell_definition.find('custom_data')
    if (custom_data != None):
        new_drug_conc = etree.SubElement(custom_data, drug + "_concentration")
        new_drug_conc.set("units", "dimensionless")
        new_drug_conc.text = str(0.0)
        new_drug_node = etree.SubElement(custom_data, drug + "_node")
        new_drug_node.set("units", "dimensionless")
        new_drug_node.text = str(0.0)

    # insert two new cell strains for the drug 
    new_cell_def_1 = etree.SubElement(cell_definitions, "cell_definition")
    new_cell_def_1.set("name", drug + "_sensitive")
    new_cell_def_1.set("ID", str(len(cell_definitions.getchildren())-1))
    new_cell_def_1.set("parent_type", "default")

    new_cell_def_2 = etree.SubElement(cell_definitions, "cell_definition")
    new_cell_def_2.set("name", drug + "_resistant")
    new_cell_def_2.set("ID", str(len(cell_definitions.getchildren())-1))
    new_cell_def_2.set("parent_type", "default")


    # get drug sensitivity parameters
    scal = get_param_from_drug_df(drug, drug_dataframe, "scal")
    xmid = get_param_from_drug_df(drug, drug_dataframe, "xmid")
    maxc = get_param_from_drug_df(drug, drug_dataframe, "maxc") 

    # insert drug in user parameters
    user_parameters = root.find('user_parameters')

    uptake_rate = etree.SubElement(user_parameters, drug + "_uptake_rate")
    uptake_rate.set("type", "double")
    uptake_rate.set("units", "")
    uptake_rate.text = str(0.0025)

    secretion_rate = etree.SubElement(user_parameters, drug + "_secretion_rate")
    secretion_rate.set("type", "double")
    secretion_rate.set("units", "fg/cell/min")
    secretion_rate.text = str(0.1)

    new_drug_scal = etree.SubElement(user_parameters, drug + "_scal")
    new_drug_scal.set("units", "dimensionless")
    new_drug_scal.set("type", "double")
    new_drug_scal.text = str(scal)
    new_drug_xmid = etree.SubElement(user_parameters, drug + "_xmid")
    new_drug_xmid.set("units", "ln(uM)")
    new_drug_xmid.set("type", "double")
    new_drug_xmid.text = str(xmid)
    new_drug_maxc = etree.SubElement(user_parameters, drug + "_maxc")
    new_drug_maxc.set("units", "uM")
    new_drug_maxc.set("type", "double")
    new_drug_maxc.text = str(maxc)

    # duration_add = etree.SubElement(user_parameters, "duration_add_" + drug)
    # duration_add.set("type", "int")
    # duration_add.set("units", "min")
    # duration_add.text = str(8000)

    # time_remove = etree.SubElement(user_parameters, "time_remove_" + drug)
    # time_remove.set("type", "int")
    # time_remove.set("units", "min")
    # time_remove.text = str(8000)

    # time_add = etree.SubElement(user_parameters, "time_add_" + drug)
    # time_add.set("type", "int")
    # time_add.set("units", "min")
    # time_add.text = str(0)

    # threshold = etree.SubElement(user_parameters, "threshold_" + drug)
    # threshold.set("type", "double")
    # threshold.set("units", "dimensionless")
    # threshold.text = str(0.14)

    # drug_rest = etree.SubElement(user_parameters, "concentration_" + drug)
    # drug_rest.set("type", "double")
    # drug_rest.set("units", "ng/mL")
    # drug_rest.text = str(0.5)

    # set the proportion of drug resistance for the drug 

    inhibition_level = etree.SubElement(user_parameters, "prop_drug_resistant_" + drug)
    inhibition_level.set("type", "double")
    inhibition_level.set("units", "dimensionless")
    inhibition_level.text = str(rest.replace("_", "."))

    # [Deprecated] set the total number of levels of drug concentration 
    # num_levels = user_parameters.find('total_concentration_levels')
    # if (num_levels == None):
    #     num_levels = etree.SubElement(user_parameters, "total_concentration_levels")
    #     num_levels.set("type", "int")
    #     num_levels.set("units", "dimensionless")
    #     num_levels.text = str(total_drug_levels)

    # # [Deprecated] set the current level of drug concentration 
    # current_level = etree.SubElement(user_parameters, "current_concentration_level_" + drug)
    # current_level.set("type", "int")
    # current_level.set("units", "dimensionless")
    # current_level.text = str(drug_level)

    # set the current drug concentration
    drug_concentration = etree.SubElement(user_parameters, "drug_concentration_" + drug)
    drug_concentration.set("type", "string")
    drug_concentration.set("units", "dimensionless")
    drug_concentration.text = str(drug_concs)
   
    # set the new bnd and cfg files 
    bnd_file = user_parameters.find('bnd_file') 
    # this path is for later when i have in the makefile saved where the files are 
    bnd_file.text = "{}/{}/{}/{}_{}.{}".format(".", "config", "boolean_network", model_name, "all_nodes", "bnd")
    cfg_file = user_parameters.find('cfg_file') 
    cfg_file.text = "{}/{}/{}/{}_{}.{}".format(".", "config", "boolean_network", model_name, "all_nodes", "cfg")

    # set the chosen simulation mode 
    simulation_mode = user_parameters.find("simulation_mode")
    if (simulation_mode == None):
        simulation_mode = etree.SubElement(user_parameters, "simulation_mode")
    if (mode == "single"):
        simulation_mode.text = "0"
    else:
        simulation_mode.text = "1"

    et = etree.ElementTree(root)
    et.write(xml_output_path, pretty_print=True)   


def add_project_to_makefile(project_name, makefile_path):
    with open(makefile_path, "r") as input_makefile:
        buf = input_makefile.readlines()

    with open(makefile_path, "w") as output_makefile:
        print("Add drug simulation " + project_name + " to makefile.")
        make_string = ""
        line_count = 11
        for line in buf:
            if  "prostate:" in line:
                # save the following lines in a string
                line_count = 0
            if line_count <= 10:
                make_string = make_string + line.replace("prostate", project_name)
                line_count += 1
            if "template_BM:" in line:
                line = make_string  + line 
            output_makefile.write(line)


def setup_drug_simulations(druglist, nodelist, bool_model_name, bool_model, project_path, rest_list, drug_concs, mode, simulation_list):


    # set variable to count the number of simulations 
    sim_count = 0
    # add nodes to network files
    add_nodes_to_network(bool_model, nodelist)
    
    translation_table = dict.fromkeys(map(ord, "[()'[] ]"), None)

    # drug_levels = list(range(1,levels+1))

    # modify nodelist and resistancelist if mode is double
    if (mode == "double"):
        drug_combinations = combinations(druglist, 2)
        rest_combinations = product(rest_list, repeat=2)
        # conc_combination = product(drug_levels, repeat=2)
        conc_combination = product(drug_concs, repeat=2)
        
        druglist = drug_combinations
        rest_combinations_list = []
        # drug_levels_list = []
        drug_concs_list = []

        for elem in rest_combinations:
            rest_combinations_list.append(elem)
        for elem in conc_combination:
            drug_concs_list.append(elem)
        rest_list = rest_combinations_list
        # drug_levels = drug_levels_list
        drug_concs = drug_concs_list 

    # create the project folder and copy the prostate project files in it
    if not os.path.exists(project_path):
        shutil.copytree(prostate_path, project_path)

    # set the output base path 
    if (mode == "single"):
        output_base_path = "output/single"
    else:
        output_base_path = "output/double"
    
    for drug in druglist:

        for rest in rest_list:

            # for drug_level in drug_levels:
            for drug_conc in drug_concs:

                for replicate in range(1, replicates + 1):

                    sim_count += 1

                    filtered_drugname = str(drug).translate(translation_table)
                    filtered_rest = str(rest).translate(translation_table)
                    # filtered_drug_level = str(drug_level).translate(translation_table)
                    filtered_drug_conc = str(drug_conc).translate(translation_table)
                    output_path = "{}/{}_{}_{}_{}_{}".format(output_base_path, bool_model_name,filtered_drug_conc.replace(",", "_"), filtered_drugname.replace(",","_"), filtered_rest.replace(",","_"), replicate)

                    # create the corresponding output folder
                    if os.path.exists(output_path):
                        shutil.rmtree(output_path)
                    os.makedirs(output_path)  

                    # modify the .xml file for the current replicate
                    xml_path = "{}/{}/{}_{}.{}".format(project_path, "config", "PhysiCell_settings", cell_line, "xml")
                    # new_xml_output_path = "{}/{}/{}_{}_{}_{}_{}.{}".format(project_path, "config", "settings", cell_line, filtered_drug_level.replace(",", "_"), filtered_drugname.replace(",","_"), filtered_rest.replace(",","_"), "xml")
                    new_xml_output_path = "{}/{}/{}_{}_{}_{}_{}_{}.{}".format(project_path, "config", "settings", cell_line, filtered_drug_conc.replace(",", "_"), filtered_drugname.replace(",","_"), filtered_rest.replace(",","_"), replicate, "xml")

                    if (type(drug) is tuple):
                        # for the tuples the first two elements of node and rest belong together
                        add_drug_to_xml(drug[0], drug_conc[0], rest[0],  xml_path, new_xml_output_path, bool_model_filename, mode, output_path)
                        add_drug_to_xml(drug[1], drug_conc[1], rest[1], new_xml_output_path, new_xml_output_path, bool_model_filename, mode, output_path)
                    else: 
                        add_drug_to_xml(drug, drug_conc, rest, xml_path, new_xml_output_path, bool_model_filename, mode, output_path)
                    xml_config_path = "{}/{}_{}_{}_{}_{}_{}.{}".format("config", "settings", cell_line, filtered_drug_conc.replace(",", "_"), filtered_drugname.replace(",","_"), filtered_rest.replace(",","_"), replicate, "xml")
                    simulation_list.append(xml_config_path) 

    # delete created folders again 
    # shutil.rmtree(project_path)
    # delete new .cfg and .bnd in the original prostate folder 
    os.remove("{}_{}.cfg".format(bool_model,"all_nodes"))
    os.remove("{}_{}.bnd".format(bool_model,"all_nodes"))

         

####################################################################
# Run drug simulations 
####################################################################

mode = args.mode
simulation_list = []
if (mode == "single" or mode == "double"):
    setup_drug_simulations(drug_list, node_list, bool_model_filename, bool_model,project_path, drug_rest_name_list, drug_concs, mode, simulation_list)
elif (mode == "both"):
    setup_drug_simulations(drug_list, node_list, bool_model_filename, bool_model, project_path, drug_rest_name_list, drug_concs, "single", simulation_list)
    setup_drug_simulations(drug_list, node_list, bool_model_filename, bool_model, project_path, drug_rest_name_list, drug_concs, "double", simulation_list)

# modify the project Makefile - rename all prostate to the project name 
project_makefile = project_path + "/Makefile"
with open(project_makefile, "r") as input_makefile:
    buf = input_makefile.readlines()

with open(project_makefile, "w") as output_makefile:
    for line in buf:
        if project in line:
            line = line.replace(project, project_name)
        output_makefile.write(line)

# add the new sample project to the Makefile in the main Physiboss folder
makefile_path = "Makefile"
add_project_to_makefile(project_name, makefile_path)

# rename the main file 
original_main_path = "{}/{}-{}.{}".format(project_path, "main", project, "cpp")
new_main_path = "{}/{}-{}.{}".format(project_path, "main", project_name, "cpp")
os.rename(original_main_path, new_main_path)

# compile the project
#subprocess.call(["make", project_name])
#subprocess.call(["make"])

# run all simulations with physiboss sequentially if cluster flag is not set
if (cluster == False):
    for xml_file in simulation_list:
        break
        #subprocess.call(["{}/{}".format(".", project_name), xml_file])
else:
    fw1 = open("./run_drug_simulations.sh", "w")
    for xml_file in simulation_list:
        fw1.write("{}/{}".format(".", project_name) + " " + xml_file + "\n")
    fw1.close()
    fw1 = open("./run_greasy.sh", "w")
    fw1.write("#!/bin/bash\n" )
    fw1.write('#SBATCH --job-name="drug_simulation"\n')
    fw1.write('#SBATCH --output=%j.out\n')
    fw1.write('#SBATCH --error=%j.err\n')
    fw1.write('#SBATCH --nodes=1\n')
    fw1.write('#SBATCH --ntasks=1\n')
    fw1.write('#SBATCH --tasks-per-node=6\n')
    fw1.write('#SBATCH --cpus-per-task=8\n')
    fw1.write('#SBATCH -t 22:00:00\n\n')
    fw1.write('/apps/GREASY/latest/INTEL/IMPI/bin/greasy run_drug_simulations.sh\n')