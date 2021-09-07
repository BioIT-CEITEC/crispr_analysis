import os
import pandas as pd
import json
from snakemake.utils import min_version


min_version("5.18.0")

GLOBAL_REF_PATH = "/mnt/references/"

# Reference processing
#

if not "adaptors" in config:
    config["adaptors"] = "null"

if not "full_name_rep" in config:
    config["full_name_rep"] = "no.rep.to.join"

# if not "contaminants" in config:
#     config['contaminants'] = ['all'] * len(config)
#
# if not "crispr_type" in config:
#     config['crispr_type'] = ['CRISPR_Brunello'] * len(config)
# if not "adapter" in config:
#     config['adapter'] = ["GGAAAGGACGAAACACCG"] * len(config)
# if not 'error_rate' in config:
#     config['error_rate'] = ['0.2'] * len(config)  # Allowed error rate of adapters
# if not 'min_overlap' in config:
#     config['min_overlap'] = [18] * len(config)  # Minimal overlap of adapters
# if not 'times' in config:
#     config['times'] = [1] * len(config)  # Maximal number of adapters to remove
# if not 'min_len' in config:
#     config['min_len'] = [20] * len(config)  # Discard length of sequences
# if not 'guide_len' in config:
#     config['guide_len'] = [20] * len(config)  # Remaining length of the guide sequence
# if not 'conditions_to_compare' in config:
#     # cfg['conditions_to_compare'] = ['all']*len(cfg)
#     config['conditions_to_compare'] = ['D0_WT:D0_dA5,D0_WT:D22_dA5,D0_WT:D7_dA5,D22_WT:D22_dA5,D7_WT:D7_dA5,D0_dA5:D22_dA5,D0_dA5:D7_dA5,D0_WT:D22_WT,D0_WT:D7_WT'] * len(config)
# if not 'top_genes' in config:
#     config['top_genes'] = [10] * len(config)  # Top genes to be shown in resulting graphs.pdf from DE part



if config["lib_ROI"] != "wgs":
    # setting reference from lib_ROI
    f = open(os.path.join(GLOBAL_REF_PATH,"reference_info","lib_ROI.json"))
    lib_ROI_dict = json.load(f)
    f.close()
    config["reference"] = [ref_name for ref_name in lib_ROI_dict.keys() if isinstance(lib_ROI_dict[ref_name],dict) and config["lib_ROI"] in lib_ROI_dict[ref_name].keys()][0]


# setting organism from reference
f = open(os.path.join(GLOBAL_REF_PATH,"reference_info","reference.json"),)
reference_dict = json.load(f)
f.close()
config["organism"] = [organism_name.lower().replace(" ","_") for organism_name in reference_dict.keys() if isinstance(reference_dict[organism_name],dict) and config["reference"] in reference_dict[organism_name].values()][0]


##### Config processing #####
# Folders
#
reference_directory = os.path.join(GLOBAL_REF_PATH,config["organism"],config["reference"])

# Samples
#
sample_tab = pd.DataFrame.from_dict(config["samples"],orient="index")

if config["lib_reverse_read_length"] == 0:
    read_pair_tags = [""]
else:
    read_pair_tags = ["_R1","_R2"]

wildcard_constraints:
    sample = "|".join(sample_tab.sample_name),
    read_pair_tag = "(_R.)?"


##### Target rules #####
rule all:
    input: "final_report.html"

##### Modules #####

include: "rules/crispr_analysis.smk"
