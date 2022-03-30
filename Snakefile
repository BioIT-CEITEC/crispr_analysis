import os
import pandas as pd
import itertools
import json
from snakemake.utils import R
from snakemake.utils import report
from os.path import split
from snakemake.utils import min_version

min_version("5.18.0")

configfile: "config.json"

GLOBAL_REF_PATH = "/mnt/references"
GLOBAL_TMPD_PATH = "/tmp/"

# Reference processing
#

if not "adaptors" in config:
    config["adaptors"] = "null"

if not "full_name_rep" in config:
    config["full_name_rep"] = "no.rep.to.join"

if not "crispr_type" in config:
    config['crispr_type'] = 'CRISPR_Brunello'
if not "adapter" in config:
    config['adapter'] = "GGAAAGGACGAAACACCG"
if not 'error_rate' in config:
    config['error_rate'] = '0.2'    # Allowed error rate of adapters
if not 'min_overlap' in config:
    config['min_overlap'] = 18   # Minimal overlap of adapters
if not 'times' in config:
    config['times'] = 1     # Maximal number of adapters to remove
if not 'min_len' in config:
    config['min_len'] = 20  # Discard length of sequences
if not 'guide_len' in config:
    config['guide_len'] = 20    # Remaining length of the guide sequence
if not 'conditions_to_compare' in config:
    config['conditions_to_compare'] = 'all'
if not 'top_genes' in config:
    config['top_genes'] = 10    # Top genes to be shown in resulting graphs.pdf from DE part
if not 'use_tag_to_pair_samples' in config:
    config['use_tag_to_pair_samples'] = False
if not 'mageck_padj_method' in config:
    config["mageck_padj_method"] = "fdr" # fdr,holm,pounds
if not 'mageck_zero_type' in config:
    config["mageck_zero_type"] = "both" # none,control,treatment,both,any
if not 'mageck_zero_value' in config:
    config["mageck_zero_value"] = 0
if not 'mageck_norm_method' in config:
    config["mageck_norm_method"] = "median" # none, median, total


# setting organism from reference
f = open(os.path.join(GLOBAL_REF_PATH,"reference_info","reference.json"),)
reference_dict = json.load(f)
f.close()
# config["organism"] = [organism_name.lower().replace(" ","_") for organism_name in reference_dict.keys() if isinstance(reference_dict[organism_name],dict) and config["reference"] in reference_dict[organism_name].values()][0]


##### Config processing #####
# Folders
#
# reference_directory = os.path.join(GLOBAL_REF_PATH,config["organism"],config["reference"])
reference_directory = os.path.join(GLOBAL_REF_PATH, "general", config["crispr_type"])

# Samples
#
sample_tab = pd.DataFrame.from_dict(config["samples"],orient="index")

if not config["is_paired"]:
    read_pair_tags = [""]
    paired = "SE"
else:
    read_pair_tags = ["_R1","_R2"]
    paired = "PE"

# if config["lib_reverse_read_length"] == 0:
#     read_pair_tags = [""]
# else:
#     read_pair_tags = ["_R1","_R2"]

wildcard_constraints:
    sample = "|".join(sample_tab.sample_name),
    read_pair_tag = "(_R.)?",
    rep = "(_rep.)?",


##### Target rules #####
rule all:
    input: "final_report.html"

##### Modules #####

include: "rules/crispr_analysis.smk"
