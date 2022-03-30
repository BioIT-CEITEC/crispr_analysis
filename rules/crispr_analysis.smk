import math
import subprocess
import json
import re
import os.path
import glob
import fnmatch
import pandas as pd
import itertools
from snakemake.utils import R
from snakemake.utils import report
from os.path import split


def DE_conditions(wildcards):
    inputs = []
    os.makedirs("DE_results", exist_ok=True)
    if 'conditions_to_compare' in config and config["conditions_to_compare"] != "":
        if config["conditions_to_compare"] == "all":
            # create all pairs of conditions
            conditions = itertools.permutations(set(sample_tab.condition),2)
        else:
            # use only specified conditions
            conditions = [x.split(":") for x in config["conditions_to_compare"].split(",")]
            
        for cond in conditions:
            if len(cond) > 2:
                raise NotImplementedError("Only comparison of two conditions is supported!")
            elif len(cond) < 2:
                raise ValueError("Not enough conditions to compare, exactly 2 conditions needed!")
            else:
                inputs.append(f"DE_results/MAGeCK/{cond[0]}_vs_{cond[1]}/{cond[0]}_vs_{cond[1]}.gene_summary.tsv")
                inputs.append(f"DE_results/edgeR/{cond[0]}_vs_{cond[1]}/{cond[0]}_vs_{cond[1]}.gene_summary.tsv")
                design = f"DE_results/edgeR/{cond[0]}_vs_{cond[1]}/{cond[0]}_vs_{cond[1]}.design_table.tsv"
                os.makedirs(os.path.dirname(design), exist_ok=True)
                sample_tab.loc[sample_tab.condition == cond[0],("sample_name","sample_name","condition","tag")].set_axis(['sample', 'name', 'condition', 'patient'], axis='columns').drop_duplicates().to_csv(design, sep="\t", index=False)
                sample_tab.loc[sample_tab.condition == cond[1],("sample_name","sample_name","condition","tag")].set_axis(['sample', 'name', 'condition', 'patient'], axis='columns').drop_duplicates().to_csv(design, sep="\t", header=False, index=False, mode="a")
    return inputs

rule final_report:
    input:  pdf = "counts/all_samples_report.pdf",
            table = "mapped/all_samples.stats.tsv",
            DEs = DE_conditions,
    output: report = "final_report.html",
    params: wdir = "./",
            name = config["task_name"],
            guide_len = config["guide_len"],
            ref = config["crispr_type"],
            top = config["top_genes"],
            cfg = json.dumps(config),
    conda:  "../wrappers/final_report/env.yaml"
    script: "../wrappers/final_report/crispr_analysis_report_template.Rmd"
    

rule DE_genes_MAGeCK:
    input:  tsv = "counts/all_samples_report.tsv",
            idx = expand("{ref_dir}/{ref}_mod.csv", ref_dir=reference_directory, ref=config["crispr_type"])[0],
    output: gene = "DE_results/MAGeCK/{c1}_vs_{c2}/{c1}_vs_{c2}.gene_summary.tsv",
            sg   = "DE_results/MAGeCK/{c1}_vs_{c2}/{c1}_vs_{c2}.sgrna_summary.tsv",
            #pdf  = "DE_results/MAGeCK/{c1}_vs_{c2}/{c1}_vs_{c2}.graphs.pdf",
    log:    run  = "DE_results/MAGeCK/{c1}_vs_{c2}/{c1}_vs_{c2}.DE_genes_MAGeCK.log",
    params: prefix = "DE_results/MAGeCK/{c1}_vs_{c2}/{c1}_vs_{c2}",
            treats = lambda ws: sample_tab.loc[sample_tab.condition == ws.c2, "sample_name"].tolist(),
            ctrls  = lambda ws: sample_tab.loc[sample_tab.condition == ws.c1, "sample_name"].tolist(),
            # script = workflow.basedir+"/../wraps/crispr_analysis/DE_results/MAGeCK/DE_genes_CRISPR_edgeR.R",
            gene = "DE_results/MAGeCK/{c1}_vs_{c2}/{c1}_vs_{c2}.gene_summary.txt",
            sg   = "DE_results/MAGeCK/{c1}_vs_{c2}/{c1}_vs_{c2}.sgrna_summary.txt",
            norm_file = "DE_results/MAGeCK/{c1}_vs_{c2}/{c1}_vs_{c2}.normalized.txt",
            report_source   = "DE_results/MAGeCK/{c1}_vs_{c2}/{c1}_vs_{c2}.report.Rmd",
            pdf_source  = "DE_results/MAGeCK/{c1}_vs_{c2}/{c1}_vs_{c2}.R",
            norm_type = config["mageck_norm_method"],
            zero_type = config["mageck_zero_type"],
            zero_value= config["mageck_zero_value"],
            adj_type =  config["mageck_padj_method"],
            top  = config["top_genes"],
            paired = config["use_tag_to_pair_samples"],
    conda:  "../wrappers/DE_genes_MAGeCK/env.yaml"
    script: "../wrappers/DE_genes_MAGeCK/script.py"
    

rule DE_genes_edgeR:
    input:  tsv = "counts/all_samples_report.tsv",
            idx = expand("{ref_dir}/{ref}_mod.csv", ref_dir=reference_directory, ref=config["crispr_type"])[0],
    output: gene = "DE_results/edgeR/{cond}/{cond}.gene_summary.tsv",
            sg   = "DE_results/edgeR/{cond}/{cond}.sgrna_summary.tsv",
            pdf  = "DE_results/edgeR/{cond}/{cond}.graphs.pdf",
    log:    "DE_results/edgeR/{cond}/{cond}.DE_genes_edgeR.log",
    params: script = workflow.basedir+"/wrappers/DE_genes_edgeR/DE_genes_CRISPR_edgeR.R",
            odir = "DE_results/edgeR/{cond}/",
            design = "DE_results/edgeR/{cond}/{cond}.design_table.tsv",
            gene = "DE_results/edgeR/{cond}/gene_summary.tsv",
            sg   = "DE_results/edgeR/{cond}/sgRNAs_summary.tsv",
            pdf  = "DE_results/edgeR/{cond}/graphs.pdf",
            top  = config["top_genes"],
            paired = config["use_tag_to_pair_samples"],
    conda:  "../wrappers/DE_genes_edgeR/env.yaml"
    script: "../wrappers/DE_genes_edgeR/script.py"

####################################
# PREPROCESS RULEs
#

# def check_contamination_input(wildcards):
#     preprocessed = "preprocessed_data"
# 
#     if config.loc[config[SAMPLE_REP] == wildcards.sample,"paired"] == "PE":
#         return {
#             'r1': ""+preprocessed+"/"+ wildcards.sample +"_R1.fastq.gz",
#             'r2': ""+preprocessed+"/"+ wildcards.sample +"_R2.fastq.gz",
#             }
#     else:
#         return {
#             'reads': ""+preprocessed+"/"+ wildcards.sample +"_SE.fastq.gz"
#             }
# 
# rule check_contamination:
#     input:  unpack(check_contamination_input),
#     output: table = "preprocessed_data/{sample}.contamination_summary.tsv",
#     log:    "logs/{sample}/check_contamination.log",
#     threads:    8 
#     resources:  mem = 30
#     params: tool = RES_DIR + "biobloomtools/BioBloomCategorizer/biobloomcategorizer",
#             prefix = "preprocessed_data/{sample}.contamination",
#             filters= lambda wildcards: config.loc[config[SAMPLE] == wildcards.sample,"contaminants"],
#     conda:  "../wrappers/fastq2bam_RNA/check_contamination/env.yaml"
#     script: "../wrappers/fastq2bam_RNA/check_contamination/script.py"

rule report_counts:
    input:  table = set(expand("counts/{sample}/{sample}_unique_inserts_tab_counts_resolved_duplicates.tsv", sample = sample_tab.sample_name)),
            idx = expand("{ref_dir}/{ref}_mod.csv", ref_dir=reference_directory, ref=config["crispr_type"])[0],
    output: pdf = "counts/all_samples_report.pdf",
            tsv = "counts/all_samples_report.tsv",
    log:    "logs/report_counts.log",
    params: script = workflow.basedir+"/wrappers/report_counts/counts_stat_CRISPR.R",
            prefix = "counts/all_samples_report",
    conda:  "../wrappers/report_counts/env.yaml"
    script: "../wrappers/report_counts/script.py"
    
    
rule resolve_dups:
    input:  counts = "counts/{sample}/{sample}_unique_inserts_tab_counts.tsv",
    output: table = "counts/{sample}/{sample}_unique_inserts_tab_counts_resolved_duplicates.tsv",
    log:    "logs/{sample}/resolve_dups.log",
    params: script = workflow.basedir+"/wrappers/resolve_dups/resolve_count_duplicates.py",
    conda:  "../wrappers/resolve_dups/env.yaml"
    script: "../wrappers/resolve_dups/script.py"
    

rule merge_hs_blast:
    input:  stats = set(expand("mapped/{sample}.stats.tsv", sample = sample_tab.sample_name)),
    output: table = "mapped/all_samples.stats.tsv",
    log:    "logs/merge_hs_blast.log",
    # conda:  "../wrappers/merge_hs_blast/env.yaml"
    script: "../wrappers/merge_hs_blast/script.py"
    

rule hs_blast_filter:
    input:  tsv = "mapped/{sample}.hsblastn.tsv",
            fa  = "preprocessed_data/{sample}.inserts.fa",
            counts = "preprocessed_data/{sample}.inserts.counts",
            idx = expand("{ref_dir}/{ref}.uniq_revcomp_ins.csv", ref_dir=reference_directory, ref=config["crispr_type"])[0],
    output: stats = "mapped/{sample}.stats.tsv",
            uct   = "counts/{sample}/{sample}_unique_inserts_tab_counts.tsv",
    log:    "logs/{sample}/hs_blast_filter.log",
    threads:    10
    params: odir = "counts/{sample}/",
            sample = "{sample}_",
            script = workflow.basedir+"/wrappers/hs_blast_filter/filter_hsblastn_counts.R",
    conda:  "../wrappers/hs_blast_filter/env.yaml"
    script: "../wrappers/hs_blast_filter/script.py"
    

rule hs_blast_alignment:
    input:  ins = "preprocessed_data/{sample}.inserts.seqs",
            idx = expand("{ref_dir}/index/{ref}.uniq_revcomp_ins.fa", ref_dir=reference_directory, ref=config["crispr_type"])[0],
    output: out = "mapped/{sample}.hsblastn.tsv",
            fa  = "preprocessed_data/{sample}.inserts.fa",
            counts = "preprocessed_data/{sample}.inserts.counts",
    log:    "logs/{sample}/hs_blast_alignment.log",
    threads: 20,
    conda:  "../wrappers/hs_blast_alignment/env.yaml"
    script: "../wrappers/hs_blast_alignment/script.py"


rule preprocess_SE:
    input:  R1 = expand("raw_fastq/{{sample}}{read_pair_tags}.fastq.gz",read_pair_tags = read_pair_tags)[0],
    output: R1 = expand("preprocessed_data/{{sample}}{read_pair_tags}.fastq.gz",read_pair_tags = read_pair_tags)[0],
            seqs = "preprocessed_data/{sample}.seqs",
            inserts = "preprocessed_data/{sample}.inserts.seqs",
            html = "preprocessed_data/{sample}.fastqc.html"
    log:    "logs/{sample}/preprocess_SE.log",
    threads:    10
    resources:  mem = 10
    params: adapter = config["adapter"],
            error_rate = config["error_rate"],
            min_overlap = config["min_overlap"],
            times = config["times"],
            min_len = config["min_len"],
            guide_len = config["guide_len"],
            prefix = "preprocessed_data/",
            qc_html_R1 = "preprocessed_data/{sample}.fastqc.html",
            qc_html_R1_tmp = "preprocessed_data/{sample}_fastqc.html",
    conda:  "../wrappers/preprocess_SE/env.yaml"
    script: "../wrappers/preprocess_SE/script.py"

