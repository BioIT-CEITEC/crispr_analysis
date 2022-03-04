######################################
# wrapper for rule: DE_genes_MAGeCK
######################################
import os
import math
import subprocess
import re
from snakemake.shell import shell

shell.executable("/bin/bash")

f = open(snakemake.log.run, 'a+')
f.write("\n##\n## RULE: DE_genes_MAGeCK \n##\n")
f.close()

version = str(subprocess.Popen("conda list 2>&1 ", shell=True, stdout=subprocess.PIPE).communicate()[0], 'utf-8')
f = open(snakemake.log.run, 'at')
f.write("## CONDA: "+version+"\n")
f.close()

add_params = ""
if snakemake.params.paired:
  add_params += " --paired"
  
# PREFIX="DE_genes_mageck/NOTCH1_d21_vs_NOTCH1_d0.default/NOTCH1_d21_vs_NOTCH1_d0"; mkdir -p $PREFIX && $(which time) --verbose mageck test -k hsblastn_filter/all_samples_report.tsv -t NOTCH1_d21_rep1,NOTCH1_d21_rep2 -c NOTCH1_d0_rep1,NOTCH1_d0_rep2 --adjust-method fdr --remove-zero both --remove-zero-threshold 0 --pdf-report --normcounts-to-file -n $PREFIX 2>&1 | tee ${PREFIX}.log
command = f"mkdir -p $(dirname {snakemake.params.prefix});"+\
          f" $(which time) --verbose mageck test"+\
          f" -k {snakemake.input.tsv}"+\
          f" -t {','.join(snakemake.params.treats)}"+\
          f" -c {','.join(snakemake.params.ctrls)}"+\
          f" --adjust-method {snakemake.params.adj_type}"+\
          f" --remove-zero {snakemake.params.zero_type}"+\
          f" --remove-zero-threshold {str(snakemake.params.zero_value)}"+\
          f" --normcounts-to-file"+\
          f" --norm-method {snakemake.params.norm_type}"+\
          f"{add_params}"+\
          f" -n {snakemake.params.prefix} >> {snakemake.log.run} 2>&1"
f = open(snakemake.log.run, 'at')
f.write("## COMMAND: "+command+"\n")
f.close()
shell(command)

# command = f"mv {snakemake.params.gene} {snakemake.output.gene} >> {snakemake.log.run} 2>&1"
# f = open(snakemake.log.run, 'at')
# f.write("## COMMAND: "+command+"\n")
# f.close()
# shell(command)

