######################################
# wrapper for rule: DE_genes_edgeR
######################################
import os
import math
import subprocess
import re
from snakemake.shell import shell

shell.executable("/bin/bash")
log_file = str(snakemake.log)

f = open(log_file, 'a+')
f.write("\n##\n## RULE: DE_genes_edgeR \n##\n")
f.close()

version = str(subprocess.Popen("conda list 2>&1 ", shell=True, stdout=subprocess.PIPE).communicate()[0], 'utf-8')
f = open(log_file, 'at')
f.write("## CONDA: "+version+"\n")
f.close()

command = f"(time Rscript {snakemake.params.script} {snakemake.input.idx} {snakemake.input.tsv} {snakemake.params.design} {snakemake.params.odir} {snakemake.params.top} {snakemake.params.paired}) >> {log_file} 2>&1"
f = open(log_file, 'at')
f.write("## COMMAND: "+command+"\n")
f.close()
shell(command)

command = f"mv {snakemake.params.gene} {snakemake.output.gene} >> {log_file} 2>&1"
f = open(log_file, 'at')
f.write("## COMMAND: "+command+"\n")
f.close()
shell(command)

command = f"mv {snakemake.params.sg} {snakemake.output.sg} >> {log_file} 2>&1"
f = open(log_file, 'at')
f.write("## COMMAND: "+command+"\n")
f.close()
shell(command)

command = f"mv {snakemake.params.pdf} {snakemake.output.pdf} >> {log_file} 2>&1"
f = open(log_file, 'at')
f.write("## COMMAND: "+command+"\n")
f.close()
shell(command)

