######################################
# wrapper for rule: hs_blast_filter
######################################
import os
import math
import subprocess
import re
from snakemake.shell import shell

shell.executable("/bin/bash")
log_file = str(snakemake.log)

f = open(log_file, 'a+')
f.write("\n##\n## RULE: hs_blast_filter \n##\n")
f.close()

version = str(subprocess.Popen("conda list 2>&1 ", shell=True, stdout=subprocess.PIPE).communicate()[0], 'utf-8')
f = open(log_file, 'at')
f.write("## CONDA: "+version+"\n")
f.close()

command = f"mkdir -p {snakemake.params.odir} >> {log_file} 2>&1"
f = open(log_file, 'at')
f.write("## COMMAND: "+command+"\n")
f.close()
shell(command)

command = f"(time Rscript {snakemake.params.script} {snakemake.input.idx} {snakemake.input.tsv} {snakemake.input.fa} {snakemake.input.counts} {snakemake.params.odir} {str(snakemake.threads)} {snakemake.params.sample} {snakemake.output.stats}) >> {log_file} 2>&1"
f = open(log_file, 'at')
f.write("## COMMAND: "+command+"\n")
f.close()
shell(command)
