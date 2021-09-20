    ######################################
# wrapper for rule: DE_genes_edgeR
######################################
import os
import math
import subprocess
import re
from snakemake.shell import shell

shell.executable("/bin/bash")

f = open(snakemake.log.run, 'a+')
f.write("\n##\n## RULE: DE_genes_edgeR \n##\n")
f.close()

version = str(subprocess.Popen("conda list 2>&1 | grep 'r-base'", shell=True, stdout=subprocess.PIPE).communicate()[0], 'utf-8')
f = open(snakemake.log.run, 'at')
f.write("## VERSION: "+version+"\n")
f.close()

command = f"Rscript {snakemake.params.script} {snakemake.input.idx} {snakemake.input.tsv} {snakemake.params.design} {snakemake.params.odir} {snakemake.params.top} >> {snakemake.log.run} 2>&1"
f = open(snakemake.log.run, 'at')
f.write("## COMMAND: "+command+"\n")
f.close()
shell(command)

command = f"mv {snakemake.params.gene} {snakemake.output.gene} >> {snakemake.log.run} 2>&1"
f = open(snakemake.log.run, 'at')
f.write("## COMMAND: "+command+"\n")
f.close()
shell(command)

command = f"mv {snakemake.params.sg} {snakemake.output.sg} >> {snakemake.log.run} 2>&1"
f = open(snakemake.log.run, 'at')
f.write("## COMMAND: "+command+"\n")
f.close()
shell(command)

command = f"mv {snakemake.params.pdf} {snakemake.output.pdf} >> {snakemake.log.run} 2>&1"
f = open(snakemake.log.run, 'at')
f.write("## COMMAND: "+command+"\n")
f.close()
shell(command)

