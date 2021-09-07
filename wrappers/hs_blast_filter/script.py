######################################
# wrapper for rule: hs_blast_filter
######################################
import subprocess
from snakemake.shell import shell
shell.executable("/bin/bash")
log_filename = str(snakemake.log)

f = open(log_filename, 'a+')
f.write("\n##\n## RULE: hs_blast_filter \n##\n")
f.close()

version = str(subprocess.Popen("conda list 2>&1 | grep 'r-base'", shell=True, stdout=subprocess.PIPE).communicate()[0], 'utf-8')
f = open(log_filename, 'at')
f.write("## VERSION: "+version+"\n")
f.close()

command = f"mkdir -p {snakemake.params.odir} >> {log_filename} 2>&1"
f = open(log_filename, 'at')
f.write("## COMMAND: "+command+"\n")
f.close()
shell(command)

command = f"(time Rscript {snakemake.params.script} {snakemake.input.idx} {snakemake.input.tsv} {snakemake.input.fa} {snakemake.input.counts} {snakemake.params.odir} {str(snakemake.threads)} {snakemake.params.sample}) >> {log_filename} 2>&1"
f = open(log_filename, 'at')
f.write("## COMMAND: "+command+"\n")
f.close()
shell(command)
