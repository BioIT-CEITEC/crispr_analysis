######################################
# wrapper for rule: report_counts
######################################
import subprocess
from snakemake.shell import shell
shell.executable("/bin/bash")
log_filename = str(snakemake.log)

f = open(log_filename, 'a+')
f.write("\n##\n## RULE: report_counts \n##\n")
f.close()

version = str(subprocess.Popen("conda list 2>&1 | grep 'r-base'", shell=True, stdout=subprocess.PIPE).communicate()[0], 'utf-8')
f = open(log_filename, 'at')
f.write("## VERSION: "+version+"\n")
f.close()

command = f"Rscript {snakemake.params.script} {snakemake.params.prefix} {snakemake.input.idx} {' '.join(snakemake.input.table)} >> {log_filename} 2>&1"
f = open(log_filename, 'at')
f.write("## COMMAND: "+command+"\n")
f.close()
shell(command)
