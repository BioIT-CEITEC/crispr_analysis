######################################
# wrapper for rule: resolve_dups
######################################
import os
import math
import subprocess
import re
from snakemake.shell import shell

shell.executable("/bin/bash")
log_file = str(snakemake.log)

f = open(log_file, 'a+')
f.write("\n##\n## RULE: resolve_dups \n##\n")
f.close()

version = str(subprocess.Popen("conda list 2>&1", shell=True, stdout=subprocess.PIPE).communicate()[0], 'utf-8')
f = open(log_file, 'at')
f.write("## CONDA: "+version+"\n")
f.close()

command = f"python {snakemake.params.script} {snakemake.input.counts} >> {log_file} 2>&1"
f = open(log_file, 'at')
f.write("## COMMAND: "+command+"\n")
f.close()
shell(command)

