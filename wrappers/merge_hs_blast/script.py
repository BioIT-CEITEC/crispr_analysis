######################################
# wrapper for rule: merge_hs_blast
######################################
import os
import math
import subprocess
import re
from snakemake.shell import shell

shell.executable("/bin/bash")
log_file = str(snakemake.log)

f = open(log_file, 'a+')
f.write("\n##\n## RULE: merge_hs_blast \n##\n")
f.close()

version = str(subprocess.Popen("conda list 2>&1 ", shell=True, stdout=subprocess.PIPE).communicate()[0], 'utf-8')
f = open(log_file, 'at')
f.write("## CONDA: "+version+"\n")
f.close()

command = "echo \"sample$(printf '\t')$(head -1 "+snakemake.input.stats[0]+")\" > "+snakemake.output.table+" && for i in "+" ".join(snakemake.input.stats)+"; do echo \"$(basename ${{i%_statistics.tsv}})$(printf '\t')$(tail -n +2 $i)\"; done >> "+snakemake.output.table+" 2>> "+log_file
f = open(log_file, 'at')
f.write("## COMMAND: "+command+"\n")
f.close()
shell(command)

