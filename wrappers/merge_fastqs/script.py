######################################
# wrapper for rule: merge_fastqs
######################################
import os
import sys
import math
import subprocess
from snakemake.shell import shell

shell.executable("/bin/bash")
log_file = str(snakemake.log)

f = open(log_file, 'a+')
f.write("\n##\n## RULE: merge_fastqs \n##\n")
f.close()

version = str(subprocess.Popen("conda list 2>&1", shell=True, stdout=subprocess.PIPE).communicate()[0], 'utf-8')
f = open(log_file, 'at')
f.write("## CONDA: "+version+"\n")
f.close()

#merge input fastq files
command = " cat " + " ".join(sorted(snakemake.input.fastq)) + " > " + snakemake.output.fastq
f = open(log_file, 'at')
f.write("## COMMAND: "+command+"\n")
f.close()
shell(command)

command = "mkdir -p `dirname "+snakemake.output.html+"`"+" >> "+log_file+" 2>&1"
f = open(log_file, 'at')
f.write("## COMMAND: "+command+"\n")
f.close()
shell(command)

command = "(time fastqc -o "+snakemake.params.prefix+" "+snakemake.params.extra+" --threads "+str(snakemake.threads)+" "+snakemake.output.fastq+") >> "+log_file+" 2>&1 "
f = open(log_file, 'at')
f.write("## COMMAND: "+command+"\n")
f.close()
shell(command)

command = ' sed -r "s:<h2>[ ]*Summary[ ]*<\/h2><ul>:&<li><b>Return to <a href=\'../final_report.html\'>start page<\/a><\/b><\/li>:" '+snakemake.params.html+' > '+snakemake.output.html
f = open(log_file, 'at')
f.write("## COMMAND: "+command+"\n")
f.close()
shell(command)

command = "rm -f "+snakemake.params.html+" >> "+log_file+" 2>&1"
f = open(log_file, 'at')
f.write("## COMMAND: "+command+"\n")
f.close()
shell(command)
