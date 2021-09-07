######################################
# wrapper for rule: resolve_dups
######################################
import subprocess
from snakemake.shell import shell
shell.executable("/bin/bash")
log_filename = str(snakemake.log)

f = open(log_filename, 'a+')
f.write("\n##\n## RULE: resolve_dups \n##\n")
f.close()

version = str(subprocess.Popen("conda list 2>&1 | grep 'python'", shell=True, stdout=subprocess.PIPE).communicate()[0], 'utf-8')
f = open(log_filename, 'at')
f.write("## VERSION: "+version+"\n")
f.close()

command = f"python {snakemake.params.script} {snakemake.input.counts} >> {log_filename} 2>&1"
f = open(log_filename, 'at')
f.write("## COMMAND: "+command+"\n")
f.close()
shell(command)


