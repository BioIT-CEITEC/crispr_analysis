######################################
# wrapper for rule: merge_hs_blast
######################################
from snakemake.shell import shell
shell.executable("/bin/bash")
log_filename = str(snakemake.log)

f = open(log_filename, 'a+')
f.write("\n##\n## RULE: merge_hs_blast \n##\n")
f.close()

# version = str(subprocess.Popen("conda list 2>&1 | grep 'r-base'", shell=True, stdout=subprocess.PIPE).communicate()[0], 'utf-8')
# f = open(snakemake.log.run, 'at')
# f.write("## VERSION: "+version+"\n")
# f.close()

command = "echo \"sample$(printf '\t')$(head -1 "+snakemake.input.stats[0]+")\" > "+snakemake.output.table+" && for i in "+" ".join(snakemake.input.stats)+"; do echo \"$(basename ${{i%_statistics.tsv}})$(printf '\t')$(tail -n +2 $i)\"; done >> "+snakemake.output.table+" 2>> "+log_filename
f = open(log_filename, 'at')
f.write("## COMMAND: "+command+"\n")
f.close()
shell(command)

