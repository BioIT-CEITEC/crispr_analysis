######################################
# wrapper for rule: merge_fastqs
######################################
import subprocess
from snakemake.shell import shell
shell.executable("/bin/bash")
log_filename = str(snakemake.log)

f = open(log_filename, 'a+')
f.write("\n##\n## RULE: merge_fastqs \n##\n")
f.close()

#merge input fastq files
command = " cat " + " ".join(sorted(snakemake.input.fastq)) + " > " + snakemake.output.fastq
f = open(log_filename, 'at')
f.write("## COMMAND: "+command+"\n")
f.close()
shell(command)

#run fastqc on the merged data

version = str(subprocess.Popen("fastqc --version 2>&1", shell=True, stdout=subprocess.PIPE).communicate()[0], 'utf-8')
f = open(log_filename, 'at')
f.write("## VERSION: "+version+"\n")
f.close()

command = "mkdir -p `dirname "+snakemake.output.html+"`"+" >> "+log_filename+" 2>&1"
f = open(log_filename, 'at')
f.write("## COMMAND: "+command+"\n")
f.close()
shell(command)

command = "fastqc -o "+snakemake.params.prefix+" "+snakemake.params.extra+" --threads "+str(snakemake.threads)+" "+snakemake.output.fastq+" >> "+log_filename+" 2>&1 "
f = open(log_filename, 'at')
f.write("## COMMAND: "+command+"\n")
f.close()
shell(command)

command = ' sed -r "s:<h2>[ ]*Summary[ ]*<\/h2><ul>:&<li><b>Return to <a href=\'../'+snakemake.params.lib_name+'.final_report.html\'>start page<\/a><\/b><\/li>:" '+snakemake.params.html+' > '+snakemake.output.html
f = open(log_filename, 'at')
f.write("## COMMAND: "+command+"\n")
f.close()
shell(command)

command = "rm -f "+snakemake.params.html+" >> "+log_filename+" 2>&1"
f = open(log_filename, 'at')
f.write("## COMMAND: "+command+"\n")
f.close()
shell(command)

# OLD STUFF:
# run:
#     shell(" {FASTQC} -o {params.prefix} {params.extra} --threads {threads} {input.reads} > {log.run} 2>&1 ")
#     shell(' sed -r "s:<h2>[ ]*Summary[ ]*<\/h2><ul>:&<li><b>Return to <a href=\'../{wildcards.run_name}.final_report.html\'>start page<\/a><\/b><\/li>:" {params.html} > {output.html} ')
#     # shell(" mv -T {params.html} {output.html} ")
