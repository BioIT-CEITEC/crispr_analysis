######################################
# wrapper for rule: preprocess_SE
######################################
import subprocess
from snakemake.shell import shell
shell.executable("/bin/bash")
log_filename = str(snakemake.log)

f = open(log_filename, 'a+')
f.write("\n##\n## RULE: preprocess_PE \n##\n")
f.close()

version = str(subprocess.Popen("conda list 2>&1 | grep cutadapt", shell=True, stdout=subprocess.PIPE).communicate()[0], 'utf-8')
f = open(log_filename, 'at')
f.write("## VERSION: cutadapt "+version+"\n")
f.close()

command = "mkdir -p $(dirname " + snakemake.output.R1 + ") >> " + log_filename + " 2>&1"
f = open(log_filename, 'at')
f.write("## COMMAND: "+command+"\n")
f.close()
shell(command)


# cutadapt -g $ADAPTER3_SEQ1_S1 --times $NUM_ADAPT_TO_REMOVE --cores=${THREADS} --trimmed-only --trim-n -e $ERROR_RATE -O $MIN_OVERLAP --minimum-length $DISC_SHORT -o ${sample%.fastq.*}_trimmed.fastq ${sample%.*} &> ${sample%.fastq.*}.cutadapt.out # -u/-U cut from both reads of a pair; -u $CUT -U $CUT
command = "(time unpigz -c " + snakemake.input.R1 + \
          " | cutadapt " + \
          " -g " + snakemake.params.adapter + \
          " --times " +str(snakemake.params.times) + \
          " --trimmed-only" + \
          " --cores " + str(snakemake.threads) + \
          " --trim-n " + \
          " -e " + str(snakemake.params.error_rate) + \
          " -O " + str(snakemake.params.min_overlap) + \
          " --minimum-length " + str(snakemake.params.min_len) + \
          " -o " + snakemake.output.R1 + \
          " - " + \
          ") >> " + log_filename + " 2>&1 "
f = open(log_filename, 'at')
f.write("## COMMAND: "+command+"\n")
f.close()
shell(command)

# awk 'NR == 2 || NR % 4 == 2' ${sample%.fastq.*}_trimmed.fastq > ${sample%.fastq.*}_trimmed.seqs
command = "(time zcat " + snakemake.output.R1 + " | awk 'NR == 2 || NR % 4 == 2' > " + snakemake.output.seqs + ") 2>> " + log_filename
f = open(log_filename, 'at')
f.write("## COMMAND: "+command+"\n")
f.close()
shell(command)

# cut -c -${GUIDE_LEN} ${sample%.fastq.*}_trimmed.seqs > ${sample%.fastq.*}_trimmed_20bpInserts.seqs
command = "(time cut -c 1-" + str(snakemake.params.guide_len) + " " + snakemake.output.seqs + " > " + snakemake.output.inserts + ") 2>> " + log_filename
f = open(log_filename, 'at')
f.write("## COMMAND: "+command+"\n")
f.close()
shell(command)

command = "mkdir -p " + snakemake.params.prefix + " >> " + log_filename + " 2>&1"
f = open(log_filename, 'at')
f.write("## COMMAND: "+command+"\n")
f.close()
shell(command)

extra = "--noextract --format fastq --nogroup"
command = "(time fastqc -o "+snakemake.params.prefix+" "+extra+" --threads "+str(snakemake.threads)+" "+snakemake.output.R1+") >> "+log_filename+" 2>&1"
f = open(log_filename, 'at')
f.write("## COMMAND: "+command+"\n")
f.close()
shell(command)

command = ' sed -r "s:<h2>[ ]*Summary[ ]*<\/h2><ul>:&<li><b>Return to <a href=\'../'+snakemake.params.lib_name+'.final_report.html\'>start page<\/a><\/b><\/li>:" '+snakemake.params.qc_html_R1_tmp+' > '+snakemake.params.qc_html_R1
f = open(log_filename, 'at')
f.write("## COMMAND: "+command+"\n")
f.close()
shell(command)

command = "rm -f "+snakemake.params.qc_html_R1_tmp+" >> "+log_filename+" 2>&1"
f = open(log_filename, 'at')
f.write("## COMMAND: "+command+"\n")
f.close()
shell(command)

# command = "rm -f "+snakemake.params.qc_html_R1_tmp.replace(".html", ".zip")+" >> "+log_filename+" 2>&1"
# f = open(log_filename, 'at')
# f.write("## COMMAND: "+command+"\n")
# f.close()
# shell(command)

# command = "mv "+snakemake.params.R1_tmp + " " + snakemake.output.R1
# f = open(log_filename, 'at')
# f.write("## COMMAND: "+command+"\n")
# f.close()
# shell(command)

