######################################
# wrapper for rule: hs_blast_alignment
######################################
import subprocess
from collections import Counter, defaultdict
from snakemake.shell import shell
shell.executable("/bin/bash")
log_filename = str(snakemake.log)

f = open(log_filename, 'a+')
f.write("\n##\n## RULE: hs_blast_alignment \n##\n")
f.close()

version = str(subprocess.Popen("conda list 2>&1 | grep 'hs-blastn'", shell=True, stdout=subprocess.PIPE).communicate()[0], 'utf-8')
f = open(log_filename, 'at')
f.write("## VERSION: "+version+"\n")
f.close()

f = open(log_filename, 'at')
f.write(f"## INFO: counting identical reads for reduction in {snakemake.input.ins}\n")
f.close()

read_file = snakemake.input.ins
reads_as_fasta = open(snakemake.output.fa, "w")
read_counts = open(snakemake.output.counts, "w")
read_list = []
with open(read_file,'r') as reads:  
    for i,line in enumerate(reads):
        read_list.append(line.strip())
    counter = Counter(read_list)
    i=1
    for key,value in counter.most_common():
        read_counts.write(key + "\t" + str(value)+"\n")
        reads_as_fasta.write(">"+str(i)+"\n")
        reads_as_fasta.write(key + "\n")
        i += 1
reads_as_fasta.close()
read_counts.close()

command = f"(time hs-blastn align -db {snakemake.input.idx} -query {snakemake.output.fa} -reward 1 -penalty -3 -gapopen 0 -gapextend 2 -word_size 12 -max_target_seqs 1 -num_threads {str(snakemake.threads)} -out {snakemake.output.out} -outfmt 6) >> {log_filename} 2>&1"
f = open(log_filename, 'at')
f.write("## COMMAND: "+command+"\n")
f.close()
shell(command)
