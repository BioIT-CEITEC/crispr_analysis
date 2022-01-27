#!/usr/bin/python3
import sys
from collections import defaultdict
def read_file(count_file):
    gene_dict = defaultdict(list)
    duplicates_line_list = []
    output_lines = []
    with open(count_file) as c_file:                  
        output_lines.append(c_file.readline())
        for line in c_file:               
            data = line.strip().split("\t")
            if not ";" in data[0]:
                if not data[0] in gene_dict.keys():
                    gene_dict[data[0]] = [int(data[3])]
                else:
                    gene_dict[data[0]].append(int(data[3]))
                output_lines.append(line)
            else:
                duplicates_line_list.append(line)
    return gene_dict,duplicates_line_list,output_lines

def resolve_duplicates(gene_dict,duplicates_line_list):
    new_lines = []
    random_genes = []
    for line in duplicates_line_list:
        data = line.strip().split("\t")
        genes = data[0].split(";")
        if len(set(genes)) != 1:
            for gene in genes:
                if gene not in gene_dict.keys():
                    if not gene in random_genes:
                        random_genes.append(gene)
    for line in duplicates_line_list:        
            gene_coverage = {}
            data = line.strip().split("\t")
            genes = data[0].split(";")
            sgRNAs = data[1].split(";")
            counts_number = int(data[3])
            if len(set(genes)) == 1:
                for i,gene in enumerate(genes):
                    new_lines.append(gene + "\t" + sgRNAs[i] + "\t" + data[2] + "\t" + str(int(counts_number/len(genes))) + "\n")
            else: 
                for gene in genes:
                    if gene not in gene_dict.keys():
                        gene_coverage[gene] = 0
                    else:
                        gene_coverage[gene] = sum(gene_dict[gene])/len(gene_dict[gene])
                cov_sum = sum(gene_coverage.values())
                if cov_sum != 0:
                    if not 0 in gene_coverage.values():
                        for gene,coverage in gene_coverage.items():            
                            gene_coverage[gene] = float((coverage/cov_sum))
                        for i,gene in enumerate(genes): 
                            new_lines.append(gene + "\t" + sgRNAs[i] + "\t" + data[2] + "\t" + str(int(counts_number*float(gene_coverage[gene]))) + "\n")
                            if len(gene_dict[gene]) <=6 :
                                gene_dict[gene].append(int(counts_number*float(gene_coverage[gene])))
                    else:
                        if cov_sum > counts_number:
                            zero_num_list = []
                            [zero_num_list.append(value) for value in gene_coverage.values() if value == 0]
                            percentage_of_zeros = len(zero_num_list)/len(gene_coverage)
                            for gene,coverage in gene_coverage.items():
                                if coverage == 0:
                                    gene_coverage[gene] = percentage_of_zeros
                                else:
                                    gene_coverage[gene] = 0.75*coverage/cov_sum
                            for i,gene in enumerate(genes): 
                                new_lines.append(gene + "\t" + sgRNAs[i] + "\t" + data[2] + "\t" + str(int(counts_number*float(gene_coverage[gene]))) + "\n")
                                if len(gene_dict[gene]) <=6 :
                                    gene_dict[gene].append(int(counts_number*float(gene_coverage[gene])))
                        else:
                            non_zero_genes = []
                            non_zero_sgRNAs = []
                            zero_genes = []
                            zero_sgRNAs = []
                            for i,gene in enumerate(genes):
                                if gene_coverage[gene] != 0:
                                    non_zero_genes.append(genes[i])
                                    non_zero_sgRNAs.append(sgRNAs[i])
                                else:
                                    zero_genes.append(genes[i])
                                    zero_sgRNAs.append(sgRNAs[i])
                            for i,gene in enumerate(non_zero_genes):
                                counts_number -= int(gene_coverage[gene])
                                new_lines.append(gene + "\t" + non_zero_sgRNAs[i] + "\t" + data[2] + "\t" + str(int(gene_coverage[gene])) + "\n")
                                if len(gene_dict[gene]) <=6 :
                                    gene_dict[gene].append(int(gene_coverage[gene]))                                
                            for i,gene in enumerate(zero_genes):
                                new_lines.append(gene + "\t" + zero_sgRNAs[i] + "\t" + data[2] + "\t" + str(int(int(counts_number)/len(zero_genes))) + "\n")
                                if not gene in gene_dict.keys():
                                    gene_dict[gene] = [int(int(counts_number)/len(zero_genes))]
                                else:
                                    gene_dict[gene].append(int(int(counts_number)/len(zero_genes)))                        
                                                                        
                else:                
                    for i,gene in enumerate(genes):
                        new_lines.append(gene + "\t" + sgRNAs[i] + "\t" + data[2] + "\t" + str(int(counts_number/len(genes))) + "\n")
    return(new_lines,random_genes)
   
def print_table(output_lines,new_lines,count_file,random_genes):
    with open(count_file.split(".")[0]+"_resolved_duplicates.tsv","w") as output_file:
        with open(count_file.split(".")[0]+"_zero_counts_after_deduplication.tsv","w") as zero_file:
            for line in output_lines:
                output_file.write(line)
            for line in new_lines:
                if int(line.strip().split("\t")[3]) == 0:
                    zero_file.write(line)
                else:    
                    output_file.write(line)
    with open(count_file.split(".")[0]+"_random_genes_list.txt","w") as random_genes_out:
        for line in random_genes:
            random_genes_out.write(line + "\n")
           
def main():
    count_file =  sys.argv[1]
    gene_dict,duplicates_line_list,output_lines = read_file(count_file)
    new_lines,random_genes = resolve_duplicates(gene_dict,duplicates_line_list)
    print_table(output_lines,new_lines,count_file,random_genes)
    
if __name__ == '__main__':
    main()
