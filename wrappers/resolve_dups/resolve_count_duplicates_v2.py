#!/usr/bin/python3
import sys
from collections import defaultdict
def read_file(count_file):
    gene_dict = defaultdict(list)
    duplicates_line_list = []
    output_lines = []
    with open(count_file) as c_file:                  
        output_lines.append(c_file.readline())
        #go over all lines
        for line in c_file:               
            data = line.strip().split("\t")
            #check whether the line contains duplicated record
            if not ";" in data[0]:
                #if line has no duplicatated record, add gene as a key to dictionary and read count as a value to dictionary
                if not data[0] in gene_dict.keys():
                    gene_dict[data[0]] = [int(data[3])]
                #if gene is already in dictionary keys, append a read count of another sgRNA for the same gene to values of dictionary  
                else:
                    gene_dict[data[0]].append(int(data[3]))
                # print the lines without duplicates directly to output
                output_lines.append(line)
            # if line contains duplicates appends it to the list of lines to resolve
            else:
                duplicates_line_list.append(line)
    return gene_dict,duplicates_line_list,output_lines

def resolve_duplicates(gene_dict,duplicates_line_list):
    new_lines = []
    random_genes = []
    #just getting records of genes which has no unique guide 
    for line in duplicates_line_list:
        data = line.strip().split("\t")
        genes = data[0].split(";")
        if len(set(genes)) != 1:
            # getting genes in line and chcek whether they are in gene-counts dictionary or not
            for gene in genes:
                if gene not in gene_dict.keys():
                    # if they are not in dictioary, it means, that they have no unique sgRNA record, and their counts are based only on the logic of algortihm under
                    if not gene in random_genes:
                        random_genes.append(gene)
    for line in duplicates_line_list:        
            gene_coverage = {}
            data = line.strip().split("\t")
            genes = data[0].split(";")
            sgRNAs = data[1].split(";")
            counts_number = int(data[3])
            # if all sgRNAs in duplicated record belongs to the same gene, their counts are divided equally
            if len(set(genes)) == 1:
                for i,gene in enumerate(genes):
                    new_lines.append(gene + "\t" + sgRNAs[i] + "\t" + data[2] + "\t" + str(int(counts_number/len(genes))) + "\n")
            # if sgRNAs in duplicater record does not belong to the same gene
            else: 
                for gene in genes:
                    #if gene is not in gene-counts dictionary(reads will be completely random according algorithm)
                    if gene not in gene_dict.keys():
                        #then the gene coverage is set to zero
                        gene_coverage[gene] = 0
                    else:
                        #otherwise gene coverage is calculated as average of other sgRNAs for particular gene, which were unique 
                        gene_coverage[gene] = sum(gene_dict[gene])/len(gene_dict[gene])
                # get the information if sum of all coverages for each gene is zero or not 
                cov_sum = sum(gene_coverage.values())
                #if sum of all coverages for each gene is not zero
                if cov_sum != 0:
                    #testing whether there is some zero coverage (from random genes)
                    # if not
                    if not 0 in gene_coverage.values():
                        #convert gene coverages into ratios between 0 and 1
                        for gene,coverage in gene_coverage.items():            
                            gene_coverage[gene] = float((coverage/cov_sum))
                        for i,gene in enumerate(genes): 
                            # append new line with unique record, containing gene,sgRNA,sequence and count derived as original count for this row multiplied by gene coverage ratio
                            new_lines.append(gene + "\t" + sgRNAs[i] + "\t" + data[2] + "\t" + str(int(counts_number*float(gene_coverage[gene]))) + "\n")
                    # if some zero in gene_coverage (comes from random genes)
                    else:
                        # to avoid to get negative read counts in the end
                        zero_num_list = []
                        #get number of zeros in coverage
                        [zero_num_list.append(value) for value in gene_coverage.values() if value == 0]
                        # calculate ratio of zeros in record according to number of all duplicated sgRNAs in this record
                        percentage_of_zeros = len(zero_num_list)/len(gene_coverage)
                        for gene,coverage in gene_coverage.items():
                            # if coverage is zero, gene coerage is assigned based on the ratio of zeros and their count
                            if coverage == 0:
                                gene_coverage[gene] = percentage_of_zeros/len(zero_num_list)
                            else:
                                # if coverage is not zero, coverage ratio is given as above on line 71, but restricted only to the ratio left after removal of parts taken by zeros 
                                gene_coverage[gene] = (1-percentage_of_zeros)*coverage/cov_sum
                        for i,gene in enumerate(genes): 
                            # append new line with unique record, containing gene,sgRNA,sequence and count derived as original count for this row multiplied by gene coverage ratio
                            new_lines.append(gene + "\t" + sgRNAs[i] + "\t" + data[2] + "\t" + str(int(counts_number*float(gene_coverage[gene]))) + "\n")                                           
                # if sum of all coverages for each gene is zero                                                        
                else:                
                    for i,gene in enumerate(genes):
                        # it means that all genes are random and divided equally
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