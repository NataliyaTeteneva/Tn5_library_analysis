#!/usr/bin/env python
# coding: utf-8


# necessary libraries
import pandas as pd
import re, os, glob
from Bio.Blast import NCBIXML

# gene positions corresponding to the W3110 reference genome (downloaded from the GenBank, entry number NC_007779.1)
W3110 = pd.read_csv("gene_positions_W3110_NC.csv")

# list and sort the files containing the 
alis = sorted(glob.glob('alis/*.xml'))
barcodes = sorted(glob.glob('barcodes/*_barcodesData.csv'))

# reads the BLASTN XML file and extracts the alignment information
def read_ali_data(file):
    #file=XML alignment file
    handle = open(file)
    records = NCBIXML.parse(handle)
    # create emtpy lists to collect read ID, length of alignment,
    # position in genome and read respectively, and the strand(sense va antisense)
    id_read = []
    ali_length = []
    genome = []
    read = []
    strand = []
    for record in records:
        for alignment in record.alignments:
            if alignment.title == 'NC_007779.1 NC_007779.1 Escherichia coli str. K-12 substr. W3110, complete sequence':
                for hsp in alignment.hsps:
                    # QC for alignment
                    if hsp.expect<0.05:
                        # collect the data
                        ali_length.append(hsp.align_length)
                        genome.append(hsp.sbjct_start)
                        read.append(hsp.query_start)
                        if hsp.strand[0]==hsp.strand[1]:
                            strand.append(1)
                        else: strand.append(-1)
                        id_read.append(re.split('\s+', record.query)[0])
    # combine all collected data into a single DataFrame
    ali = pd.DataFrame({'id_read': id_read,
                        "ali_length": ali_length,
                       "genome": genome,
                       "read": read,
                       "strand": strand}).set_index("id_read")
    return ali

# takes the alignment results (extracted from the XML file) and the barcode data (generated in the phase 1)
def analyze_files(blast, barcodes):
    # blast=alignment results generated with read_ali_data function
    # barcodes=barcode data read from the csv files

    # read the input data
    ali = read_ali_data(blast)
    barcode = pd.read_csv(barcodes).rename(columns={"id": "id_read"}).set_index("id_read")

    # inner merge using the read ID as a key
    total = pd.merge(ali, barcode, on = "id_read")
    total = total.drop_duplicates()

    # determine the position of each barcode on the genome
    # using the positions of alignment and strand data
    total.loc[(total['strand']==1) & (total['direction']==1), 'position'] = total.genome+total.begin-total.read
    total.loc[(total['strand']==1) & (total['direction']==-1), 'position'] = total.genome+total.end-total.read
    total.loc[(total['strand']==-1) & (total['direction']==1), 'position'] = total.genome+total.ali_length+total.read-total.begin
    total.loc[(total['strand']==-1) & (total['direction']==-1), 'position'] = total.genome+total.ali_length+total.read-total.end                                                                    
    # update barcode direction relative to the gene, not to the genome
    total['direction_fin'] = total.direction*total.strand

    # select only relevant columns (barcode sequence, genome position and direction)
    total = total[['seq', 'direction_fin', 'position']].reset_index(drop = True)
    total.rename(columns={"direction": "direction_fin"})
    # return DataFrame
    return total


# using the determined barcode position and the reference genome, determine the gene that was disrupted
def find_the_gene(position):
    # position = barcode position

    # to accelerate computations, already assigned positions are stored in the dictionary
    if position in genepositions:
        return genepositions[position]
    else:
        # annotate the positions close to the edges of the genome
        if (position < 190) | (position > 4646308):
            gene = "yjtD_thrL"
        else:
            # find a gene in the reference genome where desired position is between start and end
            subset_gene = W3110.query('start <= @position & end >= @position')
            # if there is no such gene then the barcode is between the genes
            if subset_gene.empty:
                # gene downstream of insertion
                down = W3110.query('start >= @position')["gene"].iloc[0]
                # gene upstream of insertion
                up = W3110.query('start <= @position')["gene"].iloc[-1]
                # insertion is annotated as intergenic 
                gene = up+'_'+down
            else:
                # insertion is annotated as within the gene
                gene = subset_gene["gene"].iloc[0]
        # update the dictionary for the future use
        genepositions.update({position: gene})
        return gene

# results from different files are temporarily stored as DataFrames in a list
results = []
print("Analyzing files..")
# run the function collecting all the data
for ali, bc in zip(alis, barcodes):
    r = analyze_files(ali, bc)
    results.append(r)
    print("Alignment", ali,  "is done", sep = " ")


print("Determining gene positions...")

# collect all the DataFrames into a single large DataFrame, drop duplicates
total = pd.concat(results).sort_values("position").reset_index(drop = True)
total = total.drop_duplicates().reset_index()

# create a helper dictionary to store the determined gene positions
genepositions = {}

# determine the genes disrupted by each barcode
total['gene'] = total['position'].apply(lambda x: find_the_gene(x))

# write csv
total.to_csv("barcode_positions.csv")

# determine the total number of unique barcode insertions
print(total.shape)



