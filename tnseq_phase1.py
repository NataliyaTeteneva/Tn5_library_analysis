#!/usr/bin/env python
# coding: utf-8


# necessary libraries
import pandas as pd
import statistics, re, os, glob, gzip
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import reverse_complement

# sequences for the script to look for

U1 = "GATGTCCACGAGGTCTCT"
U2 = "CGTACGCTGCAGGTCGAC"
Tn5 = "CTGTCTCTTATACACATCT"

# read the sequencing files placed in the same folder as the script
files = glob.glob('*.gz')
print(files)
# crop the file extensions to generate files with the same name but different formats
names = [os.path.splitext(os.path.splitext(file)[0])[0] for file in files]

# prepare folders for the helper files
os.mkdir("enrich")
os.mkdir("barcodes")

# auxilary function to count the unique barcodes
def update_dict(key, d):
    # key=key (barcode in this case)
    # d=dictionary (where keys are barcodes and values are counts)
    # if the barcode is not in the dictionary, adds it to the dictionary
    # if the barcode is in the dictionary, updates the counter
    d.setdefault(key, 0)
    d[key] += 1

# function to find a requested DNA site
def find_priming_site(read, site):
    # read=record (Biopython SeqRecord type object)
    # site=site (U1, U2 or Tn5)
    # returns the position of the discovered DNA site on the read
    
    # looking for exact match
    pos_start = read.seq.find(site)
    pos_end = pos_start+len(site)
    # if search failed try to find the reverse complement
    if pos_start==-1:
        pos_start = read.seq.find(reverse_complement(site))
        pos_end = pos_start+len(site)
        # and only then return "no found"
        if (pos_start==-1):
            return "no_found"
    position = [pos_start, pos_end]
    return position

# function to find a barcode and check its quality
def look_for_target(read):
    # read=record (Biopython SeqRecord type object)
    # finds a barcode, stores the information about it
    # and transforms it to one-entry pandas DataFrame

    # Four possible outputs of this function:
    # 1. Valid barcode (then it is stored as one-entry pandas DataFrame
    # 2. Barcode found but failed the quality check (str "bad_barcode")
    # 3. Barcode found but the remainder of the genomic DNA
    #    was too short to reliably align (str "too_short")
    # 4. No barcode found (str "no_barcode")

    # Information stored in the DataFrame:
    # seq: barcode sequence, str
    # begin: beginning of the barcode, int
    # end: end of the barcode, int
    # direction: direction of the barcode relative to the genome, 1/-1
    # id: read ID from the FASTQ file, str

    # define positions of the flanking sites (both should be present in the read)
    pos_U1 = find_priming_site(read, U1)
    pos_U2 = find_priming_site(read, U2)

    # if they are not found, try to match the shorter versions of them
    if (pos_U1=="no_found")|(pos_U2=="no_found"):
        pos_U1 = find_priming_site(read, U1[-5:])
        pos_U2 = find_priming_site(read, U2[0:5])
        # if they are still not found, there is no barcode then
        if (pos_U1=="no_found")|(pos_U2=="no_found"):
            return "no_barcode"

    # check whether the Tn5 flanking site is nearby (part of the QC)
    pos_IR = find_priming_site(read, Tn5)
    if pos_IR=="no_found":
        pos_IR = find_priming_site(read, Tn5[-5:])
        if pos_IR=="no_found":
            return("bad_barcode")

    # check the length of the remaining genomic DNA
    if pos_U2[0] < pos_IR[0]:
        remainder = len(read.seq)-pos_IR[1]
    else:
        remainder = pos_IR[0]
    if remainder < 15:
        return "too_short"

    # define direction of the barcode relative to the genome    
    pos_target = pos_U1 + pos_U2
    if pos_target[0]<pos_target[2]:
        direction=1
    else: direction=-1

    # use the position of the flanking sites U1 and U2 to define the position of the barcode
    # which is in between them
    s_pos_target = sorted(pos_target)
    
    # define the length of the barcode
    barcode_length = s_pos_target[2]-s_pos_target[1]
    
    # use fastq quality string to check the sequencing quality (>30 for each letter)
    barcode_qual = read.letter_annotations['phred_quality'][s_pos_target[1]:s_pos_target[2]]
    qual_check = all(x > 30 for x in barcode_qual)

    # if the barcode successfully passed the QC, turn it into the pandas DataFrame
    if (qual_check)&(barcode_length==20):
        bc = {'seq': str(read.seq[s_pos_target[1]:s_pos_target[2]]), 'begin': s_pos_target[1],
              'end': s_pos_target[2], 'direction': direction, 'id': read.id}
        update_dict(bc['seq'], barcodes)
        bc = pd.DataFrame.from_dict(bc, orient='index').transpose().set_index('id')
        return bc
    else:
        return "bad_barcode"

# function that reads a single fastq file

def analyze_file(infile, name):
    # infile=fastq file
    # name=just name of this file
    
    # returns 2 files:
    # 1. fasta file with all the reads that contain barcodes (stored in enrich/ folder)
    # 2. csv file with all the barcode information (stored in barcodes/ folder)

    # count of valid barcodes
    c = 0
    # count of empty reads
    n = 0
    # count of bad barcodes
    b = 0

    # list to put the sequences to the new fasta file
    outfile = []
    # list to put the barcode data to the new csv file
    bc_data = pd.DataFrame()

    # open fastq file
    with gzip.open(infile, "rt") as handle:
        print("File {} is in process".format(name))
        for record in SeqIO.parse(handle, "fastq"):
            # look for barcode in every read
            bc = look_for_target(record)

            # if there is no barcode, just record it
            if isinstance(bc, pd.DataFrame)==False:
                if bc=="no_barcode":
                    n+=1
                elif (bc=="bad_barcode")|(bc=="too_short"):
                    b+=1
            # if there is a barcode, add the corresponing read to the fasta file
            # and the barcode itself to the csv file
            elif isinstance(bc, pd.DataFrame)==True:
                outfile.append(record)
                bc_data = pd.concat([bc_data, bc])
                c+=1
                # print the length of the dataset after every 1000th barcode
                if c%1000 == 0:
                    print(len(bc_data))
                # dump the sequences to the fasta file after they reach certain size
                if len(outfile)>300000:
                    with open("enrich/"+name+"_enriched.fa", "a") as output_handle:
                        SeqIO.write(outfile, output_handle, "fasta")
                        outfile=[]
                    with open("barcodes/"+name+"_barcodesData.csv", "a") as bc_handle:
                        bc_data = bc_data.drop_duplicates()
                        bc_data.to_csv(bc_handle, header=bc_handle.tell()==0)
                        bc_data=[]
        print("Empty reads: ", n)
        print("Invalid barcodes: ", b)
        print("Valid barcodes: ", c)
        
        with open("enrich/"+name+"_enriched.fa", "a") as output_handle:
            SeqIO.write(outfile, output_handle, "fasta")

        # drop duplicates and record the barcode data as csv file
        with open("barcodes/"+name+"_barcodesData.csv", "a") as bc_handle:
            bc_data = bc_data.drop_duplicates()
            bc_data.to_csv(bc_handle, header=bc_handle.tell()==0)

# create a dictionary to collect unique barcodes
barcodes = {}

# run a for loop for all existing fastq files in the folder
for file, name in zip(files, names):
    analyze_file(file, name)
print("Unique barcodes: ", len(barcodes))

