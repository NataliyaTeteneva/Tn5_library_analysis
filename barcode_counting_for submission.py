#!/usr/bin/env python
# coding: utf-8

# necessary libraries
import gzip, glob, time, re
import pandas as pd
import numpy as np
from itertools import compress
from Bio import SeqIO
from Bio.Seq import reverse_complement
from tqdm import trange, tqdm


# LIST THE FILES

# sample information
meta = pd.read_csv("../samples.csv")
# barcodes annotation
anno = pd.read_csv("../anno_appended.csv")
anno = anno[['seq']].set_index('seq')
# sequencing results
files = glob.glob('*.gz')


# In[17]:


# FUNCTIONS

# auxilary function to count the barcodes
def update_dict(key, d):
    # key=key (barcode in this case)
    # d=dictionary (where keys are barcodes and values are counts)
    # if the barcode is not in the dictionary, adds it to the dictionary
    # if the barcode is in the dictionary, updates the counter
    d.setdefault(key, 1)
    d[key] += 1
    
# iterate througt the file to collect the barcodes and count them
def read_barcodes(file, anno, metainfo):
    # create the dictionary to keep the barcodes
    barcodes = {}
    # open a file with a progress bar to visualize
    with gzip.open(file, "rt") as handle:
        for record in SeqIO.parse(handle, "fastq"):
            # extract the barcode sequence based on known position
            r = record[23:43]
            # convert the barcode to the reverse complement
            bc = reverse_complement(r.seq)
            # check if the barcode is in the annotation
            if bc in anno.index:
                # record it
                update_dict(str(bc), barcodes)
                    
    print("Analysis is done, data preprocessing...")

    # if there are any barcodes collected:
    if len(barcodes.keys())>0:
        # convert the dictionary to the dataframe and index it
        bc = pd.DataFrame(list(barcodes.items()), columns = ["seq", "abundance"])
        bc = bc.set_index('seq')
        # merge the resulting barcode counts with the snip of sample information
        aux = pd.concat([metainfo]*bc.shape[0]).reset_index(drop = True)
        bc = pd.concat([bc.reset_index(drop = False), aux],axis=1)
        return bc
    else:
        return None


# actual analysis

# create the progress bar to visualize
for i in trange(len(files), desc='Files'):
    # extract the sample index trom the file name
    index = int(re.match(r'(sample_)([0-9]{2})(\.fastq\.gz)', files[i]).group(2))
    # extract the sample information from the sample description
    metainfo = meta[meta['primer']==index]
    metainfo = metainfo.drop(['primer', 'index'], axis=1)

    # count the barcodes in each file, add the sample information
    output = read_barcodes(files[i], anno, metainfo)
    print(output)
    
    if output is None:
        print("Not here.")
    else:
        # add the results of the computation to the table
        with open("results.csv", 'a') as f:
            output.to_csv(f, header=f.tell()==0)
print("DONE")





