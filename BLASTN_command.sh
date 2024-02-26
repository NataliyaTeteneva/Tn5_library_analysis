#!/bin/sh

# create folder to place the alignment files and descend to the folder storing the sequences
mkdir alis
cd enrich

# list the fasta files in the folder
files="*"
echo $files

for file in $files
do
	# align each file to the W3110 reference genome (downloaded from the GenBank, entry number NC_007779.1)
	# for each read, keep only the best hit
	# generate output in XML format
	nm=${file%.fa}
	echo $nm
	blastn -subject ../W3110_reference_genome.fa -outfmt 5 -max_hsps 1 -query $nm.fa -out ../alis/$nm.xml
	echo $nm "is done"
done


echo "ALL IS DONE"

