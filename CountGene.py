#!/usr/local/bin/python3

#script to count the number of occurences of genes in a fasta file generated from polyester
#Julin Maloof
#Dec 6, 2015

import sys, re, gzip

from pathlib import PurePath

from collections import Counter

if len(sys.argv)  != 2 :
    print("This script takes a fasta or fasta.gz file generated from polyester and counts gene occurences")
    print("It currently is hard-wired for Slyc and Spen, looking for the following pattern")
    print("'Sopen.+Solyc[0-9]{2}g[0-9\\.]+'")
    print("Counts are output to standard out")
    exit()

count = Counter()


if PurePath(sys.argv[1]).suffix == '.gz' :
    f = gzip.open(sys.argv[1],'rt')
else :
    f = open(sys.argv[1],'r')

for line in f:
    gene = re.search('Sopen.+Solyc[0-9]{2}g[0-9\\.]+',line)
    if gene:
        count[gene.group(0)] += 1

for gene, count in count.items():
    print(gene,count)
