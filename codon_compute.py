#!/usr/bin/env python3

#!/usr/bin/env python3

import os, gzip, itertools
import numpy as np 
import pandas as pd
from ast import literal_eval
import array as arr
from Bio.SeqUtils import GC
from Bio import  SeqIO
import warnings
from collections import defaultdict


warnings.filterwarnings("ignore", category=DeprecationWarning)

# this is code which will parse FASTA files
# define what a header looks like in FASTA format
def isheader(line):
    return line[0] == '>'

def aspairs(f):
    seq_id = ''
    sequence = ''
    for header,group in itertools.groupby(f, isheader):
        if header:
            line = next(group)
            seq_id = line[1:].split()[0]
        else:
            sequence = ''.join(line.strip() for line in group)
            yield seq_id, sequence


url1="ftp://ftp.ensemblgenomes.org/pub/bacteria/release-45/fasta/bacteria_0_collection/salmonella_enterica_subsp_enterica_serovar_typhimurium_str_lt2/cds/Salmonella_enterica_subsp_enterica_serovar_typhimurium_str_lt2.ASM694v2.cds.all.fa.gz"
url2="ftp://ftp.ensemblgenomes.org/pub/bacteria/release-45/fasta/bacteria_0_collection/mycobacterium_tuberculosis_h37rv/cds/Mycobacterium_tuberculosis_h37rv.ASM19595v2.cds.all.fa.gz"
file1="Salmonella_enterica_subsp_enterica_serovar_typhimurium_str_lt2.ASM694v2.cds.all.fa.gz"
file2="Mycobacterium_tuberculosis_h37rv.ASM19595v2.cds.all.fa.gz"

if not os.path.exists(file1):
    os.system("curl -O %s"%(url1))

if not os.path.exists(file2):
    os.system("curl -O %s"%(url2)) 

# 2 Total length of these gene sequences for each file

with gzip.open(file1,"rt") as f:
    seqs = dict(aspairs(f))
    genesum1=0
    codons1 = {}
    total_codons1 = 0
    for gl in seqs.values():
        tgl = len(gl)
        genesum1 += tgl
    for x in seqs.values():
        gc1 = GC(x)
    for A in {'A','T','C','G'}:
        for B in {'A','T','C','G'}:
            for C in {'A','T','C','G'}:
                new_codon1 = A+B+C
                codons1[new_codon1] = 0
    for seq in seqs:
        for n in range(0, len(seqs[seq]), 3):
            codons1[seqs[seq][n:n+3]] += 1
    for value in codons1:
        total_codons1 += codons1[value]
    for nk, nv in codons1.items():
        codons1[nk] = ((nv/total_codons1)*100)
    
    print("The total number of genes in species 1, Salmonella, is", len(seqs))
    print("The total length of gene sequence for species 1, Salmonella, is {} bp ({} kb)".format(genesum1, genesum1/1000), ".")
    print ("The GC content of species 1, Salmonella, is {} %" .format(gc1))
    print ("The total number of codons in species 1, Salmonella, is", total_codons1)
    

with gzip.open(file2,"rt") as f:
    seqs = dict(aspairs(f))
    genesum2=0
    codons2 = {}
    total_codons2 = 0
    for gl in seqs.values():
        tgl2 = len(gl)
        genesum2 += tgl2
    for x in seqs.values():
        gc2 = GC(x)
    for A in {'A','T','C','G'}:
        for B in {'A','T','C','G'}:
            for C in {'A','T','C','G'}:
                new_codon2 = A+B+C
                codons2[new_codon2] = 0
    for seq in seqs:
        for n in range(0, len(seqs[seq]), 3):
            codons2[seqs[seq][n:n+3]] += 1
    for value in codons2:
        total_codons2 += codons2[value]
    for nk, nv in codons2.items():
        codons2[nk] = ((nv/total_codons2)*100)
        
        
    
    # Make table (new dictionary sorted by codon)
    ds = [codons1, codons2]
    codons2_sorted = {i:codons2[i] for i in codons1.keys()}
    keys = codons1.keys()
    values = zip(codons1.values(), codons2_sorted.values())
    new_dict = dict(zip(keys, values))
    df = pd.DataFrame(new_dict).T.reset_index()
    df.columns = ['CODONS','FREQ. Sp1','FREQ Sp2']
    
    print("The total number of genes in species 2, Mycobacterium, is", len(seqs))
    print("The total length of gene sequence for species 2, Mycobacterium, is {} bp ({} kb)".format(genesum2, genesum2/1000), ".")
    print ("The GC content of species 2, Mycobacterium, is {} %" .format(gc2))
    print ("The total number of codons in species 2, Mycobacterium, is", total_codons2)
    
    
    print (df)
    
