#!/usr/bin/env python3

import os, gzip, itertools
import fastaparser as fp
import array as arr
from Bio.SeqUtils import GC
from Bio import  SeqIO

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

with gzip.open(file1,"rt") as fasta_file:
    parser = fp.Reader(fasta_file)
    gene_num = 1
    for seq in parser:
        genes1 = (str(seq.id)).count('L')
        gene_num += genes1
    print("The total number of genes in species 1, Salmonella, is", gene_num)

with gzip.open(file1,"rt") as f:
    seqs = dict(aspairs(f))
    genesum1=0
    for gl in seqs.values():
        tgl = len(gl)
        genesum1 += tgl
    for x in seqs.values():
        gc1 = GC(x)
    print("The total length of gene sequence for species 1, Salmonella, is {} bp ({} kb)".format(genesum1, genesum1/1000), ".")
    print ("The GC content of species 1, Salmonella, is {} %" .format(gc1))

    
with gzip.open(file2,"rt") as fasta_file:
    parser = fp.Reader(fasta_file)
    gene_num2 = 0
    for seq in parser:
        genes2 = (str(seq.id)).count('P')
        gene_num2 += genes2
    print("The total number of genes in species 2, Mycobacterium, is", gene_num2)
       
with gzip.open(file2,"rt") as f:
    seqs = dict(aspairs(f))
    genesum2=0
    for gl in seqs.values():
        tgl2 = len(gl)
        genesum2 += tgl2
    for x in seqs.values():
        gc2 = GC(x)
    print("The total length of gene sequence for species 2, Mycobacterium, is {} bp ({} kb)".format(genesum2, genesum2/1000), ".")
    print ("The GC content of species 2, Mycobacterium, is {} %" .format(gc2))

