gff = "Escherichia_coli_str_k_12_substr_mg1655.ASM584v2.37.gff3.gz"
fasta = "Escherichia_coli_str_k_12_substr_mg1655.ASM584v2.dna.chromosome.Chromosome.fa.gz"

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

if not os.path.exists(gff):
    os.system("wget -O ftp://ftp.ensemblgenomes.org/pub/bacteria/release-45/gff3/bacteria_0_collection/escherichia_coli_str_k_12_substr_mg1655/Escherichia_coli_str_k_12_substr_mg1655.ASM584v2.37.gff3.gz")

if not os.path.exists(fasta):
    os.system("wget -O ftp://ftp.ensemblgenomes.org/pub/bacteria/release-45/fasta/bacteria_0_collection/escherichia_coli_str_k_12_substr_mg1655/dna/Escherichia_coli_str_k_12_substr_mg1655.ASM584v2.dna.chromosome.Chromosome.fa.gz")

# 2 Count up and print out the number genes (gene feature)
# 3 Compute the total length of the genes (length is the END - START)

with gzip.open(gff,"rt") as fh:
    mygff = csv.reader(fh,delimiter="\t")
    count = 0
    gene_len = 0
    CDS_count = 0
    CDS_len = 0
    for row in mygff:
        if row[0].startswith("#"):
            continue
        if row[2] == "gene":
            count += 1
            total_gene_len =int(row[4])-int(row[3])
            gene_len += total_gene_len
        if row[2] == "CDS":
            CDS_count += 1
            total_CDS_len =int(row[4])-int(row[3])
            CDS_len +=  total_CDS_len
    print ("The number of genes is", count,".")
    print ("The number of CDS is", CDS_count,".")
    print ("Total length of genes is {} bp ({} kb)".format(gene_len,gene_len/1000), ".")
    print ("Total length of CDS is {} bp ({} kb)".format(CDS_len,CDS_len/1000), ".")
    
# 4 Use the FASTA file to compute the total length of genome (by adding up the length of each sequence in the file). 

with gzip.open(fasta,"rt") as fasta_file:
    parser = fp.Reader(fasta_file)
    for seq in parser:
        total = (len(seq.sequence_as_string()))
    print("The total length of the E. coli genome is {} bp ({} kb)".format(total, total/1000), ".")
    print("The percentage of the genome which is coding is {} %".format((CDS_len/total)*100))gff = "Escherichia_coli_str_k_12_substr_mg1655.ASM584v2.37.gff3.gz"
fasta = "Escherichia_coli_str_k_12_substr_mg1655.ASM584v2.dna.chromosome.Chromosome.fa.gz"

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

if not os.path.exists(gff):
    os.system("wget -O ftp://ftp.ensemblgenomes.org/pub/bacteria/release-45/gff3/bacteria_0_collection/escherichia_coli_str_k_12_substr_mg1655/Escherichia_coli_str_k_12_substr_mg1655.ASM584v2.37.gff3.gz")

if not os.path.exists(fasta):
    os.system("wget -O ftp://ftp.ensemblgenomes.org/pub/bacteria/release-45/fasta/bacteria_0_collection/escherichia_coli_str_k_12_substr_mg1655/dna/Escherichia_coli_str_k_12_substr_mg1655.ASM584v2.dna.chromosome.Chromosome.fa.gz")

# 2 Count up and print out the number genes (gene feature)
# 3 Compute the total length of the genes (length is the END - START)

with gzip.open(gff,"rt") as fh:
    mygff = csv.reader(fh,delimiter="\t")
    count = 0
    gene_len = 0
    CDS_count = 0
    CDS_len = 0
    for row in mygff:
        if row[0].startswith("#"):
            continue
        if row[2] == "gene":
            count += 1
            total_gene_len =int(row[4])-int(row[3])
            gene_len += total_gene_len
        if row[2] == "CDS":
            CDS_count += 1
            total_CDS_len =int(row[4])-int(row[3])
            CDS_len +=  total_CDS_len
    print ("The number of genes is", count,".")
    print ("The number of CDS is", CDS_count,".")
    print ("Total length of genes is {} bp ({} kb)".format(gene_len,gene_len/1000), ".")
    print ("Total length of CDS is {} bp ({} kb)".format(CDS_len,CDS_len/1000), ".")
    
# 4 Use the FASTA file to compute the total length of genome (by adding up the length of each sequence in the file). 

with gzip.open(fasta,"rt") as fasta_file:
    parser = fp.Reader(fasta_file)
    for seq in parser:
        total = (len(seq.sequence_as_string()))
    print("The total length of the E. coli genome is {} bp ({} kb)".format(total, total/1000), ".")
    print("The percentage of the genome which is coding is {} %".format((CDS_len/total)*100))
