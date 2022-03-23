##
# 23/07/2021
##

import re
import pandas as pd


def read_gtf(filename):
    prot_coding = []
    fp = open("prot_coding.txt", "w")
    with open(filename) as gtf:
        while True:
            line = gtf.readline().rstrip("\n")
            if len(line) == 0:
                # EOF
                break
            if not line.startswith("#"):
                # GTF body
                fields = line.split("\t")
                feature = fields[2]
                INFO = fields[8]
                gene_biotype = re.search("gene_biotype \"(.*?)\"", INFO).group(1)
                if feature == "gene" and gene_biotype == "protein_coding":
                    chromosome = "chr" + fields[0]
                    start = fields[3]
                    end = fields[4]
                    strand = fields[6]
                    gene_id = re.search("gene_id \"(.*?)\"", INFO).group(1)
                    gene_name = re.search("gene_name \"(.*?)\"", INFO).group(1)
                    prot_coding.append({
                        "gene_id": gene_id,
                        "gene_name": gene_name,
                        "chromosome": chromosome,
                        "start": int(start),
                        "end": int(end),
                        "strand": strand
                    })
                    fp.write(gene_id+"\t"+gene_name+"\t"+chromosome+"\t"+start+"\t"+end+"\t"+strand+"\n")
        # end-while
        fp.close()
        return prot_coding


def read_csv(filename):
    return pd.read_csv(filename, names=['PR', 'chr', 'start', 'end'], dtype={'start': 'int64', 'end': 'int64'})


def check_overlap(genes, df):
    overlapped_genes = []
    for gene in genes:
        for row in range(len(df)-1):
            # check chromosome
            if gene["chromosome"] == df["chr"][row]:
                # check if promoter in gene
                if gene["start"] <= df["start"][row] and gene["end"] >= df["end"][row]:
                    overlapped_genes.append(gene)
    return overlapped_genes


if __name__ == '__main__':
    prot_coding_genes = read_gtf("../Homo_sapiens.GRCh38.95.gtf")
    promoters = read_csv("promoters.csv")
    selected_genes = check_overlap(prot_coding_genes, promoters)
    print(selected_genes)
    # print_genomic_sequence todo