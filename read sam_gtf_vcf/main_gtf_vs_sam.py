##
# GTF (Gene transfer format)
# - genes/chromosomes structure inside a reference genome
# - GTF occupies about 1GB so can fit entirely in memory
##

import re
import pandas as pd


def read_gtf(filename):
    genes = {"gene_id": [], "gene_name": [], "chromosome": [], "start": [], "end": []}
    with open(filename) as gtf:
        while True:
            line = gtf.readline().rstrip("\n")
            if len(line) == 0:
                # EOF
                break
            if not line.startswith("#"):
                # GTF body
                # fields[0] = CHROMOSOME <
                # fields[1] = havana (gene source)
                # fields[2] = GENE / TRANSCRIPT / EXON / CDS (feature)
                # fields[3] = START POS <
                # fields[4] = END POS <
                # fields[5] = .
                # fields[6] = +
                # fields[7] = .
                # fields[8] = INFO
                #   gene_id "XXX"; <
                #   gene_version "n";
                #   gene_name "XXX";  <
                #   gene_source "havana";
                #   gene_biotype "XXX" [ex: protein_coding]
                #   ...
                #   transcript_name ""; [if TRANSCRIPT]
                #   ...
                #   exon_id ""; [if EXON]
                fields = line.split("\t")
                chromosome = fields[0]
                feature = fields[2]
                start = int(fields[3])
                end = int(fields[4])
                INFO = fields[8]
                gene_id = re.search("gene_id \"(.*?)\"", INFO).group(1)
                gene_name = re.search("gene_name \"(.*?)\"", INFO).group(1)
                gene_biotype = re.search("gene_biotype \"(.*?)\"", INFO).group(1)
                if feature == "gene" and gene_biotype == "protein_coding":
                    genes["gene_id"].append(gene_id)
                    genes["gene_name"].append(gene_name)
                    genes["chromosome"].append(chromosome)
                    genes["start"].append(start)
                    genes["end"].append(end)
        # end-while
        df = pd.DataFrame(genes)
        return df


def read_sam(filename, df_genes):
    gene_expressed = {}
    with open(filename) as sam:
        while True:
            line = sam.readline().rstrip("\n")
            if len(line) == 0:
                # EOF
                break
            elif not line.startswith("@"):
                # SAM body
                # fields[0] = read name
                # fields[1] = flag
                # fields[2] = chromosome
                # fields[3] = position
                # fields[4] = mapping quality
                # fields[5] = CIGAR string
                # ...
                fields = line.split("\t")
                flag = int(fields[1])
                chromosome = fields[2]
                position = int(fields[3])
                # filter out if read:
                # - is unmapped (4=0x4=0000|0000|0100 )
                # - has no secondary alignments (256=0x100=0001|0000|0000)
                if flag & (1 << (3 - 1)) or flag & (1 << (9 - 1)):
                    print("> Not mapped:", fields[0], "| chr:", chromosome, "| flag:", int(bin(flag)[2:]))
                else:
                    for row in range(0, len(df_genes.index)-1):
                        # read chromosomes matches
                        if df_genes["chromosome"][row] == chromosome:
                            # read position matches
                            if df_genes["start"][row] <= position <= df_genes["end"][row]:
                                # add to the dictionary
                                if df_genes["gene_name"][row] not in gene_expressed:
                                    # first time
                                    gene_expressed[df_genes["gene_name"][row]] = 1
                                else:
                                    # accumulate count
                                    gene_expressed[df_genes["gene_name"][row]] += 1
        # end-while
        df_expression = pd.DataFrame(gene_expressed.items(), columns=['Gene Name', 'Expression'])
        df_expression.sort_values(by=['Expression'], ascending=True, ignore_index=True, inplace=True)
        df_expression.to_csv("output_gene_expressed.csv", encoding='utf-8', index=False)
        print(df_expression)


if __name__ == '__main__':
    gtf_file = "./Homo_sapiens.GRCh38.95.gtf"
    sam_file = "./results.sam"
    protein_coding_genes = read_gtf(gtf_file)
    read_sam(sam_file, protein_coding_genes)
