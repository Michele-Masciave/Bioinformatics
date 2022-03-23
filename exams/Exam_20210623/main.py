##
# 23/06/2021
##

import re
import json
import pandas as pd


def read_gtf(filename):
    five_prime_utrs = []
    first_CDS_found = False
    # gtf file can fit entirely in memory generally (≈1GB)
    with open(filename) as gtf:
        while True:
            line = gtf.readline().rstrip("\n")
            if len(line) == 0:
                # EOF
                break
            if not line.startswith("#"):
                # GTF body
                # fields[0] = chromosome <
                # fields[1] = gene_source
                # fields[2] = feature <
                # fields[3] = start <
                # fields[4] = end <
                # fields[5] = score
                # fields[6] = strand (+,-) <
                # fields[7] = frame (0,1,2)
                # fields[8] = INFO <
                #   gene_id
                #   gene_version
                #   gene_name
                #   gene_source
                #   gene_biotype
                #   ...
                #   transcript_id
                #   ...
                #   exon_id
                #   exon_number
                #   ...
                fields = line.split("\t")
                chromosome = fields[0]
                feature = fields[2]
                start = int(fields[3])
                end = int(fields[4])
                strand = fields[6]
                INFO = fields[8]
                # from INFO get gene id and biotype
                gene_id = re.search("gene_id \"(.*?)\"", INFO).group(1)
                gene_biotype = re.search("gene_biotype \"(.*?)\"", INFO).group(1)
                # check if protein coding gene
                if gene_biotype == "protein_coding":
                    # check if gene
                    if feature == "gene":
                        # (re)set CDS indicator to false
                        first_CDS_found = False
                        # + strand -- START [[5UTR][CDS] exon1][exon2][exon3]   END
                        # - strand -- END   [exon3][exon2][exon1 [CDS][5UTR]]   START
                        if strand == "+":
                            five_prime_utrs.append({
                                "gene_id": gene_id,
                                "chromosome": "chr"+chromosome,
                                "strand": strand,
                                "start": start,
                                "end": ''
                            })
                        elif strand == "-":
                            five_prime_utrs.append({
                                "gene_id": gene_id,
                                "chromosome": "chr"+chromosome,
                                "strand": strand,
                                "start": '',
                                "end": end
                            })
                    # check if CDS
                    elif feature == "CDS" and not first_CDS_found:
                        # NOTE:
                        #  exon number is generally one, but pay
                        #  attention cause exons could not have CDS!
                        # + strand --> START [[5UTR][CDS] exon1][exon2][exon3]   END
                        # - strand <-- END   [exon3][exon2][exon1 [CDS][5UTR]]   START
                        for prime in five_prime_utrs:
                            # CDS refer to the same gene
                            if prime["gene_id"] == gene_id:
                                if strand == "+":
                                    prime["end"] = start - 1
                                elif strand == "-":
                                    prime["start"] = end + 1
                                # set CDS indicator to true
                                first_CDS_found = True
        # end-while
        return five_prime_utrs


def read_TF_as_csv(filename):
    df = pd.read_csv(filename, names=['TF', 'chr', 'start', 'end'])
    return df


if __name__ == '__main__':
    # five_primes = read_gtf("../Homo_sapiens.GRCh38.95.gtf")
    # fd = open("Exam_20210623/3_prime_utrs.txt", "w")
    # fd.write(json.dumps(five_primes, indent=1))
    # fd.close()

    json_file = open("Exam_20210623/5_prime_utrs.txt", "r")
    five_primes = json.load(json_file)

    # get transcription factors

    tf = read_TF_as_csv("Exam_20210623/TF.csv")

    for prime in five_primes:
        # Print chromosome, strand, start and end position of the 5p UTR.
        print("chr:", prime["chromosome"], "|", prime["strand"], "|", prime["start"], "->", prime["end"], end="")
        # Print the transcription factors whose coordinates overlap the 5p UTR region previously identified
        # chromosome check
        for row in range(0, len(tf.index)-1):
            if prime["chromosome"] == tf["chr"][row]:
                prime_range = set(range(prime["start"], prime["end"]))
                prime_range.add(prime["end"])
                tf_range = set(range(tf["start"][row], tf["end"][row]))
                tf_range.add(tf["end"][row])
                overlap = tf_range.intersection(prime_range)
                if len(overlap) > 0:
                    print(" |", tf["TF"][row], end="")
        print()

