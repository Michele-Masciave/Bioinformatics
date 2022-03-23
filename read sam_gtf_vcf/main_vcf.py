##
# VCF (Variant Call Format)
# - typically too big for fitting in memory (5-30 GB)
# - variations between reference (gtf) and reads (fastq) according to the SAM file
# - SAM <--- bcftools ---> VCF
##

import json

def vcf_check_snps(filename):
    SNPs = []
    INDELs = []
    with open(filename) as vcf:
        while True:
            line = vcf.readline().rstrip("\n")
            if len(line) == 0:
                # EOF
                break
            if not line.startswith("#"):
                # VCF body
                # fields[0] = CHROMOSOME<
                # fields[1] = POS<
                # fields[2] = id
                # fields[3] = REF<
                # fields[4] = ALT<
                # fields[5] = quality
                # fields[6] = info
                # fields[7] = FORMAT<
                fields = line.split("\t")
                chromosome = fields[0]
                position = fields[1]
                ref_seq = fields[3]
                alt_seq = fields[4]
                FORMAT = fields[7]
                # Looking for SNP (Single Nucleotides Polymorphism) alias mutations
                if alt_seq != "<*>" and len(ref_seq) == len(alt_seq.split(",")[0]):
                    SNPs.append({
                        "chromosome": chromosome,
                        "position": position,
                        "reference": ref_seq,
                        "mutation": alt_seq.split(",")[0]
                    })
                # Looking form INDEL alias of insertions/deletions
                elif "INDEL" in FORMAT:
                    if len(alt_seq.split(",")[0]) > len(ref_seq):
                        INDELs.append({
                            "chromosome": chromosome,
                            "position": position,
                            "reference": ref_seq,
                            "mutation": alt_seq,
                            "type": 'INSERTION'
                        })
                    else:
                        INDELs.append({
                            "chromosome": chromosome,
                            "position": position,
                            "reference": ref_seq,
                            "mutation": alt_seq,
                            "type": 'DELETION'
                        })
    # end-while
    return SNPs, INDELs


if __name__ == '__main__':
    # read VCF looking for SNPs (mutations) and INDEL (insertions/deletions)
    snps, indels = vcf_check_snps("sorted.vcf")
    json_snps = json.dumps(snps, indent=2)
    json_indels = json.dumps(indels, indent=2)
    print(json_snps)
    print(json_indels)