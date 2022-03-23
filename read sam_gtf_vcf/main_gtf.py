##
# GTF (Gene Transfer Format)
# 1GB can fit entirely in memory
##

import pandas as pd
import re
import json


def read_gtf_feature_count(filename):
    features = {"feature": [], "chromosome": [], "start": [], "end": []}
    with open(filename) as gtf:
        while True:
            line = gtf.readline().rstrip("\n")
            if len(line) == 0:
                # EOF
                break
            if not line.startswith("#"):
                # GTF body
                # fields[0] = CHROMOSOME
                # fields[1] = havana (gene_source)
                # fields[2] = GENE/TRANSCRIPT/EXON/CDS/... (feature)
                # fields[3] = START
                # fields[4] = END
                # fields[5] = .
                # fields[6] = +
                # fields[7] = .
                # fields[8] = INFO
                #   gene_id
                #   gene_version
                #   gene_name
                #   gene_source
                #   gene_biotype
                #   ...
                #   transcript_name (if transcript)
                #   ...
                #   exon_id (if exon)
                fields = line.split("\t")
                chromosome = fields[0]
                feature = fields[2]
                start = fields[3]
                end = fields[4]
                # append elements
                features["feature"].append(feature)
                features["chromosome"].append(chromosome)
                features["start"].append(start)
                features["end"].append(end)
        # end-while
        gtf_df = pd.DataFrame(features)
        df_feature_count = gtf_df['feature'].value_counts()
        return df_feature_count


def read_gtf_get_transcripts(filename, GENE_ID):
    exons = {"transcript_id": [], "exon_id": [], "exon_number": [], "exon_start": [], "exon_r_end": []}
    with open(filename) as gtf:
        while True:
            line = gtf.readline().rstrip("\n")
            if len(line) == 0:
                # EOF
                break
            if not line.startswith("#"):
                # fields[0] = CHROMOSOME
                # fields[1] = havana (gene_source)
                # fields[2] = GENE/TRANSCRIPT/EXON/CDS/etc. (feature) <
                # fields[3] = start <
                # fields[4] = end <
                # fields[5] = .
                # fields[6] = +
                # fields[7] = .
                # fields[8] = INFO
                #   gene_id <
                #   gene_version
                #   gene_name
                #   gene_source
                #   gene_biotype
                #   ...
                #   transcript_id <
                #   ...
                #   exon_id <
                #   exon_number <
                #   ...
                fields = line.split("\t")
                feature = fields[2]
                start = int(fields[3])
                end = int(fields[4])
                INFO = fields[8]
                gene_id = re.search("gene_id \"(.*?)\"", INFO).group(1)
                if gene_id == GENE_ID:
                    if feature == "transcript":
                        transcript_id = re.search("transcript_id \"(.*?)\"", INFO).group(1)
                        print(transcript_id, start, end)
                    elif feature == "exon":
                        transcript_id = re.search("transcript_id \"(.*?)\"", INFO).group(1)
                        exon_id = re.search("exon_id \"(.*?)\"", INFO).group(1)
                        exon_number = int(re.search("exon_number \"(.*?)\"", INFO).group(1))
                        exons["transcript_id"].append(transcript_id)
                        exons["exon_id"].append(exon_id)
                        exons["exon_number"].append(exon_number)
                        exons["exon_start"].append(start)
                        exons["exon_r_end"].append(end - start)
        # end-while
        df = pd.DataFrame(exons).sort_values(
            by=['transcript_id', 'exon_number'],
            ascending=True,
            ignore_index=True)
        return df


def read_gtf_get_cds(filename, GENE_ID):
    transcripts = []
    with open(filename) as gtf:
        while True:
            line = gtf.readline().rstrip("\n")
            if len(line) == 0:
                # EOF
                break
            if not line.startswith("#"):
                # fields[0] = CHROMOSOME
                # fields[1] = havana (gene_source)
                # fields[2] = feature (gene, transcript, exon, CDS, ...)
                # fields[3] = start
                # fields[4] = end
                # fields[5] = score
                # fields[6] = strand (+,-)
                # fields[7] = frame (0,1,2)
                # fields[8] = INFO
                fields = line.split("\t")
                feature = fields[2]
                strand = fields[6]
                INFO = fields[8]
                gene_id = re.search("gene_id \"(.*?)\"", INFO).group(1)
                if gene_id == GENE_ID:
                    # desired gene structure
                    if feature == "transcript":
                        transcript_id = re.search("transcript_id \"(.*?)\"", INFO).group(1)
                        start = int(fields[3])
                        end = int(fields[4])
                        transcript = {
                            "transcript_id": transcript_id,
                            "strand": strand,
                            "start": start,
                            "end": end,
                            "exons": [],
                            "five_prime_utr": [],
                            "three_prime_utr": [],
                        }
                        transcripts.append(transcript)

                    if feature == "five_prime_utr":
                        five_transcript_id = re.search("transcript_id \"(.*?)\"", INFO).group(1)
                        start = int(fields[3])
                        end = int(fields[4])
                        for transcript in transcripts:
                            if transcript["transcript_id"] == five_transcript_id:
                                transcript["five_prime_utr"].append({"start": start, "end": end})
                                break

                    if feature == "three_prime_utr":
                        five_transcript_id = re.search("transcript_id \"(.*?)\"", INFO).group(1)
                        start = int(fields[3])
                        end = int(fields[4])
                        for transcript in transcripts:
                            if transcript["transcript_id"] == five_transcript_id:
                                transcript["three_prime_utr"].append({"start": start, "end": end})
                                break

                    if feature == "exon":
                        exon_transcript_id = re.search("transcript_id \"(.*?)\"", INFO).group(1)
                        exon_id = re.search("exon_id \"(.*?)\"", INFO).group(1)
                        exon_number = int(re.search("exon_number \"(.*?)\"", INFO).group(1))
                        exon_start = int(fields[3])
                        exon_end = int(fields[4])
                        exon = {
                            "exon_id": exon_id,
                            "exon_number": exon_number,
                            "exon_start": exon_start,
                            "exon_end": exon_end,
                            "CDS": [],
                        }
                        for transcript in transcripts:
                            if transcript["transcript_id"] == exon_transcript_id:
                                transcript["exons"].append(exon)
                                break

                    if feature == "CDS":
                        CDS_transcript_id = re.search("transcript_id \"(.*?)\"", INFO).group(1)
                        exon_number = int(re.search("exon_number \"(.*?)\"", INFO).group(1))
                        cds_start = int(fields[3])
                        cds_end = int(fields[4])
                        CDS = {
                            "start": cds_start,
                            "end": cds_end,
                            "start_codon": [cds_start, cds_start+2] if strand == '+' else [cds_end-2, cds_end],
                            "end_codon": [cds_end+1, cds_end+3] if strand == '+' else [cds_start-3, cds_start-1]
                        }
                        for transcript in transcripts:
                            if transcript["transcript_id"] == CDS_transcript_id:
                                for exon in transcript["exons"]:
                                    if exon["exon_number"] == exon_number:
                                        exon["CDS"].append(CDS)
                                        break
        # end-while
        return transcripts


if __name__ == '__main__':
    with pd.option_context('display.max_rows', None, 'display.max_columns', None):
        # <feature, count>
        # print(read_gtf_feature_count("./Homo_sapiens.GRCh38.95.gtf"))
        # gene_id ENSG00000273554: <transcript_id, exons>
        # df_exons = read_gtf_get_transcripts("./Homo_sapiens.GRCh38.95.gtf", "ENSG00000273554")

        # ENSG00000273554
        # ENSG00000278633
        # ENSG00000243485
        gene = "ENSG00000284733"
        splices = read_gtf_get_cds("./Homo_sapiens.GRCh38.95.gtf", gene)
        json_obj = json.dumps(splices, indent=2)
        try:
            output_name = "gene_" + gene + "_structure"
            fp = open(output_name, "w")
            fp.write(json_obj)
            fp.close()
        finally:
            print(json_obj)
