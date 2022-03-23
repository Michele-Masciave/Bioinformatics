##
# SAM (Sequence Alignment Map)
# - alignments found between a gene structure of a reference genome and many chromosomes reads
# - 10-30 GB cannot be read entirely in memory
# - read line by line
##

import pandas as pd


def check_overlapping(df, row, next_row):
    # check same chromosome and overlapping
    return df['chromosome'][row] == df['chromosome'][next_row] and \
           df['pos'][next_row] - df['pos'][row] <= df['length'][row]


def read_sam(filename):
    with open(filename) as sam:
        reads = {"read_name": [], "chromosome": [], "pos": [], "length": []}
        overlapped_reads = []
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
                # fields[9] = SEQUENCE
                # ...
                fields = line.split("\t")
                read_name = fields[0]
                chromosome = fields[2]
                POS = int(fields[3])
                MAPQ = int(fields[4])
                length = len(fields[9])
                # read action
                if MAPQ == 60:
                    reads["read_name"].append(read_name)
                    reads["chromosome"].append(chromosome)
                    reads["pos"].append(POS)
                    reads["length"].append(length)
        # end of while
        # create dataframe
        df = pd.DataFrame(reads)
        # sort values for next step better performances
        df.sort_values(by=['chromosome', 'pos'], ascending=True, ignore_index=True, inplace=True)
        # acquire overlapped segments
        for row in range(0, len(df.index)-1):
            next_row = row + 1
            while check_overlapping(df, row, next_row) and next_row < len(df.index)-1:
                overlapped_reads.append({
                    "chromosome": df["chromosome"][row],
                    "read_name_current": df["read_name"][row],
                    "read_name_next": df["read_name"][next_row],
                    "read_pos_current": df['pos'][row],
                    "read_pos_next": df["pos"][next_row],
                    "read_overlap_length": df['pos'][next_row] - df['pos'][row]
                })
                next_row += 1
        # print overlapped segments
        for read in overlapped_reads:
            print(read["read_name_current"], read["read_name_next"], " -- #overlapped bps:", read["read_overlap_length"])


if __name__ == '__main__':
    same_file = "./results.sam"
    read_sam(same_file)
