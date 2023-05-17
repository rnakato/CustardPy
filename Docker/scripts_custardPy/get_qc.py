#!/usr/bin/env python

import argparse
from tabulate import tabulate

parser = argparse.ArgumentParser()
parser.add_argument('-p', help="Pairtools stat output")

args = parser.parse_args()

output_dict = {}
with open(args.p,'r') as f:
    for line in f:
        attrs = line.split()
        output_dict[attrs[0]] = attrs[1]

table = []
total_reads = int(output_dict["total"])
total_reads_str = format(total_reads, ",d")
table.append(["Total Read Pairs", total_reads_str, "100%"])
unmapped_reads = int(output_dict["total_unmapped"])
percent_unmapped = round(unmapped_reads*100.0/total_reads,2)
unmapped_reads = format(unmapped_reads,",d")
table.append(["Unmapped Read Pairs", unmapped_reads, f"{percent_unmapped}%"])
mapped_reads = int(output_dict["total_mapped"])
percent_mapped = round(mapped_reads*100.0/total_reads,2)
mapped_reads_str = format(mapped_reads,",d")
table.append(["Mapped Read Pairs", mapped_reads_str, f"{percent_mapped}%"])
dup_reads = int(output_dict["total_dups"])
percent_dups = round(dup_reads*100.0/total_reads,2)
dup_reads = format(dup_reads,",d")
table.append(["PCR Dup Read Pairs", dup_reads, f"{percent_dups}%"])
nodup_reads = int(output_dict["total_nodups"])
percent_nodups = round(nodup_reads*100.0/total_reads,2)
nodup_reads_str = format(nodup_reads,",d")
table.append(["No-Dup Read Pairs", nodup_reads_str, f"{percent_nodups}%"])
cis_reads = int(output_dict["cis"])
percent_cis = round(cis_reads*100.0/nodup_reads,2)
cis_reads_str = format(cis_reads,",d")
table.append(["No-Dup Cis Read Pairs", cis_reads_str, f"{percent_cis}%"])
trans_reads = int(output_dict["trans"])
percent_trans = round(trans_reads*100.0/nodup_reads,2)
trans_reads = format(trans_reads,",d")
table.append(["No-Dup Trans Read Pairs", trans_reads, f"{percent_trans}%"])
cis_gt1kb = int(output_dict["cis_1kb+"])
cis_lt1kb = cis_reads - cis_gt1kb
percent_cis_lt1kb = round(cis_lt1kb*100.0/nodup_reads,2)
percent_cis_gt1kb = round(cis_gt1kb*100.0/nodup_reads,2)
cis_gt1kb = format(cis_gt1kb,",d")
cis_lt1kb = format(cis_lt1kb,",d")
valid_read_pairs = int(output_dict["trans"]) + int(output_dict["cis_1kb+"])
percent_valid_read_pairs = round(valid_read_pairs*100.0/nodup_reads, 2)
valid_read_pairs = format(valid_read_pairs,",d")
table.append(["No-Dup Valid Read Pairs (cis >= 1kb + trans)", valid_read_pairs, f"{percent_valid_read_pairs}%"])
table.append(["No-Dup Cis Read Pairs < 1kb", cis_lt1kb, f"{percent_cis_lt1kb}%"])
table.append(["No-Dup Cis Read Pairs >= 1kb", cis_gt1kb, f"{percent_cis_gt1kb}%"])
cis_gt20kb = int(output_dict["cis_20kb+"])
percent_cis_gt20kb = round(cis_gt20kb*100.0/nodup_reads,2)
cis_gt20kb = format(cis_gt20kb,",d")
table.append(["No-Dup Cis Read Pairs >= 20kb", cis_gt20kb, f"{percent_cis_gt20kb}%"])


print(tabulate(table,  tablefmt="plain"))
