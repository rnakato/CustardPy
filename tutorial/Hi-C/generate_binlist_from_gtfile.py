#!/usr/bin/env python
# -*- coding: utf-8 -*-

import sys
import os
def extract_chromosome_length(file_path, chromosome):
    with open(file_path, 'r') as file:
        for line in file:
            parts = line.split()
            if parts[0] == chromosome:
                return int(parts[1])
    return None


usage = 'Usage: \n    {} <genometable> <chromosome> <resolution>'.format(os.path.basename(__file__))
arguments = sys.argv
if len(arguments) == 1:
    print(usage)
    exit()

gtfile = sys.argv[1]
chromosome = sys.argv[2]
resolution = int(sys.argv[3])

# 染色体22の長さを抽出
chrlen = extract_chromosome_length(gtfile, chromosome)

# ファイルに書き出す
for start in range(0, chrlen, resolution):
    end = start + resolution
    print(f"{chromosome}\t{start}\t{end}\t{start}")
