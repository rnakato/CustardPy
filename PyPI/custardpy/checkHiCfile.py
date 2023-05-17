#! /usr/bin/env python
# -*- coding: utf-8 -*-
import numpy as np
import hicstraw
import argparse

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument("hicfile", help="input file (.hic format)", type=str)

    args = parser.parse_args()

    print("Input file: " + args.hicfile)
    hic = hicstraw.HiCFile(args.hicfile)
    print("Genome ID: " + hic.getGenomeID())
    print("Available resolutions: ", end = "")
    print(hic.getResolutions())

    print("\nChromosome name and length")
    for chrom in hic.getChromosomes():
        print(chrom.name, chrom.length)
