#!/usr/bin/env python2

from biom import load_table
from biom.table import Table
import argparse
import numpy as np

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--input_file", help="location of biom table")
    parser.add_argument("-m", "--otu_map", help="location of OTU map")
    parser.add_argument("-o", "--output_file", help="location of output biom table")
    args = parser.parse_args()
    
    otu_map = open(args.otu_map, 'U')
    otu_map = otu_map.readlines()
    otu_map = [i.strip().split() for i in otu_map]
    otu_map = {i[0]: i[1:] for i in otu_map}
    
    table = load_table(args.input_file)
    new_table = np.zeros((len(otu_map), table.shape[1]))
    for i, otu in enumerate(otu_map):
        for seq in otu_map[otu]:
            new_table[i,]+=table.data(seq, axis="observation")
    new_table = Table(new_table, otu_map.keys(), table.ids())
    new_table.to_json("dada2_to_otu_table", open(args.output_file, 'w'))

if __name__ == "__main__":
    main()