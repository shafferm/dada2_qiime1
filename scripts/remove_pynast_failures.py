#!/usr/bin/env python2

from biom import load_table
import argparse


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--input_file", help="location of biom table")
    parser.add_argument("-o", "--output_file", help="location of output biom table")
    parser.add_argument("-f", "--pynast_fasta", help="location of pynast failures fasta file to be removed")
    args = parser.parse_args()

    p = open(args.pynast_fasta)
    f = p.readlines()
    headers = [f[i].strip() for i in xrange(len(f)) if i % 2 == 0]
    ids_to_toss = [i[1:] for i in headers]

    table = load_table(args.input_file)
    set_to_toss = set(table.ids(axis="observation")) & set(ids_to_toss)

    table.filter(set_to_toss, invert=True, axis="observation")
    table.to_json("remove_pynast_failures.py", open(args.output_file, 'w'))

if __name__ == "__main__":
    main()
