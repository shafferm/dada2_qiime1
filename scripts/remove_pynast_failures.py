#!/usr/bin/env python2

from biom import load_table
from skbio.io import read
import argparse


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--input_file", help="location of biom table or fasta file")
    parser.add_argument("-o", "--output_file", help="location of output biom table or fasta file")
    parser.add_argument("-f", "--pynast_fasta", help="location of pynast failures fasta file to be removed")
    args = parser.parse_args()

    ids_to_toss = set()
    for seq in read(args.pynast_fasta, format='fasta'):
        ids_to_toss.update(seq.id)

    if args.input_file.endswith(".biom"):
        table = load_table(args.input_file)
        set_to_toss = set(table.ids(axis="observation")) & ids_to_toss

        table.filter(set_to_toss, invert=True, axis="observation")
        table.to_json("remove_pynast_failures.py", open(args.output_file, 'w'))
    elif args.input_file.endswith(".fasta") or args.input_file.endswith(".fa"):
        with open(args.output_file, 'w') as f:
            for seq in read(args.input_file, format='fasta'):
                if seq.id not in ids_to_toss:
                    f.write('>%s\n%s\n' % (seq.id, str(seq)))
    else:
        raise ValueError("Input file must of type .biom, .fasta or .fa")

if __name__ == "__main__":
    main()
