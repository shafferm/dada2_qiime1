#!/usr/bin/env python2

import argparse
from dada2_qiime1.dada2_qiime import run
import os


def main(args):
    # parse arguments
    input_fastq = os.path.abspath(args.input_fastq)
    barcode_fastq = os.path.abspath(args.barcode_fastq)
    mapping = os.path.abspath(args.mapping_file)
    output_dir = os.path.abspath(args.output_directory)

    # make output directory
    os.mkdir(output_dir)
    os.chdir(output_dir)

    # run qiime with dada2 denoising
    run(input_fastq, barcode_fastq, mapping)


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--input_fastq", help="location of read 1 fastq file")
    parser.add_argument("-b", "--barcode_fastq", help="location barcode fastq file")
    parser.add_argument("-m", "--mapping_file", help="location of mapping file")
    parser.add_argument("-o", "--output_directory", help="location of output directory")
    parser.add_argument("--pick_OTUs", help="pick otus on dada2 results")
    parser.add_argument("--force", help="force overwrite of output directory if it already exists")
    args = parser.parse_args()
    main(args)
