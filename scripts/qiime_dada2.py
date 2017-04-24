#!/usr/bin/env python2

import argparse
from dada2_qiime1.dada2_qiime import run
import os
import shutil


def main(args):
    # parse arguments
    input_fastq = os.path.abspath(args.input_fastq)
    barcode_fastq = os.path.abspath(args.barcode_fastq)
    mapping = os.path.abspath(args.mapping_file)
    output_dir = os.path.abspath(args.output_directory)
    rev_comp_barcodes = args.rev_comp_mapping_barcodes
    pick_otus = args.pick_OTUs

    # make output directory
    try:
        os.mkdir(output_dir)
    except OSError as e:
        if args.force:
            shutil.rmtree(output_dir)
            os.mkdir(output_dir)
        else:
            raise e
    os.chdir(output_dir)

    # run qiime with dada2 denoising
    run(input_fastq, barcode_fastq, mapping, rev_comp_barcodes, pick_otus)


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--input_fastq", help="location of read 1 fastq file")
    parser.add_argument("-b", "--barcode_fastq", help="location barcode fastq file")
    parser.add_argument("-m", "--mapping_file", help="location of mapping file")
    parser.add_argument("-o", "--output_directory", help="location of output directory")
    parser.add_argument("--pick_OTUs", help="pick otus on dada2 results", default=False, action='store_true')
    parser.add_argument("--rev_comp_mapping_barcodes", help="Reverse complement barcodes from mapping file",
                        default=False, action='store_true')
    parser.add_argument("--force", help="force overwrite of output directory if it already exists", default=False,
                        action='store_true')
    args = parser.parse_args()
    main(args)
