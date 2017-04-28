#!/usr/bin/env python2

import argparse
from dada2_qiime1.dada2_qiime import run
import os
import shutil


def main(args):
    # parse arguments
    input_fastq = os.path.abspath(args.input_fastq)
    skip_split = args.skip_split
    if skip_split:
        barcode_fastq = None
        mapping = None
    else:
        barcode_fastq = os.path.abspath(args.barcode_fastq)
        mapping = os.path.abspath(args.mapping_file)
    output_dir = os.path.abspath(args.output_directory)
    rev_comp_barcodes = args.rev_comp_mapping_barcodes
    pick_otus = args.pick_OTUs
    similarity = args.similarity

    if similarity <= 0 or similarity > 1:
        raise ValueError("Similarity value must be between 0 and 1, %s given" % similarity)

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
    run(input_fastq, barcode_fastq, mapping, rev_comp_barcodes, pick_otus, similarity, skip_split)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument("-i", "--input_fastq", help="location of read 1 fastq file(s)")
    parser.add_argument("-b", "--barcode_fastq", help="location barcode fastq file(s)")
    parser.add_argument("-m", "--mapping_file", help="location of mapping file(s)")
    parser.add_argument("-o", "--output_directory", help="location of output directory")
    parser.add_argument("--pick_OTUs", help="pick otus on dada2 results", default=False, action='store_true')
    parser.add_argument("--similarity", help="similarity threshold for OTU picking, only used with --pick_OTUs flag",
                        default=.99, type=float)
    parser.add_argument("--rev_comp_mapping_barcodes", help="Reverse complement barcodes from mapping file",
                        default=False, action='store_true')
    parser.add_argument("--skip_split", help="skip split libraries and split sequences, give folder of fastq files to "
                                             "-i parameter", default=False, action='store_true')
    parser.add_argument("--force", help="force overwrite of output directory if it already exists", default=False,
                        action='store_true')
    args = parser.parse_args()
    main(args)
