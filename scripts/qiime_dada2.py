#!/usr/bin/env python2

import argparse
import subprocess
import os
from os import path
from rpy2 import robjects

def call_failed(cmd, error_path):
    error_msg = open(error_path).read()
    os.remove(error_path)
    raise Error("Error in %s\n" % cmd + error_msg)

def main(args):
    # make output directory
    os.mkdir(args.output_directory)
    log = open(path.join(args.output_directory, 'log.txt'), 'w')
    error_path = path.join(args.output_directory, 'error.txt')
    error = open(error_path, 'w')
    
    # qiime split_library command
    cmd = 'split_libraries_fastq.py'
    failed = subprocess.call([cmd, '-i', args.input_fastq, '-b', args.barcode_fastq, '-o',
                              path.join(args.output_directory, 'slout'), '-m', args.mapping_file, '-r', '1000', "-p",
                              "0.0000001", "-n", "1000", "-q", '0', '--store_demultiplexed_fastq'], stdout=log, stderr=error)
    if failed:
        call_failed(cmd, error_path)
    
    # qiime split_sequence_file_on_sample_ids.py command
    cmd = 'split_sequence_file_on_sample_ids.py'
    failed = subprocess.call(['split_sequence_file_on_sample_ids.py', '-i', path.join(args.output_directory,
                              'slout/seqs.fastq'), '-o', path.join(args.output_directory, 'slout_split/'), '--file_type',
                              'fastq'], stdout=log, stderr=error)
    if failed:
        call_failed(cmd, error_path)
    
    # now ruse rpy2 to run dada2.run
    r_source = robjects.r['source']
    _ = r_source('dada2_single_end_auto.R', echo=False, verbose=False)
    r_run_dada2('data/raw_data/slout_split')
    
    os.remove(error_path)
        
if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--input_fastq", help="location of read 1 fastq file")
    parser.add_argument("-b", "--barcode_fastq", help="location barcode fastq file")
    parser.add_argument("-m", "--mapping_file", help="location of mapping file")
    parser.add_argument("-o", "--output_directory", help="location of output directory")
    args = parser.parse_args()
    main(args)
