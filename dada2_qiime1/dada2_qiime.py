import subprocess
from rpy2 import robjects
from dada2_qiime1.util import get_dir
import os
from os import path
from collections import deque


class CommandCaller(object):
    """A class that calls command line tasks via subprocess"""
    def __init__(self, log_path='log.txt', error_path='error.txt'):
        """Takes in log and error paths to record stdout and stderr from called commands"""
        self.log_path = log_path
        self.log = open(log_path, 'w')
        self.error_path = error_path
        self.error = open(error_path, 'w')
        self.queue = deque()

    def add_command(self, cmd):
        """add a command to the queue"""
        self.queue.append(cmd)

    def call_command(self, cmd):
        """
        Calls commands (given as lists of strings), records stdout and stderr
        Raises an informative error if cmd fails
        """
        if subprocess.call(cmd, stdout=self.log, stderr=self.error):
            error_msg = open(self.error_path).read()
            self.exit()
            raise RuntimeError("Error in %s\n" % cmd[-1] + error_msg)

    def call_commands(self):
        """Pop off and call every command in the queue"""
        while len(self.queue) > 0:
            self.call_command(self.queue.popleft())

    def exit(self):
        """Exists """
        self.log.close()
        self.error.close()
        os.unlink(self.error_path)


def run(input_fastq, barcode_fastq, mapping_file):
    commander = CommandCaller()
    # qiime split_library command
    commander.add_command(['split_libraries_fastq.py', '-i', input_fastq, '-b', barcode_fastq, '-o', 'slout', '-m',
                           mapping_file, '-r', '1000', '-p', '0.0000001', '-n', '1000', '-q', '0',
                           '--store_demultiplexed_fastq'])
    # qiime split_sequence_file_on_sample_ids.py command
    commander.add_command(['split_sequence_file_on_sample_ids.py', '-i', 'slout/seqs.fastq', '-o', 'slout_split/',
                           '--file_type', 'fastq'])
    # call first two commands so we are ready for dada2
    commander.call_commands()

    # now use rpy2 to run dada2.run
    r_source = robjects.r['source']
    _ = r_source(path.join(get_dir(), 'dada2_single_end_auto.R'), echo=False, verbose=False)
    r_run_dada2 = robjects.r['run.dada2']
    r_run_dada2('data/raw_data/slout_split')

    # qiime assign taxonomy
    commander.add_command(['assign_taxonomy.py', '-i', 'dada2.fasta'])
    # qiime add metadata to biom
    commander.add_command(['biom add-metadata', '-o', 'dada2_w_tax.biom', '--observation-metadata-fp',
                           'uclust_assigned_taxonomy/dada2_tax_assignments.txt', '--sc-separated', 'taxonomy',
                           '--observation-header', 'OTUID,taxonomy'])
    # qiime align sequences
    commander.add_command(['align_seqs.py', '-i', 'dada2.fasta'])
    # qiime make phylogeny
    commander.add_command(['make_phylogeny.py', '-i', 'pynast_aligned/dada2_aligned.fasta', '-o', 'dada2.tre'])
    # call remove pynast failures
    commander.add_command(['remove_pynast_failures.py', '-f', 'pynast_aligned/dada2_failures.fasta', '-i',
                           'dada2_w_tax.biom', '-o', 'dada2_w_tax_no_pynast_failures.biom'])
