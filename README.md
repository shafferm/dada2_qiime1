# dada2_qiime1
Generating all required files for `core_diversity_analyses.py` using QIIME 1 and DADA2 for denoising.

## Setting up
It is recommended to run all of this in a conda environment. If you have not installed conda and bioconda then I recommend you check out their [website](https://bioconda.github.io/) and follow the instructions to both install miniconda and add the bioconda and r channels.

### Setting up in a conda environment
First create a new conda environment and activate it:

```
conda create --name dada2_qiime python=2.7 qiime rpy2 bioconductor-dada2 r-argparser
source activate dada2_qiime
``` 

Install this repository:
```
pip install git+https://github.com/shafferm/dada2_qiime1.git
```

### Setting up without conda (not recommended)
First R must be installed on your system. Then install this repository:
```
pip install git+https://github.com/shafferm/dada2_qiime1.git
```

Install R dependencies
```
install_dada2_qiime_dependencies.py
```

## Running dada2 without OTU picking
```
qiime_dada2.py -i {INSERT_PATH_TO_YOUR_READ_1} -b {INSERT_PATH_TO_YOUR_BARCODE} -m {INSERT_PATH_TO_YOUR_MAPPING_FILE} -o out
```

## Running dada2 with OTU picking
```
qiime_dada2.py -i {INSERT_PATH_TO_YOUR_READ_1} -b {INSERT_PATH_TO_YOUR_BARCODE} -m {INSERT_PATH_TO_YOUR_MAPPING_FILE} -o out --pick_OTUs
```

## Advanced
### Using dada2\_single\_end_auto.R to get DADA2 sequences as OTUs:
1. Run split libraries:
	```
	split_libraries_fastq.py -i {INSERT_PATH_TO_YOUR_READ_1} -b {INSERT_PATH_TO_YOUR_BARCODE} -o slout/ -m {INSERT_PATH_TO_YOUR_MAPPING_FILE} -r 1000 -p 0.0 -n 1000 -q 0 --rev_comp_mapping_barcodes --store_demultiplexed_fastq
	```
2. Split slout into one fastq per sample:
	```
	split_sequence_file_on_sample_ids.py -i slout/seqs.fastq -o slout_split/ --file_type fastq
	```
3. Run DADA2 to denoise samples:
	```
	Rscript dada2_single_end_auto.R --input_dir slout_split
	```

### Assigning taxonomy and getting a tree for DADA2 seqs
1. Assign taxonomy and add to biom table
	```
	assign_taxonomy.py -i dada2.fasta
	biom add-metadata -i dada2.tsv -o dada2_w_tax.biom --observation-metadata-fp uclust_assigned_taxonomy/dada2_tax_assignments.txt --sc-separated taxonomy --observation-header OTUID,taxonomy
	```
2. Align sequences, make a tree, remove pynast failues
	```
	align_seqs.py -i dada2.fasta
	filter_alignment.py -i pynast_aligned/dada2_aligned.fasta
	make_phylogeny.py -i dada2_aligned_pfiltered.fasta -o dada2.tre
	remove_pynast_failures.py -f pynast_aligned/dada2_failures.fasta -i dada2_w_tax.biom -o dada2_w_tax_no_pynast_failures.biom
	```

### Picking OTUs on output of dada2\_single\_end_auto.R
1. Align sequences and prefilter dada2 seqs before OTU picking
    ```
    align_seqs.py -i dada2.fasta
    remove_pynast_failures.py -i dada2.fasta -f pynast_aligned/dada2_failures.fasta -o dada2_no_pynast_failures.fasta
    ```

2. Pick closed reference OTUs, filter out failures and pick rep set:
	```
	pick_otus.py -i dada2_no_pynast_failures.fasta -C -m sortmerna -s .97
	filter_fasta.py -f dada2.fasta -s sortmerna_picked_otus/dada2_failures.txt -o sortmerna_picked_otus/failures.fasta
	pick_rep_set.py -i sortmerna_picked_otus/dada2_otus.txt -o sortmerna_picked_otus/rep_set.fna -f dada2.fasta
	```

3. Pick de novo OTUs and pick rep set:
	```
	pick_otus.py -i sortmerna_picked_otus/failures.fasta -s .97
	pick_rep_set.py -i uclust_picked_otus/failures_otus.txt -o uclust_picked_otus/rep_set.fna -f sortmerna_picked_otus/failures.fasta
	```

4. Concatenate rep_set.fna files from closed and denovo steps:
	```
	cat sortmerna_picked_otus/rep_set.fna uclust_picked_otus/rep_set.fna > rep_set.fna
	```

5. Build otu map and otu table:
	```
	cat sortmerna_picked_otus/dada2_otus.txt uclust_picked_otus/failures_otus.txt > final_otu_map.txt
	dada2_to_otu_table.py -i dada2.tsv -m final_otu_map.txt -o dada2_otu_table.biom
	```

6. Assign taxonomy and add to biom table:
	```
	assign_taxonomy.py -i rep_set.fna
	biom add-metadata -i dada2_otu_table.biom --observation-metadata-fp uclust_assigned_taxonomy/rep_set_tax_assignments.txt -o dada2_otu_table_w_tax.biom --sc-separated taxonomy --observation-header OTUID,taxonomy
	```

7. Align sequences, make a tree, remove pynast failues:
	```
	align_seqs.py -i rep_set.fna
	filter_alignment.py -i pynast_aligned/rep_set_aligned.fasta
	make_phylogeny.py -i rep_set_aligned_pfiltered.fasta -o rep_set.tre
	remove_pynast_failures.py -f pynast_aligned/rep_set_failures.fasta -i dada2_otu_table_w_tax.biom -o dada2_otu_table_w_tax_no_pynast_failures.biom
	```
