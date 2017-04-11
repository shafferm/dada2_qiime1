# dada2_qiime1
Generating all required files for `core_diversity_analyses.py` using QIIME 1 and DADA2 for denoising.

## Setting up
First create a new conda environment:

```
conda create --name dada2_qiime python=2.7 qiime R
``` 

Then download this repository, cd into it and then run:

```
Rscript setup_dada2.R
```

## Using dada2\_single\_end_auto.R to get DAD2 sequences as OTUs:
1. Run split libraries:
	```
	split_libraries_fastq.py -i {INSERT_PATH_TO_YOUR_READ_1} -b {INSERT_PATH_TO_YOUR_READ_1} -o slout/ -m {INSERT_PATH_TO_YOUR_MAPPING_FILE} -r 1000 -p 0.0 -n 1000 -q 0 --rev_comp_mapping_barcodes --store_demultiplexed_fastq
	```
2. Split slout into one fastq per sample:
	```
	split_sequence_file_on_sample_ids.py -i slout/seqs.fastq -o slout_split/ --file_type fastq
	```
3. Run DADA2 to denoise samples:
	```
	Rscript dada2_single_end_auto.R --input_dir slout_split
	```
4. Assign taxonomy and add to biom table
	```
	assign_taxonomy.py -i dada2.fasta
	biom add-metadata -i dada2.tsv -o dada2_w_tax.biom --observation-metadata-fp uclust_assigned_taxonomy/dada2_tax_assignments.txt --sc-separated taxonomy --observation-header OTUID,taxonomy
	```
5. Align sequences, make a tree, remove pynast failues
	```
	align_seqs.py -i dada2.fasta
	make_phylogeny.py -i pynast_aligned/dada2_aligned.fasta
	python remove_pynast_failures.py -f pynast_aligned/dada2_failures.fasta -i dada2_w_tax.biom -o dada2_w_tax_no_pynast_failures.biom
	```

NOTE: You'll need to point to `dada2_single_end_auto.R` and `remove_pynast_failures.py` directly as these are not installed and are just custom scripts I wrote that exist in this repository.
