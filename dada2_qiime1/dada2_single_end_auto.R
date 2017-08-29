library(argparser, quietly=T)
library(dada2, quietly=T)
library(ShortRead, quietly=T)
library(parallel, quietly=T)
library(RcppParallel, quietly=T)

# methods
find.first.bad <- function(i, first.under=30, ignore.bases=10) {
	qqF <- qa(i)[["perCycle"]]
	num_cycles <- max(qqF$quality$Cycle)
	mean_quals <- vector(length = num_cycles)
	for(j in seq(1:length(mean_quals))) {
	curr_quals <- qqF$quality[qqF$quality$Cycle==j,]
	# mean_quals[j] <- weighted.mean(curr_quals$Score, curr_quals$Count)
	mean_quals[j] <- median(rep(curr_quals$Score, times=curr_quals$Count))
	}
	# skips first ten bases as they are often low quality
	first_under_30 <- match(TRUE, mean_quals[ignore.bases:length(mean_quals)]<first.under, nomatch = num_cycles-ignore.bases)
	first_under_30 <- first_under_30 + ignore.bases
	return(first_under_30)
}

find.read.len <- function(i) {
  	qqF <- qa(i)[["perCycle"]]
  	return(max(qqF$quality$Cycle)-2)
}               

fastqFilter.multi <- function(i, inputs, outputs, trim.len.F, trunc.len.F) {
  	fastqFilter(inputs[i], outputs[i], maxEE=2,  rm.phix=TRUE, trimLeft=trim.len.F, truncLen=trunc.len.F, compress=TRUE, OMP=F)
}

run.dada2 <- function(path, analysis.name='dada2', tmp.dir='tmp', min.qual=30, quant=.2, threads=3, skip.len=10, keep.tmp=F) {
	# setup
	fnFs <- list.files(path)
	dir.create(tmp.dir)
	sample.names <- sapply(fnFs, function (x) {unlist(strsplit(x, '[.]'))[1]})
	print(length(sample.names))

	# determine truncation point
	first.bad.F <- unlist(mclapply(paste0(path, '/', fnFs), find.first.bad, first.under=min.qual, ignore.bases=skip.len, mc.cores=threads))
	read.len.F <- unlist(mclapply(paste0(path, '/', fnFs), find.read.len, mc.cores=threads))
	trunc.len.F <- min(c(quantile(first.bad.F, quant), read.len.F))
	print(trunc.len.F)

	# quality filter files
	fnFs.filt <- paste0(tmp.dir, '/', sample.names, ".filt.fastq.gz")
	print("before filt")
	junk <- mclapply(seq_along(fnFs), fastqFilter.multi, inputs=paste0(path, '/', fnFs), outputs=fnFs.filt, trim.len.F=5, trunc.len.F=trunc.len.F, mc.cores=threads)
    fnFs.filt <- paste0(tmp.dir, '/', list.files(tmp.dir))
    sample.names <- sapply(list.files(tmp.dir), function (x) {unlist(strsplit(x, '[.]'))[1]})
  
	# dereplicate and run dada2
    print("before drep")
	derepFs <- derepFastq(fnFs.filt)
	names(derepFs) <- sample.names
	print ("before dada2")
	dadaFs.lrn <- dada(derepFs, err=NULL, selfConsist=TRUE, multithread=threads)
	errF <- dadaFs.lrn[[1]]$err_out
	dadaFs <- dada(derepFs, err=errF, multithread=threads)
	seqtab <- makeSequenceTable(dadaFs)

	# remove chimeras
	seqtab.nochim <- removeBimeraDenovo(seqtab, multithreaded=T, verbose=T)
	print(dim(seqtab.nochim))

	# make a fasta of sequences and rename rows to sequence names
	headers <- paste0(rep("seq_", dim(seqtab.nochim)[2]), seq(1:dim(seqtab.nochim)[2]))
	seqs <- colnames(seqtab.nochim)
	colnames(seqtab.nochim) <- headers
	fasta <- paste(paste0(">", headers), seqs, sep = '\n', collapse = '\n')
	write(fasta, file = paste0(analysis.name, ".fasta"))

	# make into tsv (classic) biom table (sequences as otu names)
	seqtab.str <- capture.output(write.table(as.data.frame(t(seqtab.nochim)), sep = '\t', quote = FALSE))
	seqtab.str2 <- paste(seqtab.str, sep = '\n', collapse = "\n")
	header <- "#OTU table from DADA2\n#OTU ID\t"
	seqtab.str2 <- paste0(header, seqtab.str2, '\n')

	# this file can then be converted with biom convert to a normal biom file
	write(seqtab.str2, file = paste0(analysis.name, ".tsv"))

	# clean up
	if (!keep.tmp) {
		unlink(tmp.dir, recursive=T)
	}
}

if (!interactive()) {
	# define argparser
	p <- arg_parser("Run DADA2 starting with qiime split_sequence_file_on_sample_ids.py output.")
	p <- add_argument(p, "--input_dir", help="Directory with one fastq file per sample")
	p <- add_argument(p, "--analysis_name", help="Name to give output biom and fasta files", default="dada2")
	p <- add_argument(p, "--temp_dir", help="Directory to store temporary files", default="tmp")
	p <- add_argument(p, "--num_threads", help="Number of threads to use for multiprocessing", type="numeric", default=3)
	p <- add_argument(p, "--min_qual", help="Minimum average quality score to call a read position good", type="numeric", default=30)
	p <- add_argument(p, "--quantile", help="Quantile of lengths where mean quality meeds min_qual to use as truncation lenght [0-1]", default=.2)
	p <- add_argument(p, "--skip_len", help="Number of bases to skip before starting analysis", default=10)

	# parse args
	argv <- parse_args(p)
	analysis.name <- argv$analysis_name #change this to whatever you want your analysis to be called
	path <- argv$input_dir #path to the folder of your split libraries output (multiple fastq files)
	tmp.dir <- argv$temp_dir
	threads <- argv$num_threads
	if(is.na(threads)) {
	  threads <- defaultNumThreads() - 1
	}
	min.qual <- argv$min_qual
	quant <- argv$quantile
	skip.len = argv$skip_len

	run.dada2(path, analysis.name, tmp.dir, min.qual, quant, threads, skip.len)
}
