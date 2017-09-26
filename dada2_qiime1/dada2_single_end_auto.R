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


run.dada2 <- function(path, analysis.name='dada2', tmp.dir='tmp', min.qual=30, quant=.2, threads=3, skip.len=10, keep.tmp=F) {
	# setup
	fastq.re <- "\\.(fq|fastq)(\\.gz)?$"
	fnFs <- list.files(path, fastq.re)
	sample.names <- sub(fastq.re, "", fnFs)
	if (length(sample.names) == 0) {
		stop("No files ending in fq, fastq, fq.gz or fastq.gz found in dir")
	}
	fnFs <- list.files(path)
	dir.create(tmp.dir)

	setThreadOptions(threads)

	# determine truncation point
	first.bad.F <- unlist(mclapply(paste0(path, '/', fnFs), find.first.bad, first.under=min.qual, ignore.bases=skip.len, mc.cores=threads))
	read.len.F <- unlist(mclapply(paste0(path, '/', fnFs), find.read.len, mc.cores=threads))
	# TODO: fix read.len.F to find shortest position in distribution of read lens within a file
	trunc.len.F <- min(c(quantile(first.bad.F, quant), read.len.F))
	print(trunc.len.F)

	# quality filter files
	fnFs.filt <- paste0(tmp.dir, '/', sample.names, "_filt.fastq.gz")
	out <- filterAndTrim(fnFs, fnFs.filt, truncLen=trunc.len.F, maxN=0, maxEE=c(2,2), truncQ=2, rm.phix=TRUE,
              			 compress=TRUE, multithread=TRUE)
    fnFs.filt <- paste0(tmp.dir, '/', list.files(tmp.dir))
	sample.names <- sapply(list.files(tmp.dir), function(i) {substr(i, 1, nchar(i)-14)})

	# learn error rate
	errF <- learnErrors(fnFs.filt, multithread=TRUE)

	# dereplicate and run dada2
	derepFs <- derepFastq(fnFs.filt)
	names(derepFs) <- sample.names
	dadaFs <- dada(derepFs, err=errF, multithread=TRUE)
	seqtab <- makeSequenceTable(dadaFs)

	# remove chimeras
	seqtab.nochim <- removeBimeraDenovo(seqtab, verbose=TRUE, multithread=TRUE)
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
