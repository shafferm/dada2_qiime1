library(argparser, quietly=T)
library(dada2, quietly=T)
library(ShortRead, quietly=T)
library(parallel, quietly=T)
library(RcppParallel, quietly=T)

# define argparser
p <- arg_parser("Run DADA2 starting with qiime split_sequence_file_on_sample_ids.py output.")
p <- add_argument(p, "--input_dir", help="Directory with one fastq file per sample")
p <- add_argument(p, "--analysis_name", help="Name to give output biom and fasta files", default="dada2")
p <- add_argument(p, "--temp_dir", help="Directory to store temporary files", default="temp")
p <- add_argument(p, "--num_threads", help="Number of threads to use for multiprocessing", type="numeric", default=3)
p <- add_argument(p, "--min_qual", help="Minimum average quality score to call a read position good", type="numeric", default=30)
p <- add_argument(p, "--quantile", help="Quantile of lengths where mean quality meeds min_qual to use as truncation lenght [0-1]", default=.2)

# parse args
argv <- parse_args(p)
analysis_name <- argv$analysis_name #change this to whatever you want your analysis to be called
path <- argv$input_dir #path to the folder of your split libraries output (multiple fastq files)
temp.dir <- argv$temp_dir
threads <- argv$num_threads
if(is.na(threads)) {
  threads <- defaultNumThreads() - 1
}
min.qual <- argv$min_qual
quant <- argv$quantile

# methods
find.first.bad <- function(i, first.under=30) {
  ##TODO: Add way to filter out really bad samples
  qqF <- qa(i)[["perCycle"]]
  num_cycles <- max(qqF$quality$Cycle)
  mean_quals <- vector(length = num_cycles)
  for(j in seq(1:length(mean_quals))) {
    curr_quals <- qqF$quality[qqF$quality$Cycle==j,]
    # mean_quals[j] <- weighted.mean(curr_quals$Score, curr_quals$Count)
    mean_quals[j] <- median(rep(curr_quals$Score, times=curr_quals$Count))
  }
  # skips first ten bases as they are often low quality
  first_under_30 <- match(TRUE, mean_quals[10:length(mean_quals)]<first.under, nomatch = num_cycles-10)
  first_under_30 <- first_under_30 + 10
  return(first_under_30)
}

fastqFilter.multi <- function(i, inputs, outputs, trim.len.F, trunc.len.F) {
  fastqFilter(inputs[i], outputs[i], maxN=0, maxEE=2, truncQ=2, trimLeft=trim.len.F, truncLen=trunc.len.F, compress=TRUE)
}

# setup
fnFs <- list.files(path)
sample.names <- sapply(strsplit(fnFs, "\\.[^\\.]*$"), `[`, 1)
dir.create(temp.dir)

# determine truncation point
first.bad.F <- unlist(mclapply(paste0(path, '/', fnFs), find.first.bad, first.under=min.qual, mc.cores=threads))
trunc.len.F <- quantile(first.bad.F, .2)

# quality filter files
fnFs.filt <- paste0(temp.dir, '/', sample.names, "_filt.fastq.gz")
junk <- mclapply(seq_along(fnFs), fastqFilter.multi, inputs=paste0(path, '/', fnFs), outputs=fnFs.filt, trim.len.F=5, trunc.len.F=trunc.len.F, mc.cores=threads)

# dereplicate and run dada2
derepFs <- derepFastq(fnFs.filt, verbose=T)
names(derepFs) <- sample.names
dadaFs.lrn <- dada(derepFs, err=NULL, selfConsist=TRUE, multithread=threads)
errF <- dadaFs.lrn[[1]]$err_out
dadaFs <- dada(derepFs, err=errF, multithread=threads)
seqtab <- makeSequenceTable(dadaFs)

# remove chimeras
seqtab.nochim <- removeBimeraDenovo(seqtab, verbose=T)
print(dim(seqtab))

# make a fasta of sequences and rename rows to sequence names
headers <- paste0(rep("seq_", dim(seqtab.nochim)[2]), seq(1:dim(seqtab.nochim)[2]))
seqs <- colnames(seqtab.nochim)
colnames(seqtab.nochim) <- headers
fasta <- paste(paste0(">", headers), seqs, sep = '\n', collapse = '\n')
write(paste0(fasta, '\n'), file = paste0(analysis_name, ".fasta"))

# make into tsv (classic) biom table (sequences as otu names)
seqtab_str <- capture.output(write.table(as.data.frame(t(seqtab.nochim)), sep = '\t', quote = FALSE))
seqtab_str2 <- paste(seqtab_str, sep = '\n', collapse = "\n")
header <- "#OTU table from DADA2\n#OTU ID\t"
seqtab_str2 <- paste0(header, seqtab_str2, '\n')

# this file can then be converted with biom convert to a normal biom file
write(seqtab_str2, file = paste0(analysis_name, ".tsv"))

# clean up
unlink(temp.dir, recursive=T)