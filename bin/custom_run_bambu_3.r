#!/usr/bin/env Rscript
library(bambu)
library(Rsamtools)
library(S4Vectors)
library(BiocParallel)
args = commandArgs(trailingOnly=TRUE)
output_tag <- strsplit(grep('--tag*', args, value = TRUE), split = '=')[[1]][[2]]
ncore <- as.integer(strsplit(grep('--ncore*', args, value = TRUE), split = '=')[[1]][[2]])
genomeseq <- strsplit(grep('--fasta*', args, value = TRUE), split = '=')[[1]][[2]]
annot_gtf <- strsplit(grep('--annotation*', args, value = TRUE), split = '=')[[1]][[2]]
readlist <- args[5:length(args)]
if (!file.exists(paste0(genomeseq, ".fai"))) {
  Rsamtools::indexFa(genomeseq)
}
genomeSequence <- Rsamtools::FaFile(genomeseq)
grlist <- prepareAnnotations(annot_gtf)
cat("Processing BAM files:", readlist, "\n")
cat("Annotations loaded from:", annot_gtf, "\n")
cat("Genome sequence from:", genomeseq, "\n")
rc_dir <- "./bambu_rc"
dir.create(rc_dir, showWarnings = TRUE, recursive = TRUE)
cat("Created read class output directory:", rc_dir, "\n")
# Use serial processing
register(SerialParam())
# Run bambu with quantification-only settings
se <- bambu(reads = readlist,
            annotations = grlist,
            genome = genomeSequence,
            ncore = 1,
            lowMemory = TRUE,
            discovery = FALSE,
            rcOutDir = rc_dir,
            verbose = TRUE,
            trackReads = TRUE,
            opt.discovery = list(min.readCount = 2, min.readFractionByGene = 0.05))
cat("Bambu output summary:\n")
print(summary(se))
writeBambuOutput(se, path=output_tag)


#Kept register(SerialParam()) and ncore=1 for serial processing.
#Retained lowMemory=TRUE, discovery=FALSE, and trackReads=TRUE.
#Kept opt.discovery for filtering, even though discovery=FALSE may ignore some options.
####BUT STILL ERROR:
#Primary Error: BiocParallel with CompressedGRangesLis
#Error: BiocParallel errors
#1 remote errors, element index: 1
#first remote error:
#Error in h(simpleError(msg, call)): ... Subsetting operation on CompressedGRangesList object 'x' produces a result that is too big to be represented as a CompressedList object.
#Cause: As per GitHub issue #402, the bambu.processReadsByFile() function generates a massive CompressedGRangesList (e.g., 4.2 billion elements after unlisting) during overlap calculations (generanges[subjectHits(ov)[multiHits]]). This exceeds the memory limits of CompressedList, even in serial mode.
#Fix: Split the overlap processing into chunks, as suggested in the GitHub issue, to reduce memory usage.
