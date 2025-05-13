#!/usr/bin/env Rscript

################################################
## REQUIREMENTS                               ##
################################################
## TRANSCRIPT ISOFORM DISCOVERY AND QUANTIFICATION
## - ALIGNED READS IN BAM FILE FORMAT
## - GENOME SEQUENCE (FASTA)
## - ANNOTATION GTF FILE
## - REQUIRED PACKAGES: bambu, Rsamtools
################################################

library(bambu)
library(Rsamtools)

################################################
## PARSE COMMAND-LINE PARAMETERS              ##
################################################
args <- commandArgs(trailingOnly = TRUE)

# Extract arguments
output_tag <- strsplit(grep('--tag', args, value = TRUE), split = '=')[[1]][[2]]
ncore <- as.integer(strsplit(grep('--ncore', args, value = TRUE), split = '=')[[1]][[2]])
genomeseq <- strsplit(grep('--fasta', args, value = TRUE), split = '=')[[1]][[2]]
annot_gtf <- strsplit(grep('--annotation', args, value = TRUE), split = '=')[[1]][[2]]
readlist <- args[5:length(args)]  # Remaining positional args are BAM file paths

# Index FASTA if not already indexed
if (!file.exists(paste0(genomeseq, ".fai"))) {
  message("Indexing FASTA...")
  Rsamtools::indexFa(genomeseq)
}
genomeSequence <- Rsamtools::FaFile(genomeseq)

################################################
## LOAD ANNOTATION AND RUN BAMBU              ##
################################################

# Prepare annotations
grlist <- prepareAnnotations(annot_gtf)

# SAFETY: Convert to SimpleList if GRangesList is compressed
if (inherits(grlist, "CompressedGRangesList")) {
  grlist <- as(grlist, "SimpleList")
}

# Run bambu (with lowMemory flag for large input)
se <- bambu(
  reads = readlist,
  annotations = grlist,
  genome = genomeSequence,
  ncore = ncore,
  lowMemory = TRUE,
  discovery = FALSE,
  verbose = TRUE
)

################################################
## FILTER LOW-CONFIDENCE TRANSCRIPTS          ##
################################################

read_classes <- rowRanges(se)

# Filtering step: keep transcripts supported by at least 2 reads
if ("nReads" %in% names(mcols(read_classes))) {
  read_classes <- read_classes[mcols(read_classes)$nReads >= 2]
  se <- se[names(read_classes), ]
} else {
  warning("No 'nReads' metadata found; skipping transcript filtering.")
}

################################################
## WRITE OUTPUT                               ##
################################################

writeBambuOutput(se, path = output_tag)