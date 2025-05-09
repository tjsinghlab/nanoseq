#!/usr/bin/env Rscript

################################################
################################################
## REQUIREMENTS                               ##
################################################
################################################
## TRANSCRIPT ISOFORM DISCOVERY AND QUANTIFICATION
## - ALIGNED READS IN BAM FILE FORMAT
## - GENOME SEQUENCE
## - ANNOTATION GTF FILE
## - THE PACKAGES BELOW NEED TO BE AVAILABLE

################################################
################################################
## LOAD LIBRARY                               ##
################################################
################################################
library(bambu)
library(Rsamtools)

################################################
################################################
## PARSE COMMAND-LINE PARAMETERS              ##
################################################
################################################
args = commandArgs(trailingOnly=TRUE)
output_tag <- strsplit(grep('--tag*', args, value = TRUE), split = '=')[[1]][[2]]
ncore <- as.integer(strsplit(grep('--ncore*', args, value = TRUE), split = '=')[[1]][[2]])
genomeseq <- strsplit(grep('--fasta*', args, value = TRUE), split = '=')[[1]][[2]]
annot_gtf <- strsplit(grep('--annotation*', args, value = TRUE), split = '=')[[1]][[2]]
readlist <- args[5:length(args)]

# Index FASTA if not indexed
if (!file.exists(paste0(genomeseq, ".fai"))) {
  Rsamtools::indexFa(genomeseq)
}
genomeSequence <- Rsamtools::FaFile(genomeseq)

################################################
################################################
## RUN BAMBU                                  ##
################################################
################################################
grlist <- prepareAnnotations(annot_gtf)

# Run bambu with lowMemory option for large datasets
se <- bambu(reads = readlist, annotations = grlist, genome = genomeSequence, 
            ncore = ncore, lowMemory = TRUE, discovery = TRUE, verbose = TRUE)

# Filter low-confidence transcripts
read_classes <- rowRanges(se)
read_classes <- read_classes[mcols(read_classes) >= 2]
se <- se[names(read_classes), ]

# Write outputs
writeBambuOutput(se, path=output_tag)
