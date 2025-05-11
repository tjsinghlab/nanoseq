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
## - THE PACKAGES BELOW NEED TO BE AVAILABLE TO LOAD WHEN RUNNING R

################################################
################################################
## LOAD LIBRARIES                             ##
################################################
################################################
library(bambu)
library(Rsamtools)
library(BiocParallel)

################################################
################################################
## PARSE COMMAND-LINE PARAMETERS              ##
################################################
################################################
args = commandArgs(trailingOnly=TRUE)

output_tag     <- strsplit(grep('--tag*', args, value = TRUE), split = '=')[[1]][[2]]
ncore          <- as.numeric(strsplit(grep('--ncore*', args, value = TRUE), split = '=')[[1]][[2]])
genomeseq      <- strsplit(grep('--fasta*', args, value = TRUE), split = '=')[[1]][[2]]
annot_gtf      <- strsplit(grep('--annotation*', args, value = TRUE), split = '=')[[1]][[2]]
readlist       <- args[5:length(args)]

################################################
################################################
## PREPARE GENOME AND ANNOTATIONS             ##
################################################
################################################
if (!file.exists(paste0(genomeseq, ".fai"))) {
  Rsamtools::indexFa(genomeseq)
}
genomeSequence <- Rsamtools::FaFile(genomeseq)
grlist <- prepareAnnotations(annot_gtf)

################################################
################################################
## RUN BAMBU                                  ##
################################################
################################################
se <- bambu(
  reads = readlist,
  annotations = grlist,
  genome = genomeSequence,
  ncore = ncore,
  verbose = TRUE,
  lowMemory = TRUE,
  discovery = TRUE,
  BPPARAM = SerialParam()
)

writeBambuOutput(se, path = output_tag)
