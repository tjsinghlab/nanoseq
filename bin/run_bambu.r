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
    ## - THE PACKAGE BELOW NEEDS TO BE AVAILABLE TO LOAD WHEN RUNNING R

################################################
################################################
## LOAD LIBRARY                               ##
################################################


#!/usr/bin/env Rscript

# Load Libraries
suppressPackageStartupMessages({
  library(bambu)
  library(GenomicRanges)
  library(parallel)
  library(dplyr)
  library(Rsamtools)
})

# Helper Function for Safe Overlap Fix
split_intersection <- function(grl, geneRanges, ov, ncores = 20) {
  multiHits <- which(queryHits(ov) %in% which(countQueryHits(ov) > 1))
  if (length(multiHits) == 0) return(NULL)
  
  dfs = parallel::mclapply(split(multiHits, cut(seq_along(multiHits), 100)),
    function(multiHit) {
      rangeIntersect = GenomicRanges::intersect(
        grl[queryHits(ov)[multiHit]],
        geneRanges[subjectHits(ov)[multiHit]]
      )
      data.frame(
        queryHits = queryHits(ov)[multiHit],
        intersectWidth = sum(width(rangeIntersect)),
        subjectHits = subjectHits(ov)[multiHit]
      ) %>%
        group_by(queryHits) %>%
        summarise(
          subjectHits = subjectHits[which.max(intersectWidth)],
          intersectWidth = max(intersectWidth),
          .groups = "drop"
        )
    },
    mc.cores = ncores
  )
  do.call("rbind", dfs)
}

# Parse Command-Line Arguments
args <- commandArgs(trailingOnly=TRUE)

output_tag     <- strsplit(grep('--tag*', args, value = TRUE), split = '=')[[1]][[2]]
ncore          <- as.integer(strsplit(grep('--ncore*', args, value = TRUE), split = '=')[[1]][[2]])
genomeseq      <- strsplit(grep('--fasta*', args, value = TRUE), split = '=')[[1]][[2]]
genomeSequence <- Rsamtools::FaFile(genomeseq)
Rsamtools::indexFa(genomeseq)
annot_gtf      <- strsplit(grep('--annotation*', args, value = TRUE), split = '=')[[1]][[2]]
readlist       <- args[5:length(args)]

################################################
################################################
## RUN BAMBU                                  ##
################################################
################################################

# Run Annotations Preparation
grlist <- prepareAnnotations(annot_gtf)

# Run Bambu With Error Handling
message("Running bambu with ", length(readlist), " input files...")
se <- tryCatch({
  bambu(
    reads = readlist,
    annotations = grlist,
    genome = genomeSequence,
    ncore = ncore,
    verbose = TRUE,
    lowMemory = TRUE,
    discovery = FALSE
  )
}, error = function(e) {
  if (grepl("CompressedGRangesList", e$message)) {
    message("CompressedGRangesList error detected. Consider increasing memory or using a patched version of bambu.")
    stop(e)
  } else {
    stop(e)
  }
})

# Write Output
writeBambuOutput(se, output_tag)
