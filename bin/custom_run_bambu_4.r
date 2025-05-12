#!/usr/bin/env Rscript

# ==============================================================================
# Title: Bambu quantification script for long-read RNA-seq
# Description:
#   - Takes BAM files, genome FASTA, and a GTF annotation.
#   - Quantifies known transcripts (no discovery).
#   - Fixes BiocParallel/CompressedGRangesList errors (GitHub #402).
#   - Simplifies metadata to avoid as.vector errors.
# ==============================================================================

# --- Load required libraries ---
library(bambu)
library(Rsamtools)
library(S4Vectors)
library(BiocParallel)
library(GenomicRanges)
library(dplyr)

# --- Read command-line arguments ---
args = commandArgs(trailingOnly = TRUE)
output_tag <- strsplit(grep('--tag*', args, value = TRUE), split = '=')[[1]][[2]]
ncore <- as.integer(strsplit(grep('--ncore*', args, value = TRUE), split = '=')[[1]][[2]])
genomeseq <- strsplit(grep('--fasta*', args, value = TRUE), split = '=')[[1]][[2]]
annot_gtf <- strsplit(grep('--annotation*', args, value = TRUE), split = '=')[[1]][[2]]
readlist <- args[5:length(args)]

# --- Index the genome FASTA if missing ---
if (!file.exists(paste0(genomeseq, ".fai"))) {
  cat("Indexing genome FASTA...\n")
  Rsamtools::indexFa(genomeseq)
}
genomeSequence <- Rsamtools::FaFile(genomeseq)

# --- Load and preprocess annotations ---
cat("Loading transcript annotation from GTF: ", annot_gtf, "\n")
grlist <- prepareAnnotations(annot_gtf)

# --- Simplify metadata to avoid as.vector errors ---
cat("Simplifying GRanges metadata...\n")
grlist <- lapply(grlist, function(gr) {
  mcols(gr) <- as.data.frame(lapply(mcols(gr), function(col) {
    if (is(col, "CompressedList")) as.character(unlist(col)) else as.character(col)
  }))
  gr
})
grlist <- as(grlist, "SimpleList")

# --- Validate GRangesList ---
cat("Validating GRangesList...\n")
if (any(sapply(grlist, function(x) any(sapply(mcols(x), is, "CompressedList"))))) {
  stop("Metadata still contains CompressedList objects!")
}

# --- Show input BAM files ---
cat("Processing BAM files:\n")
print(readlist)

# --- Prepare read class output directory ---
rc_dir <- "./bambu_rc"
dir.create(rc_dir, showWarnings = TRUE, recursive = TRUE)
cat("Created read class output directory: ", rc_dir, "\n")

# --- Force serial processing ---
register(SerialParam())

# --- Custom overlap function (GitHub #402 fix) ---
split_intersection <- function(grl, geneRanges, ov) {
  multiHits <- which(queryHits(ov) %in% which(countQueryHits(ov) > 1))
  if (length(multiHits) == 0) return(data.frame(queryHits = integer(), subjectHits = integer(), intersectWidth = numeric()))
  cat("Processing", length(multiHits), "multi-hits in chunks...\n")
  dfs <- lapply(split(multiHits, cut(seq_along(multiHits), breaks = 100)), function(multiHit) {
    rangeIntersect <- GenomicRanges::intersect(grl[queryHits(ov)[multiHit]], geneRanges[subjectHits(ov)[multiHit]])
    data.frame(
      queryHits = queryHits(ov)[multiHit],
      intersectWidth = sum(width(rangeIntersect)),
      subjectHits = subjectHits(ov)[multiHit]
    ) %>%
      group_by(queryHits) %>%
      summarise(subjectHits = subjectHits[which.max(intersectWidth)], intersectWidth = max(intersectWidth))
  })
  do.call("rbind", dfs)
}

# --- Patch bambu environment (override internal function) ---
environment(bambu:::processReadsByFile)$split_intersection <- split_intersection

# --- Run bambu quantification ---
cat("Running Bambu quantification...\n")
se <- bambu(
  reads = readlist,
  annotations = grlist,
  genome = genomeSequence,
  ncore = 1,
  lowMemory = TRUE,
  discovery = FALSE,
  rcOutDir = rc_dir,
  verbose = TRUE,
  trackReads = TRUE
)

# --- Print and save output ---
cat("Bambu output summary:\n")
print(summary(se))
writeBambuOutput(se, path = output_tag)
cat("Output written to: ", output_tag, "\n")
