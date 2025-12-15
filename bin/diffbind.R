#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)

if(length(args) < 2) {
    stop("Usage: diffbind.R <sample_sheet.csv> <output_prefix>")
}

sample_sheet  <- args[1]
output_prefix <- args[2]

library(DiffBind)
library(rtracklayer)

# Read sample sheet
cdc <- dba(sampleSheet=sample_sheet)

# Count reads
cdc <- dba.count(cdc)

# Set up contrast
cdc <- dba.contrast(cdc, 
                    categories=DBA_CONDITION, 
                    minMembers = 2)

# Analyze using DESeq2
cdc <- dba.analyze(cdc, method=DBA_DESEQ2)

# Report significant peaks
cdc_sig <- dba.report(cdc, bUsePval=TRUE, th=0.01)

# Export significant peaks
export.bed(cdc_sig, paste0(output_prefix, "_significant_peaks_p001.bed"))

peaks_df <- as.data.frame(cdc_sig)

# Apply threshold matching the paper
GAIN  <- peaks_df[peaks_df$Fold > 0, ]   # Gain
LOSS  <- peaks_df[peaks_df$Fold < 0, ]  # Loss

# Export as BED files
export.bed(makeGRangesFromDataFrame(GAIN, keep.extra.columns=TRUE), paste0(output_prefix, "_gain.bed"))
export.bed(makeGRangesFromDataFrame(LOSS, keep.extra.columns=TRUE), paste0(output_prefix, "_loss.bed"))

atac_df <- peaks_df[, c(
  "seqnames",
  "start",
  "end",
  "Fold",
  "p.value",
  "FDR"
)]

# Save as CSV
write.csv(
  atac_df,
  paste0(output_prefix, "_diffbind_results.csv"),
  row.names = FALSE
)