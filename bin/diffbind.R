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