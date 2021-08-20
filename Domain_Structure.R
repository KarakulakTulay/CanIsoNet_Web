library("wiggleplotr")
library("dplyr")
library("GenomicRanges")
library("GenomicFeatures")
library("biomaRt")
library("ensembldb")
library("EnsDb.Hsapiens.v75")

normal_transcript <- commandArgs(trailingOnly = TRUE)
cancer_transcript <- commandArgs(trailingOnly = TRUE)
genename = commandArgs(trailingOnly = TRUE)
figure <- plotTranscriptsFromEnsembldb(EnsDb.Hsapiens.v75, gene_names = genename,
                                       transcript_ids = c(cancer_transcript, normal_transcript))

image <- image_write(figure, path = NULL, format = "png")


