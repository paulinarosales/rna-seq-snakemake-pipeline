log <- file(snakemake@log[[1]], open="wt")
sink(log)
sink(log, type="message")

suppressPackageStartupMessages({
        library(readr)
        library(tximeta)
        library(SummarizedExperiment)
})

# ---------- Snakemake parsing ---------- #
quant_files = snakemake@input[["quant_files"]]
sample_manifest = snakemake@input[["sample_manifest"]]
se_countsRDS = snakemake@output[["se_countsRDS"]]
gse_countsRDS = snakemake@output[["gse_countsRDS"]]
gse_countsTSV = snakemake@output[["gse_countsTSV"]]


cat("Formating input data...", sep="\n")
coldata <- read.table(sample_manifest, stringsAsFactors=FALSE, header=T, sep="\t")

coldata$names <- basename(dirname(quant_files))
coldata$files <- quant_files # CHECK CORRESPONDANCE
cat("\n")

cat("Reading metadata from...", sep="\n")
head(coldata)
cat("\n")

cat("Importing Salmon data...", sep="\n")
se <- tximeta(coldata, type="salmon")
gse <- summarizeToGene(se)

read_counts = as.data.frame(round(assays(gse)[["counts"]]))
cat("\n")

cat("Saving data...", sep="\n")
saveRDS(se, file=se_countsRDS)
saveRDS(gse, file=gse_countsRDS)
write.table(read_counts, file=gse_countsTSV, sep="\t", quote=F)
cat("DONE!", sep="\n")
cat(paste("Outputs:", se_countsRDS, gse_countsRDS, gse_countsTSV, sep="\n\t"))