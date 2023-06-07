log <- file(snakemake@log[[1]], open="wt")
sink(log)
sink(log, type="message")

suppressPackageStartupMessages({
        library(DESeq2)
})


# ---------- Snakemake parsing ---------- #
gse_countsRds = snakemake@input[["gse_countsRds"]]
ddsRds = snakemake@output[["ddsRds"]]
tcountsRData = snakemake@output[["tcountsRData"]]
reads_th = snakemake@params[["reads_th"]]


cat("Creating DESeq data set...", sep="\n")
gse <- readRDS(gse_countsRds)
dds <- DESeqDataSet(gse, design = ~Treatment)

ref_condition <- unique(dds$Treatment[dds$Control == 1])
dds$Treatment <- relevel(dds$Treatment, ref = ref_condition)
# dds <- DESeqDataSet(gse, design = ~Treatment+Batch)
cat("\n")

cat(paste("Pre-filtering counts with <", reads_th, " read counts...", sep=""), sep="\n")
# Exclude low reads
keep <- rowSums(counts(dds)) >=  reads_th
dds <- dds[keep,]
cat("\n")

cat("Running DESeq2...", sep="\n")
dds <- DESeq(dds)
cat("\n")

cat("Transforming DESeq counts with Shifted log-normal, VST and RLog...", sep="\n")
ntd <- normTransform(dds) # Shifted log
vsd <- vst(dds, blind=FALSE) # Variance Stabilizing Transformation 
rld <- rlog(dds, blind=FALSE) # Regularized-logarithm transformation
cat("\n")


cat("Saving output data...",  sep="\n")
saveRDS(dds, file=ddsRds)
saveRDS(ntd, vsd, rld, file=tcountsRData)
cat("DONE!", sep="\n")
cat(paste("Output:", ddsRds, countsTransfRds, sep="\n\t"), sep="\n")