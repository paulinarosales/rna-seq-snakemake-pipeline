log <- file(snakemake@log[[1]], open="wt")
sink(log)
sink(log, type="message")

suppressPackageStartupMessages({
        library(DESeq2)
})


# ---------- Snakemake parsing ---------- #
gse_countsRds = snakemake@input[["gse_countsRds"]]
ddsRds = snakemake@output[["ddsRds"]]
ddsCollapsedRds = snakemake@output[["ddsCollapsedRds"]]
reads_th = snakemake@params[["reads_th"]]


cat("Creating DESeq data set...", sep="\n")
gse <- readRDS(gse_countsRds)
dds <- DESeqDataSet(gse, design = ~Code)
# dds <- DESeqDataSet(gse, design = ~Treatment+Batch)
cat("\n")


cat("Raw DESeqDataSet (dds):", sep="\n")
dds
cat("\n")

cat("Collapsing samples by biological replicates...", sep="\n")
# Create label with no "Bio_rep_#" to collpase replicates
# dds$Unique_sample  <- paste(dds$Sample_type, "_", dds$Treatment, sep="")
dds$Unique_sample  <- paste(dds$Sample_type, "_", dds$Treatment, "_Bio-rep_", dds$Bio_rep, sep="") # when we only have one sample for each sample type
ddsCollapsed <- collapseReplicates(dds, dds$Unique_sample)

# Check if collapsed and individual counts correspond
for(i in length(levels(dds$Unique_sample))){
    matchLevel <- dds$Unique_sample == levels(dds$Unique_sample)[i]
    stopifnot(all(rowSums(counts(dds[,matchLevel])) == counts(ddsCollapsed[,i])))
}

cat(paste("Collapsed ", ncol(counts(dds))-ncol(counts(ddsCollapsed)), " samples."), sep="\n")
cat("\n")


cat("Pre-filtering counts data for collapsed samples...", sep="\n")
# Exclude low reads
keep <- rowSums(counts(ddsCollapsed)) >=  reads_th
ddsCollapsed <- ddsCollapsed[keep,]
cat("\n")


cat("Collapsed DESeqDataSet (ddsCollapsed):", sep="\n")
ddsCollapsed
cat("\n")


cat("Running DESeq2...", sep="\n")
ddsCollapsed <- DESeq(ddsCollapsed)

cat("Saving output data...",  sep="\n")
saveRDS(dds, file=ddsRds)
saveRDS(ddsCollapsed, file=ddsCollapsedRds)
cat("DONE!", sep="\n")
cat(paste("Output:", ddsRds, ddsCollapsedRds, sep="\n\t"), sep="\n")