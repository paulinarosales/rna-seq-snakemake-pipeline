log <- file(snakemake@log[[1]], open="wt")
sink(log)
sink(log, type="message")

suppressPackageStartupMessages({
        library(DESeq2)
        library(tibble)
})


# ---------- Snakemake parsing ---------- #
# gse_countsRds <- snakemake@input[["gse_countsRds"]]
sample_manifestTSV <- snakemake@input[["sample_manifestTSV"]]
gene_countsTSV <- snakemake@input[["gene_countsTSV"]]
ddsRds <- snakemake@output[["ddsRds"]]
tcountsRData <- snakemake@output[["tcountsRData"]]
design_formula <- as.formula(snakemake@params[["design_formula"]])
relevel_col <- snakemake@params[["relevel_columns"]]
# subset_col <- snakemake@params[["subset_col"]]
# subset_val <- snakemake@params[["subset_val"]]
reads_th <- snakemake@params[["reads_th"]]
aligner <- snakemake@params[["aligner"]]


if (aligner == "HISAT2"){
        cat("Creating DESeqDataSet from Hisat2-featureCounts count matrix...", sep="\n")
        coldata <- read.table(sample_manifestTSV, header=TRUE, sep="\t")
        rownames(coldata) <- sub("-", ".", paste(coldata$Sample_type, coldata$Treatment, "Bio.rep", coldata$Bio_rep, sep="_"))
        # subset <- rownames(coldata[coldata[[subset_col]]==subset_val,])

        cts <- read.table(gene_countsTSV, header=TRUE, row.names=1, sep="\t", comment.char = "#")
        cts <- cts[,6:ncol(cts)]
        names(cts) <- rownames(coldata)

        # coldata <- coldata[subset,]
        # cts <- cts[,subset]

        dds <- DESeqDataSetFromMatrix(countData = round(cts),
                              colData = coldata,
                              design = design_formula)

}else if (aligner == "SALMON"){
        cat("Creating DESeqDataSet from Salmon-tximeta count matrix...", sep="\n")
        gse <- readRDS(gse_countsRds)
        dds <- DESeqDataSet(gse, design = design_formula)

}else if (aligner == "SALMON-TECOUNT"){
        cat("Creating DESeqDataSet from Salmon-TEcounts matrix...", sep="\n")
        coldata <- read.table(sample_manifestTSV, header=TRUE, sep="\t")
        rownames(coldata) <- sub("-", ".", paste(coldata$Sample_type, coldata$Treatment, "Bio.rep", coldata$Bio_rep, sep="_"))
        # subset <- rownames(coldata[coldata[[subset_col]]==subset_val,])

        cts <- read.table(gene_countsTSV, header=TRUE, row.names=1, sep="\t")
        names(cts) <- rownames(coldata)

        # coldata <- coldata[subset,]
        # cts <- cts[,subset]
        dds <- DESeqDataSetFromMatrix(countData = round(cts),
                              colData = coldata,
                              design = design_formula)
}
cat("\n")

cat("Releveling DESeqDataSet to according to:", sep="\n")
for(col in relevel_col){
        ref_lvl <- as.character(unique(dds[[col]][dds$Control == 1]))
        dds[[col]] <- relevel(dds[[col]], ref = ref_lvl)
        cat(paste("\t", col, ": ", ref_lvl, sep=""), sep="\n")
}

cat("\n")

cat(paste("Pre-filtering counts with <", reads_th, " read counts...", sep=""), sep="\n")
# Exclude low reads
keep <- rowSums(counts(dds)) >=  reads_th
dds <- dds[keep,]
cat("\n")

cat("Running DESeq2...", sep="\n")
dds <- DESeq(dds)
cat("\n")
cat("DESeqDataSet (dds):", sep="\n")
dds

cat("\n")

cat("Transforming DESeq counts with Shifted log-normal, VST and RLog...", sep="\n")
ntd <- normTransform(dds) # Shifted log
vsd <- vst(dds, blind=FALSE) # Variance Stabilizing Transformation 
rld <- rlog(dds, blind=FALSE) # Regularized-logarithm transformation
cat("\n")


cat("Saving output data...",  sep="\n")
saveRDS(dds, file=ddsRds)
save(ntd, vsd, rld, file=tcountsRData)

cat("DONE!", sep="\n")
cat(paste("Output:", ddsRds, tcountsRData, sep="\n\t"), sep="\n")