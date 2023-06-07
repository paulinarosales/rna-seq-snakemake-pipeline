log <- file(snakemake@log[[1]], open="wt")
sink(log)
sink(log, type="message")

suppressPackageStartupMessages({
        library(DESeq2)
        library(biomaRt)
        library(dplyr)
        library(tibble)
        library(ashr)
})

# ---------- Snakemake parsing ---------- #
ddsRds = snakemake@input[["ddsRds"]]
ensembl_geneset = snakemake@input[["ensembl_geneset"]]

contrastRds = snakemake@output[["contrastRds"]]
degs_dir = snakemake@output[["degs_dir"]]
degs_sumTSV = snakemake@output[["degs_summaryTSV"]]
degs_freqTSV = snakemake@output[["degs_freqTSV"]]

fdr_th = as.numeric(snakemake@params[["fdr_th"]])
log2fc_th = as.numeric(snakemake@params[["log2fc_th"]])
padj_th = as.numeric(snakemake@params[["padj_th"]])

de_genes <- c()

cat("Reading input data...", sep="\n")
dds <- readRDS(ddsRds)
resultsNames(dds)

ensembl_geneset <- read.table(ensembl_geneset, header=TRUE, sep="\t")

contrast <- resultsNames(dds)[2:length(resultsNames(dds))]
cat("\n")

cat("Getting DESeq2 results for different contrasts...\n", sep="\n")
cat(paste("FDR used:", fdr_th), sep="\n")

for (con in ncol(contrasts)){
        cat(paste("Contrast:\t", contrasts[1,n], "vs." contrasts[2,n]), sep="\n")
        res <- results(dds, contrast=contrasts[,n], alpha=fdr_th)
        head(res)
        summary(res)
        cat(paste("Number of significant records padj < 0.05:", sum(res$padj < 0.05, na.rm=TRUE)), sep="\n")

        # LFC shrinkage
        cat("Performing LFC shrinkage using ashr method...", sep="\n")
        resLFC <- lfcShrink(dds, coef=con, type="ashr", quiet=TRUE)
        resLFC <- as.data.frame(resLFC) %>% rownames_to_column("ensembl_gene_id")
        resLFC <- inner_join(ensembl_geneset, resLFC, by="ensembl_gene_id")
        resLFC <- resLFC[order(resLFC$log2FoldChange, decreasing = T),] 

        cat("Identifying DEGs...", sep="\n")
        degsAll <- as.data.frame(res) %>% rownames_to_column("ensembl_gene_id_version")
        degsAll <- inner_join(ensembl_geneset, degsAll, by="ensembl_gene_id_version") 
        degsAll <- degsAll[order(degsAll$log2FoldChange, decreasing = T),] # order by fold-change
        cat(paste("Filtering DEGs with >= ", log2fc_th, " log2FC values and < ", padj_th, " p adjusted value..."), sep="\n")
        degsFilter <- filter(degsAll, abs(log2FoldChange) >= log2fc_th, padj < padj_th) # filter thresholds can be changed from config file (default log2FC= 1 padj=0.05)

        con_handle = paste(contrasts[1,n], "_vs_" contrasts[2,n], sep="")
        write.table(resLFC, file=paste(degs_dir, "/", con_handle, "_DE_all_genes_LFCshrink.tsv", sep=""), sep="\t", quote=FALSE, row.names=FALSE)
        write.table(degsAll, file=paste(degs_dir, "/", con_handle, "_DE_all_genes.tsv", sep=""), sep="\t", quote=FALSE, row.names=FALSE)
        write.table(degsFilter, file=paste(degs_dir, "/", con_handle, "_DEGs_only.tsv", sep=""), sep="\t", quote=FALSE, row.names=FALSE)

        stopifnot(nrow(degsAll) > 0 & nrow(degsFilter) > 0)

        de_genes <- c(de_genes, degsFilter$ensembl_gene_id)

        cat("\n")
}

cat("Identifying most frecuent DEGs...\n", sep="\n")
de_genes <- as.data.frame(table(de_genes))
names(de_genes) <- c("ensembl_gene_id", "Freq")
de_genes <- inner_join(ensembl_geneset[,1:2], de_genes, by="ensembl_gene_id")
de_genes <- de_genes[order(de_genes$Freq, decreasing = TRUE),]

write.table(de_genes, file=degs_freqTSV, sep="\t", quote=FALSE, row.names=FALSE)
cat("\n")

cat("Counting total DEGs for all contrasts...\n", sep="\n")
file.create(degs_sumTSV)
command <- paste("wc -l ",  degs_dir, "/*/*_DEGs_only.tsv > ", degs_sumTSV, sep="")
system(command)

degs_sum <- read.table(degs_sumTSV, header=FALSE, sep="")
degs_sum <- degs_sum[-nrow(degs_sum),]
names(degs_sum) <- c("Total_DEGs", "DEGs_only")
degs_sum$Total_DEGs <- degs_sum$Total_DEGs-1
degs_sum <- degs_sum[order(degs_sum$Total_DEGs, decreasing = FALSE),]
degs_sum$DEGs_all <- sub("_DEGs_only.tsv", "_DE_all_genes.tsv", degs_sum$DEGs_only)
degs_sum$DEGs_lfc <- sub("_DEGs_only.tsv", "_DE_all_genes_LFCshrink.tsv", degs_sum$DEGs_only)
degs_sum$Contrast <- sub("_DEGs_only.tsv", "", basename(degs_sum$DEGs_only))

write.table(degs_sum, file=degs_sumTSV, sep="\t", quote=FALSE, row.names=FALSE)
cat("\n")


cat("DONE!", sep="\n")
cat(paste("Output:", resRds, degs_dir, degs_sumTSV, sep="\n\t"), sep="\n")