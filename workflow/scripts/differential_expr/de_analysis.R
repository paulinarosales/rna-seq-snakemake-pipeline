log <- file(snakemake@log[[1]], open="wt")
sink(log)
sink(log, type="message")

suppressPackageStartupMessages({
        library(DESeq2)
        library(biomaRt)
        library(dplyr)
        library(tibble)
})

# ---------- Snakemake parsing ---------- #
ddsCollapsedRds = snakemake@input[["ddsCollapsedRds"]]
ensembl_geneset = snakemake@input[["ensembl_geneset"]]
sample_manifest = snakemake@input[["sample_manifest"]]

contrastRds = snakemake@output[["contrastRds"]]
degs_dir = snakemake@output[["degs_dir"]]

log2fc_th = as.numeric(snakemake@params[["log2fc_th"]])
padj_th = as.numeric(snakemake@params[["padj_th"]])

cat("Reading input data...", sep="\n")
ddsCollapsed <- readRDS(ddsCollapsedRds)
resultsNames(ddsCollapsed)

ensembl_geneset <- read.table(ensembl_geneset, header=TRUE, sep="\t")
sample_manifest <- read.table(sample_manifest, header=TRUE, sep="\t")

pairs <- unique(sample_manifest$Pair)
contrasts <- as.data.frame(strsplit(pairs, ""))
cat("\n")

cat("Getting DESeq2 results for different contrasts...\n", sep="\n")
con_list <- c()
for (con in ncol(contrasts)){
        cat(paste("Contrast:\t", contrasts[1,n], "vs." contrasts[2,n]), sep="\n")
        res <- results(ddsCollapsed, contrast=contrasts[,n])
        head(res)
        summary(res)
        cat(paste("Adjusted p-values < 0.05:", sum(res$padj < 0.05, na.rm=TRUE)), sep="\n")

        cat("Identifying DEGs...", sep="\n")
        # Get gene_symbols for ensembl gene IDs
        degsAll <- as.data.frame(res) %>% rownames_to_column("ensembl_gene_id_version")
        degsAll <- inner_join(ensembl_geneset, degsAll, by="ensembl_gene_id_version") 
        # rownames(degsAll) <- degsAll$hgnc_symbol
        # degsAll <- degsAll %>% column_to_rownames("hgnc_symbol") # hgnc_symbol as rownames
        degsAll <- degsAll[order(degsAll$log2FoldChange, decreasing = T),] # order by fold-change
        cat(paste("Filtering DEGs with >= ", log2fc_th, " log2FC values and < ", padj_th, " p adjusted value..."), sep="\n")
        degsFilter <- filter(degsAll, abs(log2FoldChange) >= log2fc_th, padj < padj_th) # filter threshlds can be changed from config file (default log2FC= 1 padj=0.05)

        con_handle = paste(contrasts[1,n], "_vs_" contrasts[2,n], sep="")
        write.table(degsAll, file=paste(degs_dir, "/", con_handle, "_DE_all_genes.tsv", sep=""), sep="\t", quote=FALSE, row.names=FALSE)
        write.table(degsFilter, file=paste(degs_dir, "/", con_handle, "_DEGs_only.tsv", sep=""), sep="\t", quote=FALSE, row.names=FALSE)
        # degsFilter <- degsAll[order(degsAll$log2FoldChange, decreasing = T),] # order by fold-change
        # stopifnot(nrow(degsAll) > 0 & nrow(degsFilter) > 0)
        con_list <- c(con_list, con_handle)
        cat("\n")
}

# cat("Saving output data...", sep="\n")
# saveRDS(res, file=resRds)
saveRDS(con_list, file=contrastRds)
cat("DONE!", sep="\n")
cat(paste("Output:", resRds, degs_dir, sep="\n\t"), sep="\n")