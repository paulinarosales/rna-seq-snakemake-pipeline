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
ddsRds <- snakemake@input[["ddsRds"]]
geneset <- snakemake@input[["geneset"]]

degs_dir <- snakemake@output[["degs_dir"]]
degs_sumTSV <- snakemake@output[["degs_summaryTSV"]]
degs_freqTSV <- snakemake@output[["degs_freqTSV"]]

fdr_th <- as.numeric(snakemake@params[["fdr_th"]])
log2fc_th <- as.numeric(snakemake@params[["log2fc_th"]])
padj_th <- as.numeric(snakemake@params[["padj_th"]])

de_genes <- c()

cat("Reading input data...", sep="\n")
dds <- readRDS(ddsRds)

geneset <- read.table(geneset, header=TRUE, sep="\t")

contrast <- resultsNames(dds)[2:length(resultsNames(dds))]
cat("\n")

# Summary table
degs_sum <- data.frame(Contrast = contrast, 
                        Total_DEGs = numeric(length(contrast)), 
                        DEGs_only = paste(degs_dir, "/", contrast, "/", contrast, "_DEGs_only.tsv", sep=""), 
                        DEGs_all = paste(degs_dir, "/", contrast, "/", contrast, "_DE_all_genes.tsv", sep=""), 
                        DEGs_lfc = paste(degs_dir, "/", contrast, "/", contrast, "_DE_all_genes_LFCshrink.tsv", sep=""))


cat("Getting DESeq2 results for the following contrasts:\n", sep="\n")
cat(contrast)
cat("\n")
cat(paste("FDR used:", fdr_th), sep="\n")

for (i in 1:nrow(degs_sum)){
        cat(paste("Contrast:\t", degs_sum$Contrast[i]), sep="\n")
        res <- results(dds, name=degs_sum$Contrast[i], alpha=fdr_th)
        head(res)
        summary(res)
        cat(paste("Number of significant records padj 0.05:", sum(res$padj < 0.05, na.rm=TRUE)), sep="\n")

        # LFC shrinkage
        # cat("Performing LFC shrinkage using apeglm method...", sep="\n")
        # resLFC <- lfcShrink(dds, coef=degs_sum$Contrast[i], type="apeglm", quiet=TRUE)
        cat("Performing LFC shrinkage using ashr method...", sep="\n")
        resLFC <- lfcShrink(dds, coef=degs_sum$Contrast[i], type="ashr", quiet=TRUE)
        resLFC <- as.data.frame(resLFC) %>% rownames_to_column("ensembl_gene_id")
        resLFC <- inner_join(geneset, resLFC, by="ensembl_gene_id")
        resLFC <- resLFC[order(resLFC$log2FoldChange, decreasing = T),] 

        cat("Identifying DEGs...", sep="\n")
        degsAll <- as.data.frame(res) %>% rownames_to_column("ensembl_gene_id")
        degsAll <- inner_join(geneset, degsAll, by="ensembl_gene_id") 
        degsAll <- degsAll[order(degsAll$log2FoldChange, decreasing = T),] # order by fold-change
        cat(paste("Filtering DEGs with >= ", log2fc_th, " log2FC values and < ", padj_th, " p adjusted value..."), sep="\n")
        degsFilter <- filter(degsAll, abs(log2FoldChange) >= log2fc_th, padj < padj_th) # filter thresholds can be changed from config file (default log2FC= 1 padj=0.05)

        # con_handle = paste(condition[1], "_vs_" condition[2], sep="")

        # stopifnot(nrow(degsAll) > 0 & nrow(degsFilter) > 0)

        if (!dir.exists(dirname(degs_sum$DEGs_all[i]))){
        dir.create(dirname(degs_sum$DEGs_all[i]))
        }
        
        if(nrow(degsAll) > 0){
                write.table(resLFC, file=degs_sum$DEGs_lfc[i], sep="\t", quote=FALSE, row.names=FALSE)
                write.table(degsAll, file=degs_sum$DEGs_all[i], sep="\t", quote=FALSE, row.names=FALSE)
        }else{
                cat("\tNo entries produced with results()", sep="\n")
        }
        if(nrow(degsFilter) > 0){
                degs_sum$Total_DEGs[i] <- nrow(degsFilter)
                write.table(degsFilter, file=degs_sum$DEGs_only[i], sep="\t", quote=FALSE, row.names=FALSE)
                de_genes <- c(de_genes, degsFilter$ensembl_gene_id)
        }else{  
                cat("\tNo DEGs found after filtering", sep="\n")
        }

        cat("\n")
}

degs_sum <- degs_sum[order(degs_sum$Total_DEGs, decreasing = FALSE),]

write.table(degs_sum, file=degs_sumTSV, sep="\t", quote=FALSE, row.names=FALSE)

cat("Identifying most frequent DEGs...\n", sep="\n")
de_genes <- as.data.frame(table(de_genes))
names(de_genes) <- c("ensembl_gene_id", "Freq")
# de_genes <- inner_join(geneset[,c("ensembl_gene_id", "gene_name")], de_genes, by="ensembl_gene_id")
de_genes <- de_genes[order(de_genes$Freq, decreasing = TRUE),]

write.table(de_genes, file=degs_freqTSV, sep="\t", quote=FALSE, row.names=FALSE)
cat("\n")

cat("DONE!", sep="\n")
cat(paste("Outputs:", degs_dir, degs_sumTSV, degs_freqTSV, sep="\n\t"), sep="\n")