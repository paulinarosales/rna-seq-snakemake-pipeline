log <- file(snakemake@log[[1]], open="wt")
sink(log)
sink(log, type="message")

suppressPackageStartupMessages({
        library(biomaRt)
})

# ---------- Snakemake parsing ---------- #
genesetTSV <- snakemake@output[[1]]
ensembl_params <- snakemake@params[[1]]
ensembl_dataset <- ensembl_params$ensembl_dataset
gene_symbol <- ensembl_params$gene_symbol
ensembl_version <- ensembl_params$ensembl_version



cat("Fetching data from Ensembl...", sep="\n")
mart <- useDataset(ensembl_dataset, useEnsembl(biomart="ensembl", version=ensembl_version))
ensembl <-  getBM(attributes= c("ensembl_gene_id_version", gene_symbol , "chromosome_name", 
                                "start_position", "end_position", "gene_biotype"), mart = mart)
cat("\n")
cat("Filtering non-standard chromosomes and lncRNA genes out...", sep="\n")
# Exclude empty gene symbols and lncRNA genes
ensembl <- ensembl[grep("CHR|GL|KI", ensembl$chromosome_name, invert=T),]
ensembl <- ensembl[ensembl$gene_biotype != "lncRNA",]
ensembl <- ensembl[ensembl[[gene_symbol]] != "",]
names(ensembl[[gene_symbol]]) <- "gene_symbol"
cat("\n")

cat("Saving output data...", sep="\n")
write.table(as.data.frame(ensembl), file=genesetTSV, sep="\t", quote=FALSE, row.names=FALSE)
cat("DONE!", sep="\n")
cat(paste("Output:", genesetTSV, sep="\n\t"), sep="\n")
cat("\n")