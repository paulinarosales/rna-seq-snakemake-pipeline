rule deseq_model:
    """
    Creates DESeqDataSet objectes for all raw counts (dds) and collapsed by replicate and filtered (ddsCollapsed) read counts.
    """
    input:
        gse_countsRds = 'results/read_counts/tximeta/GeneSE_counts.Rds'
    output:
        ddsRds = 'results/downstream_analysis/differential_expr/raw_dds.Rds',
        ddsCollapsedRds = 'results/downstream_analysis/differential_expr/dds_collapsed.Rds'
    log:
        'logs/deseq2/dds_object.log'
    params:
        reads_th = config['LOW_READS_THRESHOLD']
    conda:
        '../../envs/downstream_analysis/differential_expr.yaml'
    script:
        '../../scripts/differential_expr/deseq_model.R'

def _input_de_analysis(wildcards):
    return expand('resources/external/genesets/{genome}_ensembl_geneset.tsv', genome=GENOME)


rule de_analysis:
    input:
        ddsCollapsedRds = 'results/downstream_analysis/differential_expr/dds_collapsed.Rds',
        ensembl_geneset = _input_de_analysis
    output:
        degs_summaryTSV = 'results/downstream_analysis/differential_expr/DEGs/DEGs_summary.tsv',
        degs_dir =  directory('results/downstream_analysis/differential_expr/DEGs')
    log:
        'logs/deseq2/de_analysis.log'
    params:
        fdr_th = config['DESEQ2']['FDR_THRESHOLD'],
        log2fc_th = config['DESEQ2']['LOG2FC_THRESHOLD'],
        padj_th = config['DESEQ2']['PADJ_THRESHOLD']
    conda:
        '../../envs/downstream_analysis/differential_expr.yaml'
    script:
        '../../scripts/differential_expr/de_analysis.R'
