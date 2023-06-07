rule deseq_model:
    """
    Creates DESeqDataSet objectes for all raw counts (dds) and collapsed by replicate and filtered (ddsCollapsed) read counts.
    """
    input:
        gse_countsRds = 'results/read_counts/tximeta/GeneSE_counts.Rds'
    output:
        ddsRds = 'results/downstream/differential_expr/dds.Rds',
        tcountsRData = 'results/downstream/differential_expr/transformed_counts.RData'
    log:
        'logs/deseq2/dds_object.log'
    params:
        reads_th = config['LOW_READS_THRESHOLD']
    conda:
        '../../envs/downstream/differential_expr.yaml'
    script:
        '../../scripts/differential_expr/deseq_model.R'

def _input_de_analysis(wildcards):
    return expand('resources/external/genesets/{genome}_ensembl_geneset.tsv', genome=GENOME)


rule de_analysis:
    input:
        ddsRds = 'results/downstream/differential_expr/dds.Rds',
        ensembl_geneset = _input_de_analysis
    output:
        degs_dir =  directory('results/downstream/differential_expr/DEGs'),
        degs_summaryTSV = 'results/downstream/differential_expr/DEGs/DEGs_summary.tsv',
        degs_freqTSV = 'results/downstream/differential_expr/DEGs/DEGs_frequency.tsv'
    log:
        'logs/deseq2/de_analysis.log'
    params:
        fdr_th = config['DESEQ2']['FDR_THRESHOLD'],
        log2fc_th = config['DESEQ2']['LOG2FC_THRESHOLD'],
        padj_th = config['DESEQ2']['PADJ_THRESHOLD']
    conda:
        '../../envs/downstream/differential_expr.yaml'
    script:
        '../../scripts/differential_expr/de_analysis.R'