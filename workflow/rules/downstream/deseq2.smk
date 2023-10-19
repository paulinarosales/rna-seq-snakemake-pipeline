# def _input_deseq_model(wildcards):
#     _aligner = config['ALIGNER']
    
#     if _aligner == 'HISAT2':
#         cts_mtxTSV = TARGETS['feature_counts']
        
#     elif _aligner == 'SALMON-TXIMETA':
#         cts_mtxTSV = TARGETS['tximeta']

#     elif _aligner == 'SALMON-TECOUNT':
#         cts_mtxTSV = TARGETS['te_mtx']
        
#     return cts_mtxTSV



rule deseq_model:
    """
    Creates DESeqDataSet objectes for all raw counts (dds) and collapsed by replicate and filtered (ddsCollapsed) read counts.
    """
    input:
        # gse_countsRds = 'results/read_counts/tximeta/GeneSE_counts.Rds',
        sample_manifestTSV =  config['SAMPLE_MANIFEST'],
        gene_countsTSV = TARGETS['counts']
    output:
        ddsRds = 'results/downstream/differential_expr/dds.Rds',
        tcountsRData = 'results/downstream/differential_expr/transformed_counts.RData'
    log:
        'logs/deseq2/dds_object.log'
    params:
        design_formula = config['DESEQ2']['DESIGN_FORMULA'],
        relevel_columns = config['DESEQ2']['RELEVEL_COLUMNS'],
        subset_col = config['DESEQ2']['SUBSET_COLUMN'],
        subset_val = config['DESEQ2']['SUBSET_VALUE'],
        reads_th = config['DESEQ2']['LOWREADS_THRESHOLD'],
        aligner = config['ALIGNER']
    conda:
        '../../envs/downstream/differential_expr.yaml'
    script:
        '../../scripts/differential_expr/deseq_model.R'

def _input_de_analysis(wildcards):
    return expand('resources/external/gencode_{release}/{genome}.geneset.tsv', release=GENCODE_RELEASE, genome=GENOME)


rule de_analysis:
    input:
        ddsRds = 'results/downstream/differential_expr/dds.Rds',
        geneset = _input_de_analysis
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