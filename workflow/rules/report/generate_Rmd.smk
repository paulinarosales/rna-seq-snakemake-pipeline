def _input_rmd(wildcards):
    return expand('resources/external/genesets/{genome}_ensembl_geneset.tsv', genome=GENOME)



rule de_rmarkdown:
    input:
        rmd = 'workflow/scripts/report/DE_report.Rmd',
        sample_manifestTSV = 'config/sample_manifest.tsv',
        ensembl_geneset = _input_rmd,
        ddsRds = 'results/downstream_analysis/differential_expr/dds.Rds',
        degs_summaryTSV = 'results/downstream_analysis/differential_expr/DEGs/DEGs_summary.tsv',
        degs_freqTSV = 'results/downstream_analysis/differential_expr/DEGs/DEGs_frequency.tsv',
        tcountsRData = 'results/downstream_analysis/differential_expr/transformed_counts.RData'
    output:
        html = 'reports/DE_analysis.html',
        fig_dir = directory(config['FIGURE_DIR'])
    params:
        fdr_th = config['DESEQ2']['FDR_THRESHOLD'],
        log2fc_th = config['DESEQ2']['LOG2FC_THRESHOLD'],
        padj_th_deseq = config['DESEQ2']['PADJ_THRESHOLD'],
        padj_th_cluster = config['CLUSTER_PROF']['PADJ_THRESHOLD'],
        min_gss = config['CLUSTER_PROF']['MIN_GENESET_SIZE'],
        max_gss = config['CLUSTER_PROF']['MAX_GENESET_SIZE']
    log:
        'logs/summary.log'
    conda:
        'envs/Rmd.yaml'
    resources:
        mem_mb = 4000
    shell:
        """
        Rscript -e \"rmarkdown::render('{input.rmd}', output_file='../../../{output.html}', \
        params=list(sample_manifestTSV='../../../{input.sample_manifestTSV}', fig_dir='../../../{output.fig_dir}',\
        ensembl_geneset='../../../{input.ensembl_geneset}', ddsRds='../../../{input.ddsRds}', \
        degs_summaryTSV='../../../{input.degs_summaryTSV}', tcountsRData='../../../{input.tcountsRData}', \
        fdr_th='{params.fdr_th}', log2fc_th='{params.log2fc_th}', padj_th_deseq='{params.padj_th_deseq}, \
        padj_th_cluster='{params.padj_th_cluster}', min_gss='{params.min_gss}', max_gss='{params.max_gss}' \
        log_file='../../../{log}'))\"
        """
    