def _input_for_salmon_to_tximeta(wildcards):
    return  TARGETS['salmon']


rule salmon_to_tximeta:
    """
    Collapse all salmon quant.sf files in a single count matrix with Summarized Experiment metadata (SE_counts) 
    and to gene-level (GeneSE_counts)
    """
    input: 
        quant_files = _input_for_salmon_to_tximeta,
        sample_manifest = config['SAMPLE_MANIFEST']
    output:
        se_countsRDS = 'results/read_counts/tximeta/SE_counts.Rds',
        gse_countsRDS = 'results/read_counts/tximeta/GeneSE_counts.Rds',
        gse_countsTSV = 'results/read_counts/tximeta/GeneSE_counts.tsv'
    log: 
        'logs/tximeta/salmon_to_tximeta.log'
    conda:
        '../../envs/downstream_analysis/differential_expr.yaml'
    script:
        '../../scripts/read_counts/tximeta_GeneSE_counts.R'
