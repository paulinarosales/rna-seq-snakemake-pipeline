# ---------- HISAT2 ---------
def _input_gtf(wildcards):
    return f'resources/external/gencode_{GENCODE_RELEASE}/{GENOME}.annotation.gtf'

rule featurecounts_matrix:
    input:
        bam = TARGETS['samtools'],
        genome_gtf = _input_gtf
    output:
        matrix = 'results/read_counts/featureCounts_allSamples_genecounts.tsv',
        summary = 'results/read_counts/featureCounts_allSamples_genecounts.tsv.summary'
    log: 
        'logs/featureCounts/count_matrix.log'
    conda:
        '../../envs/alignment/hisat2.yaml'
    threads: 24
    params:
        count_mode = config['FEATURE_COUNTS']['COUNT_MODE'],
        ftr_type = config['FEATURE_COUNTS']['FEATURE_TYPE'],
        atr_type = config['FEATURE_COUNTS']['ATTRIBUTE_TYPE'],
        extra =config['FEATURE_COUNTS']['EXTRA']
    shell:
        """
        featureCounts -T {threads} -p --countReadPairs {params.count_mode}\
        -t {params.ftr_type} -g {params.atr_type} {params.extra} \
        -a {input.genome_gtf} -o {output.matrix} {input.bam} 2> {log} 
        """
        
        
# ---------- Salmon ---------        

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
        '../../envs/downstream/differential_expr.yaml'
    script:
        '../../scripts/read_counts/tximeta_GeneSE_counts.R'
