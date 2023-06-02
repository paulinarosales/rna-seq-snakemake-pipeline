def _biomart_handles(wildcards):
    if GENCODE_REALEASE.startswith('M'):
        _organism = 'mmusculus_gene_ensembl'
        _geneset = 'mgi_symbol'
    else:
        _organism = 'hsapiens_gene_ensembl'
        _geneset = 'hgnc_symbol'
    return  {'ensembl_dataset': _organism, 'gene_symbol': _geneset, 'ensembl_version': config['ENSEMBL_VERSION']}

rule ensembl_geneset:
    output:
        'resources/external/genesets/{genome}_ensembl_geneset.tsv'
    log:
        'logs/external_data/{genome}_ensembl_geneset.log'
    params:
       _biomart_handles
    conda:
        '../../envs/downstream_analysis/differential_expr.yaml'
    script:
        '../../scripts/external_data/ensembl_geneset.R'


