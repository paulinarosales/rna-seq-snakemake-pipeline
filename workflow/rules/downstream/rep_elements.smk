
def _input_gtf(wildcards):
    genome = f'resources/external/gencode_{GENCODE_RELEASE}/{GENOME}.annotation.gtf'
    TE = f'resources/external/gencode_{GENCODE_RELEASE}/{GENOME}_Ensembl_rmsk_TE.gtf'
    return  {'genome_gtf': genome, 'TE_gtf': TE}

def _input_for_tecount_mtx(wildcards):
    return  TARGETS['te_count']

rule TEcount:
    input:
        unpack(_input_gtf),
        bam = 'results/bam_files/{sample_type}_{treatment}_Bio-rep_{bio_rep}/{sample_type}_{treatment}_Bio-rep_{bio_rep}.salmon_map.bam'
        # bam = 'results/bam_files/{sample_type}_{treatment}_Bio-rep_{bio_rep}/{sample_type}_{treatment}_Bio-rep_{bio_rep}.hisat2_align.sorted.filtered.bam'
    output:
        outdir = directory('results/read_counts/te_count/{sample_type}_{treatment}_Bio-rep_{bio_rep}'),
        cts = 'results/read_counts/te_count/{sample_type}_{treatment}_Bio-rep_{bio_rep}/{sample_type}_{treatment}_Bio-rep_{bio_rep}_TEcounts.cntTable'
    params:
        basename = '{sample_type}_{treatment}_Bio-rep_{bio_rep}_TEcounts',
        strnd = config['TE_COUNT']['STRANDED'],
        mode = config['TE_COUNT']['MODE'],
        iteration = config['TE_COUNT']['ITERATION'],
        extra = config['TE_COUNT']['EXTRA']
    resources:
        mem_mb = 24000,
        tmpdir = 'tmp'
    threads: 16
    conda:
        '../../envs/downstream/rep_elements.yaml'
    log:
        'logs/te_count/{sample_type}_{treatment}_Bio-rep_{bio_rep}_TEcounts.log'
    shell:
        """
            mkdir -p {output.outdir} && \
            TEcount --sortByPos --format BAM -b {input.bam} \
            --stranded {params.strnd} --mode {params.mode} -i {params.iteration} {params.extra} \
            --GTF {input.genome_gtf} --TE {input.TE_gtf} \
            --project {params.basename} --outdir {output.outdir} 2> {log}
        """

rule TEcount_mtx:
    input:
        cts_files = _input_for_tecount_mtx
    output:
        bothTSV = 'results/read_counts/te_count/allSamples/TEcount_allSamples_geneTEcounts.tsv',
        genectsTSV = 'results/read_counts/te_count/allSamples/TEcount_allSamples_genecounts.tsv',
        tectsTSV = 'results/read_counts/te_count/allSamples/TEcount_allSamples_TEcounts.tsv'
    threads: 4
    log:
        'logs/te_count/allSamples_TEcount_mtx.log'
    conda:
        '../../envs/downstream/differential_expr.yaml'
    script:
        '../../scripts/read_counts/collapse_TEcounts.R'

    # shell:
    #     """
    #         awk 'NR==FNR {{h[$2] = $2; next}} {{print $1,$2,h[$2]}}' {input} > {output}
    #     """
# $ awk 'NR==FNR{array[$1]=$3FS$4FS$5FS$6FS$7FS$8; next} { print $0,array[$1]}' file2 file1
