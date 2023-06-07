def _params_trim_galore(wildcards):

    sample_type = wildcards['sample_type']
    treatment = wildcards['treatment']
    bio_rep = wildcards['bio_rep']

    _sequencer = SAMPLES.loc[( 
                          sample_type, 
                          treatment,
                          bio_rep), 
                        'Sequencer']

    qual_th = config['TRIM_GALORE']['QUALITY']

    if _sequencer == 'HiSeq4000':
        quality = f'--quality {qual_th}'
        
    elif _sequencer == 'NovaSeq':
        quality = f'--2colour {qual_th}'

    elif _sequencer == 'NextSeq500':
        quality = f'--2colour {qual_th}'

    return quality

rule trim_galore:
    input:
        fq_1 = 'resources/fastq_seq/merged/{sample_type}_{treatment}_Bio-rep_{bio_rep}/{sample_type}_{treatment}_Bio-rep_{bio_rep}_R1.fq.gz',
        fq_2 = 'resources/fastq_seq/merged/{sample_type}_{treatment}_Bio-rep_{bio_rep}/{sample_type}_{treatment}_Bio-rep_{bio_rep}_R2.fq.gz'
    output:
        out_dir = directory('results/fastq/trimmed/{sample_type}_{treatment}_Bio-rep_{bio_rep}'),
        fastqc_dir = directory('results/quality_control/trimmed/{sample_type}_{treatment}_Bio-rep_{bio_rep}'),
        report1 = 'results/fastq/trimmed/{sample_type}_{treatment}_Bio-rep_{bio_rep}/{sample_type}_{treatment}_Bio-rep_{bio_rep}_R1.fq.gz_trimming_report.txt',
        report2 = 'results/fastq/trimmed/{sample_type}_{treatment}_Bio-rep_{bio_rep}/{sample_type}_{treatment}_Bio-rep_{bio_rep}_R2.fq.gz_trimming_report.txt',
        fq1 = 'results/fastq/trimmed/{sample_type}_{treatment}_Bio-rep_{bio_rep}/{sample_type}_{treatment}_Bio-rep_{bio_rep}_val_1.fq.gz',
        fq2 = 'results/fastq/trimmed/{sample_type}_{treatment}_Bio-rep_{bio_rep}/{sample_type}_{treatment}_Bio-rep_{bio_rep}_val_2.fq.gz'
    threads: 15
    resources:
        mem_mb = 5000
    conda:
        '../../envs/alignment/trim-galore.yaml'
    log: 
        'logs/fastq/trim_galore/{sample_type}_{treatment}_Bio-rep_{bio_rep}.log'
    params:
        basename = '{sample_type}_{treatment}_Bio-rep_{bio_rep}',
        quality = _params_trim_galore,
        extra = config['TRIM_GALORE']['EXTRA']
    shell:
        """
            mkdir -p {output.fastqc_dir} && \
            trim_galore {params.quality} {params.extra} -j {threads} \
            --basename {params.basename} \
            --gzip --paired --fastqc_args "--outdir {output.fastqc_dir}"\
            -o {output.out_dir} {input} 2> {log}
        """
        