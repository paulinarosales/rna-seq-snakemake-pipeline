# def _params_trim_galore(wildcards):

#     sample_type = wildcards['sample_type']
#     treatment = wildcards['treatment']
#     bio_rep = wildcards['bio_rep']

#     _sequencer = SAMPLES.loc[( 
#                           sample_type, 
#                           treatment,
#                           bio_rep), 
#                         'Sequencer']
#
#     if _sequencer == 'HiSeq4000':
#         params = config['TRIM_GALORE']['HISEQ']
#     if _sequencer == 'NovaSeq':
#         params = config['TRIM_GALORE']['NOVASEQ']
#     if _sequencer == 'NextSeq500':
#         params = config['TRIM_GALORE']['NEXTSEQ']


#     return params

# def _params_trim_galore (wildcards):
#     for sample_type, treatment, bio_rep in SAMPLES.index:
#         basename = f'{sample_type}_{treatment}_Bio-rep_{bio_rep}'
#         out_dir = str('results/fastq/trimmed' / f'{sample_type}_{treatment}_Bio-rep_{bio_rep}')
#         fastqc_dir = str('results/quality_control' / f'{sample_type}_{treatment}_Bio-rep_{bio_rep}')
#     return {'basename':basename, 'out_dir' : out_dir, 'fastqc_dir': fastqc_dir}

# def _params_trim_galore (wildcards):
#     basename = f'{wildcards.sample_type}_{wildcards.treatment}_Bio-rep_{wildcards.bio_rep}'
#     out_dir = f'results/fastq/trimmed/{wildcards.sample_type}_{wildcards.treatment}_Bio-rep_{wildcards.bio_rep}'
#     fastqc_dir = f'results/quality_control/{wildcards.sample_type}_{wildcards.treatment}_Bio-rep_{wildcards.bio_rep}'
#     return {'basename': basename, 'out_dir' : out_dir, 'fastqc_dir': fastqc_dir}


rule trim_galore:
    input:
        # fq_1 = str(FASTQ_DIR / f'{sample_type}_{treatment}_Bio-rep_{bio_rep}_R1.fastq.gz'),
        # fq_2 =  str(FASTQ_DIR / f'{sample_type}_{treatment}_Bio-rep_{bio_rep}_R2.fastq.gz')
        fq_1 = 'resources/raw_data/{sample_type}_{treatment}_Bio-rep_{bio_rep}/{sample_type}_{treatment}_Bio-rep_{bio_rep}_R1.fq.gz',
        fq_2 = 'resources/raw_data/{sample_type}_{treatment}_Bio-rep_{bio_rep}/{sample_type}_{treatment}_Bio-rep_{bio_rep}_R2.fq.gz'
    output:
        # str(TEMP_DIR / 'fastq' / 'trimmed' / f'{sample_type}_{treatment}_Bio-rep_{bio_rep}_R1.fq.gz_trimming_report.txt'),
        # str(TEMP_DIR / 'fastq' / 'trimmed' / f'{sample_type}_{treatment}_Bio-rep_{bio_rep}_R2.fq.gz_trimming_report.txt'),
        # fq1 = temp(str(TEMP_DIR / 'fastq' / 'trimmed' / f'{Sample_type}_{treatment}_Bio-rep_{bio_rep}_R1_val_1.fq.gz')),
        # fq2 = temp(str(TEMP_DIR / 'fastq' / 'trimmed' / f'{Sample_type}_{treatment}_Bio-rep_{bio_rep}_R2_val_2.fq.gz')),
        # out_dir = directory(str(TEMP_DIR / 'fastq' / 'trimmed' / f'{Sample_type}_{treatment}_Bio-rep_{bio_rep}')),
        # fastqc_dir = directory(str(QC_DIR / 'fastqc' / 'trimmed' / f'{Sample_type}_{treatment}_Bio-rep_{bio_rep}')),
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
        # str(LOG_DIR / 'trim_galore' / f'{sample_type}_{treatment}_Bio-rep_{bio_rep}.log')
    params:
#        length = config['trim_galore']['length'],
        # extra = _params_trim_galore
        # out_dir = directory('results/fastq/trimmed/{sample_type}_{treatment}_Bio-rep_{bio_rep}'),
        # fastqc_dir = directory('results/quality_control/{sample_type}_{treatment}_Bio-rep_{bio_rep}'),
        basename = '{sample_type}_{treatment}_Bio-rep_{bio_rep}',
        # out_dir = _params_trim_galore['out_dir'],
        # fastqc_dir = _params_trim_galore['fastqc_dir'],
        # unpack(_params_trim_galore),
        extra = '--quality 20'  # SOLVE PARAMS PARSING
    shell:
        """
            mkdir -p {output.fastqc_dir} && \
            trim_galore {params.extra} -j {threads} \
            --basename {params.basename} \
            --gzip --paired --fastqc_args "--outdir {output.fastqc_dir}"\
            -o {output.out_dir} {input} 2> {log}
        """
        
        # """
        #     mkdir -p {params}[0]['fastqc_dir'] && \
        #     trim_galore {params.extra} -j {threads} \
        #     --basename {params[basename]} \
        #     --gzip --paired --fastqc_args "--outdir {params}[0][fastqc_dir]"\
        #     -o {params}[0]['out_dir'] {input} 2> {log}
        # """
        