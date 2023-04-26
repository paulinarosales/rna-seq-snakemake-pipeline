rule hisat2_build:
    input:
        'resources/external/gencode_{realease}/{genome}_genome.fa'
    output:
        multiext('resources/external/index/hisat2/gencode_{realease}/{genome}', '.1.ht2', '.2.ht2', '.3.ht2', '.4.ht2', '.5.ht2', '.6.ht2', '.7.ht2', '.8.ht2')
    log:
        'logs/hisat2/gencode_{realease}_{genome}_index.log'
    params:
        basename = directory('resources/external/index/hisat2/gencode_{realease}/{genome}')
    conda:
        '../../envs/alignment/hisat2.yaml'
    threads: 24
    resources:
        mem_mb = 30000
    shell:
        'hisat2-build -p {threads} {input} {params.basename} > {log} 2>&1'

# OPTIMIZE
def _input_for_hisat2_align(wildcards):
    return expand('resources/external/index/hisat2/gencode_{realease}/{genome}.8.ht2', realease=GENCODE_REALEASE, genome=GENOME)

def _params_for_hisat2_align(wildcards):
    return expand('resources/external/index/hisat2/gencode_{realease}/{genome}', realease=GENCODE_REALEASE, genome=GENOME)


rule hisat2_align:
    input:
        index = _input_for_hisat2_align,
        # fq1 = temp(str(TEMP_DIR / 'fastq' / 'trimmed' / '{sample_type}_{treatment}_Bio-rep_{bio_rep}' / '{sample_type}_{treatment}_Bio-rep_{bio_rep}_R1_val_1.fq.gz')),
        # fq2 = temp(str(TEMP_DIR / 'fastq' / 'trimmed' / '{sample_type}_{treatment}_Bio-rep_{bio_rep}' / '{sample_type}_{treatment}_Bio-rep_{bio_rep}_R2_val_2.fq.gz'))
        fq1 = 'results/fastq/trimmed/{sample_type}_{treatment}_Bio-rep_{bio_rep}/{sample_type}_{treatment}_Bio-rep_{bio_rep}_val_1.fq.gz',
        fq2 = 'results/fastq/trimmed/{sample_type}_{treatment}_Bio-rep_{bio_rep}/{sample_type}_{treatment}_Bio-rep_{bio_rep}_val_2.fq.gz'
    output:
        # str(RESULTS_DIR / 'bam_files' / '{sample_type}_{treatment}_Bio-rep_{bio_rep}.bam')
        'results/bam_files/{sample_type}_{treatment}_Bio-rep_{bio_rep}/{sample_type}_{treatment}_Bio-rep_{bio_rep}.hisat2_align.bam'
    log:
        # str(LOG_DIR / 'hisat2' / '{sample_type}_{treatment}_Bio-rep_{bio_rep}_align.log')
        'logs/hisat2/{sample_type}_{treatment}_Bio-rep_{bio_rep}_align.log'
    conda:
        '../../envs/alignment/hisat2.yaml'
    threads: 24
    resources:
        mem_mb = 16000,
        time = '1-00:00:00'
    params:
        index_basename = _params_for_hisat2_align,
        # index_dir = 'data/genome_index/hisat2/{genome}/{genome}',     # basename for indexes
        km = config['HISAT2']['KM']
    shell:
        """
        hisat2 -p {threads} -k {params.km} -x {params.index_basename} \
        -1 {input.fq1} -2 {input.fq2} 2> {log} | \
        samtools view -Sbh -o {output}
        """
