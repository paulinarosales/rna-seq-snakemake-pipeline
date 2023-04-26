rule samtools_sort:
    # input:
    #      str(RESULTS_DIR / 'bam_files' / '{sample_type}_{treatment}_Bio-rep_{bio_rep}.bam')
    # output:
    #     str(RESULTS_DIR / 'bam_files' / '{sample_type}_{treatment}_Bio-rep_{bio_rep}.sorted.bam')
    # log:
    #     str(LOG_DIR / 'samtools' / '{sample_type}_{treatment}_Bio-rep_{bio_rep}_sort.log')
    input:
        'results/bam_files/{sample_type}_{treatment}_Bio-rep_{bio_rep}/{sample_type}_{treatment}_Bio-rep_{bio_rep}.hisat2_align.bam'
    output:
        'results/bam_files/{sample_type}_{treatment}_Bio-rep_{bio_rep}/{sample_type}_{treatment}_Bio-rep_{bio_rep}.hisat2_align.sorted.bam'
    log:
        'logs/samtools/{sample_type}_{treatment}_Bio-rep_{bio_rep}/{sample_type}_{treatment}_Bio-rep_{bio_rep}_sort.log'
    conda:
        '../../envs/quality_control/samtools.yaml'
    resources:
        mem_mb = 24000
    threads: 16
    shell:
        """
            samtools sort -@ {threads} -o {output} {input} 2> {log}
        """

rule samtools_filter:
    input:
        'results/bam_files/{sample_type}_{treatment}_Bio-rep_{bio_rep}/{sample_type}_{treatment}_Bio-rep_{bio_rep}.hisat2_align.sorted.bam'
    output:
        'results/bam_files/{sample_type}_{treatment}_Bio-rep_{bio_rep}/{sample_type}_{treatment}_Bio-rep_{bio_rep}.hisat2_align.sorted.filtered.bam'
#        temp('results/samtools/{genome}/{label}.filtered.bam')
    log:
        'logs/samtools/{sample_type}_{treatment}_Bio-rep_{bio_rep}/{sample_type}_{treatment}_Bio-rep_{bio_rep}_filter.log'
    conda:
        '../../envs/quality_control/samtools.yaml'
    resources:
        mem_mb = 24000
    threads: 16
    params:
        min_qual = config['SAMTOOLS']['MIN_MAPQ']
    shell:
        """
            samtools view -b -q {params.min_qual} -@ {threads} -o {output} {input} 2> {log}
        """

rule samtools_index:
    input:
        'results/bam_files/{sample_type}_{treatment}_Bio-rep_{bio_rep}/{sample_type}_{treatment}_Bio-rep_{bio_rep}.hisat2_align.sorted.filtered.bam'
    output:
        'results/bam_files/{sample_type}_{treatment}_Bio-rep_{bio_rep}/{sample_type}_{treatment}_Bio-rep_{bio_rep}.hisat2_align.filtered.bam.bai'
    log: 
        'logs/samtools/{sample_type}_{treatment}_Bio-rep_{bio_rep}/{sample_type}_{treatment}_Bio-rep_{bio_rep}_index.log'
    conda:
        '../../envs/quality_control/samtools.yaml'
    resources:
        mem_mb = 24000
    threads: 16
    shell:
        """
            samtools index -@ {threads} {input} {output} 2> {log}
        """