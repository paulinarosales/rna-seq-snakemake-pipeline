def _input_merge_fq(wildcards):
    sample_type = wildcards.sample_type
    treatment = wildcards.treatment
    bio_rep = wildcards.bio_rep
    
    identifier = SAMPLES.loc[(sample_type, treatment, bio_rep), 'Identifier']
    handle = SAMPLES.loc[(sample_type, treatment, bio_rep), 'Fastq_handle']

    l_r1 = [f'resources/fastq_seq/raw/{identifier}_S{handle}_L001_R1_001.fastq.gz',  f'resources/fastq_seq/raw/{identifier}_S{handle}_L002_R1_001.fastq.gz']
    l_r2 = [f'resources/fastq_seq/raw/{identifier}_S{handle}_L001_R2_001.fastq.gz',  f'resources/fastq_seq/raw/{identifier}_S{handle}_L002_R2_001.fastq.gz']

    return {'lanes_r1': l_r1, 'lanes_r2': l_r2}



rule merge_fq_lanes:
    input:
        unpack(_input_merge_fq)
    output:
        outdir = directory('resources/fastq_seq/merged/{sample_type}_{treatment}_Bio-rep_{bio_rep}'),
        merged_fq1 = 'resources/fastq_seq/merged/{sample_type}_{treatment}_Bio-rep_{bio_rep}/{sample_type}_{treatment}_Bio-rep_{bio_rep}_R1.fq.gz',
        merged_fq2 = 'resources/fastq_seq/merged/{sample_type}_{treatment}_Bio-rep_{bio_rep}/{sample_type}_{treatment}_Bio-rep_{bio_rep}_R2.fq.gz'
    threads: 6
    shell:
        """
            mkdir -p {output.outdir} && \
            zcat {input.lanes_r1} | gzip > {output.merged_fq1} && \
            zcat {input.lanes_r2} | gzip > {output.merged_fq2}
        """

