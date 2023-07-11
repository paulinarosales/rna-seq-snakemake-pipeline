
rule salmon_decoy:
    input:
        # genome_fa = rules.get_fasta.output,
        # transcriptome = rules.get_transcriptome.output
        genome_fa = 'resources/external/gencode_{release}/{genome}_genome.fa',
        transcriptome_fa = 'resources/external/gencode_{release}/{genome}_transcriptome.fa'
    output:
        # decoys = str(GENOME_INDEX_DIR / 'salmon' / '{genome}_decoys.txt'),
        # gentrome = str(GENOME_INDEX_DIR / 'salmon' / '{genome}_gentrome.fa')
        decoys = 'resources/external/index/salmon/gencode_{release}/{genome}_decoys.txt',
        gentrome = 'resources/external/index/salmon/gencode_{release}/{genome}_gentrome.fa'
    log:
        'logs/salmon/gencode_{release}_{genome}_decoy.log'
    threads: 4
    shell:
        # Gathering decoy sequences names
        # Sed command works as follow:
        # -n       = do not print all lines
        # s/ .*//g = Remove anything after spaces. (remove comments)
        # s/>//p  = Remove '>' character at the begining of sequence names. Print names
        # Building big gentrome file
        """
            sed -n 's/ .*//g;s/>//p' {input.genome_fa} > {output.decoys} | \
            cat {input.transcriptome_fa} {input.genome_fa} > {output.gentrome} 2> {log}
        """

rule salmon_index:
    input:
        decoys = 'resources/external/index/salmon/gencode_{release}/{genome}_decoys.txt',
        gentrome = 'resources/external/index/salmon/gencode_{release}/{genome}_gentrome.fa'
    output:
        # str(GENOME_INDEX_DIR / salmon / '{genome}_transcript_index' / 'pos.bin')
        # out_dir = str(GENOME_INDEX_DIR / 'salmon' / '{genome}_transcript_index')
        'resources/external/index/salmon/gencode_{release}/{genome}_transcript_index/pos.bin',
        out_dir = directory('resources/external/index/salmon/gencode_{release}/{genome}_transcript_index')
    log:
        'logs/salmon/gencode_{release}_{genome}_transcript_index.log'
    threads: 16
    resources:
        mem_mb = 28000
    params:
        extra = '--gencode'
    conda:
        '../../envs/alignment/salmon.yaml'
    shell:
        """
            salmon index --transcripts {input.gentrome} --decoys {input.decoys} \
            --index {output.out_dir} -p {threads} {params.extra} > {log} 2>&1
        """

def _input_for_salmon_quant(wildcards):
    index = expand('resources/external/index/salmon/gencode_{release}/{genome}_transcript_index/pos.bin', release=GENCODE_RELEASE, genome=GENOME)
    fq1 = f'results/fastq/trimmed/{wildcards.sample_type}_{wildcards.treatment}_Bio-rep_{wildcards.bio_rep}/{wildcards.sample_type}_{wildcards.treatment}_Bio-rep_{wildcards.bio_rep}_val_1.fq.gz'
    fq2 =  f'results/fastq/trimmed/{wildcards.sample_type}_{wildcards.treatment}_Bio-rep_{wildcards.bio_rep}/{wildcards.sample_type}_{wildcards.treatment}_Bio-rep_{wildcards.bio_rep}_val_2.fq.gz'
    return {'index': index, 'fq1': fq1, 'fq2': fq2}
def _param_for_salmon_quant(wildcards):
    return expand('resources/external/index/salmon/gencode_{release}/{genome}_transcript_index', release=GENCODE_RELEASE, genome=GENOME)



rule salmon_quant:
    input:
        unpack(_input_for_salmon_quant)
        # fq1 = 'results/fastq/trimmed/{sample_type}_{treatment}_Bio-rep_{bio_rep}/{sample_type}_{treatment}_Bio-rep_{bio_rep}_R1_val_1.fq.gz',
        # fq2 = 'results/fastq/trimmed/{sample_type}_{treatment}_Bio-rep_{bio_rep}/{sample_type}_{treatment}_Bio-rep_{bio_rep}_R2_val_2.fq.gz'
    output:
        out_dir = directory('results/read_counts/salmon/{sample_type}_{treatment}_Bio-rep_{bio_rep}'),
        quant = 'results/read_counts/salmon/{sample_type}_{treatment}_Bio-rep_{bio_rep}/quant.sf'
        # lib = str(RESULTS_DIR / read_counts / salmon / '{sample_type}_{treatment}_Bio-rep_{bio_rep}' / 'lib_format_counts.json')
    log:
        'logs/salmon/{sample_type}_{treatment}_Bio-rep_{bio_rep}_quant.log'
    params:
        index_dir = _param_for_salmon_quant,
        libtype = 'A',      # automatic detection of lib
        extra = '-- gcBias'       # https://salmon.readthedocs.io/en/latest/salmon.html
    threads: 24
    resources:
        mem_mb = 28000,
        time = '1-00:00:00'
    conda:
        '../../envs/alignment/salmon.yaml'
    shell:
        """
            salmon quant -i {params.index_dir} -l {params.libtype} \
            -1 {input.fq1} -2 {input.fq2} \
            -p {threads} {params.extra} --validateMappings -o {output.out_dir} > {log} 2>&1
        """