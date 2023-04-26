
rule get_unaligned_reads:
    """
    Extracts a set of common un-aligned sequences
    """
    input:
        # Add here the output of hisat2 - make sure the bams contain unaligned reads
        bam = 'results/bam_files/{sample_type}_{treatment}_Bio-rep_{bio_rep}/{sample_type}_{treatment}_Bio-rep_{bio_rep}.hisat2_align.sorted.filtered.bam',
        bai = 'results/bam_files/{sample_type}_{treatment}_Bio-rep_{bio_rep}/{sample_type}_{treatment}_Bio-rep_{bio_rep}.hisat2_align.filtered.bam.bai'
    output:
        # Placeholder - you can change that :)
        'results/quality_control/unaligned_reads_hisat2/{sample_type}_{treatment}_Bio-rep_{bio_rep}/{sample_type}_{treatment}_Bio-rep_{bio_rep}.unaligned_sequences.txt'
    params:
        # Consider to parameterize these params in the config
        n = config['UNALIGNED_READS']['HEAD_N'],
        fraction = config['UNALIGNED_READS']['FRACTION']
    threads:
        1
    conda:
        '../../envs/quality_control/samtools.yaml'
    shell:
        """
            set +o pipefail; samtools view -f4 -s {params.fraction} {input.bam} | \
             cut -f 10 | sort | uniq -c | sort -nr | head -n {params.n} > {output}
        """


rule unaligned_sequences_fasta:
    input:
        rules.get_unaligned_reads.output
    output:
        # Placeholder - you can change that :)
         'results/quality_control/unaligned_reads_hisat2/{sample_type}_{treatment}_Bio-rep_{bio_rep}/{sample_type}_{treatment}_Bio-rep_{bio_rep}.raw.unaligned_sequences.fa'
    threads:
        1
    script:
        # This script used to be a run statement, which was stupid in my
        # pipeline since this would only allow for local execution on the entry
        # node. I shortly parsed it into a proper python script - can you make
        # sure that the sample name in the output is porper? Since I added the
        # os.path.basename(path) call just now
        '../../scripts/quality_control/parse_unaligned_reads.py'


rule blast_unaligned_sequences:
    input:
        unalignedFA = rules.unaligned_sequences_fasta.output
    output:
        blastTSV = 'results/quality_control/unaligned_reads_hisat2/{sample_type}_{treatment}_Bio-rep_{bio_rep}/{sample_type}_{treatment}_Bio-rep_{bio_rep}.blast.unaligned_sequences.tsv'
    threads:
        1
    resources:
        http=1
    conda:
        # In the past blast sometimes dropped an error
        '../../envs/quality_control/unaligned_reads.yaml'
    script:
        '../../scripts/quality_control/postprocess_unaligned_reads.py'


# def _input_unaligned_sequence_table (wildcards):
#         sample_ids = SAMPLE_PAIRING.index
#         sequences = [str(QUALITY_CONTROL_DIRECTORY / 'unaligned_reads_bowtie2_alignment' / f'{condition}_{target}_{type}_bio_rep_{bio_rep}' / f'{condition}_{target}_{type}_bio_rep_{bio_rep}.unaligned_sequences.tsv') for condition, target, type, bio_rep in sample_ids]
#         return {'subdata' : sequences}


def _input_unaligned_sequence_table (wildcards):
        sequences =  TARGETS['unaligned_reads']
        return {'subdata' : sequences}


rule unaligned_sequence_table:
    input:
        unpack(_input_unaligned_sequence_table)
    output:
        'results/quality_control/unaligned_reads_hisat2/unaligned_sequences.tsv'
    threads:
        1
    conda:
        '../../envs/quality_control/unaligned_reads.yaml'
    script:
        '../../scripts/quality_control/unaligned_sequence_table.py'