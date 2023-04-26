def _params_get_seq(wildcards):
  
    if wildcards.realease.startswith('M'):
        organism = 'mouse'
    else:
        organism = 'human'

    if wildcards.realease == "M1":
        link_handle = 'genome'
    else:
        link_handle = 'primary_assembly.genome'

    return dict(organism=organism, link_handle=link_handle)


rule get_genome:
    output:
        'resources/external/gencode_{realease}/{genome}_genome.fa.gz'
    params:
        _params_get_seq
    log:
        'logs/external_data/get_gencode_{realease}_{genome}_genome.log'
    threads: 1
    resources:
        http = 1
    shell:
        """
            wget --quiet http://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_{params[0][organism]}/release_{wildcards.realease}/{wildcards.genome}.{params[0][link_handle]}.fa.gz -O  {output}
        """


rule unzip_genome:
    input:
        rules.get_genome.output
    output:
        temp('resources/external/gencode_{realease}/{genome}_genome.fa')
    shell:
        """
             gzip -dc {input} > {output}
        """

rule get_transcriptome:
    output:
        'resources/external/gencode_{realease}/{genome}_transcriptome.fa.gz'
    params:
        _params_get_seq
    log:
        'logs/external_data/get_gencode_{realease}_{genome}_transcriptome.log'
    threads: 1
    resources:
        http = 1
    shell:
        """
            wget --quiet http://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_{params[0][organism]}/release_{wildcards.realease}/gencode.v{wildcards.realease}.transcripts.fa.gz -O  {output}
        """
        # For M1 genome use:
        # """
        # wget --quiet http://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_{params[0][organism]}/release_{wildcards.realease}/gencode.v{wildcards.realease}.pc_transcripts.fa.gz | \
        # gunzip > {output} > {log} 2>&1
        # """


rule unzip_transcriptome:
    input:
        rules.get_transcriptome.output
    output:
        temp('resources/external/gencode_{realease}/{genome}_transcriptome.fa')
    shell:
        """
             gzip -dc {input} > {output}
        """