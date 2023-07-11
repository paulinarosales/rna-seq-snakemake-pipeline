def _params_get_seq(wildcards):
  
    if wildcards.release.startswith('M'):
        organism = 'mouse'
    else:
        organism = 'human'

    if wildcards.release == "M1":
        link_handle = 'genome'
    else:
        link_handle = 'primary_assembly.genome'

    return dict(organism=organism, link_handle=link_handle)


rule get_genome:
    output:
        'resources/external/gencode_{release}/{genome}_genome.fa.gz'
    params:
        _params_get_seq
    log:
        'logs/external_data/get_gencode_{release}_{genome}_genome.log'
    threads: 1
    resources:
        http = 1
    shell:
        """
            wget --quiet http://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_{params[0][organism]}/release_{wildcards.release}/{wildcards.genome}.{params[0][link_handle]}.fa.gz -O  {output}
        """


rule unzip_genome:
    input:
        rules.get_genome.output
    output:
        temp('resources/external/gencode_{release}/{genome}_genome.fa')
    shell:
        """
             gzip -dc {input} > {output}
        """

rule get_transcriptome:
    output:
        'resources/external/gencode_{release}/{genome}_transcriptome.fa.gz'
    params:
        _params_get_seq
    log:
        'logs/external_data/get_gencode_{release}_{genome}_transcriptome.log'
    threads: 1
    resources:
        http = 1
    shell:
        """
            wget --quiet http://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_{params[0][organism]}/release_{wildcards.release}/gencode.v{wildcards.release}.transcripts.fa.gz -O  {output}
        """
        # For M1 genome use:
        # """
        # wget --quiet http://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_{params[0][organism]}/release_{wildcards.release}/gencode.v{wildcards.release}.pc_transcripts.fa.gz | \
        # gunzip > {output} > {log} 2>&1
        # """


rule unzip_transcriptome:
    input:
        rules.get_transcriptome.output
    output:
        temp('resources/external/gencode_{release}/{genome}_transcriptome.fa')
    shell:
        """
             gzip -dc {input} > {output}
        """

rule get_annotation:
    output:
        'resources/external/gencode_{release}/{genome}.annotation.gtf.gz'
    params:
        organism = _params_get_seq
    log:
        'logs/external_data/get_gencode_{release}_{genome}_annotation.log'
    threads: 1
    resources:
        http = 1
    shell:
        """
            wget --quiet http://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_{params.organism}/release_{wildcards.release}/gencode.v{wildcards.release}.annotation.gtf.gz -O  {output} 2> {log}
        """

rule unzip_annotation:
    input:
        rules.get_annotation.output
    output:
        temp('resources/external/gencode_{release}/{genome}.annotation.gtf')
    shell:
        """
             gzip -dc {input} > {output}
        """

rule unzip_te:
    input:
        'resources/external/gencode_{release}/{genome}_Ensembl_rmsk_TE.gtf.gz'
    output:
        temp('resources/external/gencode_{release}/{genome}_Ensembl_rmsk_TE.gtf')
    shell:
        """
             gzip -dc {input} > {output}
        """