# RNA-Seq analysis standard pipeline for paired-end data

<blockquote class="callout warning">
  <h4>⚠️ Warning</h4>
  <p>To use this pipeline use the <b>Fork</b> option, on the top right of the repo.</p>
</blockquote>

## Description
This pipeline uses the [Snakemake](https://snakemake.readthedocs.io/en/stable/) workflow manager for PE RNA sequencing data processing.

## Workflow
This is the current workflow for the snakemake rules.

![rulegraph](./snakeflow.svg)


## Getting started
1. Fork the pipeline to create your own project repo.

2. Create the `config/sample_manifest.tsv` using the following structure:

    1. **Identifier:** Facility ID 
    2. **Sample_type:** Tissue, cell type, etc.
    3. **Treatment**: Condition (treated/untreated, fed/starved, etc.)
    4. **Bio_rep:** Biological replicate number
    5. **Target_genome:** ENCODE target genome for sequencing (*supported options:* M25, M32, 19, 38)
    6. **Sequencer**: Sequencer to define —2colour parameter for the trimming (*supported options:* HiSeq4000, NovaSeq, NextSeq500)
    7. **Code:** Single factor defined by important characteristics for 1-1 comparisons for downstream analysis (A, B, C…)
    8. **Pair:** Two factor combination indicating which samples should be compared using the code (AB, AC, BC…)
    9. **Fastq_handle:** Particular handle useful for raw `.fastq.gz` files selection using the name from the facility (number, extension, etc.)


        > Sample manifest example available [here](./config/sample_manifest_example.tsv).
<br/><br/>

3. Copy raw sequencing data. For this there are 2 `input` options (data with 2 sequencing lanes, unique lane/merged files).
    1. **For `fastq.gz` files from 2 sequencing lanes.** Copy _all_ your `fastq.gz` to `resources/fastq_seq/raw`, with no additional subfolder and keeping the original facility names (see bellow for the usual facility filenames). Your folder should look like this:

        ```    
        resources/fastq_seq/raw
            ├── {Identifier}_S{Fastq_handle}_L001_R1_001.fastq.gz
            ├── {Identifier}_S{Fastq_handle}_L002_R1_001.fastq.gz
            ├── {Identifier}_S{Fastq_handle}_L001_R2_001.fastq.gz
            └── {Identifier}_S{Fastq_handle}_L002_R2_001.fastq.gz
        ```
    
    As a result you should have 4 `fastq.gz` files per sample (2 reads with 2 lanes each).
   

    <blockquote class="callout warning">
    <h4>⚠️ Warning</h4>
    <p> Note that the <code>{wildcards}</code> correspond to the columns on the <code>config/sample_manifest.tsv </code>.</p>
    </blockquote>

    2. **For `fastq.gz` files from merged lanes or unique.** Manually copy your `fastq.gz` to `resources/fastq_seq/merged`, with a subfolder for each sample, named using the following structure:


        ```    
        resources/fastq_seq/merged
            ├── {Sample_type}_{Treatment}_Bio-rep_{Bio_rep}
            │   ├── {Sample_type}_{Treatment}_Bio-rep_{Bio_rep}_R1.fq.gz
            │   └── {Sample_type}_{Treatment}_Bio-rep_{Bio_rep}_R2.fq.gz
            │
            └── {Sample_type}_{Treatment}_Bio-rep_{Bio_rep}
                ├── {Sample_type}_{Treatment}_Bio-rep_{Bio_rep}_R1.fq.gz
                └── {Sample_type}_{Treatment}_Bio-rep_{Bio_rep}_R2.fq.gz
         ```


    As a result you should have 2 `fq.gz` files per sample (2 reads each).

    <blockquote class="callout warning">
    <h4>⚠️ Warning</h4>
    <p> Note that the <code>{wildcards}</code> correspond to the columns on the  <code>config/sample_manifest.tsv </code>. And the extension is  <code>fq.gz </code> instead of  <code>fastq.gz </code>.</p>
    </blockquote>

4. Change processing parameters on your `config/config.yaml`.

5. Run the pipeline. 

    Your pipeline should look like this:

    ```    
    <YOUR PROJECT FOLDER>
        ├── config
        │   ├── config.yaml
        │   └── sample_manifest.tsv
        ├── profiles
        |   └── profile
        ├── resources
        |   └── fastq_seq
        |       └── raw or merged
        └── workflow
    ```
    
    Use the following command on your project folder:

    ```
        snakemake --profile profiles/<YOUR_PROFILE> --use-conda -j<N_JOBS>
    ```

<br/><br/>

> For more options checkout the [--help](https://snakemake.readthedocs.io/en/stable/executing/cli.html)

    
## Important notes

- This pipeline only works with _paired-end_ Illumnia sequencing reads <!-- and trimming depending on the specified `Sequencer` in a color-chemistry aware mode. -->


- The original pipeline has a `merge_fq_lanes` rule adapted for particular `input fastq.gz filenames`, for more information check **step 3 on Getting started.**


- The default configurations for this workflow are suitable to run snakemake on the `Helmholtz-Munich HPC` and using conda to deal with software dependencies. If one wishes to run the pipeline on a different computing platform, the profiles need to be adapted accordingly.


- The pipeline only supports GENCODE genomes as reference, choose the right `GENCODE Realease` as value for the `Target_genome` entry in your `config/sample_manifest.tsv`. These are the equivalent realeases for GENCODE and UCSC databases:


    | Organism | GENCODE Realease | GENCODE Genome | UCSC Genome |
    | ----- | ---- | ----- | ---- |
    | Mouse | M25 | GRCm38 | mm10 |
    | Mouse | M32 (Latest) | GRCm38 | mm39 |
    | Human | 19 | GRCh37 | h19 |
    | Human | 38 | GRCh38 | h38 |



## Authors and acknowledgment

Paulina Rosales-Becerra, Kevin Brokers, Saulius Lukauskas & Robert Schneider

## Contact
paulina.rosales@helmholtz-muenchen.de

## Project status
Actual version works just fine :)
