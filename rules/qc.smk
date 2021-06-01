
# ------------------------------------------------------------------------------------------------
# Perform read quality control and evaluation

rule qc:
    input:
        expand("data/nanoplot_raw/{sample}/NanoStats.txt", sample = SAMPLES),
        expand("data/nanofilt/{sample}_nanofilt.fastq.gz", sample = SAMPLES),
        expand("data/nanoplot_filtered/{sample}/NanoStats.txt", sample = SAMPLES)
    output:
        touch("finished_QC")

rule nanoplot_raw:
    input:
         reads = config["LONG_READ_DIR"] + "/{sample}.fastq.gz"
    output:
          "data/nanoplot_raw/{sample}/NanoStats.txt"
    conda:
         "envs/nanoplot.yaml"
    shell:
        """
        mkdir -p data/nanoplot_raw/{wildcards.sample}
        NanoPlot -o data/nanoplot_raw/{wildcards.sample} --fastq {input.reads}
        """

rule nanofilt:
    input:
         reads = config["LONG_READ_DIR"] + "/{sample}.fastq.gz"
    output:
        "data/nanofilt/{sample}_nanofilt.fastq.gz"
    conda:
         "envs/nanofilt.yaml"
    params:
        length = 200,
        headcrop = 25,
        quality = 7,
        readtype = "1D"
    shell:
        """
        gunzip -c {input.reads} | NanoFilt --length {params.length} --headcrop {params.headcrop} \
        --quality {params.quality} --readtype {params.readtype} | 
        gzip > data/nanofilt/{wildcards.sample}_nanofilt.fastq.gz
        """

rule nanoplot_filtered:
    input:
         reads = config["LONG_READ_DIR"] + "/{sample}.fastq.gz"
    output:
          "data/nanoplot_filtered/{sample}/NanoStats.txt"
    conda:
         "envs/nanoplot.yaml"
    shell:
        """
        mkdir -p data/nanoplot_filtered/{wildcards.sample}
        NanoPlot -o data/nanoplot_filtered/{wildcards.sample} --fastq {input.reads}
        """

# Filter reads against a reference
rule filter_reads_reference:
    input:
        reads = "data/nanofilt/{sample}_nanofilt.fastq.gz",
        reference_filter = config["REFERENCE_FILTER"]
    output:
        "data/reference_filtered_reads/{sample}_RF.fastq.gz"
    conda:
        "envs/coverm.yaml"
    threads:
         config["MAX_THREADS"]
    shell:
        """
        mkdir -p data/reference_filtered_reads && \
        minimap2 -ax map-ont -t {threads} {input.reference_filter} {input.reads} | samtools fastq -n -f 4 - \
        | gzip > {output}
        """
