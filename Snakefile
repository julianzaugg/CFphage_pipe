import sys
import os
import re
import tempfile
import glob

import pandas as pd

configfile: "config.yaml"

onsuccess:
    print("Workflow finished, no error")

onerror:
    print("An error occurred")

onstart:
    # sample_sheet = config["sample_sheet"].strip().split()
    long_reads_dir = config["LONG_READ_DIR"].strip().split()
    # short_reads_1 = config["short_reads_1"].strip().split()
    # short_reads_2 = config["short_reads_2"].strip().split()
    # if sample_sheet == "none":
        # sys.exit("Provide sample_sheet")
    if long_reads_dir == "none":
        sys.exit("Need to specify long reads directory")

# Load sample_sheet
# sample_df=pd.read_csv (config["sample_sheet"], comment="#", skip_blank_lines=True, sep="\t", index_col=0)

# Get list of samples
# SAMPLES=list(sample_df.index)
LONG_READ_DIR = config["LONG_READ_DIR"]
SHORT_READ_DIR = config["SHORT_READ_DIR"]

# The list of samples to be processed
SAMPLES = glob.glob(f"{LONG_READ_DIR}**/*.fastq.gz",
                    recursive=True)
print(SAMPLES)
SAMPLES = [sample.replace(f"{LONG_READ_DIR}/","").replace(".fastq.gz","")
           for sample
           in SAMPLES]
print(SAMPLES)

# rule pycoqc:
#     input:
#         guppy_summary = config["guppy_summary_file"]
#     output:
#         "data/pycoqc/{sample}_read_qc.html"
#     conda:
#         "envs/pycoqc.yaml"

rule nanoplot_raw_all:
    input:
         expand("data/nanoplot_raw/{sample}/done", sample = SAMPLES)

rule nanoplot_raw:
    input:
         reads = config["LONG_READ_DIR"] + "/{sample}.fastq.gz"
    output:
          "data/nanoplot_raw/{sample}/done"
    conda:
         "envs/nanoplot.yaml"
    shell:
        "mkdir -p data/nanoplot_raw/{wildcards.sample} && " \
        "NanoPlot -o data/nanoplot_raw/{wildcards.sample} --fastq {input.reads} && " \
        "touch {output}"


rule nanofilt_all:
    input:
         expand("data/nanofilt/{sample}_nanofilt.fastq.gz", sample = SAMPLES)

rule nanofilt:
    input:
         reads = config["LONG_READ_DIR"] + "/{sample}.fastq.gz"
    output:
        "data/nanofilt/{sample}_nanofilt.fastq.gz"
    conda:
         "envs/nanofilt.yaml"
    shell:
        "gunzip -c {input.reads} | NanoFilt --length 200 --headcrop 25 --quality 7 --readtype 1D | " \
        "gzip > data/{wildcards.sample}_nanofilt.fastq.gz"

rule nanoplot_filtered_all:
    input:
         expand("data/nanoplot_filtered/{sample}/done", sample = SAMPLES)

rule nanoplot_filtered:
    input:
         reads = config["LONG_READ_DIR"] + "/{sample}.fastq.gz"
    output:
          "data/nanoplot_filtered/{sample}/done"
    conda:
         "envs/nanoplot.yaml"
    shell:
        "mkdir -p data/nanoplot_filtered/{wildcards.sample} && " \
        "NanoPlot -o data/nanoplot_filtered/{wildcards.sample} --fastq {input.reads} && " \
        "touch {output}"

rule flye_all:
    input:
         expand("data/nanofilt/{sample}_nanofilt.fastq.gz", sample = SAMPLES)

rule flye:
    input:
        reads = "data/nanofilt/{sample}_nanofilt.fastq.gz"
    output:
        "data/flye/{sample}/{sample}.flye.fasta"
    conda:
         "envs/flye.yaml"
    params:
        genome_size = config["GENOME_SIZE"],
        flye_parameters = config["FLYE_PARAMS"]
    threads:
        config["MAX_THREADS"]
    shell:
         "mkdir -p data/flye/{wildcards.sample} &&" \
         "flye --nano-raw {input.reads} -o data/flye/{wildcards.sample} " \
         "-g {genome_size} -t {threads} {params.flye_parameters} && " \
         "cp data/flye/{wildcards.sample}/assembly.fasta data/flye/{wildcards.sample}/wildcards.sample.flye.fasta"