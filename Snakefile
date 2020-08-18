import sys
import os
import re
import tempfile
import glob

import pandas as pd

configfile: "config.yaml"


# Get list of samples
LONG_READ_DIR = config["LONG_READ_DIR"]
SHORT_READ_DIR = config["SHORT_READ_DIR"]

# Get list of specified assemblers
ASSEMBLERS = config["ASSEMBLERS"].strip().split(",")

onsuccess:
    print("Workflow finished, no error")

onerror:
    print("An error occurred")

onstart:
    # sample_sheet = config["sample_sheet"].strip().split()
    long_reads_dir = config["LONG_READ_DIR"].strip().split()
    # if sample_sheet == "none":
        # sys.exit("Provide sample_sheet")
    if long_reads_dir == "none":
        sys.exit("Need to specify long reads directory")
    valid_assemblers = ["flye", "canu", "miniasm","wtdbg2", "raven"]
    for assembler in ASSEMBLERS:
        if assembler not in valid_assemblers:
            sys.exit(f'The specified assembler \'{assembler}\' must be one of: ' \
                     f'{", ".join(str(x) for x in valid_assemblers)}')

# rule test:
    # run: print("test")

# Load sample_sheet
# sample_df=pd.read_csv (config["sample_sheet"], comment="#", skip_blank_lines=True, sep="\t", index_col=0)


# The list of samples to be processed
# SAMPLES=list(sample_df.index)
SAMPLES = glob.glob(f"{LONG_READ_DIR}**/*.fastq.gz",
                    recursive=True)
SAMPLES = [sample.replace(f"{LONG_READ_DIR}/","").replace(".fastq.gz","")
           for sample
           in SAMPLES]

# ------------------------------------------------------------------------------------------------
# Read quality control and evaluation
rule qc:
    input:
        expand("data/nanoplot_raw/{sample}/NanoStats.txt", sample = SAMPLES),
        expand("data/nanofilt/{sample}_nanofilt.fastq.gz", sample = SAMPLES),
        expand("data/nanoplot_filtered/{sample}/NanoStats.txt", sample = SAMPLES)
        # report= "reports/QC_report.html",
    output:
        temp(touch("finished_QC"))

# rule pycoqc:
#     input:
#         guppy_summary = config["guppy_summary_file"]
#     output:
#         "data/pycoqc/{sample}_read_qc.html"
#     conda:
#         "envs/pycoqc.yaml"

# rule nanoplot_raw_all:
    # input:
         # expand("data/nanoplot_raw/{sample}/NanoStats.txt", sample = SAMPLES)

rule nanoplot_raw:
    input:
         reads = config["LONG_READ_DIR"] + "/{sample}.fastq.gz"
    output:
          "data/nanoplot_raw/{sample}/NanoStats.txt"
    conda:
         "envs/nanoplot.yaml"
    shell:
        "mkdir -p data/nanoplot_raw/{wildcards.sample} && " \
        "NanoPlot -o data/nanoplot_raw/{wildcards.sample} --fastq {input.reads}"

# rule nanofilt_all:
    # input:
    #      expand("data/nanofilt/{sample}_nanofilt.fastq.gz", sample = SAMPLES)

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

# rule nanoplot_filtered_all:
#     input:
#          expand("data/nanoplot_filtered/{sample}/NanoStats.txt", sample = SAMPLES)

rule nanoplot_filtered:
    input:
         reads = config["LONG_READ_DIR"] + "/{sample}.fastq.gz"
    output:
          "data/nanoplot_filtered/{sample}/NanoStats.txt"
    conda:
         "envs/nanoplot.yaml"
    shell:
        "mkdir -p data/nanoplot_filtered/{wildcards.sample} && " \
        "NanoPlot -o data/nanoplot_filtered/{wildcards.sample} --fastq {input.reads}"

# ------------------------------------------------------------------------------------------------
# Assembly
rule assembly:
    input:
        expand("data/assembly/{sample}/{assembler}/finished_assembly", sample = SAMPLES, assembler = ASSEMBLERS),
        "finished_QC"
    output:
        temp(touch("finished_assembly"))


# rule flye_all:
#     input:
#          expand("data/assembly/{sample}/flye/{sample}.flye.fasta", sample = SAMPLES)

rule flye:
    input:
        reads = "data/nanofilt/{sample}_nanofilt.fastq.gz"
    output:
        "data/assembly/{sample}/flye/{sample}.flye.fasta",
        temp(touch("data/assembly/{sample}/flye/finished_assembly"))
    conda:
         "envs/flye.yaml"
    params:
        genome_size = config["GENOME_SIZE"],
        flye_parameters = config["FLYE_PARAMS"]
    threads:
        config["MAX_THREADS"]
    shell:
         """
         mkdir -p data/assembly/{wildcards.sample}/flye
         flye --nano-raw {input.reads} -o data/{wildcards.sample}/flye \
         -g {params.genome_size} -t {threads} {params.flye_parameters}
         cp data/assembly/{wildcards.sample}/flye/assembly.fasta \
         data/assembly/{wildcards.sample}/flye/{wildcards.sample}.flye.fasta
         """

rule canu:
    input:
        reads = "data/nanofilt/{sample}_nanofilt.fastq.gz"
    output:
        "data/assembly/{sample}/canu/{sample}.canu.fasta",
        temp(touch("data/assembly/{sample}/canu/finished_assembly"))
    params:
        genome_size = config["GENOME_SIZE"],
        min_read_length=500,
        min_overlap_length=200,
        min_input_coverage=0,
        stop_on_low_coverage=0,
        max_memory=config["MAX_MEMORY"]
    conda:
         "envs/canu.yaml"
    threads:
        config["MAX_THREADS"]
    shell:
        """
        mkdir -p data/assembly/{wildcards.sample}/canu
        canu -p {wildcards.sample} -d data/assembly/{wildcards.sample}/canu/ \
        -nanopore {input.reads} genomeSize={params.genome_size} -minReadLength={params.min_read_length} \
        -minOverlapLength={params.min_overlap_length} -minInputCoverage {params.min_input_coverage} \
        -stopOnLowCoverage={params.stop_on_low_coverage} -maxMemory={params.max_memory} -maxThreads={threads}
        cp data/assembly/{wildcards.sample}/canu/assembly.fasta \
        data/assembly/{wildcards.sample}/canu/{wildcards.sample}.canu.fasta
        """

# ------------------------------------------------------------------------------------------------
# Polishing of assemblies

rule polishing:
    input:
        expand("data/polishing/{sample}/medaka/{sample}.{assembler}.medaka.fasta",
               sample = SAMPLES, assembler = ASSEMBLERS),
        "finished_assembly"
    output:
        temp(touch("finished_polishing"))

# rule polish_racon_all:
#     input:
#          expand("data/polishing/racon/{sample}/{sample}.{assembler}.racon.fasta",
#                 sample = SAMPLES)

rule racon_polish:
    input:
         reads = config["LONG_READ_DIR"] + "/{sample}.fastq.gz",
         assembly = "data/assembly/{sample}/{assembler}/{sample}.{assembler}.fasta"
    output:
          "data/polishing/{sample}/racon/{sample}.{assembler}.racon.fasta"
    conda:
         "envs/racon.yaml"
    threads:
        config["MAX_THREADS"]
    message:
        "Polishing {input.assembly} with Racon"
    params:
        rounds = 3
    script:
        "scripts/racon_polish.py"

rule medaka_polish:
    input:
         reads = config["LONG_READ_DIR"] + "/{sample}.fastq.gz",
         assembly = "data/polishing/{sample}/racon/{sample}.{assembler}.racon.fasta"
    output:
          "data/polishing/{sample}/medaka/{sample}.{assembler}.medaka.fasta"
    conda:
         "envs/medaka.yaml"
    threads:
        config["MAX_THREADS"]
    message:
        "Polishing {input.assembly} with Medaka"
    params:
        guppy_model = config["GUPPY_MODEL"]
    shell:
        "medaka_consensus -d {input.assembly} -i {input.reads} -t {threads} "
        "-o data/polishing/{wildcards.sample}/medaka -m {params.guppy_model} && " \
        "cp data/polishing/{wildcards.sample}/medaka/consensus.fasta " \
        "data/polishing/{wildcards.sample}/medaka/{wildcards.sample}.medaka.fasta"


# ------------------------------------------------------------------------------------------------