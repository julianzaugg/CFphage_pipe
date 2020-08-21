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
# Perform read quality control and evaluation

rule qc:
    input:
        expand("data/nanoplot_raw/{sample}/NanoStats.txt", sample = SAMPLES),
        expand("data/nanofilt/{sample}_nanofilt.fastq.gz", sample = SAMPLES),
        expand("data/nanoplot_filtered/{sample}/NanoStats.txt", sample = SAMPLES)
    output:
        temp(touch("finished_QC"))

# rule pycoqc:
#     input:
#         guppy_summary = config["guppy_summary_file"]
#     output:
#         "data/pycoqc/{sample}_read_qc.html"
#     conda:
#         "envs/pycoqc.yaml"

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

# ------------------------------------------------------------------------------------------------
# Assemble reads

rule assemble:
    input:
        expand("data/assembly/{sample}/{assembler}/{sample}.{assembler}.fasta",
               sample = SAMPLES, assembler = ASSEMBLERS),
        "finished_QC"
    output:
        temp(touch("finished_assembly"))

rule flye:
    input:
        reads = "data/nanofilt/{sample}_nanofilt.fastq.gz"
    output:
        "data/assembly/{sample}/flye/{sample}.flye.fasta"
    conda:
         "envs/flye.yaml"
    params:
        genome_size = config["GENOME_SIZE"],
        flye_parameters = config["FLYE_PARAMS"]
    threads:
        config["MAX_THREADS"]
    message:
        "Assembling {wildcards.sample} with flye"
    shell:
         """
         mkdir -p data/assembly/{wildcards.sample}/flye
         flye --nano-raw {input.reads} -o data/assembly/{wildcards.sample}/flye \
         -g {params.genome_size} -t {threads} {params.flye_parameters}
         cp data/assembly/{wildcards.sample}/flye/assembly.fasta \
         data/assembly/{wildcards.sample}/flye/{wildcards.sample}.flye.fasta
         """

rule canu:
    input:
        reads = "data/nanofilt/{sample}_nanofilt.fastq.gz"
    output:
        "data/assembly/{sample}/canu/{sample}.canu.fasta"
    params:
        genome_size = config["GENOME_SIZE"],
        min_read_length=1000,
        min_overlap_length=500,
        min_input_coverage=0,
        stop_on_low_coverage=0,
        corrected_error_rate=0.105,
        use_grid="false",
        max_memory=config["MAX_MEMORY"]
    conda:
         "envs/canu.yaml"
    threads:
        config["MAX_THREADS"]
    message:
        "Assembling {wildcards.sample} with canu"
    shell:
        """
        mkdir -p data/assembly/{wildcards.sample}/canu
        canu -p {wildcards.sample} -d data/assembly/{wildcards.sample}/canu/ \
        -nanopore-raw {input.reads} genomeSize={params.genome_size} -maxMemory={params.max_memory} -maxThreads={threads} \
        -corThreads={threads} -useGrid={params.use_grid} correctedErrorRate={params.corrected_error_rate} \ 
        -minReadLength={params.min_read_length} -minOverlapLength={params.min_overlap_length} \
        -minInputCoverage {params.min_input_coverage} -stopOnLowCoverage={params.stop_on_low_coverage} -fast
        cp data/assembly/{wildcards.sample}/canu/assembly.fasta \
        data/assembly/{wildcards.sample}/canu/{wildcards.sample}.canu.fasta
        """

rule raven:
    input:
        reads = "data/nanofilt/{sample}_nanofilt.fastq.gz"
    output:
        "data/assembly/{sample}/raven/{sample}.raven.fasta"
    conda:
         "envs/raven.yaml"
    params:
        polishing_rounds = config["RACON_ROUNDS"]
    threads:
        config["MAX_THREADS"]
    message:
        "Assembling {wildcards.sample} with raven"
    shell:
        """
        mkdir -p data/assembly/{wildcards.sample}/raven
        raven --polishing-rounds {params.polishing_rounds} --threads {threads} {input.reads} > {output}
        """

rule wtdbg2: # also known as redbean
    input:
        reads = "data/nanofilt/{sample}_nanofilt.fastq.gz"
    output:
        "data/assembly/{sample}/wtdbg2/{sample}.wtdbg2.fasta"
    conda:
         "envs/wtdbg2.yaml"
    params:
        genome_size = config["GENOME_SIZE"]
    threads:
        config["MAX_THREADS"]
    message:
        "Assembling {wildcards.sample} with wtdbg2"
    shell:
        """
        mkdir -p data/assembly/{wildcards.sample}/wtdbg2
        wtdbg2 -x ont -g {params.genome_size} -t {threads} -i {input.reads} -f \
        -o data/assembly/{wildcards.sample}/wtdbg2/{wildcards.sample}.wtdbg2
        wtpoa-cns -t {threads} \
        -i data/assembly/{wildcards.sample}/wtdbg2/{wildcards.sample}.wtdbg2.ctg.lay.gz -fo {output}
        """

rule miniasm:
    input:
        reads = "data/nanofilt/{sample}_nanofilt.fastq.gz"
    output:
        "data/assembly/{sample}/miniasm/{sample}.miniasm.fasta"
    conda:
         "envs/miniasm.yaml"
    params:
        genome_size = config["GENOME_SIZE"]
    threads:
        config["MAX_THREADS"]
    message:
        "Assembling {wildcards.sample} with miniasm"
    shell:
        """
        mkdir -p data/assembly/{wildcards.sample}/miniasm
        minimap2 -t {threads} -x ava-ont {input.reads} {input.reads} \
        > data/assembly/{wildcards.sample}/miniasm/{wildcards.sample}.reads.paf.gz
        miniasm -f {input.reads} data/assembly/{wildcards.sample}/miniasm/{wildcards.sample}.reads.paf.gz \
        > data/assembly/{wildcards.sample}/miniasm/{wildcards.sample}.gfa
        awk '$1 ~/S/ {{print ">"$2"\\n"$3}}' data/assembly/{wildcards.sample}/miniasm/{wildcards.sample}.gfa > {output}
        """
# ------------------------------------------------------------------------------------------------
# Polish assemblies

rule polish:
    input:
        expand("data/polishing/{sample}/medaka/{assembler}/{sample}.{assembler}.medaka.fasta",
               sample = SAMPLES, assembler = ASSEMBLERS),
        "finished_assembly"
    output:
        temp(touch("finished_polishing"))

rule racon_polish:
    input:
         reads = config["LONG_READ_DIR"] + "/{sample}.fastq.gz",
         assembly = "data/assembly/{sample}/{assembler}/{sample}.{assembler}.fasta"
    output:
          "data/polishing/{sample}/racon/{assembler}/{sample}.{assembler}.racon.fasta"
    conda:
         "envs/racon.yaml"
    threads:
        config["MAX_THREADS"]
    message:
        "Polishing {input.assembly} with Racon"
    params:
        rounds = 4,
        match = 8,
        mismatch = -6,
        gap = -8,
        window_length = 500
    script:
        "scripts/racon_polish.py"

rule medaka_polish:
    input:
         reads = config["LONG_READ_DIR"] + "/{sample}.fastq.gz",
         assembly = "data/polishing/{sample}/racon/{assembler}/{sample}.{assembler}.racon.fasta"
    output:
          "data/polishing/{sample}/medaka/{assembler}/{sample}.{assembler}.medaka.fasta"
    conda:
         "envs/medaka.yaml"
    threads:
        config["MAX_THREADS"]
    message:
        "Polishing {input.assembly} with Medaka"
    params:
        guppy_model = config["GUPPY_MODEL"]
    shell:
        """
        mkdir -p data/polishing/{wildcards.sample}/medaka/{wildcards.assembler}
        medaka_consensus -d {input.assembly} -i {input.reads} -t {threads} \
        -o data/polishing/{wildcards.sample}/medaka/{wildcards.assembler} -m {params.guppy_model}
        
        cp data/polishing/{wildcards.sample}/medaka/{wildcards.assembler}/consensus.fasta \
        data/polishing/{wildcards.sample}/medaka/{wildcards.assembler}/{wildcards.sample}.{wildcards.assembler}.medaka.fasta
        """

# ------------------------------------------------------------------------------------------------
# Coverage
#   - Of provided reference genomes
#   - Of assemblies


# For all samples map to all reference genomes
# rule coverage_references:
    # input:
        # expand("data/coverage/references_coverage",
               # sample = SAMPLES, assembler = ASSEMBLERS)

# For each sample (reads) map to all assemblies produced (medaka) for that sample
# {sample}_mean_coverage.tsv
# {sample}_read_counts.tsv
# {sample}_relative_abundance.tsv (genome)
# {sample}_covered_fraction.tsv
# rule coverage_assemblies:
# ------------------------------------------------------------------------------------------------

# TODO
#          checkm, gtdbtk, busco?
#          coverm/to_{reference}
#          coverm/to_assembly
#          kaiju (profile reads)
#          compile read stats
#          virsorter2, vibrant, checkv (separate snakemake?)
#          plasmid identification (PlasFlow/gplas or PlasClass)
#          log failed assemblys
