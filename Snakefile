from snakemake.utils import logger, min_version
min_version("6.0")

import sys
import os
import re
import tempfile
import glob
import subprocess
import pandas as pd


sys.path.append(
    os.path.join(os.path.dirname(os.path.abspath(workflow.snakefile)), "scripts")
)

configfile: "config.yaml"

# Get list of samples
LONG_READ_DIR = config["LONG_READ_DIR"]
SHORT_READ_DIR = config["SHORT_READ_DIR"]

# Get list of specified assemblers
ASSEMBLERS = config["ASSEMBLERS"].strip().split(",")
VIRAL_TOOLS = config["VIRAL_TOOLS"].strip().split(",")


onsuccess:
    print("Workflow finished, no error")

onerror:
    print("An error occurred")

onstart:
    long_reads_dir = config["LONG_READ_DIR"].strip().split()
    if long_reads_dir == "none":
        sys.exit("Need to specify long reads directory")

    valid_assemblers = ["flye", "metaflye","canu", "miniasm","wtdbg2", "raven"]
    for assembler in ASSEMBLERS:
        if assembler not in valid_assemblers:
            sys.exit(f'The specified assembler \'{assembler}\' must be one of: ' \
                     f'{", ".join(str(x) for x in valid_assemblers)}')


# Load sample_sheet
# sample_df=pd.read_csv (config["sample_sheet"], comment="#", skip_blank_lines=True, sep="\t", index_col=0)

# The list of samples to be processed
# SAMPLES=list(sample_df.index)
SAMPLES = glob.glob(f"{LONG_READ_DIR}**/*.fastq.gz",
                    recursive=True)
SAMPLES = [sample.replace(f"{LONG_READ_DIR}/","").replace(".fastq.gz","")
           for sample
           in SAMPLES]
# SAMPLES = glob_wildcards("%s/{id}.fastq.gz".format(LONG_READ_DIR))

# List of reference genome files
REFERENCE_GENOMES = glob.glob(f"{config['REFERENCE_GENOMES_DIR']}**/*.fasta",
                    recursive=True)
REFERENCE_GENOMES = [reference_genome.replace(f"{config['REFERENCE_GENOMES_DIR']}/","").replace(".fasta","")
           for reference_genome
           in REFERENCE_GENOMES]

ABSOLUTE_DATA_PATH = os.getcwd()
SNAKE_PATH= workflow.basedir

module qc:
    snakefile: "rules/qc.smk"
    config: config

module assemble:
    snakefile: "rules/assemble.smk"
    config: config

module assembly_annotate:
    snakefile: "rules/assembly_annotate.smk"
    config: config


module viral_predict:
    snakefile: "rules/viral_predict.smk"
    config: config

module utils:
    snakefile: "rules/utils.smk"
    config: config


use rule * from qc as *
use rule * from assembly as *
use rule * from assembly_annotate as *
use rule * from viral_predict as *
use rule * from utils as *


rule qc:
    input:
        expand("data/nanoplot_raw/{sample}/NanoStats.txt", sample = SAMPLES),
        expand("data/nanofilt/{sample}_nanofilt.fastq.gz", sample = SAMPLES),
        expand("data/nanoplot_filtered/{sample}/NanoStats.txt", sample = SAMPLES)
    output:
        touch("finished_QC")

rule assemble:
    input:
        expand("data/assembly/{sample}/{assembler}/{sample}.{assembler}.fasta",
               sample = SAMPLES, assembler = ASSEMBLERS),
        "finished_QC"
    output:
        touch("finished_assembly")

rule viral_annotate:
    input:
        "finished_viral_clustering",
        "data/viral_annotation/prodigal/done",
        "data/viral_annotation/abricate/done"
    output:
        touch("finished_viral_annotation")


rule assembly_annotate:
    input:
        expand("data/assembly_annotation/prodigal/{sample}/{assembler}/{sample}.{assembler}.faa",
        sample = SAMPLES, assembler = ASSEMBLERS),
        expand("data/assembly_annotation/abricate/{sample}/{assembler}/{sample}.{assembler}.{db}.tsv",
           sample = SAMPLES, assembler = ASSEMBLERS, db = ["card", "vfdb"]),
        "finished_polishing"
    output:
        touch("finished_assembly_annotation")

# ------------------------------------------------------------------------------------------------
# Coverage
#   - Of provided reference genomes
#   - Of assemblies

rule coverage_reference_genomes_all:
    input:
        expand("data/coverage/reference_genomes/{reference_genome}/{reference_genome}_coverage_table.tsv",
               reference_genome = REFERENCE_GENOMES)

def get_coverm_reference_params(wildcards):
    """
    Return coverM and other coverage calculation parameters for a reference genome.
    Will use parameters from config file if provided, else default.
    """
    param_dict = {
        "min_read_percent_identity" : 0.9,
        "min_read_aligned_percent" : 0.75,
        "multiple_genomes" : False
    }
    if wildcards in config["REFERENCE_GENOMES_PARAMS"]["COVERM_PARAMS"]:
        param_dict.update(config["REFERENCE_GENOMES_PARAMS"]["COVERM_PARAMS"][wildcards])

    return(param_dict)

rule coverage_reference_genomes:
    input:
        reference_fasta = config["REFERENCE_GENOMES_DIR"] + "/{reference_genome}.fasta",
        qc_complete = "finished_QC"
    output:
        "data/coverage/reference_genomes/{reference_genome}/{reference_genome}_coverage_table.tsv"
    conda:
        "envs/coverm.yaml"
    params:
        read_files = "data/nanofilt/*.fastq.gz",
        reference_coverm_parameters_dict = lambda wildcards: get_coverm_reference_params(wildcards.reference_genome)
    threads:
        config["MAX_THREADS"]
    message:
        "Mapping reads to {input.reference_fasta} with CoverM"
    script:
        "scripts/get_coverage_reference.py"


# Calculate coverage of all assemblies generated for a sample
# Combine all the assemblies into single, one-entry-per-genome FASTA file?
# data/final_assemblies/{sample}_assemblies.fasta
# data/final_assemblies/{sample}.{assembler}.medaka.fasta
# rule coverage_assemblies_all:
#     input:
#         expand("data/coverage/{sample}/{sample}_coverage_table.tsv", sample = SAMPLES),
#
# rule coverage_assemblies:
#     input:
#         reads = "data/nanofilt/{sample}_nanofilt.fastq.gz",
#         assembly = "data/polishing/{sample}/medaka/{assembler}/{sample}.{assembler}.medaka.fasta"
#     output:
#         "data/coverage/{sample}/{sample}_coverage_table.tsv"
#     threads:
#         config["MAX_THREADS"]
#     message:
#         "Mapping reads with CoverM to assemblies generated from {wildcards.sample}"


# ------------------------------------------------------------------------------------------------
# Assembly evaluation

# CheckM is used to evaluate the quality of the assemblies
# rule checkm:
#     input:
#         assemblies="data/polishing/{sample}/medaka/{assembler}/{sample}.{assembler}"
#     output:
#         "data/checkm.out"
#     priority: 1
#     conda:
#         "envs/checkm.yaml"
#     threads:
#         config["max_threads"]
#     shell:
#         'checkm lineage_wf -t {threads} -x fa data/das_tool_bins/das_tool_DASTool_bins data/checkm > data/checkm.out'
# rule polish:
#     input:
#         expand("data/polishing/{sample}/medaka/{assembler}/{sample}.{assembler}.medaka.fasta",
#                sample = SAMPLES, assembler = ASSEMBLERS),
#         "finished_assembly"
#     output:
#         temp(touch("finished_polishing"))

# rule gtdbtk:

# ------------------------------------------------------------------------------------------------
# TODO
#          Fix handling of failed assembly/polishing
#          checkm, gtdbtk, busco?
#          Unassembled contigs > polish (viral not making it into assembly)
#          coverm/to_assembly
#          kraken / bracken (or alterantive) to profile qc reads. Can map to GTDB though slow.
#          compile read stats
#          seeker - input needs to be > 200bp
#          a) virsorter2 / vibrant (phage identication) DONE
#          b) checkv  (identify likely true phage) DONE
#          c) fastANI + MCL (dereplication of phage) DONE
#          plasmid identification (PlasFlow/gplas or PlasClass)
#          log failed assemblies
#          min contig size for assemblies, 2000bp?
#          assembly annotation : prodigal >
#              abricate (VFDB, CARD, etc.), amrfinderplus
#          viral annotation : prodigal >
#              abricate (VFDB, CARD, etc.), BLAST against img/vr (majority tax)

