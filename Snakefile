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

# Get list of specified viral tools
VIRAL_TOOLS_ASSEMBLY = config["VIRAL_TOOLS"]["ASSEMBLY"].strip().split(",")
VIRAL_TOOLS_READS = config["VIRAL_TOOLS"]["READS"].strip().split(",")

onsuccess:
    print("Workflow finished, no error")

onerror:
    print("An error occurred")

onstart:
    long_reads_dir = config["LONG_READ_DIR"].strip().split()
    if long_reads_dir == "none":
        sys.exit("Need to specify long reads directory")

    valid_assemblers = ["flye", "metaflye","canu", "miniasm","wtdbg2", "raven"]
    valid_viral_tools_assembly = ["virsorter", "vibrant"]
    valid_viral_tools_reads = ["seeker", "virsorter", "vibrant"]
    for assembler in ASSEMBLERS:
        if assembler not in valid_assemblers:
            sys.exit(f'The specified assembler \'{assembler}\' must be one of: ' \
                     f'{", ".join(str(x) for x in valid_assemblers)}')
    for viral_tool in VIRAL_TOOLS_ASSEMBLY:
        if viral_tool not in valid_viral_tools_assembly:
            sys.exit(f'The specified viral tool \'{viral_tool}\' to be run on assemblies must be one of: ' \
                     f'{", ".join(str(x) for x in valid_viral_tools_assembly)}')
    for viral_tool in VIRAL_TOOLS_READS:
        if viral_tool not in valid_viral_tools_reads:
            sys.exit(f'The specified viral tool \'{viral_tool}\' to be run on assembly filtered reads must be one of: ' \
                     f'{", ".join(str(x) for x in valid_viral_tools_assembly)}')


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

include: "rules/qc.smk"
include: "rules/assemble.smk"
include: "rules/assembly_annotate.smk"
include: "rules/viral_assembly_predict.smk"
include: "rules/viral_reads_predict.smk"
include: "rules/viral_annotate.smk"
include: "rules/utils.smk"

rule run_full_pipeline:
    input:
        "data/viral_summary.tsv"
        "finished_viral_annotation",
        "finished_viral_lineage",
        "finished_assembly",
    output:
        touch("finished_full_pipeline")



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


# rule gtdbtk:

# ------------------------------------------------------------------------------------------------
# TODO
#          Fix handling of failed assembly/polishing
#          handle short read - hybrid assembly and/or polishing
#          plasmid identification (PlasFlow/gplas or PlasClass)
#          log failed assemblies
#          annotation : viral and assembly
#              amrfinderplus
#              dram?
#          Generate final summary table(s) and results
#              Summary read stats
#              Assembly stats (N50 etc), checkm, gtdbtk, GUNC
#              Assembly parsnp tree (all assemblies? Apply filtering criteria?)
#              Viral sequences:
#                  quality, cluster membership and cluster representative T/F (note not all sequences are clustered)
#                  Number of genes,
#                  majority lineage (IMG/VR), Abricate/AMRFinderPlus,
#

