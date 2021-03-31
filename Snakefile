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
VIRAL_TOOLS = config["VIRAL_TOOLS"].strip().split(",")


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

    valid_assemblers = ["flye", "metaflye","canu", "miniasm","wtdbg2", "raven"]
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
# SAMPLES = glob_wildcards("%s/{id}.fastq.gz".format(LONG_READ_DIR))

# List of reference genome files
REFERENCE_GENOMES = glob.glob(f"{config['REFERENCE_GENOMES_DIR']}**/*.fasta",
                    recursive=True)
REFERENCE_GENOMES = [reference_genome.replace(f"{config['REFERENCE_GENOMES_DIR']}/","").replace(".fasta","")
           for reference_genome
           in REFERENCE_GENOMES]

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
        > {output}
        """

# ------------------------------------------------------------------------------------------------

rule reads2fasta:
    input:
        reads = "data/nanofilt/{sample}_nanofilt.fastq.gz"
    output:
        "data/read_fastas/{sample}.fasta"
    conda:
        "envs/seqtk.yaml"
    shell:
        """
        mkdir -p data/read_fastas
        seqk seq -a {input.reads} > {output}
        """

# ------------------------------------------------------------------------------------------------
# Run viral tools on reads

rule viral_reads_predict:
    input:
        expand("data/viral_reads_predict/{sample}/virsorter/done", sample = SAMPLES),
        "finished_QC"
    output:
        touch("finished_viral_reads_predict")

rule virsorter_reads:
    input:
        reads = "data/nanofilt/{sample}_nanofilt.fastq.gz"
    params:
        virsorter_database = config["VIRSORTER"]["DATABASE_DIR"],
        virsorter_min_length = config["VIRSORTER"]["MIN_LENGTH"]
    output:
        "data/viral_reads_predict/{sample}/virsorter/done"
    conda:
        "envs/virsorter.yaml"
    threads:
        config["MAX_THREADS"]
    shell:
        """
        mkdir -p data/viral_reads_predict/{wildcards.sample} && \
        virsorter run \
        --rm-tmpdir \
        --seqfile {input.reads} \
        --working-dir data/viral_reads_predict/{wildcards.sample}/virsorter \
        --db-dir {params.virsorter_database} \
        --min-length {params.virsorter_min_length} \
        --jobs {threads} \
        all
        touch data/viral_reads_predict/{wildcards.sample}/virsorter/done
        """

# ------------------------------------------------------------------------------------------------
# Assemble reads

rule assemble:
    input:
        expand("data/assembly/{sample}/{assembler}/{sample}.{assembler}.fasta",
               sample = SAMPLES, assembler = ASSEMBLERS),
        "finished_QC"
    output:
        touch("finished_assembly")

rule flye:
    input:
        reads = "data/nanofilt/{sample}_nanofilt.fastq.gz"
    output:
        "data/assembly/{sample}/flye/{sample}.flye.fasta"
    conda:
         "envs/flye.yaml"
    params:
        genome_size = config["GENOME_SIZE"],
        flye_parameters = config["ASSEMBLER_PARAMS"]["FLYE_PARAMS"]
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
         cp data/assembly/{wildcards.sample}/flye/assembly_graph.gfa \
         data/assembly/{wildcards.sample}/flye/{wildcards.sample}.flye.gfa
         """

rule metaflye:
    input:
        reads = "data/nanofilt/{sample}_nanofilt.fastq.gz"
    output:
        "data/assembly/{sample}/metaflye/{sample}.metaflye.fasta"
    conda:
         "envs/flye.yaml"
    params:
        genome_size = config["GENOME_SIZE"],
        metaflye_parameters = config["ASSEMBLER_PARAMS"]["METAFLYE_PARAMS"]
    threads:
        config["MAX_THREADS"]
    message:
        "Assembling {wildcards.sample} with metaflye"
    shell:
         """
         mkdir -p data/assembly/{wildcards.sample}/metaflye
         flye --meta --nano-raw {input.reads} -o data/assembly/{wildcards.sample}/metaflye \
         -g {params.genome_size} -t {threads} {params.metaflye_parameters}
         cp data/assembly/{wildcards.sample}/metaflye/assembly.fasta \
         data/assembly/{wildcards.sample}/metaflye/{wildcards.sample}.metaflye.fasta
         cp data/assembly/{wildcards.sample}/metaflye/assembly_graph.gfa \
         data/assembly/{wildcards.sample}/metaflye/{wildcards.sample}.metaflye.gfa
         """

rule canu:
    input:
        reads = "data/nanofilt/{sample}_nanofilt.fastq.gz"
    output:
        "data/assembly/{sample}/canu/{sample}.canu.fasta"
    params:
        genome_size = config["GENOME_SIZE"],
        min_read_length=2000,
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
        canu -fast -p {wildcards.sample} -d data/assembly/{wildcards.sample}/canu/ \
        -nanopore {input.reads} genomeSize={params.genome_size} -maxMemory={params.max_memory} -maxThreads={threads} \
        -corThreads={threads} -useGrid={params.use_grid} -correctedErrorRate={params.corrected_error_rate} \
        -minReadLength={params.min_read_length} -minOverlapLength={params.min_overlap_length} \
        -minInputCoverage={params.min_input_coverage} -stopOnLowCoverage={params.stop_on_low_coverage}
        cp data/assembly/{wildcards.sample}/canu/{wildcards.sample}.contigs.fasta \
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
        raven \
        --polishing-rounds {params.polishing_rounds} \
        --graphical-fragment-assembly data/assembly/{wildcards.sample}/raven/{wildcards.sample}.raven.gfa \
        --threads {threads} {input.reads} > {output}
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
        > data/assembly/{wildcards.sample}/miniasm/{wildcards.sample}.miniasm.gfa
        awk '$1 ~/S/ {{print ">"$2"\\n"$3}}' data/assembly/{wildcards.sample}/miniasm/{wildcards.sample}.miniasm.gfa > {output}
        """
# ------------------------------------------------------------------------------------------------
# Polish assemblies
# TODO - add/replace with HYPO or Nextpolish

rule polish:
    input:
        expand("data/polishing/{sample}/medaka/{assembler}/{sample}.{assembler}.medaka.fasta",
               sample = SAMPLES, assembler = ASSEMBLERS),
        "finished_assembly"
    output:
        touch("finished_polishing")

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
# Run viral tools on polished assemblies

rule viral_assembly_predict:
    input:
        expand("data/viral_assembly_predict/{sample}/{viral_predict_tool}/{assembler}/done",
            sample = SAMPLES,
            assembler = ASSEMBLERS,
            viral_predict_tool = VIRAL_TOOLS),
        "data/checkv/done",
        "finished_polishing"
    output:
        touch("finished_viral_assembly_predict")

rule virsorter_assembly:
    input:
        assembly = "data/polishing/{sample}/medaka/{assembler}/{sample}.{assembler}.medaka.fasta"
    params:
        virsorter_database = config["VIRSORTER"]["DATABASE_DIR"],
        virsorter_min_length = config["VIRSORTER"]["MIN_LENGTH"]
    output:
        "data/viral_assembly_predict/{sample}/virsorter/{assembler}/done"
    message:
        "Running virsorter on {input.assembly}"
    conda:
        "envs/virsorter.yaml"
    threads:
        config["MAX_THREADS"]
    shell:
        """
        virsorter_sample_assembler_base_path="data/viral_assembly_predict/{wildcards.sample}/virsorter/{wildcards.assembler}"
        
        mkdir -p data/viral_assembly_predict/{wildcards.sample}/virsorter/{wildcards.assembler} && \
        virsorter run \
        --rm-tmpdir \
        --seqfile {input.assembly} \
        --working-dir $virsorter_sample_assembler_base_path \
        --db-dir {params.virsorter_database} \
        --min-length {params.virsorter_min_length} \
        --jobs {threads} \
        all && \
        touch $virsorter_sample_assembler_base_path/done
        
        if [[ -f $virsorter_sample_assembler_base_path/final-viral-combined.fa ]]; then
            cp $virsorter_sample_assembler_base_path/final-viral-combined.fa \
            $virsorter_sample_assembler_base_path/{wildcards.sample}.{wildcards.assembler}.virsorter.fasta
            sed -i "s/>/>{wildcards.sample}__{wildcards.assembler}__virsorter____/g" \
            $virsorter_sample_assembler_base_path/{wildcards.sample}.{wildcards.assembler}.virsorter.fasta
        fi 
        """

def collect_viral_outputs(wildcards):
    files = expand("data/viral_assembly_predict/{sample}/{viral_predict_tool}/{assembler}/"
                   "{sample}.{assembler}.{viral_predict_tool}.fasta",
        sample = SAMPLES,
        assembler = ASSEMBLERS,
        viral_predict_tool = VIRAL_TOOLS)
    files = [file for file in files if os.path.isfile(file)]
    return files

# Run checkV on predicted viral sequences
rule checkv_assembly:
    input:
        viral_tool_output = collect_viral_outputs
    output:
        checkv_selected = "data/checkv/checkv_selected.fasta",
        done = "data/checkv/done"
    message:
        "Running checkv"
    params:
        checkv_db = config["CHECKV_DB"]
    conda:
        "envs/checkv.yaml"
    threads:
        config["MAX_THREADS"]
    shell:
        """
        mkdir -p data/checkv
        if [[ -f {output.checkv_selected} ]];then
            rm data/checkv/checkv_selected.fasta
        fi
        
        cat {input.viral_tool_output} \
        > data/checkv/all_samples_viral_sequences.fasta
        
        checkv end_to_end \
        -d {params.checkv_db} \
        -t {threads} \
        --restart \
        data/checkv/all_samples_viral_sequences.fasta \
        data/checkv/
        
        # Combine viruses.fna and proviruses.fna, assume they exist
        cat data/checkv/viruses.fna data/checkv/proviruses.fna > data/checkv/checkv_all.fasta
        
        # Grab all Medium and High quality, and Complete, viral genomes and write to fasta file 
        while read contig; do
        grep -h $contig data/checkv/checkv_all.fasta -A 1 >> {output.checkv_selected}
        done < <(awk -F "\t" '$8~/(Complete|[Medium,High]-quality)$/{{print $1}}' data/checkv/quality_summary.tsv)
        touch {output.done}
        """

rule viral_cluster:
    input:
        "data/viral_clustering/fastani/fastani_viral.tsv",
        "data/viral_clustering/mcl/mcl_viral_clusters.tsv"
    output:
        touch("finished_viral_clustering")

# Run FastANI on selected viral sequences
rule fastani_viral:
    input:
        checkv_selected = "data/checkv/checkv_selected.fasta"
    output:
        "data/viral_clustering/fastani/fastani_viral.tsv"
    message:
        "Clustering selected viral sequences"
    conda:
        "envs/fastani.yaml"
    params:
        frag_len = 600,
        min_fraction = 0.85
    threads:
        config["MAX_THREADS"]
    script:
        "scripts/viral_fastani.py"

# Takes the raw output from FastANI and calculates average for each bidirectional pair of genomes
rule fastani_average:
    input:
        "data/viral_clustering/fastani/fastani_viral.tsv"
    output:
        fastani_viral_ani95_mcl = "data/viral_clustering/fastani/fastani_viral_ani95_mcl.tsv",
        done = "data/viral_clustering/fastani/done"
    conda:
        "envs/fastani_average.yaml"
    shell:
        """
        Rscript --vanilla scripts/fastani_average.R {input} \
        data/viral_clustering/fastani/fastani_viral_mean.tsv \
        {output.fastani_viral_ani95_mcl}
        touch {output.done}
        """

rule mcl:
    input:
        "data/viral_clustering/fastani/fastani_viral_ani95_mcl.tsv"
    output:
        mcl_viral_clusters = "data/viral_clustering/mcl/mcl_viral_clusters.tsv",
        done = "data/viral_clustering/mcl/done"
    params:
        inflation_factor = 3.5
    conda:
        "envs/mcl.yaml"
    shell:
        """
        BASE_DIR="data/viral_clustering/mcl"
        mkdir -p $BASE_DIR
        
        cd $BASE_DIR
        mcxload -abc $1 --stream-mirror -write-tab viral.tab -o viral.mci
        mcl viral.mci -I {params.inflation_factor}
        mcxdump -icl out.viral.mci.I -tabr viral.tab -o dump.viral.mci.I

        sed 's/\t/,/g' dump.viral.mci.I | cat -n | sed -e 's/^[ \t]*//' | sed 's/^/cluster_/g' \
        > {output.mcl_viral_clusters}
        touch {output.done}
        """

# ------------------------------------------------------------------------------------------------
# Circularise assemblies (if possible)
rule circularise:
    input:
        expand("data/circlator/{sample}/{assembler}/{sample}.{assembler}.final.fasta",
               sample = SAMPLES, assembler = ASSEMBLERS)
    output:
        temp(touch("finished_circularising"))

rule circlator:
    input:
        assembly = "data/polishing/{sample}/medaka/{assembler}/{sample}.{assembler}.medaka.fasta",
        reads = config["LONG_READ_DIR"] + "/{sample}.fastq.gz"
    output:
        "data/circlator/{sample}/{assembler}/{sample}.{assembler}.final.fasta"
    threads:
        config["MAX_THREADS"]
    conda:
        "envs/circlator.yaml"
    message:
        "Attempting to circularise {input.assembly} with Circlator"
    shell:
        """
        circlator all \
        --threads {threads} \
        --verbose \
        --b2r_min_read_length 4000 \
        --b2r_length_cutoff 30000 \
        --b2r_discard_unmapped \
        --bwa_opts "-x ont2d" \
        --data_type nanopore-raw \
        {input.assembly} {input.reads} data/circlator/temp
         
        if [ -f data/circlator/temp/06.fixstart.fasta ]
        then
            cp data/circlator/temp/06.fixstart.fasta {output}
        else
            cp {input.assembly} {output}
        fi
        mv data/circlator/temp/* data/circlator/{wildcards.sample}/{wildcards.assembler}/
        rm -r data/circlator/temp
        """

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
#          checkm, gtdbtk, busco?
#          Unassembled contigs > polish (viral not making it into assembly)
#          coverm/to_assembly
#          kraken / bracken to profile qc reads. Can map to GTDB though slow.
#          compile read stats
#          a) virsorter2 / vibrant (phage identication)
#          b) checkv  (identify likely true phage)
#          c) fastANI + MCL (dereplication of phage)
#          plasmid identification (PlasFlow/gplas or PlasClass)
#          log failed assemblies
#          min contig size for assemblies, 2000bp?

