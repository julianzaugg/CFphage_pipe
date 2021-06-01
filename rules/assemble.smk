

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
         -g {params.genome_size} -t {threads} {params.flye_parameters} && \
         cp data/assembly/{wildcards.sample}/flye/assembly.fasta \
         {output} && \
         cp data/assembly/{wildcards.sample}/flye/assembly_graph.gfa \
         data/assembly/{wildcards.sample}/flye/{wildcards.sample}.flye.gfa \
         || touch {output}
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
         -g {params.genome_size} -t {threads} {params.metaflye_parameters} && \
         cp data/assembly/{wildcards.sample}/metaflye/assembly.fasta \
         {output} && \
         cp data/assembly/{wildcards.sample}/metaflye/assembly_graph.gfa \
         data/assembly/{wildcards.sample}/metaflye/{wildcards.sample}.metaflye.gfa \
         || touch {output}
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
        {output} \
        || touch {output}
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
        --threads {threads} {input.reads} > {output} \
        || touch {output}
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
        -o data/assembly/{wildcards.sample}/wtdbg2/{wildcards.sample}.wtdbg2 && \
        wtpoa-cns -t {threads} \
        -i data/assembly/{wildcards.sample}/wtdbg2/{wildcards.sample}.wtdbg2.ctg.lay.gz -fo {output} \
        || touch {output}
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
        awk '$1 ~/S/ {{print ">"$2"\\n"$3}}' data/assembly/{wildcards.sample}/miniasm/{wildcards.sample}.miniasm.gfa \
        > {output} \
        || touch {output}
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
        -o data/polishing/{wildcards.sample}/medaka/{wildcards.assembler} -m {params.guppy_model} \
        || touch data/polishing/{wildcards.sample}/medaka/{wildcards.assembler}/consensus.fasta
        
        cp data/polishing/{wildcards.sample}/medaka/{wildcards.assembler}/consensus.fasta {output}
        """

# ------------------------------------------------------------------------------------------------
# Circularise assemblies (if possible)
rule circularise:
    input:
        expand("data/circlator/{sample}/{assembler}/{sample}.{assembler}.final.fasta",
            sample=SAMPLES,assembler=ASSEMBLERS)
    output:
        temp(touch("finished_circularising"))

rule circlator:
    input:
        assembly="data/polishing/{sample}/medaka/{assembler}/{sample}.{assembler}.medaka.fasta",
        reads=config["LONG_READ_DIR"] + "/{sample}.fastq.gz"
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
