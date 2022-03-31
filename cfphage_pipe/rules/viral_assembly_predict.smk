# ------------------------------------------------------------------------------------------------
# Run viral tools on polished assemblies

rule viral_assembly_predict:
    input:
        expand("data/viral_predict/assembly/{sample}/{viral_predict_tool}/{assembler}/done",
            sample = SAMPLES,
            assembler = ASSEMBLERS,
            viral_predict_tool = VIRAL_TOOLS_ASSEMBLY),
        "finished_polishing"
    priority: 1
    output:
        touch("finished_viral_assembly_predict")

rule virsorter_assembly:
    input:
        assembly = "data/polishing/{sample}/medaka/{assembler}/{sample}.{assembler}.medaka.fasta"
    params:
        virsorter_database = config["VIRSORTER"]["DATABASE_DIR"],
        virsorter_min_length = config["VIRSORTER"]["MIN_LENGTH"]
    output:
        touch("data/viral_predict/assembly/{sample}/virsorter/{assembler}/done")
    message:
        "Running virsorter on {input.assembly}"
    conda:
        "../envs/virsorter.yaml"
    threads:
        config["MAX_THREADS"]
    shell:
        """
        VIRSORTER_DIR="data/viral_predict/assembly/{wildcards.sample}/virsorter/{wildcards.assembler}"
        
        mkdir -p $VIRSORTER_DIR && \
        if [ -s {input.assembly} ]; then
            # --high-confidence-only
            virsorter run \
            --rm-tmpdir \
            --seqfile {input.assembly} \
            --working-dir $VIRSORTER_DIR \
            --db-dir {params.virsorter_database} \
            --min-length {params.virsorter_min_length} \
            --jobs {threads} \
            all 
        fi 
        
        if [[ -f $VIRSORTER_DIR/final-viral-combined.fa ]]; then
            cp $VIRSORTER_DIR/final-viral-combined.fa \
            $VIRSORTER_DIR/{wildcards.sample}.{wildcards.assembler}.virsorter.fasta
            sed -i "s/>/>{wildcards.sample}__{wildcards.assembler}__virsorter____/g" \
            $VIRSORTER_DIR/{wildcards.sample}.{wildcards.assembler}.virsorter.fasta
        fi 
        """


rule vibrant_assembly:
    input:
        assembly="data/polishing/{sample}/medaka/{assembler}/{sample}.{assembler}.medaka.fasta"
    output:
        touch("data/viral_predict/assembly/{sample}/vibrant/{assembler}/done")
    params:
        vibrant_database = config["VIBRANT"]["DATABASE_DIR"],
        vibrant_min_length=config["VIBRANT"]["MIN_LENGTH"]
    conda:
        "../envs/vibrant.yaml"
    message:
        "Running vibrant on {input.assembly}"
    threads:
        config["MAX_THREADS"]
    shell:
        """
        VIBRANT_DIR="data/viral_predict/assembly/{wildcards.sample}/vibrant/{wildcards.assembler}"
        mkdir -p $VIBRANT_DIR
        
        if [ -s {input.assembly} ]; then
            seqkit seq -m {params.vibrant_min_length} {input.assembly} \
            > $VIBRANT_DIR/length_filtered.fasta
            
            if [ -s $VIBRANT_DIR/length_filtered.fasta ]; then
                VIBRANT_run.py \
                -i $VIBRANT_DIR/length_filtered.fasta \
                -folder $VIBRANT_DIR \
                -d {params.vibrant_database} \
                -t {threads}
            fi
            
            if [ -f $VIBRANT_DIR/VIBRANT_length_filtered/VIBRANT_phages_length_filtered/length_filtered.phages_combined.fna ]; then
                cp $VIBRANT_DIR/VIBRANT_length_filtered/VIBRANT_phages_length_filtered/length_filtered.phages_combined.fna \
                $VIBRANT_DIR/{wildcards.sample}.{wildcards.assembler}.vibrant.fasta
            else
                touch $VIBRANT_DIR/{wildcards.sample}.{wildcards.assembler}.vibrant.fasta
            fi
            if [ -s $VIBRANT_DIR/{wildcards.sample}.{wildcards.assembler}.vibrant.fasta ]; then
                sed -i "s/>/>{wildcards.sample}__{wildcards.assembler}__vibrant____/g" \
                $VIBRANT_DIR/{wildcards.sample}.{wildcards.assembler}.vibrant.fasta
            fi

            rm $VIBRANT_DIR/length_filtered.fasta
        fi
        """

