
# ------------------------------------------------------------------------------------------------
# Run viral tools on raw reads
rule viral_raw_reads_predict:
    input:
        expand("data/viral_reads_predict.smk/{sample}/virsorter/done", sample = SAMPLES),
        "finished_QC"
    output:
        touch("finished_viral_raw_reads_predict")

rule virsorter_raw_reads:
    input:
        reads = "data/nanofilt/{sample}_nanofilt.fastq.gz"
    params:
        virsorter_database = config["VIRSORTER"]["DATABASE_DIR"],
        virsorter_min_length = config["VIRSORTER"]["MIN_LENGTH"]
    output:
        "data/viral_reads_predict.smk/{sample}/virsorter/done"
    conda:
        "../envs/virsorter.yaml"
    threads:
        config["MAX_THREADS"]
    shell:
        """
        mkdir -p data/viral_reads_predict.smk/{wildcards.sample} && \
        virsorter run \
        --rm-tmpdir \
        --seqfile {input.reads} \
        --working-dir data/viral_reads_predict.smk/{wildcards.sample}/virsorter \
        --db-dir {params.virsorter_database} \
        --min-length {params.virsorter_min_length} \
        --jobs {threads} \
        all
        touch data/viral_reads_predict.smk/{wildcards.sample}/virsorter/done
        """

# ------------------------------------------------------------------------------------------------
# Run viral tools on assembly filtered reads

rule viral_assembly_filtered_reads_predict:
    input:
        expand("data/viral_predict/assembly_filtered_reads/{sample}/{viral_predict_tool}/{assembler}/done",
            sample = SAMPLES,
            assembler = ASSEMBLERS,
            viral_predict_tool = VIRAL_TOOLS_READS),
        "data/assembly_filtered_reads/{sample}_{assembler}_AF.fasta",
        "finished_polishing"
    output:
        touch("finished_viral_assembly_filtered_reads_predict")

rule virsorter_assembly_filtered_reads:
    input:
        filtered_fasta = "data/assembly_filtered_reads/{sample}_{assembler}_AF.fasta"
    params:
        virsorter_database=config["VIRSORTER"]["DATABASE_DIR"],
        virsorter_min_length=config["VIRSORTER"]["MIN_LENGTH"]
    output:
        touch("data/viral_predict/assembly_filtered_reads/{sample}/virsorter/{assembler}/done")
    message:
        "Running virsorter on {input.assembly}"
    conda:
        "../envs/virsorter.yaml"
    threads:
        config["MAX_THREADS"]
    shell:
        """
        VIRSORTER_DIR="data/viral_predict/assembly_filtered_reads/{wildcards.sample}/virsorter/{wildcards.assembler}"

        mkdir -p $VIRSORTER_DIR && \
        if [ -s {input.filtered_fasta} ]; then
            # --high-confidence-only
            virsorter run \
            --rm-tmpdir \
            --seqfile {input.filtered_fasta} \
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

rule seeker_assembly_filtered_reads:
    input:
        filtered_fasta = "data/assembly_filtered_reads/{sample}_{assembler}_AF.fasta"
    output:
        touch("data/viral_predict/assembly_filtered_reads/{sample}/seeker/{assembler}/done")
    params:
        seeker_min_length = config["SEEKER"]["MIN_LENGTH"]
    conda:
        "../envs/seeker.yaml"
    message:
        "Running seeker on {input.filtered_fasta}"
    threads:
        config["MAX_THREADS"]
    shell:
        """
        SEEKER_DIR="data/viral_predict/assembly_filtered_reads/{wildcards.sample}/seeker/{wildcards.assembler}"
        mkdir -p $SEEKER_DIR
        if [ -s {input.filtered_fasta} ]; then
            seqkit seq -m {params.seeker_min_length} {input.assembly} > $SEEKER_DIR/length_filtered.fasta

            predict-metagenome $SEEKER_DIR/length_filtered.fasta \
            > $SEEKER_DIR/seeker_scores.tsv

            if grep -q --max-count 1 "Phage" $SEEKER_DIR/seeker_scores.tsv; then
                awk -F "\t" '{{if($2=="Phage" && $3 >= 0.5 )print $1 }}' \
                $SEEKER_DIR/seeker_scores.tsv \
                | seqkit grep --by-name --pattern-file - $SEEKER_DIR/length_filtered.fasta \
                > $SEEKER_DIR/{wildcards.sample}.{wildcards.assembler}.seeker.fasta

                sed -i "s/>/>{wildcards.sample}__{wildcards.assembler}__seeker____/g" \
                $SEEKER_DIR/{wildcards.sample}.{wildcards.assembler}.seeker.fasta
            else
                touch $SEEKER_DIR/{wildcards.sample}.{wildcards.assembler}.seeker.fasta
            fi
        fi
        """

rule vibrant_assembly_filtered_reads:
    input:
        filtered_fasta = "data/assembly_filtered_reads/{sample}_{assembler}_AF.fasta"
    output:
        touch("data/viral_predict/assembly_filtered_reads/{sample}/vibrant/{assembler}/done")
    params:
        vibrant_database=config["VIBRANT"]["DATABASE_DIR"],
        vibrant_min_length=config["VIBRANT"]["MIN_LENGTH"]
    conda:
        "../envs/vibrant.yaml"
    message:
        "Running vibrant on {input.filtered_fasta}"
    threads:
        config["MAX_THREADS"]
    shell:
        """
        VIBRANT_DIR="data/viral_predict/assembly_filtered_reads/{wildcards.sample}/vibrant/{wildcards.assembler}"
        mkdir -p $VIBRANT_DIR

        if [ -s {input.assembly} ]; then
            seqkit seq -m {params.vibrant_min_length} {input.assembly} \
            > $VIBRANT_DIR/length_filtered.fasta

            VIBRANT_run.py \
            -i $VIBRANT_DIR/length_filtered.fasta \
            -folder $VIBRANT_DIR \
            -d {params.vibrant_database} \
            -t {threads}

            cp $VIBRANT_DIR/VIBRANT_length_filtered/VIBRANT_phages_length_filtered/length_filtered.phages_combined.fna \
            $VIBRANT_DIR/{wildcards.sample}.{wildcards.assembler}.vibrant.fasta

            if [ -s $VIBRANT_DIR/{wildcards.sample}.{wildcards.assembler}.vibrant.fasta ]; then
                sed -i "s/>/>{wildcards.sample}__{wildcards.assembler}__vibrant____/g" \
                $VIBRANT_DIR/{wildcards.sample}.{wildcards.assembler}.vibrant.fasta
            fi

            rm $VIBRANT_DIR/length_filtered.fasta
        fi
        """
