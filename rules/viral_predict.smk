# ------------------------------------------------------------------------------------------------
# Run viral tools on reads

# TODO instead of rule order, separate handling to a 'run_viral_prediction' rule and 'run_viral_filtering' rules
ruleorder: virsorter_assembly > checkv_assembly

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
        "../envs/virsorter.yaml"
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
        touch("data/viral_assembly_predict/{sample}/virsorter/{assembler}/done")
    message:
        "Running virsorter on {input.assembly}"
    conda:
        "../envs/virsorter.yaml"
    threads:
        config["MAX_THREADS"]
    shell:
        """
        VIRSORTER_DIR="data/viral_assembly_predict/{wildcards.sample}/virsorter/{wildcards.assembler}"
        
        mkdir -p $VIRSORTER_DIR && \
        if [ -s {input.assembly} ]; then
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

rule seeker_assembly:
    input:
        assembly = "data/polishing/{sample}/medaka/{assembler}/{sample}.{assembler}.medaka.fasta"
    output:
        touch("data/viral_assembly_predict/{sample}/seeker/{assembler}/done")
    params:
        seeker_min_length = config["SEEKER"]["MIN_LENGTH"]
    conda:
        "../envs/seeker.yaml"
    message:
        "Running seeker on {input.assembly}"
    threads:
        config["MAX_THREADS"]
    shell:
        """
        SEEKER_DIR="data/viral_assembly_predict/{wildcards.sample}/seeker/{wildcards.assembler}"
        mkdir -p $SEEKER_DIR
        if [ -s {input.assembly} ]; then
            seqkit seq -m {params.seeker_min_length} {input.assembly} > $SEEKER_DIR/length_filtered.fasta
            predict-metagenome $SEEKER_DIR/length_filtered.fasta \
            > $SEEKER_DIR/seeker_scores.tsv
            
            awk -F "\t" '{{if($2=="Phage" && $3 >= 0.5 )print $1 }}' \
            $SEEKER_DIR/seeker_scores.tsv \
            | seqkit grep --by-name --pattern-file - $SEEKER_DIR/length_filtered.fasta \
            > $SEEKER_DIR/{wildcards.sample}.{wildcards.assembler}.seeker.fasta

            if [ -s $SEEKER_DIR/{wildcards.sample}.{wildcards.assembler}.seeker.fasta ]; then
                sed -i "s/>/>{wildcards.sample}__{wildcards.assembler}__seeker____/g" \
                $SEEKER_DIR/{wildcards.sample}.{wildcards.assembler}.seeker.fasta
            fi
        fi
        """

rule vibrant_assembly:
    input:
        assembly="data/polishing/{sample}/medaka/{assembler}/{sample}.{assembler}.medaka.fasta"
    output:
        touch("data/viral_assembly_predict/{sample}/vibrant/{assembler}/done")
    params:
        vibrant_database = config["VIBRANT"]["DATABASE_DIR"],
        vibrant_min_length=config["VIBRANT"]["MIN_LENGTH"]
    conda:
        "../envs/vibrant.yaml"
    message:
        "Running vobrant on {input.assembly}"
    threads:
        config["MAX_THREADS"]
    shell:
        """
        VIBRANT_DIR="data/viral_assembly_predict/{wildcards.sample}/vibrant/{wildcards.assembler}"
        mkdir -p $VIBRANT_DIR
        
        if [ -s {input.assembly} ]; then
            seqkit seq -m {params.vibrant_min_length} {input.assembly} \
            > $VIBRANT_DIR/length_filtered.fasta
            
            VIBRANT_run.py \
            -i $VIBRANT_DIR/length_filtered.fasta \
            -folder $VIBRANT_DIR \
            -d {params.vibrant_database} \
            -t {threads}
            
            cat $VIBRANT_DIR/VIBRANT*/VIBRANT_phages_*/*phages_combined.txt \
            | seqkit grep --by-name --pattern-file - $VIBRANT_DIR/length_filtered.fasta \
            > $VIBRANT_DIR/{wildcards.sample}.{wildcards.assembler}.vibrant.fasta
            
            if [ -s $VIBRANT_DIR/{wildcards.sample}.{wildcards.assembler}.vibrant.fasta ]; then
                sed -i "/^>/ s/$/___vibrant/g" $VIBRANT_DIR/{wildcards.sample}.{wildcards.assembler}.vibrant.fasta
            fi
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
        "../envs/checkv.yaml"
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

# ------------------------------------------------------------------------------------------------
# Cluster and dereplicate viral sequences

rule viral_cluster:
    input:
        "data/viral_clustering/fastani/fastani_viral.tsv",
        "data/viral_clustering/mcl/mcl_viral_clusters.tsv",
        "data/viral_clustering/mcl/cluster_representative_sequences/done",
        "finished_viral_assembly_predict"
        # dynamic("data/viral_clustering/mcl/cluster_representative_sequences/cluster_{id}_rep.fasta"),
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
        "../envs/fastani.yaml"
    params:
        frag_len = 600,
        min_fraction = 0.85
    threads:
        config["MAX_THREADS"]
    script:
        "../scripts/viral_fastani.py"

# Takes the raw output from FastANI and calculates average for each bidirectional pair of genomes
rule fastani_average:
    input:
        "data/viral_clustering/fastani/fastani_viral.tsv"
    output:
        fastani_viral_mean = "data/viral_clustering/fastani/fastani_viral_mean.tsv",
        fastani_viral_ani95_mcl = "data/viral_clustering/fastani/fastani_viral_ani95_mcl.tsv",
        done = "data/viral_clustering/fastani/done"
    conda:
        "../envs/fastani_average.yaml"
    shell:
        f"""
        Rscript --vanilla {SNAKE_PATH}/scripts/fastani_avg.R \
        {ABSOLUTE_DATA_PATH}/{{input}} \
        {ABSOLUTE_DATA_PATH}/{{output.fastani_viral_mean}} \
        {ABSOLUTE_DATA_PATH}/{{output.fastani_viral_ani95_mcl}}
        touch {{output.done}}
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
        "../envs/mcl.yaml"
    shell:
        """
        BASE_DIR="data/viral_clustering/mcl"
        mkdir -p $BASE_DIR
        
        mcxload -abc {input} --stream-mirror -write-tab $BASE_DIR/viral.tab -o $BASE_DIR/viral.mci
        mcl $BASE_DIR/viral.mci -I {params.inflation_factor} -o $BASE_DIR/out.viral.mci.I
        mcxdump -icl $BASE_DIR/out.viral.mci.I -tabr $BASE_DIR/viral.tab -o $BASE_DIR/dump.viral.mci.I

        sed 's/\t/,/g' $BASE_DIR/dump.viral.mci.I | cat -n | sed -e 's/^[ \t]*//' | sed 's/^/cluster_/g' \
        > {output.mcl_viral_clusters}
        touch {output.done}
        """


rule get_cluster_representatives:
    input:
        mcl_viral_clusters = "data/viral_clustering/mcl/mcl_viral_clusters.tsv",
    output:
        touch("data/viral_clustering/mcl/cluster_representative_sequences/done")
        #dynamic("data/viral_clustering/mcl/cluster_representative_sequences/cluster_{id}_rep.fasta")
    shell:
        """
         if [[ -f data/viral_clustering/mcl/clusters_all_members_lengths.tsv ]];then
             rm data/viral_clustering/mcl/clusters_all_members_lengths.tsv
         fi
         while read cluster_info; do
             cluster=$(echo "$cluster_info" | awk -F "\t" '{{print $1}}')
             cluster_members=$(echo "$cluster_info" | awk -F "\t" '{{print $2}}')
             while read member;do
                 member_length=$(cat "${{member}}.fasta" | tail -n +2 | wc -c)
                 echo -e "${{cluster}}\t${{member}}\t${{member_length}}" \
                 >> data/viral_clustering/mcl/clusters_all_members_lengths.tsv
             done < <(echo "$cluster_members" | tr , "\n")
         done < {input.mcl_viral_clusters}

         cat data/viral_clustering/mcl/clusters_all_members_lengths.tsv | sort -k 1,1 -r -k 3,3 \
         > data/viral_clustering/mcl/temp
         mv data/viral_clustering/mcl/temp data/viral_clustering/mcl/clusters_all_members_lengths.tsv

         awk '! a[$1]++' data/viral_clustering/mcl/clusters_all_members_lengths.tsv \
         > data/viral_clustering/mcl/cluster_representatives_lengths.tsv

         mkdir -p data/viral_clustering/mcl/cluster_representative_sequences
         while read rep_entry;do
             cluster=$(echo "$rep_entry" | awk -F "\t" '{{print $1}}')
             member_name=$(echo "$rep_entry" | awk -F "\t" '{{print $2}}') # Assumes name is full path!!
             member_filename=${{member_name}}.fasta
             cp $member_filename data/viral_clustering/mcl/cluster_representative_sequences/${{cluster}}_rep.fasta
         done < data/viral_clustering/mcl/cluster_representatives_lengths.tsv
         """

# ------------------------------------------------------------------------------------------------
# Viral annotation

rule viral_annotate:
    input:
        "finished_viral_clustering",
        "data/viral_annotation/prodigal/done",
        "data/viral_annotation/abricate/done"
    output:
        touch("finished_viral_annotation")

# Run prodigal on viruses. Assumes at least one CDS will be present
rule viral_prodigal:
    input:
        "data/viral_clustering/mcl/cluster_representative_sequences/done"
    output:
        touch("data/viral_annotation/prodigal/done")
    message:
        "Running prodigal on representative viral sequences"
    params:
        procedure = "meta"
    conda:
        "../envs/prodigal.yaml"
    shell:
        """
        for rep_sequence_file in data/viral_clustering/mcl/cluster_representative_sequences/cluster_*_rep.fasta; do
            name=$(basename $rep_sequence_file .fasta)
            OUTDIR=data/viral_annotation/prodigal/$name
            mkdir -p $OUTDIR
            prodigal -i $rep_sequence_file \
            -a $OUTDIR/$name.faa \
            -d $OUTDIR/$name.fna \
            -p {params.procedure} \
            -f gff \
            > $OUTDIR/$name.gff        
        done
        """

rule viral_abricate:
    input:
        "data/viral_clustering/mcl/cluster_representative_sequences/done"
    output:
        touch("data/viral_annotation/abricate/done")
    conda:
        "../envs/abricate.yaml"
    message:
        "Running abricate on representative viral sequences"
    shell:
        """
        declare -a databases=("card" "vfdb")
        ABRICATE_DIR="data/viral_annotation/abricate"
        mkdir -p $ABRICATE_DIR
        for db in "${{databases[@]}}";do
            abricate -db $db data/viral_clustering/mcl/cluster_representative_sequences/*.fasta | sed "s/.fasta//g" \
            > $ABRICATE_DIR/viral_abricate_${{db}}.tsv
        done
        """
