# Annotate viral sequences

rule viral_annotate:
    input:
        "finished_viral_clustering",
        "finished_viral_assembly_predict",
        "finished_viral_assembly_filtered_reads_predict",
        "data/viral_predict/all_samples_viral_sequences.fasta",
        "data/viral_annotation/checkv/done",
        "data/viral_annotation/prodigal/done",
        "data/viral_annotation/abricate/done"
    output:
        touch("finished_viral_annotation")


# TODO remove
# rule viral_assembly_filter:
#     input:
#         expand("data/viral_assembly_predict/{sample}/{viral_predict_tool}/{assembler}/done",
#             sample=SAMPLES,
#             assembler=ASSEMBLERS,
#             viral_predict_tool=VIRAL_TOOLS),
#         "data/checkv_assembly/done",
#         "finished_viral_assembly_predict"
#     output:
#         touch("finished_viral_assembly_filter")

# def collect_assembly_viral_outputs(wildcards):
#     files = expand("data/viral_assembly_predict/{sample}/{viral_predict_tool}/{assembler}/"
#                    "{sample}.{assembler}.{viral_predict_tool}.fasta",
#         sample=SAMPLES,
#         assembler=ASSEMBLERS,
#         viral_predict_tool=VIRAL_TOOLS)
#     files = [file for file in files if os.path.isfile(file)]
#     return files
#
# rule collect_viral_sequences_assembly:
#     input:
#         viral_tool_output=collect_assembly_viral_outputs
#     output:
#         "data/viral_predict/assembly/all_samples_assembly_viral_sequences.fasta"
#     shell:
#         """
#         cat {input.viral_tool_output} \
#         > data/viral_predict/assembly/all_samples_assembly_viral_sequences.fasta
#         """


def collect_viral_outputs(wildcards):
    assembly_files = expand("data/viral_predict/assembly/{sample}/{viral_predict_tool}/{assembler}/"
                   "{sample}.{assembler}.{viral_predict_tool}.fasta",
        sample = SAMPLES,
        assembler = ASSEMBLERS,
        viral_predict_tool = VIRAL_TOOLS_ASSEMBLY)

    assembly_filtered_files = expand("data/viral_predict/assembly_filtered_reads/{sample}/{viral_predict_tool}/{assembler}/"
                   "{sample}.{assembler}.{viral_predict_tool}.fasta",
        sample = SAMPLES,
        assembler = ASSEMBLERS,
        viral_predict_tool = VIRAL_TOOLS_READS)
    files = [file for file in assembly_files + assembly_filtered_files if os.path.isfile(file)]
    return files

rule collect_viral_sequences:
    input:
        viral_tool_output = collect_viral_outputs
    output:
        "data/viral_predict/all_samples_viral_sequences.fasta"
    shell:
        """        
        cat {input.viral_tool_output} \
        > {output}
        """

# Run checkV on predicted viral sequences
rule checkv:
    input:
        all_viral_sequences="data/viral_predict/all_samples_viral_sequences.fasta"
    output:
        touch("data/viral_annotation/checkv/done")
    message:
        "Running checkv"
    params:
        checkv_db=config["CHECKV_DB"]
    conda:
        "../envs/checkv.yaml"
    threads:
        config["MAX_THREADS"]
    shell:
        """
        mkdir -p data/viral_annotation/checkv
        # if [[ -f data/viral_annotation/checkv/checkv_selected.fasta ]];then
        #     rm data/viral_annotation/checkv/checkv_selected.fasta
        # fi

        checkv end_to_end \
        -d {params.checkv_db} \
        -t {threads} \
        --restart \
        {input.all_viral_sequences} \
        data/viral_annotation/checkv/

        # Combine viruses.fna and proviruses.fna, assume they exist
        cat data/viral_annotation/checkv/viruses.fna data/viral_annotation/checkv/proviruses.fna \
        > data/viral_annotation/checkv/checkv_viruses.fasta

        # Grab all Medium and High quality, and Complete, viral genomes and write to fasta file 
        # while read contig; do
        #     grep -h $contig data/viral_annotation/checkv/checkv_viruses.fasta -A 1 >> data/viral_annotation/checkv/checkv_selected_viruses.fasta
        # done < <(awk -F "\t" '$8~/(Complete|[Medium,High]-quality)$/{{print $1}}' data/viral_annotation/checkv/quality_summary.tsv)
        """

# ------------------------------------------------------------------------------------------------
# Cluster and dereplicate predicted viral sequences

rule viral_cluster:
    input:
        "data/viral_clustering/fastani/fastani_viral.tsv",
        "data/viral_clustering/mcl/mcl_viral_clusters.tsv",
        "data/viral_clustering/mcl/cluster_representative_sequences/done",
        "data/viral_predict/all_samples_viral_sequences.fasta"
    # dynamic("data/viral_clustering/mcl/cluster_representative_sequences/cluster_{id}_rep.fasta"),
    output:
        touch("finished_viral_clustering")

# Run FastANI on viral sequences
rule fastani_viral:
    input:
        # checkv_selected = "data/viral_annotation/checkv/checkv_selected.fasta"
        all_viral_sequences="data/viral_predict/all_samples_viral_sequences.fasta"
    output:
        "data/viral_clustering/fastani/fastani_viral.tsv"
    message:
        "Clustering selected viral sequences"
    conda:
        "../envs/fastani.yaml"
    params:
        frag_len=600,
        min_fraction=0.85
    threads:
        config["MAX_THREADS"]
    script:
        "../scripts/viral_fastani.py"

# Takes the raw output from FastANI and calculates average for each bidirectional pair of genomes
rule fastani_average:
    input:
        "data/viral_clustering/fastani/fastani_viral.tsv"
    output:
        fastani_viral_mean="data/viral_clustering/fastani/fastani_viral_mean.tsv",
        fastani_viral_ani95_mcl="data/viral_clustering/fastani/fastani_viral_ani95_mcl.tsv",
        done="data/viral_clustering/fastani/done"
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
        mcl_viral_clusters="data/viral_clustering/mcl/mcl_viral_clusters.tsv",
        done="data/viral_clustering/mcl/done"
    params:
        inflation_factor=3.5
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

# TODO output file listing : cluster rep(TRUE/FALSE) seq_name
rule get_cluster_representatives:
    input:
        mcl_viral_clusters="data/viral_clustering/mcl/mcl_viral_clusters.tsv",
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

# Run prodigal on viruses. Assumes at least one CDS will be present
rule viral_prodigal:
    input:
        all_viral_sequences = "data/viral_predict/all_samples_viral_sequences.fasta"
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
        # for rep_sequence_file in data/viral_clustering/mcl/cluster_representative_sequences/cluster_*_rep.fasta; do
        name=$(basename {input.all_viral_sequences} .fasta)
        OUTDIR=data/viral_annotation/prodigal/$name
        mkdir -p $OUTDIR
        prodigal -i {input.all_viral_sequences} \
        -a $OUTDIR/$name.faa \
        -d $OUTDIR/$name.fna \
        -p {params.procedure} \
        -f gff \
        > $OUTDIR/$name.gff        
        # done
        """

rule viral_abricate:
    input:
        all_viral_sequences = "data/viral_predict/all_samples_viral_sequences.fasta"
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
            abricate -db $db {input.all_viral_sequences} | sed "s/.fasta//g" \
            > $ABRICATE_DIR/viral_abricate_${{db}}.tsv
        done
        """
