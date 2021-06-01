# ------------------------------------------------------------------------------------------------
# Assembly annotation

rule assembly_annotate:
    input:
        expand("data/assembly_annotation/prodigal/{sample}/{assembler}/{sample}.{assembler}.faa",
        sample = SAMPLES, assembler = ASSEMBLERS),
        expand("data/assembly_annotation/abricate/{sample}/{assembler}/{sample}.{assembler}.{db}.tsv",
           sample = SAMPLES, assembler = ASSEMBLERS, db = ["card", "vfdb"]),
        "finished_polishing"
    output:
        touch("finished_assembly_annotation")


# Run prodigal on polished assemblies.
rule assembly_prodigal:
    input:
        assembly = "data/polishing/{sample}/medaka/{assembler}/{sample}.{assembler}.medaka.fasta"
    output:
        "data/assembly_annotation/prodigal/{sample}/{assembler}/{sample}.{assembler}.faa"
    message:
        "Running prodigal on polished assemblies"
    params:
        procedure = "meta"
    conda:
        "../envs/prodigal.yaml"
    shell:
        """
        name=$(basename {input.assembly} .medaka.fasta)
        OUTDIR=data/assembly_annotation/prodigal/{wildcards.sample}/{wildcards.assembler}
        mkdir -p $OUTDIR
        prodigal -i {input.assembly} \
        -a $OUTDIR/$name.faa \
        -d $OUTDIR/$name.fna \
        -p {params.procedure} \
        -f gff \
        > $OUTDIR/$name.gff        
        """

# Run abricate on polished assemblies
rule assembly_abricate:
    input:
        assembly = "data/polishing/{sample}/medaka/{assembler}/{sample}.{assembler}.medaka.fasta"
    output:
        "data/assembly_annotation/abricate/{sample}/{assembler}/{sample}.{assembler}.faa"
    conda:
        "../envs/abricate.yaml"
    message:
        "Running abricate on polished assemblies"
    shell:
        """
        declare -a databases=("card" "vfdb")
        ABRICATE_DIR="data/assembly_annotation/abricate/{wildcards.sample}/{wildcards.assembler}"
        mkdir -p $ABRICATE_DIR
        for db in "${{databases[@]}}";do
            abricate -db $db {input.assembly} | sed "s/.medaka.fasta//g" \
            > $ABRICATE_DIR/{wildcards.sample}.{wildcards.assembler}.${{db}}.tsv
        done
        """
