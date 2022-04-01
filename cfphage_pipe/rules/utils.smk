


rule reads2fasta:
    input:
        reads = "data/nanofilt/{sample}_nanofilt.fastq.gz"
    output:
        "data/read_fastas/{sample}.fasta"
    conda:
        "../envs/seqtk.yaml"
    shell:
        """
        mkdir -p data/read_fastas
        seqtk seq -a {input.reads} > {output}
        """

