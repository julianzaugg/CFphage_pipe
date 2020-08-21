import sys
import os
import subprocess
import shutil


if not os.path.exists(snakemake.input.reference_fasta):
    sys.exit("No reference genome fasta found")


out = f"data/coverage/{snakemake.wildcards.reference_genome}"
try:
    os.makedirs(out)
except FileExistsError:
    pass

# {snakemake.params.read_files}
# {snakemake.params.reference_coverm_parameters_dict}

# if snakemake.params.reference_coverm_parameters_dict["multiple_genomes"] == False:
#     echo {params.reference_parameters['min_read_percent_identity']}
#     mkdir -p data/coverage/{wildcards.reference_genome}
#
#     coverm make --reference {input.reference_fasta} --threads {snakemake.threads} \
#     --output-directory data/coverage/{wildcards.reference_genome}/ \
#     --single {params.read_files} --mapper minimap2-ont
#
#     coverm genome --bam-files data/coverage/{wildcards.reference_genome}/*.bam --threads {threads}
#     --methods relative_abundance --single-genome {input.reference_fasta} \
#     > data/coverage/{wildcards.reference_genome}/{wildcards.reference_genome}_rel_abundance_table.tsv
# else: