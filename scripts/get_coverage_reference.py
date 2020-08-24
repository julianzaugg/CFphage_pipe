import sys
import os
import subprocess


if not os.path.exists(snakemake.input.reference_fasta):
    sys.exit("No reference genome fasta found")

out = f"data/coverage/{snakemake.wildcards.reference_genome}"
try:
    os.makedirs(out)
except FileExistsError:
    pass

MIN_READ_IDENTITY_PERCENT = snakemake.params.reference_coverm_parameters_dict["min_read_percent_identity"]
MIN_READ_ALIGNED_PERCENT = snakemake.params.reference_coverm_parameters_dict["min_read_aligned_percent"]

if snakemake.params.reference_coverm_parameters_dict["multiple_genomes"] == False:
    subprocess.Popen(f"coverm make --reference {snakemake.input.reference_fasta} --threads {snakemake.threads} " \
                     f"--output-directory data/coverage/{snakemake.wildcards.reference_genome} " \
                     f"--single {snakemake.params.read_files} --mapper minimap2-ont")

    subprocess.Popen(
        f"coverm genome --bam-files data/coverage/{snakemake.wildcards.reference_genome}/*.bam " \
        f"--threads {snakemake.threads} --methods relative_abundance " \
        f"--min-read-percent-identity {MIN_READ_IDENTITY_PERCENT} --min-read-aligned-percent {MIN_READ_ALIGNED_PERCENT}"
        f"--min-covered-fraction 0.0 --discard-unmapped " \
        f"--single-genome {snakemake.input.reference_fasta} " \
        f"> data/coverage/{snakemake.wildcards.reference_genome}/{wildcards.reference_genome}_rel_abundance_table.tsv")

    subprocess.Popen(
        f"coverm genome --bam-files data/coverage/{snakemake.wildcards.reference_genome}/*.bam " \
        f"--threads {snakemake.threads} --methods mean " \
        f"--min-read-percent-identity {MIN_READ_IDENTITY_PERCENT} --min-read-aligned-percent {MIN_READ_ALIGNED_PERCENT}"
        f"--min-covered-fraction 0.0 --discard-unmapped " \
        f"--single-genome {snakemake.input.reference_fasta} " \
        f"> data/coverage/{snakemake.wildcards.reference_genome}/{snakemake.wildcards.reference_genome}_coverage_table.tsv")

    subprocess.Popen(
        f"coverm genome --bam-files data/coverage/{snakemake.wildcards.reference_genome}/*.bam " \
        f"--threads {snakemake.threads} --methods count " \
        f"--min-read-percent-identity {MIN_READ_IDENTITY_PERCENT} --min-read-aligned-percent {MIN_READ_ALIGNED_PERCENT}"
        f"--min-covered-fraction 0.0 --discard-unmapped " \
        f"--single-genome {snakemake.input.reference_fasta} " \
        f"> data/coverage/{snakemake.wildcards.reference_genome}/{snakemake.wildcards.reference_genome}_count_table.tsv")

# else: