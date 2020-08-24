import sys
import os
import subprocess
import pathlib


if not os.path.exists(snakemake.input.reference_fasta):
    sys.exit("No reference genome fasta found")

out = f"data/coverage/{snakemake.wildcards.reference_genome}"
out_filtered = f"data/coverage/{snakemake.wildcards.reference_genome}/filtered"
pathlib.Path(out_filtered).mkdir(parents=True, exist_ok=True)

MIN_READ_IDENTITY_PERCENT = snakemake.params.reference_coverm_parameters_dict["min_read_percent_identity"]
MIN_READ_ALIGNED_PERCENT = snakemake.params.reference_coverm_parameters_dict["min_read_aligned_percent"]

# Single entry in the reference fasta
if snakemake.params.reference_coverm_parameters_dict["multiple_genomes"] == False:
    subprocess.Popen(
        f"""
        coverm make --reference {snakemake.input.reference_fasta} --threads {snakemake.threads} \
        --output-directory data/coverage/{snakemake.wildcards.reference_genome} \
        --single {snakemake.params.read_files} \
        --mapper minimap2-ont
        """,
        shell=True).wait()

    for method in ["relative_abundance", "mean", "count", "covered_fraction"]:
        min_covered_fraction_param = "--min-covered-fraction 0.0"
        if method == "relative_abundance":
            table_out = f"{snakemake.wildcards.reference_genome}_relative_abundance_table"
        elif method == "mean":
            table_out = f"{snakemake.wildcards.reference_genome}_coverage_table"
        elif method == "count":
            table_out = f"{snakemake.wildcards.reference_genome}_count_table"
        elif method == "count":
            table_out = f"{snakemake.wildcards.reference_genome}_covered_fraction_table"

        subprocess.Popen(
            f"""
            coverm genome --bam-files {out}/*.bam \
            {min_covered_fraction_param} \
            --threads {snakemake.threads} --methods {method} \
            --single-genome \
            > data/coverage/{snakemake.wildcards.reference_genome}/{table_out}.tsv
            """,
            shell=True).wait()

        for bam_file in glob.glob(f"{out}/*.bam"):
            out_bam_filename = pathlib.Path(bam_file).name.replace(".bam", "_filtered.bam")
            subprocess.Popen(
                f"""
                    coverm filter -b {bam_file} -o {out_filtered}/{out_bam_filename} \
                    --min-read-percent-identity {MIN_READ_IDENTITY_PERCENT} \
                    --min-read-aligned-percent {MIN_READ_ALIGNED_PERCENT} \
                    --threads {snakemake.threads}
                """,
                shell=True).wait()

        #             --min-read-percent-identity {MIN_READ_IDENTITY_PERCENT} \
        #             --min-read-aligned-percent {MIN_READ_ALIGNED_PERCENT} \
        #             --discard-unmapped \
        subprocess.Popen(
            f"""
            coverm genome --bam-files {out_filtered}/*.bam \
            {min_covered_fraction_param} \
            --threads {snakemake.threads} --methods {method} \
            --single-genome \
            > {out_filtered}/{table_out}_filtered.tsv
            """,
            shell=True).wait()