import sys
import os
import subprocess
import shutil

if not os.path.exists(snakemake.input.reads):
    sys.exit("No read file found")

out_dir = "data/polishing/{}/racon/{}".format(snakemake.wildcards.sample,
                                          snakemake.wildcards.assembler)

try:
    os.makedirs(out_dir)
except FileExistsError:
    pass

for round in range(snakemake.params.rounds):
    # If this is the first polishing round, generate first polished assembly
    if round == 0:
        if not os.path.isfile(f"{out_dir}/{snakemake.wildcards.sample}.pol.{round}.fasta"):
            # Map reads
            subprocess.Popen(
                f"""
                minimap2 -t {snakemake.threads} -x map-ont {snakemake.input.assembly} {snakemake.input.reads} \
                -o {out_dir}/{snakemake.wildcards.sample}.{snakemake.wildcards.assembler}.pol.{round}.paf
                """,
                shell=True).wait()
            # Polish
            subprocess.Popen(
                f"""
                racon --include-unpolished -m {snakemake.params.match} -x {snakemake.params.mismatch} \
                -g {snakemake.params.gap} -w {snakemake.params.window_length} -t {snakemake.threads} \
                {snakemake.input.reads} {out_dir}/{snakemake.wildcards.sample}.{snakemake.wildcards.assembler}.pol.{round}.paf \
                {snakemake.input.assembly} > {out_dir}/{snakemake.wildcards.sample}.{snakemake.wildcards.assembler}.pol.{round}.fasta
                """,
                shell=True).wait()
    else:
        prev_round = round - 1

        if not os.path.isfile("{}/{}.pol.{}.fasta".format(out_dir, snakemake.wildcards.sample, round)):
            # Map reads onto the assembly from the previous round
            subprocess.Popen(
                f"""
                minimap2 -t {snakemake.threads} \
                -x map-ont {out_dir}/{snakemake.wildcards.sample}.{snakemake.wildcards.assembler}.pol.{prev_round}.fasta \
                {snakemake.input.reads} \
                -o {out_dir}/{snakemake.wildcards.sample}.{snakemake.wildcards.assembler}.pol.{round}.paf
                """,
                shell=True).wait()
            # Polish the assembly from the previous round
            subprocess.Popen(
                f"""
                racon --include-unpolished -m {snakemake.params.match} -x {snakemake.params.mismatch} \
                -g {snakemake.params.gap} -w {snakemake.params.window_length} -t {snakemake.threads} \
                {snakemake.input.reads} {out_dir}/{snakemake.wildcards.sample}.{snakemake.wildcards.assembler}.pol.{round}.paf \
                {out_dir}/{snakemake.wildcards.sample}.{snakemake.wildcards.assembler}.pol.{prev_round}.fasta \
                > {out_dir}/{snakemake.wildcards.sample}.{snakemake.wildcards.assembler}.pol.{round}.fasta
                """,
                shell=True).wait()
    reference = f"{out_dir}/{snakemake.wildcards.sample}.{snakemake.wildcards.assembler}.pol.{round}.fasta"
shutil.copy2(reference, snakemake.output[0])

