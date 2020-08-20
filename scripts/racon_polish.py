import sys
import os
import subprocess
import shutil

if not os.path.exists(snakemake.input.reads):
    sys.exit("No read file found")

out = "data/polishing/{}/racon/{}".format(snakemake.wildcards.sample,
                                          snakemake.wildcards.assembler)

try:
    os.makedirs(out)
except FileExistsError:
    pass

for round in range(snakemake.params.rounds):
    # If this is the first polishing round, generate first polished assembly
    if round == 0:
        if not os.path.isfile("{}/{}.pol.{}.fasta".format(out, snakemake.wildcards.sample, round)):
            # Map reads
            subprocess.Popen("minimap2 -t {} -x map-ont {} {} -o {}/{}.{}.pol.{}.paf".format(
                snakemake.threads,
                snakemake.input.assembly,
                snakemake.input.reads,
                out, snakemake.wildcards.sample,snakemake.wildcards.assembler,round),
                shell=True).wait()
            # Polish
            subprocess.Popen("racon ‐‐include-unpolished -m {} -x {} -g {} -w {} -t {} {} {}/{}.{}.pol.{}.paf {} " \
                             "> {}/{}.{}.pol.{}.fasta".format(
                snakemake.params.match,
                snakemake.params.mismatch,
                snakemake.params.gap,
                snakemake.params.window_length,
                snakemake.threads,
                snakemake.input.reads,
                out,snakemake.wildcards.sample,snakemake.wildcards.assembler,round,
                snakemake.input.assembly,
                out, snakemake.wildcards.sample,snakemake.wildcards.assembler, round
            ), shell=True).wait()
    else:
        prev_round = round - 1

        if not os.path.isfile("{}/{}.pol.{}.fasta".format(out, snakemake.wildcards.sample, round)):
            # Map reads onto the assembly from the previous round
            subprocess.Popen("minimap2 -t {} -x map-ont {}/{}.{}.pol.{}.fasta {} -o {}/{}.{}.pol.{}.paf".format(
                snakemake.threads,
                out,snakemake.wildcards.sample,snakemake.wildcards.assembler,prev_round,
                snakemake.input.reads,
                out,snakemake.wildcards.sample,snakemake.wildcards.assembler,round),
                shell=True).wait()
            # Polish the assembly from the previous round
            subprocess.Popen("racon ‐‐include-unpolished -m {} -x {} -g {} -w {} -t {} {} {}/{}.{}.pol.{}.paf " \
                             "{}/{}.{}.pol.{}.fasta > {}/{}.{}.pol.{}.fasta".format(
                snakemake.params.match,
                snakemake.params.mismatch,
                snakemake.params.gap,
                snakemake.params.window_length,
                snakemake.threads,
                snakemake.input.reads,
                out,snakemake.wildcards.sample,snakemake.wildcards.assembler,round,
                out,snakemake.wildcards.sample,snakemake.wildcards.assembler,prev_round,
                out, snakemake.wildcards.sample,snakemake.wildcards.assembler, round
            ), shell=True).wait()
    reference = "{}/{}.{}.pol.{}.fasta".format(out, snakemake.wildcards.sample,snakemake.wildcards.assembler,round)
shutil.copy2(reference, snakemake.output[0])

