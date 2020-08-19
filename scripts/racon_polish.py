import sys
import os
import subprocess
import shutil

if not os.path.exists(snakemake.input.reads):
    sys.exit("No read file found")

out = "data/polishing/{}/racon".format(snakemake.wildcards.sample)

try:
    os.makedirs(out)
except FileExistsError:
    pass

for round in range(snakemake.params.rounds):
    # If this is the first polishing round, generate first polished assembly
    if round == 0:
        # Map reads
        subprocess.Popen("minimap2 -t {} -x map-ont {} {} -o {}/{}.pol.{}.paf".format(snakemake.threads,
                                                                        snakemake.input.assembly,
                                                                        snakemake.input.reads,
                                                                        out,
                                                                        snakemake.wildcards.sample,
                                                                                  round),
                         shell=True).wait()
        # Polish
        subprocess.Popen("racon -m 8 -x -6 -g -8 -w 500 -t {} {} {}/{}.pol.{}.paf {} > {}/{}.pol.{}.fasta".format(
            snakemake.threads,
            snakemake.input.reads,
            out,snakemake.wildcards.sample,round,
            snakemake.input.assembly,
            out, snakemake.wildcards.sample, round + 1
        ), shell=True).wait()
    else:
        prev_round = round - 1
        # Map reads onto the assembly from the previous round
        subprocess.Popen("minimap2 -t {} -x map-ont {}/{}.pol.{}.fasta {} -o {}/{}.pol.{}.paf".format(snakemake.threads,
                                                                        out,snakemake.wildcards.sample,prev_round,
                                                                        snakemake.input.reads,
                                                                        out,snakemake.wildcards.sample,round),
                         shell=True).wait()
        # Polish the assembly from the previous round
        subprocess.Popen("racon -m 8 -x -6 -g -8 -w 500 -t {} {} {}/{}.pol.{}.paf {}/{}.pol.{}.fasta > " \
                         "{}/{}.pol.{}.fasta".format(
            snakemake.threads,
            snakemake.input.reads,
            out,snakemake.wildcards.sample,prev_round,
            out,snakemake.wildcards.sample,prev_round,
            out, snakemake.wildcards.sample, round
        ), shell=True).wait()
    reference = "{}/{}.pol.{}.fasta".format(out, snakemake.wildcards.sample, round)
shutil.copy2(reference, snakemake.output)

