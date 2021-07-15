import sys
import subprocess
import pathlib
import glob
import os
import re

import Bio.Seq
import numpy as np

import pandas as pd
from Bio import SeqIO
from Bio.SeqIO.FastaIO import SimpleFastaParser

viral_summary_dir = "data/viral_summary"

try:
    os.makedirs(viral_summary_dir)
except FileExistsError:
    pass

# -------------------------------------------------------------------------------------------------------------------
# Load viral sequences from FASTA file.
# with open("data/viral_predict/all_samples_viral_sequences.fasta") as fasta_file:  # Will close handle cleanly
#     seq_ids = []
#     seq_lengths = []
#     for title, sequence in SimpleFastaParser(fasta_file):
#         seq_ids.append(title.split(None, 1)[0])  # First word is ID
#         seq_lengths.append(len(sequence))

viral_sequences = list(SeqIO.parse("data/viral_predict/all_samples_viral_sequences.fasta", "fasta"))
seq_ids = [s.name for s in viral_sequences]

# Determine the Sample, assembly tool and viral prediction tool from the name and create summary table
def _split_seqid(seqid):
    sequence_details = re.split("_{4}", seqid)[0].split("__")
    derived_from = sequence_details[0]
    sample = sequence_details[1]
    assembly_tool = sequence_details[2]
    viral_tool = sequence_details[3]
    return(seqid,derived_from,sample,assembly_tool,viral_tool)

summary_table = pd.DataFrame(map(_split_seqid, seq_ids),
                             columns = ["Sequence_ID", "Derived_from", "Sample", "Assembly_tool", "Viral_tool"])
# -------------------------------------------------------------------------------------------------------------------
# Load lineage results.
# There will be sequences with no lineage result. This is because they had no genes.
lineage_results = pd.read_csv("data/viral_annotation/imgvr_lineage/sequence_resolved_lineages.tsv", sep = "\t")
summary_table = summary_table.merge(lineage_results,  how='left', left_on = "Sequence_ID", right_on = "Sequence")

# Write out a new sequence file and annotate with the consensus lineage as description
with open(f"{viral_summary_dir}/viral_sequences.fasta", 'w') as fh:
    for seq_record in viral_sequences:
        taxonomy=summary_table.loc[summary_table['Sequence_ID'] == seq_record.name]["Majority_taxonomy"]
        if pd.isna(taxonomy).values[0] == False:
             seq_record.description = taxonomy.values[0]
        SeqIO.write(seq_record, fh, 'fasta-2line')

# -------------------------------------------------------------------------------------------------------------------
# Load checkv results
checkv_results = pd.read_csv("data/viral_annotation/checkv/quality_summary.tsv", sep = "\t")
checkv_results = checkv_results.add_prefix("checkv_")
summary_table = summary_table.merge(checkv_results,  how='left', left_on = "Sequence_ID", right_on = "checkv_contig_id")

# -------------------------------------------------------------------------------------------------------------------
# Load cluster results
cluster_representive_lengths = pd.read_csv("data/viral_clustering/mcl/cluster_representatives_lengths.tsv", sep = "\t",
                                           names = ["Cluster", "Sequence_ID", "Length"],usecols = [0,1])
cluster_members_lengths = pd.read_csv("data/viral_clustering/mcl/clusters_all_members_lengths.tsv", sep = "\t",
                                      names = ["Cluster", "Sequence_ID", "Length"], usecols = [0,1])
cluster_representive_lengths["Sequence_ID"] = cluster_representive_lengths["Sequence_ID"].apply(os.path.basename)
cluster_members_lengths["Sequence_ID"] = cluster_members_lengths["Sequence_ID"].apply(os.path.basename)

cluster_members_lengths['Cluster_representative'] = \
    np.where(cluster_members_lengths.Sequence_ID.isin(cluster_representive_lengths.Sequence_ID), 'Yes', 'No')

summary_table = summary_table.merge(cluster_members_lengths,
                                    how='left', left_on = "Sequence_ID", right_on = "Sequence_ID")
summary_table["Cluster"] = summary_table["Cluster"].fillna("Singleton")
summary_table["Cluster_representative"] = summary_table["Cluster_representative"].fillna("Singleton")

# -------------------------------------------------------------------------------------------------------------------
# Write summary table
summary_table.to_csv(f"{viral_summary_dir}/viral_summary.tsv",sep = "\t")
summary_table.to_excel(f"{viral_summary_dir}/viral_summary.xlsx")

# TODO - summary table for all proteins?
# amrfinder_plus_N_hits
# abricate_card_N_hits
# abricate_vfdb_N_hits
