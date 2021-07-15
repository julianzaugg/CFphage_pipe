import sys
import subprocess
import pathlib
import glob
import os
import re

import pandas as pd
from Bio import SeqIO
from Bio.SeqIO.FastaIO import SimpleFastaParser

viral_summary_dir = "data/viral_summary"

try:
    os.makedirs(viral_summary_dir)
except FileExistsError:
    pass

# Load viral sequences from FASTA file. Just need the names and lengths
# viral_sequences = list(SeqIO.parse("data/viral_predict/all_samples_viral_sequences.fasta"))
with open("data/viral_predict/all_samples_viral_sequences.fasta") as fasta_file:  # Will close handle cleanly
    seq_ids = []
    seq_lengths = []
    for title, sequence in SimpleFastaParser(fasta_file):
        seq_ids.append(title.split(None, 1)[0])  # First word is ID
        seq_lengths.append(len(sequence))


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
