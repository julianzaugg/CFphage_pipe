import sys
import subprocess
import pathlib
import glob
import os

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

# Create summary table from IDs and lengths
summary_table = pd.DataFrame({"Seq_ID" : seq_ids, "Seq_length" : seq_lengths})

# Determine the Sample, assembly tool and viral prediction tool from the name

# Get list of names
# viral_ids = [seq_record.name for seq_record in viral_sequences]

# Make dictionary, add each sequence
# from seqID, resolve Batch_ID, Barcode_ID, Assembly_method, viral_tool
# data = {'row_1': [3, 2, 1, 0], 'row_2': ['a', 'b', 'c', 'd']}
# >>> pd.DataFrame.from_dict(data, orient='index')

# Predicted_from_AF_reads (yes/no)
# Split on "____" and
# >MN200521_batch5_barcode01__metaflye__virsorter____contig_95||0_partial
# >MN200521_batch5_barcode01__metaflye__virsorter____4d63269c-53e4-455c-a732-19bf0959b096||full