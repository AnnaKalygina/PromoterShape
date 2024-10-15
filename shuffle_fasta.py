#!/usr/bin/env python3

from Bio import SeqIO
import random
import argparse

parser = argparse.ArgumentParser()
parser.add_argument("input_dir")
parser.add_argument("output_dir")
args = parser.parse_args()

input_fasta = args.input_dir  
output_fasta = args.output_dir

def shuffle_sequence(seq_record):
    seq = list(seq_record.seq) 
    random.shuffle(seq)
    seq_record.seq = "".join(seq)
    return seq_record

sequences = list(SeqIO.parse(input_fasta, "fasta"))
shuffled_sequences = [shuffle_sequence(seq_record) for seq_record in sequences]

with open(output_fasta, "w") as output_handle:
    SeqIO.write(shuffled_sequences, output_handle, "fasta")

print(f"Nucleotide-shuffled sequences written to {output_fasta}")