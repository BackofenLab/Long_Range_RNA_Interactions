# Translate start of CDS (first 100 bases) to proteins

# Step 1: Read ids and seq5s from parameter_table.csv
# Step 2: Translate seq5 to protein codings
# Step 3: Align protein coding?

import pandas as pd
from Bio.Seq import Seq

def make_fasta(l, output_name):
    """
    Create a FASTA file out of a given list of sequences
    and write it into a file.
    The list should have the form [(seq_name1, part_5_1, part_3_1), 
                                   (seq_name2, part_5_2, part_3_2),...]
    """
    with open(output_name, "w") as f:
        for i in l:
            if len(i[1]) > 0:
                f.write(f">{i[0]}\n")
                f.write(f"{i[1]}\n")
                f.write(f"\n")

def codons_to_amino_acids(s):
    seq = Seq(s[:-(len(s)%3)]) # cut off trailing bases (eg: len 100 => cut of last base)
    acids = str(seq.translate())
    return acids


def cds_to_proteins(param_table_file, output_path, extra_bases_roi):
    params = pd.read_csv(param_table_file)
    output = pd.DataFrame()
    acid_sequences = []
    for index, row in params.iterrows():
        row["id"]
        UTR5len = row["UTR5len"]
        seq5 = row["seq5"][UTR5len:UTR5len+extra_bases_roi]
        #seq5 = row["seq5"][-200:-extra_bases_roi] # Same thing
        
        acid_sequences.append((f"{row['class']}-{row['id']}", codons_to_amino_acids(seq5)))
    make_fasta(acid_sequences, output_path)