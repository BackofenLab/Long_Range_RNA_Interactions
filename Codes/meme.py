import pandas as pd
import os
import math

def make_fasta(l, output_name):
    """
    Create a FASTA file out of a given list of sequences
    and write it into a file.
    The list should have the form [(seq_name1, seq1), (seq_name2, seq2),...]
    """
    with open(output_name, "w") as f:
        for i in l:
            if len(i[1]) > 0:
                f.write(f">{i[0]}\n")
                f.write(f"{i[1]}\n")
                f.write(f"\n")


def get_meme_sequences(inta_df_path, output_path):
    """Get the sequences from the following 4 sites:
      1. 5': -40 to +0 from CDS start
      2. 5': +20 to +60 in CDS
      3. 3': -25 to +5 at CM hit start
      4. 3': +5 to +35 from CM hit start
    inta_df_path (str): Path to the dataframe resulting from IntaRNA
    output_path (str): Output directory for the extracted sequences
    """
    site1 = (-40, 0)  # Position from CDS start
    site2 = (20, 60)  # Position from CDS start
    site3 = (-25, 5)  # Position from CMHit start
    site4 = (5, 35)   # Position from CMHit start
    os.makedirs(output_path, exist_ok=True)
    list_site_1 = []
    list_site_2 = []
    list_site_3 = []
    list_site_4 = []
    list_all_5 = []
    list_all_3 = []
    inta_df = pd.read_csv(inta_df_path)
    for index, row in inta_df.iterrows():
        if not "cm_hit_f" in row:
            raise Exception("Dataframe provided does not contain CM-search hits")
        elif math.isnan(row["cm_hit_f"]):
            continue
        seq5 = row["seq5"]
        seq3 = row["seq3"]
        CDS_start = row["UTR5len"]
        CMhit_start = int(row["cm_hit_f"]) - row["UTR3len"]
        part1 = seq5[CDS_start+site1[0]:CDS_start+site1[1]]
        part2 = seq5[CDS_start+site2[0]:CDS_start+site2[1]]
        part3 = seq3[CMhit_start+site3[0]:CMhit_start+site3[1]]
        part4 = seq3[CMhit_start+site4[0]:CMhit_start+site4[1]]
        list_site_1.append((f"{row['class']}-{row['id']}", part1))
        list_site_2.append((f"{row['class']}-{row['id']}", part2))
        list_site_3.append((f"{row['class']}-{row['id']}", part3))
        list_site_4.append((f"{row['class']}-{row['id']}", part4))
        list_all_5.append((f"{row['class']}-{row['id']}", seq5))
        list_all_3.append((f"{row['class']}-{row['id']}", seq3))
    make_fasta(list_site_1, f"{output_path}/site_1.fa")
    make_fasta(list_site_2, f"{output_path}/site_2.fa")
    make_fasta(list_site_3, f"{output_path}/site_3.fa")
    make_fasta(list_site_4, f"{output_path}/site_4.fa")
    make_fasta(list_all_5, f"{output_path}/site_all_5.fa")
    make_fasta(list_all_3, f"{output_path}/site_all_3.fa")

