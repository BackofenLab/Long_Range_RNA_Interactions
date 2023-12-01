# 1. Get tuples from IntaRNA interactions
# 2. Take the middle of these interactions and extract sequences +-15 bases in both directions 

import pandas as pd
import ast
import os
from collections import defaultdict
import math

def make_fasta(l, output_name):
    with open(output_name, "w") as f:
        for i in l:
            f.write(f">{i[0]}\n")
            f.write(f"{i[1]}\n")
            f.write(f"\n")


def get_meme_sequences(param_df_path, inta_df_path, output_path):
    """Get the sequences from the following 4 sites:
      1. 5': -30 to +10 from CDS start
      2. 5': +20 to +60 in CDS
      3. 3': -25 to +5 at CM hit start
      4. 3': +5 to +35 from CM hit start
    param_df_path (str): Path to the parameter dataframe
    inta_df_path (str): Path to the dataframe resulting from IntaRNA
    """
    os.makedirs(output_path, exist_ok=True)
    list_site_1 = []
    list_site_2 = []
    list_site_3 = []
    list_site_4 = []
    param_df = pd.read_csv(param_df_path)
    inta_df = pd.read_csv(inta_df_path)
    merged_df = pd.merge(param_df, inta_df, on=["id"], how="inner")
    for index, row in merged_df.iterrows():
        if not "cm_hit_f" in row:
            raise Exception("Dataframe provided does not contain CM-search hits")
        elif math.isnan(row["cm_hit_f"]):
            continue
        seq5 = row["seq5"]
        seq3 = row["seq3"]
        CDS_start = row["UTR5len_x"]
        CMhit_start = int(row["cm_hit_f"]) - row["UTR3len_x"]
        part1 = seq5[CDS_start-30:CDS_start+10]
        part2 = seq5[CDS_start+20:CDS_start+60]
        part3 = seq3[CMhit_start-25:CMhit_start+5]
        part4 = seq3[CMhit_start+5:CMhit_start+35]
        list_site_1.append((f"{row['class_x']}-{row['id']}", part1))
        list_site_2.append((f"{row['class_x']}-{row['id']}", part2))
        list_site_3.append((f"{row['class_x']}-{row['id']}", part3))
        list_site_4.append((f"{row['class_x']}-{row['id']}", part4))
    make_fasta(list_site_1, f"{output_path}/site_1.fa")
    make_fasta(list_site_2, f"{output_path}/site_2.fa")
    make_fasta(list_site_3, f"{output_path}/site_3.fa")
    make_fasta(list_site_4, f"{output_path}/site_4.fa")


def old_get_meme_sequences(param_df_path, inta_df_path, output_path, extra_bases):
    """Get the sequences of the main interaction.
    param_df_path (str): Path to the parameter dataframe
    inta_df_path (str): Path to the dataframe resulting from IntaRNA
    extra_bases (int): Extra bases that were used for IntaRNA. 
                       Neccessary because they need to be added to the 3' ranges
    """
    os.makedirs(output_path, exist_ok=True)
    dict5 = defaultdict(list)
    dict3 = defaultdict(list)
    param_df = pd.read_csv(param_df_path)
    inta_df = pd.read_csv(inta_df_path)
    merged_df = pd.merge(param_df, inta_df, on=["id"], how="inner")
    for index, row in merged_df.iterrows():
        t_range = ast.literal_eval(row["t_inter_range"])
        q_range = ast.literal_eval(row["q_inter_range"])
        UTR5len = row["UTR5len_x"]
        center_5 = round((t_range[1]+UTR5len + t_range[0]+UTR5len-1)/2)
        center_3 = round((q_range[0]+extra_bases-1 + q_range[1]+extra_bases)/2)
        dict5[row["class_x"]].append((f"{row['class_x']}-{row['id']}", row["seq5"][center_5 - 15: center_5 + 15]))
        dict3[row["class_x"]].append((f"{row['class_x']}-{row['id']}", row["seq3"][center_3 - 15: center_3 + 15]))
    for i5 in dict5:
        make_fasta(dict5[i5], f"{output_path}/{i5}_5.fa")
    for i3 in dict3:
        make_fasta(dict3[i3], f"{output_path}/{i3}_3.fa")
    
        #raise