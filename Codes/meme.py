# 1. Get tuples from IntaRNA interactions
# 2. Take the middle of these interactions and extract sequences +-15 bases in both directions 

import pandas as pd
import ast
import os
from collections import defaultdict

def make_fasta(l, output_name):
    with open(output_name, "w") as f:
        for i in l:
            f.write(f">{i[0]}\n")
            f.write(f"{i[1]}\n")
            f.write(f"\n")


def get_meme_sequences(param_df_path, output_path, inta_df_path, extra_bases):
    """Get the sequences neccessary for the alignment
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
    merged_df = pd.merge(param_df, inta_df, on=["id", "UTR5len"], how="inner")
    for index, row in merged_df.iterrows():
        t_range = ast.literal_eval(row["t_inter_range"])
        q_range = ast.literal_eval(row["q_inter_range"])
        UTR5len = row["UTR5len"]
        center_5 = round((t_range[1]+UTR5len + t_range[0]+UTR5len-1)/2)
        center_3 = round((q_range[0]+extra_bases-1 + q_range[1]+extra_bases)/2)
        dict5[row["class_x"]].append((row["id"], row["seq5"][center_5 - 15: center_5 + 15]))
        dict3[row["class_x"]].append((row["id"], row["seq3"][center_3 - 15: center_3 + 15]))
    for i5 in dict5:
        make_fasta(dict5[i5], f"{output_path}/{i5}_5.fa")
    for i3 in dict3:
        make_fasta(dict3[i3], f"{output_path}/{i3}_3.fa")
    
        #raise