# 1. Get tuples from IntaRNA interactions
# 2. Extract the bases of those ranges from the original sequences 

import pandas as pd
import ast

def get_alignment_sequences(param_df_path, inta_df_path, extra_bases):
    """Get the sequences neccessary for the alignment
    param_df_path (str): Path to the parameter dataframe
    inta_df_path (str): Path to the dataframe resulting from IntaRNA
    extra_bases (int): Extra bases that were used for IntaRNA. 
                       Neccessary because they need to be added to the 3' ranges
    """
    param_df = pd.read_csv(param_df_path)
    inta_df = pd.read_csv(inta_df_path)
    merged_df = pd.merge(param_df, inta_df, on=["id", "UTR5len"], how="inner")
    for index, row in merged_df.iterrows():
        t_range = ast.literal_eval(row["t_inter_range"])
        q_range = ast.literal_eval(row["q_inter_range"])
        UTR5len = row["UTR5len"]
        int_seq_5 = row["seq5"][t_range[0]+UTR5len-1:t_range[1]+UTR5len]
        int_seq_3 = row["seq3"][q_range[0]+extra_bases-1:q_range[1]+extra_bases]
        #print(row["id"])
        #print(row["t_inter_range"])
        #print(row["q_inter_range"])
        #print(int_seq_5)
        #print(int_seq_3)