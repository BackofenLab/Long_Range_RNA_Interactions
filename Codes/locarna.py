import pandas as pd
import os
import math
import subprocess
import collections
from Codes.locarna_help import run_mlocarna, run_rnaalifold, run_ps_to_pdf


def make_locarna_fasta(l, output_name, CDS_left, CDS_right):
    """
    Create a FASTA file out of a given list of sequences
    and write it into a file.
    The list should have the form [(seq_name1, part_5_1, part_3_1), 
                                   (seq_name2, part_5_2, part_3_2),...]
    """
    with open(output_name, "w") as f:
        for i in l:
            if len(i[1]) > 0 and len(i[2]) > 0:
                f.write(f">{i[0]}\n")
                f.write(f"{i[1]}NNNNNNN{i[2]}\n")
                f.write(f"{(len(i[1]))*'<'}xxxxxxx{len(i[2])*'>'} #S\n")
                f.write(f"{CDS_left*'.'}AAA{(CDS_right-3)*'.'}BBBBBBB{len(i[2])*'.'} #1\n")
                f.write(f"{CDS_left*'.'}123{(CDS_right-3)*'.'}1234567{len(i[2])*'.'} #2\n")
                f.write(f"\n")





def main_locarna(param_df_path, cm_path, output_path,
                 CDS_left, CDS_right, CMHit_left, CMHit_right):
    """Run locARNA with RNAalifold.
    param_df_path (str): Path to the parameter dataframe
    cm_path (str): Path to the dataframe resulting from CMSearch
    output_path (int): Output directory for the extracted sequences
    """
    os.makedirs(output_path, exist_ok=True)
    seq_dir = collections.defaultdict(list)
    param_df = pd.read_csv(param_df_path)
    cm_results = pd.read_csv(cm_path)  ## Take the CM-hit starts from previous CM dataframe
    param_df["cm_hit_f"] = pd.Series(cm_results["cm_hit_f"])  ## Insert them into our dataframe

    for index, row in param_df.iterrows():
        if not "cm_hit_f" in row:
            raise Exception("Dataframe provided does not contain CM-search hits")
        elif math.isnan(row["cm_hit_f"]):
            continue
        seq5 = row["seq5"]
        seq3 = row["seq3"]
        CDS_start = row["UTR5len"]
        CMhit_start = int(row["cm_hit_f"]) - row["UTR3len"]
        part5 = seq5[CDS_start-CDS_left:CDS_start+CDS_right]
        part3 = seq3[CMhit_start-CMHit_left:CMhit_start+CMHit_right]
        seq_dir["all"].append((f"{row['class']}-{row['id']}", part5, part3))
        if row['class'] != "ISFV":
            seq_dir[row['class']].append((f"{row['class']}-{row['id']}", part5, part3))
        else: # Separate cISFV and dISFV
            group_name = row['type'][:-1]
            seq_dir[group_name].append((f"{group_name}-{row['id']}", part5, part3))
        #if row['class'] == "TBFV" or row['type'][:-1] == "dISFV": 
    # dISFV + TBFV alignment
    seq_dir["dISFV+TBFV"] = seq_dir["dISFV"] + seq_dir["TBFV"]
    seq_dir["MBFV+dISFV"] = seq_dir["MBFV"] + seq_dir["dISFV"]
    seq_dir["MBFV+TBFV"] = seq_dir["MBFV"] + seq_dir["TBFV"]
    for seq_class in seq_dir:
        make_locarna_fasta(seq_dir[seq_class], f"{output_path}/locARNA_{seq_class}_input.fa", CDS_left, CDS_right)
        run_mlocarna(f"{output_path}/locARNA_{seq_class}_input.fa", f"{output_path}/{seq_class}")
        run_rnaalifold(f"{output_path}/{seq_class}/results")

        run_ps_to_pdf(f"{output_path}/{seq_class}/results/alirna.ps", f"{output_path}/{seq_class}_alirna.pdf")
        run_ps_to_pdf(f"{output_path}/{seq_class}/results/aln.ps", f"{output_path}/{seq_class}_aln.pdf")