import os
import pandas as pd
import collections
import math
from Codes.locarna_help import run_mlocarna, run_rnaalifold, run_ps_to_pdf, run_mrri


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
                f.write(f".{(len(i[1])-1)*'<'}xxxxxxx{len(i[2])*'>'} #S\n")
                f.write(f"{CDS_left*'.'}AAA{(CDS_right-3)*'.'}BBBBBBB{len(i[2])*'.'} #1\n")
                f.write(f"{CDS_left*'.'}123{(CDS_right-3)*'.'}1234567{len(i[2])*'.'} #2\n")
                f.write(f"{i[3]} #FS\n")
                f.write(f"\n")
    
def get_mrri_file(parameter_table_file, static_param_path, extra_bases, extra_bases_roi, 
                  mrri_file_output, raw_mrri_output):
    params = pd.read_csv(parameter_table_file)
    additional_args = []
    output = pd.DataFrame()
    hybridDPs = []
    start1 = []
    end1 = []
    start2 = []
    end2 = []
    hybridDPs = []
    with open(raw_mrri_output, "w") as f:
        for index, row in params.iterrows():
            print(f"MRRI: {row['id']}")
            f.write(f"{row['id']} :\n")
            f.write(f"{'#'*(len(row['id']) + 2)}\n")
            stdout, d = run_mrri(row["seq5"], row["seq3"], row['id'],
                                        extra_bases, extra_bases_roi, static_param_path)
            f.write(f"{stdout}\n")
            hybridDPs.append(d["hybridDP"])
            start1.append(d['start1'])
            end1.append(d['end1'])
            start2.append(d['start2'])
            end2.append(d['end2'])
    output = params
    output["hybridDP"] = hybridDPs
    output["start1"] = start1
    output["end1"] = end1
    output["start2"] = start2
    output["end2"] = end2
    output.to_csv(mrri_file_output, index=False)

def main_loc_with_mrri(mrri_file_path, cm_path,
                       parameter_table_file, output_path):
    CDS_left = 30
    CDS_right = 70
    CMHit_left = 30
    CMHit_right = 30
    os.makedirs(output_path, exist_ok=True)

    mrri_file = pd.read_csv(mrri_file_path)
    cm_results = pd.read_csv(cm_path)

    merged_df = pd.merge(mrri_file, cm_results, on=["id"], how="inner")
    seq_dir = collections.defaultdict(list)
    for index, row in merged_df.iterrows():
        if not "cm_hit_f" in row:
            raise Exception("Dataframe provided does not contain CM-search hits")
        elif math.isnan(row["cm_hit_f"]):
            continue
        seq5 = row["seq5"]
        seq3 = row["seq3"]
        CDS_start = row["UTR5len_x"]
        CMhit_start = int(row["cm_hit_f"]) - row["UTR3len_x"]
        part5 = seq5[CDS_start-CDS_left:CDS_start+CDS_right]
        part3 = seq3[CMhit_start-CMHit_left:CMhit_start+CMHit_right]

        hybridDP_split = row['hybridDP'].split("&")
        FS_seq_5 = ("."*row["start1"] + hybridDP_split[0] + "."*(len(seq5)-row["end1"]))[CDS_start-CDS_left:CDS_start+CDS_right]
        FS_seq_3 = ("."*row["start2"] + hybridDP_split[1] + "."*(len(seq3)-row["end2"]))[CMhit_start-CMHit_left:CMhit_start+CMHit_right]
        FS_seq = FS_seq_5 + "NNNNNNN" + FS_seq_3 

        seq_dir["all"].append((f"{row['class_x']}-{row['id']}", part5, part3, FS_seq))
        if row['class_x'] != "ISFV":
            seq_dir[row['class_x']].append((f"{row['class_x']}-{row['id']}", part5, part3, FS_seq))
        else: # Separate cISFV and dISFV
            group_name = row['type_x'][:-1]
            seq_dir[group_name].append((f"{group_name}-{row['id']}", part5, part3, FS_seq))
        if row['class_x'] == "TBFV" or row['type_x'][:-1] == "dISFV": 
            # dISFV + TBFV alignment
            seq_dir["dISFV+TBFV"].append((f"{'dISFV+TBFV'}-{row['id']}", part5, part3, FS_seq))
    for seq_class in seq_dir:
        make_locarna_fasta(seq_dir[seq_class], f"{output_path}/locARNA_{seq_class}_input.fa", CDS_left, CDS_right)
        #raise
        run_mlocarna(f"{output_path}/locARNA_{seq_class}_input.fa", f"{output_path}/{seq_class}")
        run_rnaalifold(f"{output_path}/{seq_class}/results")
        
        run_ps_to_pdf(f"{output_path}/{seq_class}/results/alirna.ps", f"{output_path}/{seq_class}_alirna.pdf")
        run_ps_to_pdf(f"{output_path}/{seq_class}/results/aln.ps", f"{output_path}/{seq_class}_aln.pdf")