import os
import pandas as pd
import collections
import math
import ast
from Codes.locarna_help import (run_mlocarna, run_rnaalifold, run_ps_to_pdf, 
                                hacked_MRRI_main, find_all, integrate_into)
import string


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
                f.write(f"{i[3]} #FS\n")
                f.write(f"\n")
    
    
def main_mrri(parameter_table_file, static_param_path, extra_bases, extra_bases_roi, 
              mrri_file_output, raw_mrri_output, param_mode):
    params = pd.read_csv(parameter_table_file)
    additional_args = []
    output = pd.DataFrame()
    hybridDPs = []
    t_opt_ranges = [] ## Interaction regions of the main interaction of all sequences
    q_opt_ranges = []
    t_constraint_ranges = []  ## All constrained interactions of all sequences
    q_constraint_ranges = []
    hybridDPs = []
    with open(raw_mrri_output, "w") as f_raw:
        for index, row in params.iterrows():
            print(f"MRRI: {row['id']}")
            f_raw.write(f"{row['id']} :\n")
            f_raw.write(f"{'#'*(len(row['id']) + 2)}\n")
            main_interaction, c_interactions = hacked_MRRI_main(row["seq5"], row["seq3"], static_param_path, param_mode)
            f_raw.write(f"{main_interaction}\n")
            f_raw.write(f"{c_interactions}\n")
            hybridDPs.append(main_interaction["hybridDP"].replace("(","A").replace(")","a"))
            t_opt_ranges.append((int(main_interaction['start1']), int(main_interaction['end1'])))
            q_opt_ranges.append((int(main_interaction['start2']), int(main_interaction['end2'])))

            c_inter_ts = [] ## Other constrained interactions of a sequence
            c_inter_qs = []
            for i in range(0, len(c_interactions)):
                interaction = c_interactions[i]
                hybDP_0, hybDP_1 = interaction["hybridDP"].split("&")
                hybDP_0 = hybDP_0.replace("(", string.ascii_uppercase[i+1])
                hybDP_1 = hybDP_1.replace(")", string.ascii_lowercase[i+1])
                c_inter_ts.append((interaction["start1"], interaction["end1"], hybDP_0))
                c_inter_qs.append((interaction["start2"], interaction["end2"], hybDP_1))
            t_constraint_ranges.append(c_inter_ts)
            q_constraint_ranges.append(c_inter_qs)
    output = params
    output["hybridDP"] = hybridDPs
    output["t_inter_range"] = t_opt_ranges
    output["q_inter_range"] = q_opt_ranges
    output["constrained_predictions_t"] = t_constraint_ranges
    output["constrained_predictions_q"] = q_constraint_ranges
    output.to_csv(mrri_file_output, index=False)

def main_loc_with_mrri(mrri_file_path, cm_path,
                       parameter_table_file, output_path,
                       CDS_left, CDS_right, CMHit_left, CMHit_right):
    os.makedirs(output_path, exist_ok=True)

    mrri_df = pd.read_csv(mrri_file_path)

    cm_results = pd.read_csv(cm_path)  ## Take the CM-hit starts from previous CM dataframe
    mrri_df["cm_hit_f"] = pd.Series(cm_results["cm_hit_f"])  ## Insert them into our dataframe

    seq_dir = collections.defaultdict(list)
    for index, row in mrri_df.iterrows():
        if not "cm_hit_f" in row:
            raise Exception("Dataframe provided does not contain CM-search hits")
        elif math.isnan(row["cm_hit_f"]): ## Skip sequences without CMHits
            continue
        start1, end1 = ast.literal_eval(row["t_inter_range"])
        start2, end2 = ast.literal_eval(row["q_inter_range"])
        constrained_t = ast.literal_eval(row["constrained_predictions_t"])
        constrained_q = ast.literal_eval(row["constrained_predictions_q"])
        seq5 = row["seq5"]
        seq3 = row["seq3"]
        CDS_start = row["UTR5len"]
        CMhit_start = int(row["cm_hit_f"]) - row["UTR3len"]
        part5 = seq5[CDS_start-CDS_left:CDS_start+CDS_right]
        #part3 = seq3[CMhit_start-CMHit_left:CMhit_start+CMHit_right]
        part3 = seq3[len(seq3)-150:len(seq3)-50] #####

        hybridDP_split = row['hybridDP'].split("&")
        
        FS_seq_5 = integrate_into("."*(len(seq5)-100), "."*(start1+len(seq5)-200) + hybridDP_split[0] + "."*(99-end1) ,"A")
        for i in range(len(constrained_t)):
            FS_seq_5 = integrate_into(FS_seq_5, "."*(int(constrained_t[i][0])+len(seq5)-200) + constrained_t[i][2] + "."*(99-int(constrained_t[i][1])), string.ascii_uppercase[i+1])
        FS_seq_5 = FS_seq_5[CDS_start-CDS_left:CDS_start+CDS_right]
        FS_seq_3 = integrate_into("."*(len(seq3)), "."*(start2+200) + hybridDP_split[1] + "."*(len(seq3)-(end2+200)),"a")
        for i in range(len(constrained_q)):
            FS_seq_3 = integrate_into(FS_seq_3, "."*(int(constrained_q[i][0])+200) + constrained_q[i][2] + "."*(len(seq3)-(int(constrained_q[i][1])+200)), string.ascii_lowercase[i+1])

        #FS_seq_3 = FS_seq_3[CMhit_start-CMHit_left:CMhit_start+CMHit_right])
        FS_seq_3 = FS_seq_3[len(FS_seq_3)-140:len(FS_seq_3)-49] #####
        FS_seq = FS_seq_5 + "NNNNNNN" + FS_seq_3

        seq_dir["all"].append((f"{row['class']}-{row['id']}", part5, part3, FS_seq))
        if row['class'] != "ISFV":
            seq_dir[row['class']].append((f"{row['class']}-{row['id']}", part5, part3, FS_seq))
        else: # Separate cISFV and dISFV
            group_name = row['type'][:-1]
            seq_dir[group_name].append((f"{group_name}-{row['id']}", part5, part3, FS_seq))
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