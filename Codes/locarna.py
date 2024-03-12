import os
import pandas as pd
import collections
import math
import ast
from Codes.locarna_help import (run_mlocarna, run_rnaalifold, run_ps_to_pdf, 
                                find_all, integrate_into, cm_compare,
                                get_special_combined_constraint)
import string


def make_locarna_fasta(l, output_name, skip_FS=False):
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
                f.write(f"{i[1]}NNNNNNN{i[2]}\n") ## part5NNNNNNNpart3
                f.write(f"{i[3]} #S\n")           ## cons_S
                f.write(f"{i[4]} #1\n")           ## cons_1
                f.write(f"{i[5]} #2\n")           ## cons_2
                if not skip_FS:
                    f.write(f"{i[6]} #FS\n")          ## cons_FS
                f.write(f"\n")


def main_locarna(input_df_file_path, cm_path, output_path,
                 CDS_left, CDS_right, 
                 CMHit_left, CMHit_right, 
                 use_carna=False, skip_FS=False, mode=1):
    """Run locARNA with RNAalifold.
    param_df_path (str): Path to the parameter dataframe
    cm_path (str): Path to the dataframe resulting from CMSearch
    output_path (int): Output directory for the extracted sequences
    """
    os.makedirs(output_path, exist_ok=True)

    input_df = pd.read_csv(input_df_file_path)
    cm_results = pd.read_csv(cm_path)  ## Take the CM-hits from previous CM dataframe
    input_df["cm_hit_f"] = pd.Series(cm_results["cm_hit_f"])  ## Insert them into our dataframe
    input_df["cm_hit_t"] = pd.Series(cm_results["cm_hit_t"])
    input_df["cm_hit_src"] = pd.Series(cm_results["cm_hit_src"])
    input_df["align_cons_3SL"] = pd.Series(cm_results["align_cons_3SL"])

    seq_dir = collections.defaultdict(list)
    for index, row in input_df.iterrows():
        if not "cm_hit_f" in row:
            raise Exception("Dataframe provided does not contain CM-search hits")
        elif math.isnan(row["cm_hit_f"]): ## Skip sequences without CMHits
            continue
        ranges_t = []
        ranges_q = []
        ranges_t += ast.literal_eval(row["predictions_t"])
        ranges_q += ast.literal_eval(row["predictions_q"])
        seq5 = row["seq5"]
        seq3 = row["seq3"]
        cm_hit_f = int(row["cm_hit_f"])
        cm_hit_t = int(row["cm_hit_t"])
        CDS_start = row["UTR5len"]
        #CMhit_start = cm_hit_f - row["UTR3len"]
        align_cons_3SL = row["align_cons_3SL"]
        
        ## How much to cut of from the 5' and 3' sequences and their constraints
        cutoff_5_left = CDS_start-CDS_left
        cutoff_5_right = CDS_start+CDS_right
        cutoff_3_left = len(seq3)-141
        cutoff_3_right = None ##### Old: len(seq3)-49
        
        ## Cut sequences:
        part5 = seq5[cutoff_5_left:cutoff_5_right]
        part3 = seq3[cutoff_3_left:cutoff_3_right]

        ## #FS Constraints:
        cons_FS_5 = "."*(len(seq5)-100)        
        cons_FS_3 = "."*(len(seq3))
        for i in range(len(ranges_t)): # ranges_t[i][0] = Start_index, [3] = hybridDP
            cons_FS_5 = integrate_into(cons_FS_5, "."*(int(ranges_t[i][0])+len(seq5)-200) + ranges_t[i][3], string.ascii_uppercase[i])
        for i in range(len(ranges_q)): # ranges_q[i][0] = Start_index, [3] = hybridDP
            cons_FS_3 = integrate_into(cons_FS_3, "."*(int(ranges_q[i][0])+199) + ranges_q[i][3], string.ascii_lowercase[i])
        cons_FS_5 = cons_FS_5[cutoff_5_left:cutoff_5_right]
        cons_FS_3 = cons_FS_3[cutoff_3_left:cutoff_3_right]
        cons_FS = cons_FS_5 + "xxxxxxx" + cons_FS_3
        
        ## #S constraint:
        if mode == 0:   ## None:    .......xxxxxxx.......
            cons_S = f"{'.'*len(part5)}xxxxxxx{'.'*len(part3)}"
        elif mode == 1: ## Classic: <<<<<<<xxxxxxx>>>>>>>
            cons_S = f"{'<'*len(part5)}xxxxxxx{'>'*len(part3)}"
        elif mode == 2: ## everything blocked ("x"), except where the cmhit is (there ".")
            covar_3SL = "x"*(cm_hit_f+199) + "."*(cm_hit_t - cm_hit_f + 1)
            covar_3SL += "x"*max(0, len(seq3) - len(covar_3SL)) ## Extend remaining dots
            covar_3SL = covar_3SL[:len(seq3)] ## One specific sequence has cm_t LARGER than the sequence itself so..
            covar_3SL = covar_3SL[cutoff_3_left:cutoff_3_right]
            cons_S = f"{(len(part5))*'x'}xxxxxxx{covar_3SL}"
            skip_FS = True  ## Enforce #FS to be skipped since #S is otherwise ignored
        elif mode == 3: ## #S constraint using interactions (using the #FS sequence)
            # Only take C and c if they are nested with A and B
            cons_S = str(cons_FS)
            first_A = cons_S.find("A")
            first_B = cons_S.find("B")
            first_C = cons_S.find("C")
            first_a = cons_S.find("a")
            first_b = cons_S.find("b")
            first_c = cons_S.find("c")
            cons_S = cons_S.replace("A", "(").replace("a", ")").replace("B", "(").replace("b", ")")
            if ((first_C < first_A and first_C < first_B and first_c > first_a and first_c > first_b) or # CAB bac and CBA abc
               (first_C < first_A and first_C > first_B and first_c > first_a and first_c < first_b) or  # BCA acb
               (first_C > first_A and first_C < first_B and first_c < first_a and first_c > first_b) or  # ACB bca
               (first_C > first_A and first_C > first_B and first_c < first_a and first_c < first_b)):   # ABC cba and BAC cab
                cons_S = cons_S.replace("C", "(").replace("c", ")")
            else: # C is not nested with A and B so it gets ignored
                cons_S = cons_S.replace("C", ".").replace("c", ".")
            #cons_S = cons_S.replace(".", "x")
            cons_FS = cons_S
            cons_S = f"{'.'*len(part5)}xxxxxxx{'.'*len(part3)}"
            skip_FS = False ## Enforce #FS to be actually used in this mode
            cons_FS = get_special_combined_constraint(f"{part5}NNNNNNN{part3}", cons_FS)
        else:
            raise ValueError("Invalid mode for #S constraint")
        ## Constraint 1/2:
        cons_1 = f"{CDS_left*'.'}AAA{(CDS_right-3)*'.'}BBBBBBB{len(part3)*'.'}"
        cons_2 = f"{CDS_left*'.'}123{(CDS_right-3)*'.'}1234567{len(part3)*'.'}"
        if not ((cons_FS.count("(") == cons_FS.count(")"))
           and (cons_FS.count("A") == cons_FS.count("a"))
           and (cons_FS.count("B") == cons_FS.count("b"))
           and (cons_FS.count("C") == cons_FS.count("c"))):
           print(f"FS constraint unbalanced for {row['id']}")

        # seq_dir["all"].append((f"{row['class']}-{row['id']}", part5, part3, cons_S, cons_1, cons_2, cons_FS))
        if row['class'] != "ISFV":
            seq_dir[row['class']].append([f"{row['class']}-{row['id']}", part5, part3, cons_S, cons_1, cons_2, cons_FS])
        else: # Separate cISFV and dISFV
            group_name = row['type'][:-1]
            seq_dir[group_name].append([f"{group_name}-{row['id']}", part5, part3, cons_S, cons_1, cons_2, cons_FS])
        if row['id'] == "NC_009942.1": ###### TEMPORARY FOR TESTING
            seq_dir[f"single_{row['virus']}"].append([f"{row['virus']}-{row['id']}", part5, part3, cons_S, cons_1, cons_2, cons_FS])
            seq_dir[f"single_{row['virus']}"].append([f"{row['virus']}-{row['id']}a", part5, part3, cons_S, cons_1, cons_2, cons_FS])   
    # Combined alignments
    seq_dir["dISFV+TBFV"] = seq_dir["dISFV"] + seq_dir["TBFV"]
    seq_dir["MBFV+dISFV"] = seq_dir["MBFV"] + seq_dir["dISFV"]
    seq_dir["MBFV+TBFV"] = seq_dir["MBFV"] + seq_dir["TBFV"]

    for seq_class in seq_dir:
        fasta_file_name = f"{output_path}/locARNA_{seq_class}_input.fa"
        make_locarna_fasta(seq_dir[seq_class], fasta_file_name, skip_FS=skip_FS)
        run_mlocarna(fasta_file_name, f"{output_path}/{seq_class}", use_carna)
        run_rnaalifold(f"{output_path}/{seq_class}/results", seq_dir[seq_class], mode=mode, locARNA_input=fasta_file_name)
        run_ps_to_pdf(f"{output_path}/{seq_class}/results/alirna.ps", f"{output_path}/{seq_class}_alirna.pdf")
        run_ps_to_pdf(f"{output_path}/{seq_class}/results/aln.ps", f"{output_path}/{seq_class}_aln.pdf")