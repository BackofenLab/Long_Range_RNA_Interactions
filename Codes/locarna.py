import pandas as pd
import os
import math
import subprocess
import collections

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
                f.write(f"\n")


def run_mlocarna(input_fasta, output_dir):
    """Run mlocarna on a given fasta file.
    
    input_fasta (str): Filepath to a fasta file to apply locarna on
    output_dir (str): Filepath of the resulting CM file
    """
    print(f"mlocarna: {input_fasta} `=> {output_dir}")
    cmd = ["mlocarna", input_fasta,
           #"--indel=-50", # Webserver parameter
           #"--indel-opening=-750", # Webserver parameter
           "--width=300",
           "--tgtdir", output_dir
           ]
    p = subprocess.Popen(cmd, stdout=subprocess.PIPE)
    p.wait()

def run_rnaalifold(input_dir):
    print(f"RNAalifold : {input_dir}/result.aln => {input_dir}/alirna.ps+aln.ps")
    input_file = f"{input_dir}/result.aln"
    with open(input_file, "r") as f:
        for line in f.readlines():
            if line.startswith("#A1"):
                anchor_seq = line.split(" ")[-1]
                s1, s2 = anchor_seq.split("BBBBBBB")
                constraint = b"."+b"<"*(len(s1)-1) + b"xxxxxxx" + b">"*(len(s2)-1)
                break
    cmd = ["RNAalifold", input_file,
           "--aln", "--color", "-C"
           ]
    p = subprocess.Popen(cmd, stdout=subprocess.PIPE, stdin=subprocess.PIPE, stderr=subprocess.PIPE)
    p.communicate(input=constraint)
    p.wait()
    os.rename("alirna.ps", f"{input_dir}/alirna.ps")
    os.rename("aln.ps", f"{input_dir}/aln.ps")

def run_ps_to_pdf(ps_file, output):
    print(f"ps2pdf: {ps_file} `=> {output}")
    cmd = ["ps2pdf", "-dEPSCrop", ps_file, output]
    p = subprocess.Popen(cmd, stdout=subprocess.PIPE)
    p.wait()
    
    
def get_locarna_sequences(param_df_path, inta_df_path, output_path):
    """Get the sequences neccessary for the locARNA alignment
    param_df_path (str): Path to the parameter dataframe
    inta_df_path (str): Path to the dataframe resulting from IntaRNA
    output_path (int): Output directory for the extracted sequences
    """
    CDS_left = 30
    CDS_right = 70
    CMHit_left = 30
    CMHit_right = 30
    os.makedirs(output_path, exist_ok=True)
    seq_dir = collections.defaultdict(list)
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
        part5 = seq5[CDS_start-CDS_left:CDS_start+CDS_right]
        part3 = seq3[CMhit_start-CMHit_left:CMhit_start+CMHit_right]
        seq_dir["all"].append((f"{row['class_x']}-{row['id']}", part5, part3))
        if row['class_x'] != "ISFV":
            seq_dir[row['class_x']].append((f"{row['class_x']}-{row['id']}", part5, part3))
        else: # Separate cISFV and dISFV
            group_name = row['type_x'][:-1]
            seq_dir[group_name].append((f"{group_name}-{row['id']}", part5, part3))
    for seq_class in seq_dir:
        make_locarna_fasta(seq_dir[seq_class], f"{output_path}/locARNA_{seq_class}_input.fa", CDS_left, CDS_right)
        run_mlocarna(f"{output_path}/locARNA_{seq_class}_input.fa", f"{output_path}/{seq_class}")
        run_rnaalifold(f"{output_path}/{seq_class}/results")
        
        run_ps_to_pdf(f"{output_path}/{seq_class}/results/alirna.ps", f"{output_path}/{seq_class}_alirna.pdf")
        run_ps_to_pdf(f"{output_path}/{seq_class}/results/aln.ps", f"{output_path}/{seq_class}_aln.pdf")