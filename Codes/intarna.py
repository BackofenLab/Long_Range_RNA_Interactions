import pandas as pd
import math
import subprocess
from Codes.intarna_help import hacked_MRRI_main
import string

def execute_IntaRNA(UTR5pCDS, UTR3pCDS, ID,
                    extra_bases, extra_bases_roi, static_param_path,
                    additional_args=[]):
    """Starts a subprocess for IntaRNA with given parameters.
    
    UTR5pCDS (str): 5'UTR+ the extra bases from the CDS
    UTR3pCDS (str): 3'UTR+ the extra bases from the CDS
    ID (str): ID code of the current sequence ("NC_XXXXX.X")
    static_param_path (str): Filepath for the static parameter file to be used
    additional_args (list): Optional additional arguments, must have the form:
                            ["--par1", str(val1), "--par2",...]
    """
    # roi_difference cuts off the difference between
    # the extra CDS bases and the extra CDS bases as region of interest
    roi_difference = extra_bases-extra_bases_roi # = 100
    tidxpos0 = -(len(UTR5pCDS)-extra_bases) # len of UTR5 alone
    tregion_s = tidxpos0 # = len(UTR5)
    tregion_e = roi_difference # = 100
    qidxpos0 = -extra_bases # = -200
    qregion_s = -roi_difference # = -100
    qregion_e = len(UTR3pCDS)-extra_bases # len of UTR3 alone
    cmd = ["IntaRNA", "-t", UTR5pCDS, "-q", UTR3pCDS,
           "--tidxpos0", str(tidxpos0), ## Transition UTR5-CDS
           "--qidxpos0", str(qidxpos0),  ## Transition CDS-UTR3
           "--tregion", f"{str(tregion_s)}-{str(tregion_e)}", ## UTR5start to end+100CDS
           "--qregion", f"{str(qregion_s)}-{str(qregion_e)}", ## 100CDS to UTR3end
           "--tId", f"{ID}.5UTR",
           "--qId", f"{ID}.3UTR",
           "--parameterFile", static_param_path,
           "--outMode", "C"
           ]
    cmd += additional_args
    #print(" ".join(cmd))
    p = subprocess.Popen(cmd, stdout=subprocess.PIPE)
    return p


def read_output(p, f):
    """Reads output of the IntaRNA process,
    writes the output in a currently open file f
    and processes the data to return the relevant parts.
    
    p (subprocess): Process created by execute_IntaRNA
    f (open file): File for the raw output
    """
    t_starts = []
    t_ends = []
    q_starts = []
    q_ends = []
    energy = []
    hybDPs = []
    p_output = list(p.stdout)
    f.write(p_output[0].decode("utf-8")) # header
    for out in p_output[1:]:
        opt = out.decode("utf-8")
        f.write(opt)
        tid, ts, te, qid, qs, qe, seqDP, hybDP, e = opt.split(";")
        t_starts.append(int(ts))
        t_ends.append(int(te))
        q_starts.append(int(qs))
        q_ends.append(int(qe))
        energy.append(float(e.rstrip("\n")))
        hybDPs.append(hybDP)
    return t_starts, t_ends, q_starts, q_ends, energy, hybDPs


def main_intarna(static_param_path, extra_bases, extra_bases_roi,
                 parameter_table_file, output_path, raw_intarna_output_path,
                 outNumber):
    """
    Main IntaRNA process.
    Iterates over the parameter file, runs IntaRNA for each row and 
    returns the outputs as a single dataframe.
    
    static_param_path (str): Filepath for the file that will be fed to IntaRNA directly
    extra_bases (int): #Extra bases of the CDS to include in IntaRNA input
    extra_bases_roi (int): #Extra bases that of the CDS that IntaRNA is allowed to find interactions in
    parameter_table_file (str): Filepath to parameter file for most variable IntaRNA parameters
    output_path (str): Where to save the output file
    raw_intarna_output_path (str): Where to save the raw IntaRNA output
    outNumber (int): Sets how many subopts are allowed (N-1)
    
    Output:
    output (df): Resulting dataframe with the IntaRNA results
    """
    params = pd.read_csv(parameter_table_file)
    output = pd.DataFrame()
    t_ranges = []  ## All constrained interactions of all sequences
    q_ranges = []
    with open(raw_intarna_output_path, "w") as f:
        for index, row in params.iterrows():
            f.write(f"{row['id']} :\n")
            f.write(f"{'#'*(len(row['id']) + 2)}\n")


            p = execute_IntaRNA(row["seq5"], row["seq3"], row['id'],
                                extra_bases, extra_bases_roi, static_param_path,
                                 ["--outNumber", str(outNumber), "--outOverlap", "T"])
            t_starts, t_ends, q_starts, q_ends, energy, hybDPs = read_output(p, f)
            p.wait()
            if outNumber > 1:
                p = execute_IntaRNA(row["seq5"], row["seq3"], row['id'],
                                    extra_bases, extra_bases_roi, static_param_path,
                                    ["--outNumber", str(outNumber), "--outOverlap", "Q"])
                t_starts_2, t_ends_2, q_starts_2, q_ends_2, energy_2, hybDPs_2 = read_output(p, f)
                t_starts += t_starts_2[1:]
                t_ends += t_ends_2[1:]
                q_starts += q_starts_2[1:]
                q_ends += q_ends_2[1:]
                energy += energy_2[1:]
                hybDPs += hybDPs_2[1:]
            f.write(f"{'#'*(len(row['id']) + 2)}\n")
            
            inter_ts = []
            inter_qs = []
            for i in range(len(t_starts)):
                inter_ts.append((t_starts[i], t_ends[i], energy[i], hybDPs[i]))
                inter_qs.append((q_starts[i], q_ends[i], energy[i], hybDPs[i]))
            t_ranges.append(inter_ts)
            q_ranges.append(inter_qs)
    output = params
    output = params
    output["predictions_t"] = t_ranges
    output["predictions_q"] = q_ranges
    output.to_csv(output_path, index=False)
    return output


def main_mrri(parameter_table_file, static_param_path, extra_bases, extra_bases_roi, 
              mrri_file_output, raw_mrri_output, param_mode):
    """
    param_mode (int): Decides region of interest for MRRI
                      1 - Whole sequence from 5' to 100 into CDS and last 100 to 3' end
                      2 - Limited to small area around 5'UTR/CDS transition and short area in 3' UTR
    """
    params = pd.read_csv(parameter_table_file)
    output = pd.DataFrame()
    t_ranges = []  ## All constrained interactions of all sequences
    q_ranges = []
    with open(raw_mrri_output, "w") as f_raw:
        for index, row in params.iterrows():
            #if row["id"] != "NC_003690.1":
            #    continue
            print(f"MRRI: {row['id']}")
            f_raw.write(f"{row['id']} :\n")
            f_raw.write(f"{'#'*(len(row['id']) + 2)}\n")
            interactions = hacked_MRRI_main(row["seq5"], row["seq3"], static_param_path, param_mode)
            f_raw.write(f"{interactions}\n")

            inter_ts = []
            inter_qs = []
            for i in range(0, len(interactions)):
                interaction = interactions[i]
                #print(interaction)
                hybDP_0, hybDP_1 = interaction["hybridDP"].split("&")
                hybDP_0 = hybDP_0.replace("(", string.ascii_uppercase[i])
                hybDP_1 = hybDP_1.replace(")", string.ascii_lowercase[i])
                s1 = int(interaction["start1"])
                e1 = int(interaction["end1"])
                if s1 >= 0:
                    s1 = str(s1-1)
                if e1 >= 0:
                    e1 = str(e1-1)
                inter_ts.append((s1, e1, interaction['E'], hybDP_0))
                inter_qs.append((interaction["start2"], interaction["end2"], interaction['E'], hybDP_1))
            #raise
            t_ranges.append(inter_ts)
            q_ranges.append(inter_qs)
    #raise
    output = params
    output["predictions_t"] = t_ranges
    output["predictions_q"] = q_ranges
    output.to_csv(mrri_file_output, index=False)