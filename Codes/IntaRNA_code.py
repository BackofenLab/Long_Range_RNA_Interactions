import pandas as pd
import math
import subprocess
import re


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
    # extra_region cuts off the difference between
    # the extra CDS bases and the extra CDS bases as region of interest
    extra_region = extra_bases-extra_bases_roi # = 100
    tidxpos0 = -(len(UTR5pCDS)-extra_bases) # len of UTR5 alone
    tregion_s = tidxpos0 # = len(UTR5)
    tregion_e = extra_region # = 100
    qidxpos0 = -extra_bases # = -200
    qregion_s = -extra_region # = -100
    qregion_e = len(UTR3pCDS)-extra_bases # len of UTR3 alone
    cmd = ["IntaRNA", "-t", UTR5pCDS, "-q", UTR3pCDS,
           "--tidxpos0", str(tidxpos0), ## Transition UTR5-CDS
           "--qidxpos0", str(qidxpos0),  ## Transition CDS-UTR3
           "--tregion", f"{str(tregion_s)}-{str(tregion_e)}", ## UTR5start to end+100CDS
           "--qregion", f"{str(qregion_s)}-{str(qregion_e)}", ## 100CDS to UTR3end
           "--tId", f"{ID}.5UTR",
           "--qId", f"{ID}.3UTR",
           "--parameterFile", static_param_path,
           "--outMode", "C" # This seems incredibly useful but I opt to not
           # use it for now to preserve the readable textfile!
           # (and because I already went through the effort to write code
           # to interpret the normal output..
           ]
    cmd += additional_args
    #print(cmd)
    #raise
    p = subprocess.Popen(cmd, stdout=subprocess.PIPE)
    return p


def read_output(p, f):
    """Reads output of the IntaRNA process,
    writes the output in a currently open file f
    and processes the data to return the relevant parts.
    """
    t_inter_range = []
    q_inter_range = []
    energy = []
    p_output = list(p.stdout)
    f.write(p_output[0].decode("utf-8")) # header
    for out in p_output[1:]:
        opt = out.decode("utf-8")
        f.write(opt)
        tid, ts, te, qid, qs, qe, seqDP, hybDP, e = opt.split(";")
        t_inter_range.append((int(ts), int(te)))
        q_inter_range.append((int(qs), int(qe)))
        energy.append(float(e.rstrip("\n")))
    opt_t = t_inter_range[0]
    opt_q = q_inter_range[0]
    opt_energy = energy[0]
    subopt_t = t_inter_range[1:]
    subopt_q = q_inter_range[1:]
    subopt_e = energy[1:]
    return opt_t, opt_q, opt_energy, subopt_t, subopt_q, subopt_e


def main_intarna(database_path, static_param_path, extra_bases, extra_bases_roi,
                 parameter_table_file, output_path, raw_intarna_output_path,
                 outNumber = 2):
    params = pd.read_csv(parameter_table_file)
    additional_args = []
    output = pd.DataFrame()
    t_opt_ranges = []
    q_opt_ranges = []
    opt_energies = []
    suboptts = []
    suboptqs = []
    suboptes = []
    with open(raw_intarna_output_path, "w") as f:
        for index, row in params.iterrows():
            f.write(f"{row['id']} :\n")
            f.write(f"{'#'*(len(row['id']) + 2)}\n")
            p = execute_IntaRNA(row["seq5"], row["seq3"], row['id'],
                                extra_bases, extra_bases_roi, static_param_path,
                                 ["--outNumber", str(outNumber), "--outOverlap", "T"])
            opt_t, opt_q, opt_energy, subopt1_t, subopt1_q, subopt1_e = read_output(p, f)
            p.wait()
            if outNumber > 1:
                p = execute_IntaRNA(row["seq5"], row["seq3"], row['id'],
                                    extra_bases, extra_bases_roi, static_param_path,
                                    ["--outNumber", str(outNumber), "--outOverlap", "Q"])
                _, _, _, subopt2_t, subopt2_q, subopt2_e = read_output(p, f)
                p.wait()
            else:
                subopt2_t, subopt2_q, subopt2_e = (), (), ()
            f.write(f"{'#'*(len(row['id']) + 2)}\n")
            t_opt_ranges.append(opt_t)      ## This is fairly suboptimal right now.
            q_opt_ranges.append(opt_q)      ## Might need to rework if used on large dbs
            opt_energies.append(opt_energy) ##
            suboptts.append((subopt1_t, subopt2_t))
            suboptqs.append((subopt1_q, subopt2_q))
            suboptes.append((subopt1_e, subopt2_e))
    output["id"] = params["id"]
    output["class"] = params["class"]
    output["UTR5len"] = params["UTR5len"]
    output["CDSlen"] = params["CDSlen"]
    output["UTR3len"] = params["UTR3len"]
    output["t_inter_range"] = t_opt_ranges
    output["q_inter_range"] = q_opt_ranges
    output["energy"] = opt_energies
    output["suboptts"] = suboptts
    output["suboptqs"] = suboptqs
    output["suboptes"] = suboptes
    output.to_csv(output_path, index=False)
    return output