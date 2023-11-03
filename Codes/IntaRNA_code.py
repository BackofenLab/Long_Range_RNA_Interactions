import pandas as pd
import math
import subprocess
import re


qidxpos0 = 0



def execute_IntaRNA(UTR5pCDS, UTR3pCDS, static_params):
    """Starts a subprocess for IntaRNA with given parameters.
    UTR5pCDS (str): 5'UTR+ the extra bases from the CDS
    UTR3pCDS (str): 3'UTR+ the extra bases from the CDS
    static_params (df): Dataframe of predefined static parameters
    """
    # extra_region cuts off the difference between
    # the extra CDS bases and the extra CDS bases as region of interest
    extra_region = -static_params["extra_bases"][0]+static_params["extra_bases_roi"][0]
    cmd = ["IntaRNA", "-t", UTR5pCDS, "-q", UTR3pCDS,
           "--energyVRNA", str(static_params["energyVRNA"][0]),
           "--intLenMax", str(static_params["intLenMax"][0]),
           "--tregion", f"1-{str(len(UTR5pCDS)+extra_region)}",
           "--qregion", f"1-{str(len(UTR3pCDS)+extra_region)}",
           "--seedBP", str(static_params["seedlength"][0]),
           #"--outMode", "C" # This seems incredibly useful but I opt to not
           # use it for now to preserve the readable textfile!
           # (and because I already went through the effort to write code
           # to interpret the normal output..
           
           ]
    p = subprocess.Popen(cmd, stdout=subprocess.PIPE)
    return p


def read_output(p, f):
    """Reads output of the IntaRNA process,
    writes the output in a currently open file f
    and processes the data to return the relevant parts.
    """
    linecounter = 0
    t_inter_range = (0,0)
    q_inter_range = (0,0)
    energy = 0###
    for line in p.stdout:
        dline = line.decode("utf-8")
        f.write(dline)
        if linecounter == 2:
            l = re.findall(r'\d+', dline)
            t_inter_range = (int(l[0]), int(l[1]))
        elif linecounter == 10:
            l = re.findall(r'\d+', dline)
            q_inter_range = (int(l[0]), int(l[1]))
        elif linecounter == 13:
            energy = float(re.findall(r'[-]\d+[.]\d+', dline)[0])
        linecounter += 1
    return t_inter_range, q_inter_range, energy




def main_intarna(database_path, static_params, params,
                 output_path, raw_intarna_output_path):
    static_params = pd.read_csv(static_params)
    params = pd.read_csv(params)
    output = pd.DataFrame()
    f = open(raw_intarna_output_path, "w")
    t_ranges = []
    q_ranges = []
    energies = []
    for index, row in params.iterrows():
            f.write(f"{row['id']} :\n")
            f.write(f"{'#'*(len(row['id']) + 2)}\n")
            p = execute_IntaRNA(row["seq5"], row["seq3"], static_params)
            t_inter_range, q_inter_range, energy = read_output(p, f)
            p.wait()
            f.write(f"{'#'*(len(row['id']) + 2)}\n")
            t_ranges.append(t_inter_range) ## This is fairly suboptimal right now.
            q_ranges.append(q_inter_range) ## Might need to rework if used on large dbs
            energies.append(energy)        ##
    output["5UTRend"] = params["5UTRend"]
    output["CDSend"] = params["CDSend"]
    output["3UTRend"] = params["3UTRend"]
    output["t_inter_range"] = t_ranges
    output["q_inter_range"] = q_ranges
    output["energy"] = energies
    output.to_csv(output_path, index=False)
    return output


if __name__ == "__main__":
    pass
