import pandas as pd
import os
import math
from Bio import SeqIO
import subprocess
#import seaborn as sns
import matplotlib.pyplot as plt
#import matplotlib as mpl##??
from matplotlib.collections import LineCollection

database_path = "Flavivirus_NCBI/Flavivirus_RefSeq_20220621"

qidxpos0 = 0
energyVRNA = "rna_andronescu2007.par"
intLenMax = 20
qintLenMax = tintLenMax = intLenMax
extra_bases = 200
extra_bases_roi = 100 ## Region of Interest
seedlength = 7


def read_out_bed6_data(file):
    f = open(file, "r")
    UTR5end, CDSend, UTR3end = None, None, None
    lines = f.readlines()
    if len(lines) == 3:
        UTR5end = lines[0].split("\t")[2]
        CDSend = lines[1].split("\t")[2]
        UTR3end = lines[2].split("\t")[2]
    return UTR5end, CDSend, UTR3end


def create_parameter_table():
    d = {}
    for root, dirs, files in os.walk(database_path):
        for file in files:
            if file.endswith(".bed6"):
                _, _, vclass, vtype, virus = root.split("/")
                split_name = file.split(".")[:-1]
                vid = ".".join(split_name)
                UTR5end, CDSend, UTR3end = read_out_bed6_data(f"{root}/{file}")
                if UTR5end is None:
                    continue
                d[vid] = [vclass, vtype, virus, UTR5end, CDSend, UTR3end,
                           energyVRNA, qintLenMax, tintLenMax, qidxpos0]

    columns = ["class", "type", "virus", "5UTRend", "CDSend", "3UTRend",
               "energyVRNA", "qintLenMax", "tintLenMax", "qidxpos0"]
    df = pd.DataFrame(d)
    df.index = pd.Index(columns, name="id")
    df = df.T
    df.to_csv("parameter_table.csv", index_label="id")


def execute_IntaRNA(UTR5, UTR3, row):
    cmd = ["IntaRNA", "-t", UTR5, "-q", UTR3,
           "--energyVRNA", row["energyVRNA"],
           "--tintLenMax", str(row["tintLenMax"]),
           "--qintLenMax", str(row["qintLenMax"]),
           "--tregion", f"1-{str(len(UTR5)-extra_bases+extra_bases_roi)}",
           "--qregion", f"1-{str(len(UTR3)-extra_bases+extra_bases_roi)}",
           "--seedBP", str(seedlength)
           ]
    p = subprocess.Popen(cmd, stdout=subprocess.PIPE)
    return p


def get_ranges(p, f):
    linecounter = 0
    t_inter_range = (0,0)
    q_inter_range = (0,0)
    for line in p.stdout:
        dline = line.decode("utf-8")
        f.write(dline)
        if linecounter == 2:
            split_line = dline.strip(" ").split(" ")
            t_inter_range = (int(split_line[0]), int(split_line[-1].rstrip("\n")))
        elif linecounter == 10:
            split_line = dline.strip(" ").split(" ")
            q_inter_range = (int(split_line[0]), int(split_line[-1].rstrip("\n")))
        linecounter += 1
    return t_inter_range, q_inter_range


def drawplot(trange, qrange, UTR5end, CDSend, UTR3end, extra_bases):
    fig, ax = plt.subplots(1)
    tlines = [((0,0), (trange[0],0)),((trange[0],0),(trange[1],0)), ((trange[1],0),(UTR5end+extra_bases,0))]
    qlines = [((CDSend-extra_bases,0), (qrange[1],0)),((qrange[1],0),(qrange[0],0)), ((qrange[0],0),(UTR3end,0))]
    clines = [((0,1), (UTR5end,1)),((UTR5end,1),(CDSend,1)), ((CDSend,1),(UTR3end,1))]
    colors1 = ["c", "r", "c"]
    colors2 = ["c", "grey", "c"]
    colored_lines = LineCollection(tlines, colors=colors1, linewidths=(2,))
    ax.add_collection(colored_lines)
    colored_lines2 = LineCollection(qlines, colors=colors1, linewidths=(2,))
    ax.add_collection(colored_lines2)
    colored_lines3 = LineCollection(clines, colors=colors2, linewidths=(2,))
    ax.add_collection(colored_lines3)
    ax.autoscale_view()
    plt.show()


def main():
    df = pd.read_csv("parameter_table.csv")
    f = open("test.txt", "w")
    for index, row in df.iterrows():
        if not math.isnan(row["3UTRend"]):
##            if row['virus'] != "JEV":
##                continue
            file = f"{database_path}/{row['class']}/{row['type']}/{row['virus']}/{row['id']}.fa"
            DNA = SeqIO.read(file, "fasta")
            UTR5 = str(DNA.seq[:int(row["5UTRend"])+extra_bases])
            UTR3 = str(DNA.seq[-int(row["3UTRend"]-row["CDSend"])-extra_bases:])
##            print(UTR5)
##            print("----")
##            print(UTR3)
            f.write(f"{row['id']} :\n")
            f.write(f"{'#'*(len(row['id']) + 2)}\n")
            p = execute_IntaRNA(UTR5, UTR3, row)
            t_inter_range, q_inter_range= get_ranges(p, f)
            p.wait()
            f.write(f"{'#'*(len(row['id']) + 2)}\n")

##create_parameter_table()
main()
#drawplot((45,51),(866,860), 55, 900, 1000, extra_bases)
