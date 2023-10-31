import pandas as pd
import os
import math
from Bio import SeqIO
import subprocess
#import seaborn as sns
import matplotlib.pyplot as plt
#import matplotlib as mpl##??
from matplotlib.collections import LineCollection
import re

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


def read_output(p, f):
    linecounter = 0
    t_inter_range = (0,0)
    q_inter_range = (0,0)
    for line in p.stdout:
        dline = line.decode("utf-8")
        f.write(dline)
        if linecounter == 2:
##            split_line = dline.strip(" ").split(" ")
            l = re.findall(r'\d+', dline)
            t_inter_range = (int(l[0]), int(l[1]))
##            t_inter_range = (int(split_line[0]), int(split_line[-1].rstrip("\n")))
        elif linecounter == 10:
##            split_line = dline.strip(" ").split(" ")
            l = re.findall(r'\d+', dline)
            q_inter_range = (int(l[0]), int(l[1]))
##            q_inter_range = (int(split_line[0]), int(split_line[-1].rstrip("\n")))
        elif linecounter == 13:
            energy = float(re.findall(r'[-]\d+[.]\d+', dline)[0])
        linecounter += 1
    return t_inter_range, q_inter_range, energy


def draw_lineplot(sequence_range, interaction_range, UTRrange, nr, ax):
    lines = [((sequence_range[0], nr), (sequence_range[1], nr)),
             ((UTRrange[0], nr), (UTRrange[1], nr)),
             ((interaction_range[0], nr), (interaction_range[1], nr)),
             ]
    colors = ["grey", "c", "r"]
    colored_lines = LineCollection(lines, colors=colors, linewidths=(2,))
    ax.add_collection(colored_lines)


def draw_energy_histo(energies):
    plt.hist(energies, bins=int(len(energies)/2))
    plt.xticks(range(int(min(energies)),int(max(energies))))
    plt.xlabel("Interaction energy kcal/mol")
    plt.ylabel("Amount of interactions")
    plt.suptitle("Interaction Energy Histogram")
    plt.savefig("energy_histo.png")


def main():
    df = pd.read_csv("parameter_table.csv")
    f = open("test.txt", "w")
    t_ranges = []
    q_ranges = []
    energies = []
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
            t_inter_range, q_inter_range, energy = read_output(p, f)
            p.wait()
            f.write(f"{'#'*(len(row['id']) + 2)}\n")
            t_ranges.append(t_inter_range) ## This is fairly suboptimal right now.
            q_ranges.append(q_inter_range) ## Might need to rework if used on large dbs
            energies.append(energy)        ##
    draw_energy_histo(energies)
    df["t_inter_range"] = t_ranges
    df["q_inter_range"] = q_ranges
    df["energy"] = energies
    df.to_csv("parameter_table.csv", index=False)

#create_parameter_table()
main()
##fig, ax = plt.subplots(1)
##draw_lineplot((0,500),(51,60), (0, 300), 0, ax)
##draw_lineplot((0,500),(51,60), (0, 55), 1, ax)
##draw_lineplot((0,500),(51,60), (300, 500), 2, ax)
##ax.autoscale_view()
##plt.show()
##df = pd.read_csv("parameter_table.csv")
##draw_energy_histo(df["energy"])
