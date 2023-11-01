import pandas as pd
import os
import math
from Bio import SeqIO
import subprocess
#import seaborn as sns
import matplotlib.pyplot as plt
from matplotlib.collections import LineCollection
import re
import ast

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
    if len(lines) == 3: ## Filter files without UTR data
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
                d[vid] = [vclass, vtype, virus, UTR5end, CDSend, UTR3end, extra_bases,
                           energyVRNA, qintLenMax, tintLenMax, qidxpos0]

    columns = ["class", "type", "virus", "5UTRend", "CDSend", "3UTRend",
               "extra_bases", "energyVRNA", "qintLenMax", "tintLenMax", "qidxpos0"]
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
            l = re.findall(r'\d+', dline)
            t_inter_range = (int(l[0]), int(l[1]))
        elif linecounter == 10:
            l = re.findall(r'\d+', dline)
            q_inter_range = (int(l[0]), int(l[1]))
        elif linecounter == 13:
            energy = float(re.findall(r'[-]\d+[.]\d+', dline)[0])
        linecounter += 1
    return t_inter_range, q_inter_range, energy


def draw_lineplots(df):
    ## Need 2 plots:
    ## 1: UTR5+some_CDS+t_inter
    ## 2: some_CDS+UTR3+q_inter
    ## inter can be either on UTR or CDS (or maybe overlapping both)
    row_dict = {}###
    #plt.figure(figsize=(12.8, 9.6))
    fig, (side5, side3) = plt.subplots(2, figsize=(16, 9))
    for index, row in df.iterrows():
        ext_end5 = row["5UTRend"]+row["extra_bases"]
        ext_end3 = (row["3UTRend"]-row["CDSend"])+row["extra_bases"]
        t_tuple = ast.literal_eval(row["t_inter_range"])
        q_tuple = ast.literal_eval(row["q_inter_range"])
        line5 = [((0, index), (ext_end5, index)), ## grey
                ((0, index), (row["5UTRend"], index)), ## blue
                ((t_tuple[0], index),
                 (t_tuple[1], index)) ## red
                 ]
        line3 = [((0, index), (-ext_end3, index)), ## grey
                 ((0, index), (-(row["3UTRend"]-row["CDSend"]), index)), ## blue
                 #((-q_tuple[0], index),
                 #(-q_tuple[1], index)) ## red
                 ((-ext_end3+q_tuple[0], index), ## This works (might have small index issue +-1?)
                                                 ## but is incredibly ugly and hacky.
                                                 ## Definitely needs to be reworked after the indexing is cleared up
                 (-ext_end3+q_tuple[1], index)) ## red
                 ]
        colors = ["grey", "c", "r"]
        side5.add_collection(LineCollection(line5, colors=colors, linewidths=(2,)))
        side3.add_collection(LineCollection(line3, colors=colors, linewidths=(2,)))
    side5.set_xlabel("Bases after 5' end")
    side3.set_xlabel("Bases before 3' end")
    side5.set_ylabel("Index")
    side3.set_ylabel("Index")
    side5.autoscale_view()
    side3.autoscale_view()
    plt.suptitle("Interaction Lineplot")
    plt.savefig("interaction_lineplot.png")


def draw_energy_histo(df):
    energies = df["energy"]
    plt.figure(figsize=(12.8, 9.6))
    plt.hist(energies, bins=int(len(energies)/2))
    plt.xticks(range(int(min(energies)),int(max(energies))))
    plt.xlabel("Interaction energy kcal/mol")
    plt.ylabel("Amount of interactions")
    plt.suptitle("Interaction Energy Histogram")
    plt.savefig("energy_histo.png")


def main():
    df = pd.read_csv("parameter_table.csv")
    f = open("IntaRNA_Output.txt", "w")
    t_ranges = []
    q_ranges = []
    energies = []
    for index, row in df.iterrows():
        if not math.isnan(row["3UTRend"]): ## Should be unneccessary now but just in case
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
    df["t_inter_range"] = t_ranges
    df["q_inter_range"] = q_ranges
    df["energy"] = energies
    df.to_csv("parameter_table.csv", index=False)
    return df


if __name__ == "__main__":
    create_parameter_table()
    df = main()
    df = pd.read_csv("parameter_table.csv")
    draw_lineplots(df)
    draw_energy_histo(df)
