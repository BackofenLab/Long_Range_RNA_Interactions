import pandas as pd
import os
import math
from Bio import SeqIO
import subprocess

database_path = "Flavivirus_NCBI/Flavivirus_RefSeq_20220621"

qidxpos0 = 0
energyVRNA = "rna_andronescu2007.par"
intLenMax = 20
qintLenMax = tintLenMax = intLenMax
extra_bases = 200
extra_bases_roi = 100 ## Region of Interest

def read_out_bed6_data(file):
    f = open(file, "r")
    UTR5end, CDSend, UTR3end = None, None, None
    lines = f.readlines()
    if len(lines) is 3:
        UTR5end = lines[0].split("\t")[2]
        CDSend = lines[1].split("\t")[2]
        UTR3end = lines[2].split("\t")[2]
    return UTR5end, CDSend, UTR3end

def create_parameter_table():
    d = {}
    for root, dirs, files in os.walk(database_path):
        #print(files)
        for file in files:
            if file.endswith(".bed6"):
                _, _, family, virus, gene_name = root.split("/")
                split_name = file.split(".")[:-1]
                gene = ".".join(split_name)
                UTR5end, CDSend, UTR3end = read_out_bed6_data(f"{root}/{file}")
                if UTR5end is None:
                    continue
                d[gene] = [family, virus, gene_name, UTR5end, CDSend, UTR3end,
                           energyVRNA, qintLenMax, tintLenMax, qidxpos0]

    columns = ["family", "virus", "gene_name", "5UTRend", "CDSend", "3UTRend",
               "energyVRNA", "qintLenMax", "tintLenMax", "qidxpos0"]
    df = pd.DataFrame(d)
    df.index = pd.Index(columns, name="gene")
    df = df.T
    df.to_csv("parameter_table.csv", index_label="name")

def execute_IntaRNA(UTR5, UTR3, row):
    cmd = ["IntaRNA", "-t", UTR5, "-q", UTR3,
           "--energyVRNA", row["energyVRNA"],
           "--tintLenMax", str(row["tintLenMax"]),
           "--qintLenMax", str(row["qintLenMax"]),
           "--tregion", f"0-{str(len(UTR5)-extra_bases+extra_bases_roi)}",
           "--qregion", f"0-{str(len(UTR3)-extra_bases+extra_bases_roi)}"
           ]
    p = subprocess.Popen(cmd, stdout=subprocess.PIPE)
    return p

def main():
    df = pd.read_csv("parameter_table.csv")
    f = open("test.txt", "w")
    for index, row in df.iterrows():
        if not math.isnan(row["3UTRend"]):
            if row['gene_name'] != "JEV":
                continue
            file = f"{database_path}/{row['family']}/{row['virus']}/{row['gene_name']}/{row['name']}.fa"
            DNA = SeqIO.read(file, "fasta")
            UTR5 = str(DNA.seq[:int(row["5UTRend"])+extra_bases])
            UTR3 = str(DNA.seq[-int(row["3UTRend"]-row["CDSend"])-extra_bases:])
##            print(UTR5)
##            print("----")
##            print(UTR3)
            f.write(f"{row['name']} :\n")
            f.write(f"{'#'*(len(row['name']) + 2)}\n")
            p = execute_IntaRNA(UTR5, UTR3, row)
            for line in p.stdout:
                f.write(line.decode("utf-8"))
            p.wait()
            f.write(f"{'#'*(len(row['name']) + 2)}\n")

#create_parameter_table()
main()
