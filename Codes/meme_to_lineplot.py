import seaborn as sns
import matplotlib.pyplot as plt
from matplotlib.collections import LineCollection
from matplotlib.lines import Line2D
import ast
import pandas as pd
import numpy as np
import math
from collections import defaultdict
from Codes.evaluation import draw_lineplots


def meme_to_lineplot(df, extra_bases_roi, meme_output):
    """Improvised function to extract relevant information from meme output text files
    to add them to a lineplot similar to the one from evaluation.py.
    Noteworth variables:
    mega_dict(dict): A dictionary of virus classes where each entry is a
                     dictionary of the 4 sites and their respective motifs
                     locations for each of the sequences as tuples.
                     {"class": {"site": {"id" : (motif_start, motif_end)}}}
    """
    virus_classes = ("MBFV", "TBFV", "ISFV", "NKV")
    mega_dict = defaultdict(dict)
    for site in ["site_1", "site_2", "site_3", "site_4"]:
        d = defaultdict(dict)
        f = open(f"{meme_output}/MEME Outputs/{site}/meme.txt", "r")
        trigger = False
        for l in f.readlines():
            if l.startswith("Sequence name             Start   P-value                    Site"):
                trigger = True
            elif trigger and l.startswith(virus_classes):
                l_split = l.split()
                l_split0_split = l_split[0].split("-")
                motif_start = int(l_split[1])
                motif_end = motif_start + len(l_split[4])
                id = l_split0_split[1]
                # if l_split[2].split("-")[1] > THRESHOLD: ## has form of "8.92e-11"
                d[l_split0_split[0]][id] = (motif_start, motif_end)
            elif trigger and l.startswith("--------------"):
                break
        for virus_class in virus_classes:
            mega_dict[virus_class][site] = d[virus_class]
    for virus_class in virus_classes:
        meme_sites = mega_dict[virus_class]
        draw_lineplots(df.loc[df["class"] == virus_class],
                       extra_bases_roi,
                       f"{meme_output}/MEME_{virus_class}.png", subopt_mode=True, meme_sites=meme_sites) #virus_class ? 
 
    
def glam2_to_lineplot(df, extra_bases_roi, meme_output):
    """
    Improvised function to extract relevant information from glam2 output text files
    to add them to a lineplot similar to the one from evaluation.py.
    Note: This will literally only work with the super specific output file by the MEME Suite webserver.
    Noteworth variables:
    mega_dict(dict): A dictionary of virus classes where each entry is a
                     dictionary of the 4 sites and their respective motifs
                     locations for each of the sequences as tuples.
                     {"class": {"site": {"id" : (motif_start, motif_end)}}}
    """
    from Bio import SeqIO ### ...
    virus_classes = ("MBFV", "TBFV", "ISFV", "NKV")
    mega_dict = defaultdict(dict)
    for site in ["site_1", "site_2", "site_3", "site_4"]:
        seq_dict = {str(rec.seq) : rec.id for rec in SeqIO.parse(f"{meme_output}/{site}.fa", "fasta")} ### ...
        d = defaultdict(dict)
        f = open(f"{meme_output}/GLAM2 Outputs/{site}/glam2.txt", "r")
        trigger = False
        for l in f.readlines():
            if l.startswith("                ************************"):
                trigger = True
            elif trigger and l.startswith(virus_classes):
                l_split = l.split()
                motif_start = int(l_split[1])
                sequence = l_split[2].replace(".", "").upper()
                motif_end = int(l_split[3])
                #print(seq_dict)
                for seq in seq_dict.keys(): ### AAAAAAA
                    if sequence in seq:
                        virus_class, id = seq_dict[seq].split("-")
                d[virus_class][id] = (motif_start, motif_end)
            elif trigger and l.startswith(" "):
                break
        for virus_class in virus_classes:
            mega_dict[virus_class][site] = d[virus_class]
    for virus_class in virus_classes:
        meme_sites = mega_dict[virus_class]
        draw_lineplots(df.loc[df["class"] == virus_class],
                       extra_bases_roi,
                       f"{meme_output}/GLAM2_{virus_class}.png", subopt_mode=True, meme_sites=meme_sites)