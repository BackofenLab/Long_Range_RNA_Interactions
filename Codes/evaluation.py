import seaborn as sns
import matplotlib.pyplot as plt
from matplotlib.collections import LineCollection
from matplotlib.lines import Line2D
import ast
import pandas as pd
import numpy as np
import math
from collections import defaultdict


main_colours = {"UTR": "gray", "interaction" : "red", "subopt": "orange"}
CDS_colours = defaultdict(lambda:"black")
for k, v in [("ISFV","m"), ("cISFVG","m"), ("dISFVG","plum"),("MBFV","blue"),("NKV","green"),("TBFV","c")]:
    CDS_colours[k] = v
linewidths = {"UTR": 4,
              "CDS": 4,
              "CMHit": 8,
              "Subopt": 4,
              "Interaction": 4
              }

def draw_lineplots(df, extra_bases_roi, output):
    """
    Draw an interaction lineplot 
    showcasing interaction position relative to UTR/CDS.
    
    df (df): the dataframe outputted by main_intarna
    extra_bases_roi (int): Amount of extra bases from the CDS used on each side
    output (str): Filepath where the plot should be saved
    """
    ## Need 3 plots:
    ## 1: UTR5+some_CDS+t_inter (aligned to UTR5/CDS transition)
    ## 2: some_CDS+UTR3+q_inter (aligned to end of UTR3)
    plt.style.use("seaborn-darkgrid")
    fig, (side5, side3) = plt.subplots(2, figsize=(16, 20))

    for index, row in df.iterrows():
        t_tuple = ast.literal_eval(row["t_inter_range"])
        q_tuple = ast.literal_eval(row["q_inter_range"])
        line5sub = []
        line3sub = []
        ## Plot CDS:
        ## Differentiate between cISFVG and dISFVG
        cds_colour = CDS_colours[row["type"]] if row["class"] == "ISFV" else  CDS_colours[row["class"]]
        side5.plot((0, extra_bases_roi), (index, index), 
                   linewidth=linewidths["CDS"], color=cds_colour)
        side3.plot((0 - row["UTR3len"], -extra_bases_roi - row["UTR3len"]), (index, index), 
                   linewidth=linewidths["CDS"], color=cds_colour)
        ## Plot UTR:
        side5.plot((0, -row["UTR5len"]), (index, index), 
                   linewidth=linewidths["UTR"], color=main_colours["UTR"])
        side3.plot((0 - row["UTR3len"], 0), (index, index), 
                   linewidth=linewidths["UTR"], color=main_colours["UTR"])
        ## Plot CM hits:
        if "cm_hit_f" in row and not math.isnan(row["cm_hit_f"]):
            side3.plot((row["cm_hit_f"] - row["UTR3len"], row["cm_hit_t"] - row["UTR3len"]), (index, index), 
                       linewidth=linewidths["CMHit"], color=CDS_colours[row["cm_hit_src"]], alpha=0.5)
        ## Plot Subopt Interactions:
        for suboptt in ast.literal_eval(row["suboptts"]): ## 5' Subopt stuff..
            if suboptt:
                for subopt in suboptt:
                    side5.plot((subopt[0], subopt[1]), (index, index), 
                               linewidth=linewidths["Subopt"], color=main_colours["subopt"])
        for suboptq in ast.literal_eval(row["suboptqs"]): ## 3' Subopt stuff..
            if suboptq:
                for subopt in suboptq:
                    side3.plot((subopt[0] - row["UTR3len"], subopt[1] - row["UTR3len"]), (index, index), 
                               linewidth=linewidths["Subopt"], color=main_colours["subopt"])
        ## Plot Main Interactions:
        side5.plot((t_tuple[0], t_tuple[1]), (index, index), 
                   linewidth=linewidths["Interaction"], color=main_colours["interaction"])
        side3.plot((q_tuple[0] - row["UTR3len"], q_tuple[1] - row["UTR3len"]), (index, index), 
                   linewidth=linewidths["Interaction"], color=main_colours["interaction"])
    ## Annotations:
    side5.set_xlabel("Distance from 5'UTR-CDS transition", fontsize=16)
    side3.set_xlabel("Distance from 3'UTR end", fontsize=16)
    side5.set_ylabel("Index", fontsize=16)
    side3.set_ylabel("Index", fontsize=16)
    side5.yaxis.set_label_position("right")
    side5.yaxis.tick_right()
    side3.yaxis.set_label_position("right")
    side3.yaxis.tick_right()

    side5.set_yticks(np.arange(len(df["id"])))
    side5.set_yticklabels(list(df["id"]))
    side3.set_yticks(np.arange(len(df["id"])))
    side3.set_yticklabels(list(df["id"]))
    side5.autoscale_view()
    side3.autoscale_view()
    legend_elements = [Line2D([0], [0], color=CDS_colours["cISFVG"], lw=8, label='cISFV-CDS'),
                       Line2D([0], [0], color=CDS_colours["dISFVG"], lw=8, label='dISFV-CDS'),
                       Line2D([0], [0], color=CDS_colours["MBFV"], lw=8, label='MBFV-CDS'),
                       Line2D([0], [0], color=CDS_colours["NKV"], lw=8, label='NKV-CDS'),
                       Line2D([0], [0], color=CDS_colours["TBFV"], lw=8, label='TBFV-CDS'),
                       Line2D([0], [0], color='gray', lw=8, label='UTR'),
                       Line2D([0], [0], color='r', lw=8, label='Interaction'),
                       Line2D([0], [0], color='orange', lw=8, label='Subopt'),
                       Line2D([0], [0], color="blue", lw=8, label='CM Hit', alpha=0.5),]
    side5.legend(handles=legend_elements[:-1],loc="upper left", prop={"size": 16})
    side3.legend(handles=legend_elements,loc="upper left", prop={"size": 16})
    plt.suptitle("Interaction Lineplot", fontsize=48)
    plt.savefig(output)
    plt.close()


def draw_energy_histo(df, output):
    """Plots a simple histogram of the of the MFEs of a dataframe.
    
    df (df): the dataframe outputted by main_intarna
    output (str): Filepath where the plot should be saved
    """
    energies = df["energy"]
    plt.figure(figsize=(12.8, 9.6))
    sns.set_style("darkgrid")
    sns.histplot(df, x="energy", hue="class", multiple="stack", bins=int(len(energies)/2))
    plt.xticks(range(int(min(energies)),int(max(energies))))
    plt.xlabel("Interaction energy kcal/mol")
    plt.ylabel("Amount of interactions")
    plt.suptitle("Interaction Energy Histogram")
    plt.savefig(output)
    plt.close()


def draw_energy_histo_subopt(df, output):
    """Plots a simple energy histogram for the subopt MFEs of a dataframe.
    Similarly to draw_energy_histo()
    
    df (df): the dataframe outputted by main_intarna
    output (str): Filepath where the plot should be saved
    """
    subopt_df = {}
    classlist = []
    energylist = []
    for index, row in df.iterrows():
        for el in ast.literal_eval(row["suboptes"]):
            for e in el: ## Not a fan of this
                classlist.append(row["class"])
                energylist.append(e)
    edf = pd.DataFrame({"class": classlist, "energy": energylist})
    sns.histplot(edf, x="energy", hue="class", multiple="stack")
    plt.savefig(f"{output.rstrip('.png')}_subopt.png")
    plt.close()
