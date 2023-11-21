import seaborn as sns
import matplotlib.pyplot as plt
from matplotlib.collections import LineCollection
from matplotlib.lines import Line2D
import ast
import pandas as pd
import numpy as np
import math


main_colours = {"UTR": "gray", "interaction" : "red", "subopt": "orange"}
CDS_colours = {"ISFV":"c","MBFV":"blue","NKV":"green","TBFV":"m"}
CMS_hit_colours = {"ISFV":"lightcyan","MBFV":"cornflowerblue","NKV":"limegreen","TBFV":"plum"}

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
        for suboptt in ast.literal_eval(row["suboptts"]): ## Subopt stuff..
            if suboptt:
                for subopt in suboptt:
                    line5sub.append(((subopt[0], index),
                              (subopt[1], index)))
        subopt_len5 = len(line5sub)
        for suboptq in ast.literal_eval(row["suboptqs"]): ## Subopt stuff..
            if suboptq:
                for subopt in suboptq:
                    line3sub.append(((subopt[0] - row["UTR3len"], index),
                                     (subopt[1] - row["UTR3len"], index))) 
        subopt_len3 = len(line3sub)
        line5 = [((0, index), (extra_bases_roi, index)), ## CDS
                 ((0, index), (-row["UTR5len"], index))] ## 5'UTR
        line3 = [((0 - row["UTR3len"], index), (-extra_bases_roi - row["UTR3len"], index)), ## CDS
                 ((0 - row["UTR3len"], index), (0, index))] ## 3'UTR
        if math.isnan(row["cm_hit_f"]):
            cmsearch_hit = 0
        else:
            cmsearch_hit = 1
            line3 += [((row["cm_hit_f"] - row["UTR3len"], index),
                       (row["cm_hit_t"] - row["UTR3len"], index))] ## CM Hit
        line5 += line5sub ## Subopt
        line3 += line3sub ## Subopt
        line5 += [((t_tuple[0], index),
                  (t_tuple[1], index))] ## Main Interaction
        line3 += [((q_tuple[0] - row["UTR3len"], index),
                 (q_tuple[1] - row["UTR3len"], index))] ## Main Interaction
        colours5 = [CDS_colours[row["class"]]] + [main_colours["UTR"]] + [main_colours["subopt"]]*subopt_len5 + [main_colours["interaction"]]
        if cmsearch_hit:
            colours3 = [CDS_colours[row["class"]]] + [main_colours["UTR"]] + [CMS_hit_colours[row["cm_hit_src"]]] + [main_colours["subopt"]]*subopt_len3 + [main_colours["interaction"]]
        else:
            colours3 = [CDS_colours[row["class"]]] + [main_colours["UTR"]] + [main_colours["subopt"]]*subopt_len3 + [main_colours["interaction"]]
        linewidths5 = [4] + [4] + [4]*subopt_len5 + [4]
        linewidths3 = [4] + [4] + [8]*cmsearch_hit + [4]*subopt_len3 + [4]
        side5.add_collection(LineCollection(line5, colors=colours5, linewidths=linewidths5))
        side3.add_collection(LineCollection(line3, colors=colours3, linewidths=linewidths3))
    side5.set_xlabel("Distance from 5'UTR-CDS transition", fontsize=16)
    side3.set_xlabel("Distance from 3'UTR end", fontsize=16)
    side5.set_ylabel("Index", fontsize=16)
    side3.set_ylabel("Index", fontsize=16)
    side5.yaxis.set_label_position("right")
    side5.yaxis.tick_right()
    side3.yaxis.set_label_position("right")
    side3.yaxis.tick_right()

    #print(list(df["id"]))
    side5.set_yticks(np.arange(len(df["id"])))
    side5.set_yticklabels(list(df["id"]))
    side3.set_yticks(np.arange(len(df["id"])))
    side3.set_yticklabels(list(df["id"]))
    side5.autoscale_view()
    side3.autoscale_view()
    legend_elements = [Line2D([0], [0], color=CDS_colours["ISFV"], lw=8, label='ISFV-CDS'),
                       Line2D([0], [0], color=CDS_colours["MBFV"], lw=8, label='MBFV-CDS'),
                       Line2D([0], [0], color=CDS_colours["NKV"], lw=8, label='NKV-CDS'),
                       Line2D([0], [0], color=CDS_colours["TBFV"], lw=8, label='TBFV-CDS'),
                       Line2D([0], [0], color='gray', lw=8, label='UTR'),
                       Line2D([0], [0], color='r', lw=8, label='Interaction'),
                       Line2D([0], [0], color='orange', lw=8, label='Subopt'),
                       Line2D([0], [0], color=CMS_hit_colours["ISFV"], lw=8, label='CM Hit ISFV'),
                       Line2D([0], [0], color=CMS_hit_colours["MBFV"], lw=8, label='CM Hit MBFV'),
                       Line2D([0], [0], color=CMS_hit_colours["NKV"], lw=8, label='CM Hit NKV'),
                       Line2D([0], [0], color=CMS_hit_colours["TBFV"], lw=8, label='CM Hit TBFV'),]
    side5.legend(handles=legend_elements[:-4],loc="upper left", prop={"size": 16})
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
