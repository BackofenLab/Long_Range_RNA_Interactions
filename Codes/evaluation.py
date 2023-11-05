import seaborn as sns
import matplotlib.pyplot as plt
from matplotlib.collections import LineCollection
from matplotlib.lines import Line2D
import ast


def draw_lineplots(df, extra_bases_roi, output):
    """
    Draw an interaction lineplot 
    showcasing interaction position relative to UTR/CDS.
    df (df): the dataframe outputted by main_intarna
    extra_bases_roi (int): Amount of extra bases from the CDS used on each side
    output (str): Filepath where the plot should be saved
    """
    ## Need 2 plots:
    ## 1: UTR5+some_CDS+t_inter
    ## 2: some_CDS+UTR3+q_inter
    ## inter can be either on UTR or CDS (or maybe overlapping both)
    plt.style.use("seaborn-darkgrid")
    fig, (side5, side3) = plt.subplots(2, figsize=(16, 9))
    for index, row in df.iterrows():
        #if row["id"] != "NC_001437.1":
        #    continue
        t_tuple = ast.literal_eval(row["t_inter_range"])
        q_tuple = ast.literal_eval(row["q_inter_range"])
        line5 = [((0, index), (extra_bases_roi, index)), ## grey
                ((0, index), (-row["UTR5len"], index)), ## blue
                ((t_tuple[0], index),
                 (t_tuple[1], index)) ## red
                 ]
        line3 = [((0, index), (-extra_bases_roi, index)),#((0, index), (-ext_end3, index)), ## grey
                 ((0, index), ((row["UTR3len"]), index)), ## blue
                 ((q_tuple[0], index),
                 (q_tuple[1], index)) ## red
                 ]
        CDS_colors = {"ISFV":"y","MBFV":"b","NKV":"g","TBFV":"m"}
        colors = [CDS_colors[row["class"]], "c", "r"]
        side5.add_collection(LineCollection(line5, colors=colors, linewidths=(2,)))
        side3.add_collection(LineCollection(line3, colors=colors, linewidths=(2,)))
    side5.set_xlabel("Distance from 5'UTR-CDS transition")
    side3.set_xlabel("Distance from CDS-3'UTR transition")
    side5.set_ylabel("Index")
    side3.set_ylabel("Index")
    side5.autoscale_view()
    side3.autoscale_view()
    legend_elements = [Line2D([0], [0], color=CDS_colors["ISFV"], lw=4, label='ISFV-CDS'),
                       Line2D([0], [0], color=CDS_colors["MBFV"], lw=4, label='MBFV-CDS'),
                       Line2D([0], [0], color=CDS_colors["NKV"], lw=4, label='NKV-CDS'),
                       Line2D([0], [0], color=CDS_colors["TBFV"], lw=4, label='TBFV-CDS'),
                       Line2D([0], [0], color='c', lw=4, label='UTR'),
                       Line2D([0], [0], color='r', lw=4, label='Interaction')]
    side5.legend(handles=legend_elements,loc="upper left")
    side3.legend(handles=legend_elements,loc="upper right")
    plt.suptitle("Interaction Lineplot")
    plt.savefig(output)


def draw_energy_histo(df, output):
    """Plots a simple histogram of the of the MFEs of a dataframe.
    df (df): the dataframe outputted by main_intarna
    output (str): Filepath where the plot should be saved
    """
    energies = df["energy"]
    plt.figure(figsize=(12.8, 9.6))
    #plt.hist(energies, bins=int(len(energies)/2))
    sns.set_style("darkgrid")
    sns.histplot(df, x="energy", hue="class", multiple="stack", bins=int(len(energies)/2))
    plt.xticks(range(int(min(energies)),int(max(energies))))
    plt.xlabel("Interaction energy kcal/mol")
    plt.ylabel("Amount of interactions")
    plt.suptitle("Interaction Energy Histogram")
    plt.savefig(output)
