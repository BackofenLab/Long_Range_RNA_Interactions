#import seaborn as sns
import matplotlib.pyplot as plt
from matplotlib.collections import LineCollection
from matplotlib.lines import Line2D
import ast


def draw_lineplots(df, extra_bases, output):
    ## Need 2 plots:
    ## 1: UTR5+some_CDS+t_inter
    ## 2: some_CDS+UTR3+q_inter
    ## inter can be either on UTR or CDS (or maybe overlapping both)
    row_dict = {}###
    #plt.figure(figsize=(12.8, 9.6))
    fig, (side5, side3) = plt.subplots(2, figsize=(16, 9))
    for index, row in df.iterrows():
        ext_end5 = row["5UTRend"]+extra_bases
        ext_end3 = (row["3UTRend"]-row["CDSend"])+extra_bases
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
    legend_elements = [Line2D([0], [0], color='grey', lw=4, label='CDS'),
                       Line2D([0], [0], color='c', lw=4, label='UTR'),
                       Line2D([0], [0], color='r', lw=4, label='Interaction')]
    side5.legend(handles=legend_elements,loc="upper right")
    side3.legend(handles=legend_elements,loc="upper left")
    plt.suptitle("Interaction Lineplot")
    plt.savefig(output)


def draw_energy_histo(df, output):
    """Plots a simple histogram of the of the MFEs of a dataframe.
    """
    energies = df["energy"]
    plt.figure(figsize=(12.8, 9.6))
    plt.hist(energies, bins=int(len(energies)/2))
    plt.xticks(range(int(min(energies)),int(max(energies))))
    plt.xlabel("Interaction energy kcal/mol")
    plt.ylabel("Amount of interactions")
    plt.suptitle("Interaction Energy Histogram")
    plt.savefig(output)
