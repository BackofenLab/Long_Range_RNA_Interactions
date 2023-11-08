import seaborn as sns
import matplotlib.pyplot as plt
from matplotlib.collections import LineCollection
from matplotlib.lines import Line2D
import ast
import pandas as pd


def draw_lineplots(df, extra_bases_roi, output):
    """
    Draw an interaction lineplot 
    showcasing interaction position relative to UTR/CDS.
    df (df): the dataframe outputted by main_intarna
    extra_bases_roi (int): Amount of extra bases from the CDS used on each side
    output (str): Filepath where the plot should be saved
    """
    ## Need 3 plots:
    ## 1: UTR5+some_CDS+t_inter
    ## 2: some_CDS+UTR3+q_inter
    ## 3: Like 2 but aligned to the genome end (end of UTR3)
    ## inter can be either on UTR or CDS (or maybe overlapping both)
    plt.style.use("seaborn-darkgrid")
    fig, (side5, side3, side3_re) = plt.subplots(3, figsize=(16, 15))

    for index, row in df.iterrows():
        #if row["id"] != "NC_001437.1":
        #    continue
        #row["suboptqs"]
        t_tuple = ast.literal_eval(row["t_inter_range"])
        q_tuple = ast.literal_eval(row["q_inter_range"])
        line5sub = []
        line3sub = []
        for suboptt in ast.literal_eval(row["suboptts"]): ## Subopt stuff..
            if suboptt:
                for subopt in suboptt:
                    line5sub.append(((subopt[0], index),
                              (subopt[1], index))) ## orange?
        subopt_len5 = len(line5sub)
        for suboptq in ast.literal_eval(row["suboptqs"]): ## Subopt stuff..
            if suboptq:
                for subopt in suboptq:
                    line3sub.append(((subopt[0], index),
                                     (subopt[1], index))) ## orange?
        subopt_len3 = len(line3sub)
        line5 = [((0, index), (extra_bases_roi, index)), ## grey
                 ((0, index), (-row["UTR5len"], index))] ## blue
        line5 += line5sub
        line5 += [((t_tuple[0], index),
                  (t_tuple[1], index))] ## red
        line3 = [((0, index), (-extra_bases_roi, index)),#((0, index), (-ext_end3, index)), ## grey
                 ((0, index), ((row["UTR3len"]), index))] ## blue
        line3 += line3sub
        line3 += [((q_tuple[0], index),
                 (q_tuple[1], index))] ## red
        lin3_realigned = []
        for range_tuple in line3:
            lin3_realigned.append(((range_tuple[0][0] - row["UTR3len"], index), (range_tuple[1][0] - row["UTR3len"], index)))
        CDS_colors = {"ISFV":"c","MBFV":"b","NKV":"g","TBFV":"m"}
        colors5 = [CDS_colors[row["class"]], "gray"] + ["orange"]*subopt_len5 + ["r"]
        colors3 = [CDS_colors[row["class"]], "gray"] + ["orange"]*subopt_len3 + ["r"]
        side5.add_collection(LineCollection(line5, colors=colors5, linewidths=(2,)))
        side3.add_collection(LineCollection(line3, colors=colors3, linewidths=(2,)))
        side3_re.add_collection(LineCollection(lin3_realigned, colors=colors3, linewidths=(2,)))
    side5.set_xlabel("Distance from 5'UTR-CDS transition")
    side3.set_xlabel("Distance from CDS-3'UTR transition")
    side3_re.set_xlabel("Distance from 3'UTR end")
    side5.set_ylabel("Index")
    side3.set_ylabel("Index")
    side3_re.set_ylabel("Index")
    side5.autoscale_view()
    side3.autoscale_view()
    side3_re.autoscale_view()
    legend_elements = [Line2D([0], [0], color=CDS_colors["ISFV"], lw=4, label='ISFV-CDS'),
                       Line2D([0], [0], color=CDS_colors["MBFV"], lw=4, label='MBFV-CDS'),
                       Line2D([0], [0], color=CDS_colors["NKV"], lw=4, label='NKV-CDS'),
                       Line2D([0], [0], color=CDS_colors["TBFV"], lw=4, label='TBFV-CDS'),
                       Line2D([0], [0], color='gray', lw=4, label='UTR'),
                       Line2D([0], [0], color='r', lw=4, label='Interaction'),
                       Line2D([0], [0], color='orange', lw=4, label='Subopt')]
    side5.legend(handles=legend_elements,loc="upper left")
    side3.legend(handles=legend_elements,loc="upper right")
    side3_re.legend(handles=legend_elements,loc="upper left")
    plt.suptitle("Interaction Lineplot")
    plt.savefig(output)
    plt.close()


def draw_energy_histo(df, output):
    """Plots a simple histogram of the of the MFEs of a dataframe.
    df (df): the dataframe outputted by main_intarna
    output (str): Filepath where the plot should be saved
    """
    #print(df["suboptes"]) # These need to go into the histogram..somehow or just a second?
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
    plt.close()


def draw_energy_histo_subopt(df, output):
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
