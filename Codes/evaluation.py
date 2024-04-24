import seaborn as sns
import matplotlib.pyplot as plt
from matplotlib.collections import LineCollection
from matplotlib.lines import Line2D
from matplotlib.patches import Rectangle
import ast
import pandas as pd
import numpy as np
import math
from collections import defaultdict


main_colours = {"UTR": "gray", "main_interaction" : "red", "subopt_interaction": "orange", "interactions": ["red", "orange", "yellow"], "MEME": "red"}
CDS_colours = defaultdict(lambda:"black")
for k, v in [("ISFV","m"), ("cISFVG","m"), ("dISFVG","plum"),("MBFV","blue"),("NKV","green"),("TBFV","c")]:
    CDS_colours[k] = v
linewidths = {"UTR": 4,
              "CDS": 4,
              "CMHit": 8,
              "Subopt": 4,
              "Interaction": 4,
              "MEME": 16,
              }

def draw_lineplots(df, extra_bases_roi, output, subopt_mode=False, draw_roi_box=False, meme_sites={}):
    """
    Draw an interaction lineplot 
    showcasing interaction position relative to UTR/CDS.
    
    df (df): the dataframe outputted by main_intarna
    extra_bases_roi (int): Amount of extra bases from the CDS used on each side
    output (str): Filepath where the plot should be saved
    subopt_mode (bool): If True: First interaction is assumed to be main, rest is subopt
                        If False: First is first colour, second is second colour etc
    meme_sites (dict): Optional dictionary of meme sites. generated by meme_to_lineplot.
                       Note: Give the class specific dict here not the mega dict!
                       (e.g.: mega_dict[virus_class])
    """
    ## Need 2 plots:
    ## 1: UTR5+some_CDS+t_inter (aligned to UTR5/CDS transition)
    ## 2: some_CDS+UTR3+q_inter (aligned to end of UTR3)
    no_subopts = False
    no_predictions = False
    no_cm_hits = False
    appeared_virus_types = {}
    plt.style.use("seaborn-v0_8-darkgrid")
    fig, (side5, side3) = plt.subplots(1, 2, figsize=(16, 9)) # Formerly (24, 20)

    for index, row in df.iterrows():
        appeared_virus_types[row["type"] if row["class"] == "ISFV" else row["class"]] = 1
        line5sub = []
        line3sub = []
        ## Plot CDS:
        ## Differentiate between cISFVG and dISFVG
        cds_colour = CDS_colours[row["type"]] if row["class"] == "ISFV" else CDS_colours[row["class"]]
        side5.plot((0, extra_bases_roi), (index, index), 
                   linewidth=linewidths["CDS"], color=cds_colour, solid_capstyle="butt")
        side3.plot((0 - row["UTR3len"], -extra_bases_roi - row["UTR3len"]), (index, index), 
                   linewidth=linewidths["CDS"], color=cds_colour, solid_capstyle="butt")
        ## Plot UTR:
        side5.plot((0, -row["UTR5len"]), (index, index), 
                   linewidth=linewidths["UTR"], color=main_colours["UTR"], solid_capstyle="butt")
        side3.plot((0 - row["UTR3len"], 0), (index, index), 
                   linewidth=linewidths["UTR"], color=main_colours["UTR"], solid_capstyle="butt")

        ## Plot CM hits:
        if "cm_hit_f" in row:
            if not math.isnan(row["cm_hit_f"]):
                side3.plot((row["cm_hit_f"] - row["UTR3len"], row["cm_hit_t"] - row["UTR3len"]), (index, index), 
                           linewidth=linewidths["CMHit"], color=CDS_colours[row["cm_hit_src"]], alpha=0.5, solid_capstyle="butt")
        else:
            no_cm_hits = True
        ## Plot Interaction Predictions:
        if "predictions_t" in row:
            counter = 0
            predictions_t = ast.literal_eval(row["predictions_t"])
            for prediction_t in predictions_t: ## 5' Interaction predictions
                if prediction_t:
                    if subopt_mode:
                        if counter == 0:
                            color = main_colours["main_interaction"]
                        else:
                            color = main_colours["subopt_interaction"]
                    else:
                        color = main_colours["interactions"][counter]
                    side5.plot((int(prediction_t[0]), int(prediction_t[1])), (index, index), 
                                   linewidth=linewidths["Interaction"], color=color, solid_capstyle="butt", zorder=(len(predictions_t)-counter)*10)
                    counter += 1
            counter = 0
            predictions_q = ast.literal_eval(row["predictions_q"])
            for prediction_q in predictions_q: ## 3' Interaction predictions
                if prediction_q:
                    if subopt_mode:
                        if counter == 0:
                            color = main_colours["main_interaction"]
                        else:
                            color = main_colours["subopt_interaction"]
                    else:
                        color = main_colours["interactions"][counter]
                    side3.plot((int(prediction_q[0]) - row["UTR3len"], int(prediction_q[1]) - row["UTR3len"]), (index, index), 
                                   linewidth=linewidths["Interaction"], color=color, solid_capstyle="butt", zorder=(len(predictions_t)-counter)*10)
                    counter += 1
        else:
            no_predictions = True

    ## Plot MEME sites:
        if meme_sites:
            ## Plot Site 1:
            if row["id"] in meme_sites["site_1"]:
                site1 = meme_sites["site_1"][row["id"]]
                site1_area_start = -40
                site1_tuple = (site1_area_start + site1[0], site1_area_start + site1[1]) 
                side5.plot(site1_tuple, (index, index),  linewidth=linewidths["MEME"],
                           color=main_colours["MEME"],  alpha=0.3,  solid_capstyle="butt", zorder=5000)
            ## Plot Site 2:
            if row["id"] in meme_sites["site_2"]:
                site2 = meme_sites["site_2"][row["id"]]
                site2_area_start = 20
                site2_tuple = (site2_area_start + site2[0], site2_area_start + site2[1])
                side5.plot(site2_tuple, (index, index),  linewidth=linewidths["MEME"],
                           color=main_colours["MEME"],  alpha=0.3,  solid_capstyle="butt", zorder=5000)
            ## Plot Site 3:
            if row["id"] in meme_sites["site_3"]:
                CMhit_start = row["cm_hit_f"] - row["UTR3len"]
                site3 = meme_sites["site_3"][row["id"]]
                site3_area_start = CMhit_start-25
                site3_tuple = (site3_area_start + site3[0], site3_area_start + site3[1])
                side3.plot(site3_tuple, (index, index),  linewidth=linewidths["MEME"],
                           color=main_colours["MEME"],  alpha=0.3,  solid_capstyle="butt", zorder=5000)
            ## Plot Site 4:
            if row["id"] in meme_sites["site_4"]:
                CMhit_start = row["cm_hit_f"] - row["UTR3len"]
                site4 = meme_sites["site_4"][row["id"]]
                site4_area_start = CMhit_start+5
                site4_tuple = (site4_area_start + site4[0], site4_area_start + site4[1])
                side3.plot(site4_tuple, (index, index),  linewidth=linewidths["MEME"],
                           color=main_colours["MEME"],  alpha=0.3, solid_capstyle="butt", zorder=5000)

    ## Draw Box marking limited region of interest:
    if draw_roi_box:
        rect5 = Rectangle((-40,-1), 110, len(df)+1, color='red', fc = 'none', zorder=2, lw = 2)
        rect3 = Rectangle((-140,-1), 90, len(df)+1, color='red', fc = 'none', zorder=2, lw = 2)
        side5.add_patch(rect5)
        side3.add_patch(rect3)
    ## Annotations:
    side5.set_xlabel("Distance from 5'UTR-CDS transition", fontsize=16)
    side3.set_xlabel("Distance from 3'UTR end", fontsize=16)
    side5.set_ylabel("Accession number", fontsize=16)
    side3.set_ylabel("Accession number", fontsize=16)
    side5.yaxis.set_label_position("right")
    side5.yaxis.tick_right()
    side3.yaxis.set_label_position("right")
    side3.yaxis.tick_right()

    side5.set_yticks(np.arange(len(df["id"])))
    side5.set_yticklabels(list(df["id"]+"_"+df["virus"]))
    side3.set_yticks(np.arange(len(df["id"])))
    side3.set_yticklabels(list(df["id"]+"_"+df["virus"]))
    side5.autoscale_view()
    side3.autoscale_view()
    ## Legend:
    legend_elements = []
    for virus_class in appeared_virus_types.keys():
        legend_elements.append(Line2D([0], [0], color=CDS_colours[virus_class], lw=8, label=f'{virus_class}-CDS'))
    legend_elements += [Line2D([0], [0], color='gray', lw=8, label='UTR'),
                        Line2D([0], [0], color='r', lw=8, label='Interaction')]
    if subopt_mode:
        legend_elements.append(Line2D([0], [0], color='orange', lw=8, label='Subopt'))
    elif not no_predictions:
        legend_elements.append(Line2D([0], [0], color='orange', lw=8, label='Constrained Interaction 1'))
        legend_elements.append(Line2D([0], [0], color='yellow', lw=8, label='Constrained Interaction 2'))
    if meme_sites:
         legend_elements.append(Line2D([0], [0], color='r', alpha=0.3, lw=linewidths["MEME"], label='MEME Hit'))
    side5.legend(handles=legend_elements,loc="upper left", prop={"size": 16})
    if not no_cm_hits:
        legend_elements.append(Line2D([0], [0], color="blue", lw=8, label='CM Hit', alpha=0.5))
    side3.legend(handles=legend_elements,loc="upper left", prop={"size": 16})
    plt.suptitle("Interaction Lineplot", fontsize=48)
    plt.savefig(output, bbox_inches='tight')
    
    # Save plots separately too:
    #extent5 = side5.get_tightbbox().transformed(fig.dpi_scale_trans.inverted())
    #extent3 = side3.get_tightbbox().transformed(fig.dpi_scale_trans.inverted())
    #o_split = output.split(".",1)
    #print(extent5[0], extent3)
    #plt.savefig(f"{o_split[0]}_5.{o_split[1]}", bbox_inches=extent5)
    #plt.savefig(f"{o_split[0]}_3.{o_split[1]}", bbox_inches=extent3) # .expanded(1.2, 1.2)

    plt.close()