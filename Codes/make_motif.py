import logomaker
import matplotlib.pyplot as plt
from matplotlib.patches import ConnectionPatch
import pandas as pd
import numpy as np
from collections import defaultdict
import os
import argparse


def fasta_to_df(fasta, reverse_sequences=False, drop_minus=False):
    valid_characters = ("A","C","G","T","U","-")
    count_dict = defaultdict(lambda: defaultdict(int))
    with open(fasta, "r") as f:
        lines = [line.rstrip() for line in f]
        for line in lines:
            if line.startswith(";interaction"):
                interaction = line.split(" ")[-1]
            elif not line.startswith(valid_characters):
                continue
            else:
                if reverse_sequences:
                    line = "".join(reversed(line))
                    interaction = "".join(reversed(interaction))
                for i in range(0, len(line)):
                    count_dict[i][line[i]] += 1
    df = pd.DataFrame(count_dict).T         # Transpose
    df = df.fillna(0)                       # Replace NaN with 0
    result = df.div(df.sum(axis=1), axis=0) # Divide by the sum of each row
    if drop_minus:
        result = result.drop("-", axis=1, errors='ignore') # Remove "-" from the results
    return result, interaction


def add_interactions(motif, motif2, inter_left, inter_right, color="cornflowerblue",
                     vertical=False, dynamic_alpha=True):
    alpha = alpha2 = 1
    linewidth = 8
    alt_bracket_color = "firebrick"
    if vertical:
        arc_rad = 0
        height_left = -0.15
        height_right = 1.05
    else:
        arc_rad = -0.2
        height_left = height_right = 1.02
        inter_right = "".join(reversed(inter_right))
    left_brackets_1 = inter_left.count("(")
    right_brackets_1 = inter_right.count(")")
    left_brackets_2 = inter_left.count("[")
    right_brackets_2 = inter_right.count("]")
    if left_brackets_1 != right_brackets_1:
        raise(ValueError(f"'(' Interactions do not contain the same amount of brackets ({left_brackets_1} vs {right_brackets_1})"))
    if left_brackets_2 != right_brackets_2:
        raise(ValueError(f"'[' Interactions do not contain the same amount of brackets ({left_brackets_2} vs {right_brackets_2})"))
    if not (len(inter_left) == len(motif.df)) or not (len(inter_right) == len(motif2.df)):
        raise(ValueError("Interactions need to have the same lengths as their sequences"))
    l_pos = 0
    r_pos = 0
    inter_count = 0
    while l_pos < len(motif.df) and r_pos < len(motif2.df):
        if vertical:
            xyB_pos = r_pos
        else:
            xyB_pos = len(motif2.df)-r_pos-1
        #print(l_pos, r_pos, inter_left[l_pos], inter_right[r_pos])
        if inter_left[l_pos] != "(" and inter_left[l_pos] != "[":
            l_pos += 1
        elif inter_right[r_pos] != ")" and inter_right[r_pos] != "]":
            r_pos += 1
        elif inter_left[l_pos] == "(" and inter_right[r_pos] == ")":
            inter_count += 1
            tmp_color = color
            if inter_count % 5 == 0:
                tmp_color = "black"
            con = ConnectionPatch(xyA=(l_pos, height_left), xyB=(xyB_pos, height_right), coordsA=motif.ax.transData, coordsB=motif2.ax.transData,
                  axesA=motif, axesB=motif2, connectionstyle=f"arc3,rad={arc_rad}",color=tmp_color, alpha=alpha, linewidth=linewidth)
            motif2.ax.add_artist(con)
            if dynamic_alpha:
                alpha -= (1/(left_brackets_1+3))
            l_pos += 1
            r_pos += 1
        elif inter_left[l_pos] == "[" and inter_right[r_pos] == "]":
            inter_count += 1
            tmp_color = alt_bracket_color
            if inter_count % 5 == 0:
                tmp_color = "black"
            con = ConnectionPatch(xyA=(l_pos, height_left), xyB=(xyB_pos, height_right), coordsA=motif.ax.transData, coordsB=motif2.ax.transData,
                  axesA=motif, axesB=motif2, connectionstyle=f"arc3,rad={arc_rad}", color=tmp_color, alpha=alpha2, linewidth=linewidth)
            motif2.ax.add_artist(con)
            if dynamic_alpha:
                alpha2 -= (1/(left_brackets_2+3))
            l_pos += 1
            r_pos += 1
        else:
            raise(ValueError("No crossing between '(' and '[' allowed."))


def create_motif(input_fasta_left, input_fasta_right, output_file="", 
                 drop_minus=True, color="cornflowerblue", vertical=False, dynamic_alpha=True,
                 plot_xaxis=False, plot_yaxis=False):
    fasta_df_left, inter_left = fasta_to_df(input_fasta_left, drop_minus=drop_minus)
    if vertical:
        plt.figure(figsize=(10, 5))
        plots_vertical = 2
        plots_horizontally = 1
        reverse_right_sequences = True
    else:
        plt.figure(figsize=(30, 5))
        plots_vertical = 1
        plots_horizontally = 2
        reverse_right_sequences = False
    ## Generate first motif
    plt.subplot(plots_vertical, plots_horizontally, 1)
    motif = logomaker.Logo(fasta_df_left,
                           #font_name='Arial Rounded MT Bold',
                           color_scheme="classic", ax=plt.gca())
    ## Generate second motif
    fasta_df_right, inter_right = fasta_to_df(input_fasta_right, 
                                 reverse_sequences=reverse_right_sequences,
                                 drop_minus=drop_minus)
    plt.subplot(plots_vertical, plots_horizontally, 2)
    motif2 = logomaker.Logo(fasta_df_right,
                            #font_name='Arial Rounded MT Bold',
                            color_scheme="classic", ax=plt.gca())
    ## Connect the motifs with their interactions
    if inter_left and inter_right:
        add_interactions(motif, motif2, inter_left, inter_right, color=color, vertical=vertical, dynamic_alpha=dynamic_alpha)
    # style using Logo methods
    motif.style_spines(visible=False)
    motif.style_spines(spines=['left', 'bottom'], visible=True)
    motif.style_xticks(rotation=90, fmt='%d', anchor=0)
    motif2.style_spines(visible=False)
    motif2.style_spines(spines=['left', 'bottom'], visible=True)
    motif2.style_xticks(rotation=90, fmt='%d', anchor=0)

    # style using Axes methods
    #motif.ax.set_ylabel("Amount", labelpad=-1)
    motif.ax.xaxis.set_ticks_position('none')
    motif2.ax.xaxis.set_ticks_position('none')
    motif.ax.xaxis.set_tick_params(pad=-1)
    motif2.ax.xaxis.set_tick_params(pad=-1)
    #motif.ax.tick_params(axis="x", labelsize=24)
    #motif2.ax.tick_params(axis="x", labelsize=24)
    #motif.ax.set_yticks([0.2, 0.4, 0.6, 0.8, 1])
    #motif2.ax.set_yticks([0.2, 0.4, 0.6, 0.8, 1])
    if vertical:
        plt.subplots_adjust(hspace=0.5, bottom=0.05, top=0.99, left=0.05, right=0.99)
    else:
        #plt.tight_layout()
        plt.subplots_adjust(wspace=0.05, bottom=0.05, top=0.4, left=0.01, right=0.99)
    if not plot_xaxis:
        motif.ax.get_xaxis().set_visible(False)
        motif2.ax.get_xaxis().set_visible(False)
        plt.subplots_adjust(bottom=0.01)
    if not plot_yaxis:
        motif.ax.get_yaxis().set_visible(False)
        motif2.ax.get_yaxis().set_visible(False)
        plt.subplots_adjust(left=0.01)

    # style and show figure
    if not output_file:
        output_file = input_fasta_left.replace("_left.fa", ".png")
    plt.savefig(output_file)
    plt.close()


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("-p", "--path", dest="directory",
                        type=str, help="Directory with fasta files", required=True)
    args = vars(parser.parse_args())
    for root, dirs, files in os.walk(args["directory"]):
        for file in files:
            if file.endswith("_left.fa"):
                valid_file_left = f"{root}/{file}"
                valid_file_right = valid_file_left.replace("_left.fa", "_right.fa")
                color = valid_file_left.split("_")[1]
                for vertical in [True, False]:
                    if vertical:
                        orientation = "vertical"
                    else:
                        orientation = "horizontal"
                    output_name = valid_file_left.replace("_left.fa", f"_{orientation}.png")
                    create_motif(valid_file_left, valid_file_right, 
                                 output_file=output_name, color=color, 
                                 vertical=vertical, dynamic_alpha=False)
