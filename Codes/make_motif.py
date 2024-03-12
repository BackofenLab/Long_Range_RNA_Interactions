import logomaker
import matplotlib.pyplot as plt
import pandas as pd
from collections import defaultdict
import os
import argparse


def fasta_to_df(fasta, drop_minus=False):
    valid_characters = ("A","C","G","T","U","-")
    count_dict = defaultdict(lambda: defaultdict(int))
    with open(fasta, "r") as f:
        lines = [line.rstrip() for line in f]
        for line in lines:
            if not line.startswith(valid_characters):
                continue
            else:
                for i in range(0, len(line)):
                    count_dict[i][line[i]] += 1
    df = pd.DataFrame(count_dict).T         # Transpose
    df = df.fillna(0)                       # Replace NaN with 0
    result = df.div(df.sum(axis=1), axis=0) # Divide by the sum of each row
    if drop_minus:
        result = result.drop("-", axis=1, errors='ignore') # Remove "-" from the results
    return result

def create_motif(input_fasta, output_file):
    fasta_df = fasta_to_df(input_fasta, drop_minus=True)
    motif = logomaker.Logo(fasta_df,
                           #font_name='Arial Rounded MT Bold',
                           color_scheme="classic")

    # style using Logo methods
    motif.style_spines(visible=False)
    motif.style_spines(spines=['left', 'bottom'], visible=True)
    motif.style_xticks(rotation=90, fmt='%d', anchor=0)

    # style using Axes methods
    motif.ax.set_ylabel("Amount", labelpad=-1)
    motif.ax.xaxis.set_ticks_position('none')
    motif.ax.xaxis.set_tick_params(pad=-1)
    motif.ax.set_yticks([0.2, 0.4, 0.6, 0.8, 1])

    # style and show figure
    motif.fig.savefig(output_file)
    #motif.fig.show()

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("-p", "--path", dest="directory",
                        type=str, help="Directory with fasta files", required=True)
    args = vars(parser.parse_args())
    for root, dirs, files in os.walk(args["directory"]):
        for file in files:
            if file.endswith(".fa"):
                valid_file = f"{root}/{file}"
                output_name = f"{valid_file.split('.')[0]}.png"
                create_motif(valid_file, output_name)
