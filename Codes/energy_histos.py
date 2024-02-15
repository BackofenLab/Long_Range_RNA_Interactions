import ast
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

def plot_energy_histos(IntaRNA_df, output_dir):
    classlist = []
    energylist = []
    grouplist = []
    interaction_number_list = []
    for index, row in IntaRNA_df.iterrows():
        constrained_predictions_t = ast.literal_eval(row["constrained_predictions_t"])
        constrained_predictions_q = ast.literal_eval(row["constrained_predictions_q"])
        for i in range(0, len(constrained_predictions_t)):
            interaction_t = constrained_predictions_t[i]
            classlist.append(row["class"])
            energylist.append(float(interaction_t[2]))
            interaction_number_list.append(i)
            if int(interaction_t[0]) <= 0:
                grouplist.append(0) ## Add to group 1 (5' pos <= 0)
            else:
                grouplist.append(1) ## Add to group 2 (5' pos > 0)
    edf = pd.DataFrame({"vclass": classlist, "interaction_number": interaction_number_list, "group": grouplist, "energy": energylist})
    for vclass in edf.vclass.unique():
        for pos_group in [0, 1]:
            filtered_edf = edf.loc[(edf["vclass"]==vclass) & (edf["group"]==pos_group)]
            #print(filtered_edf)
            #raise
            sns.histplot(filtered_edf, x="energy", hue="interaction_number", multiple="stack", palette={0:"red", 1:"orange", 2:"yellow"})
            plt.savefig(f"{output_dir}/energy_{vclass}_{pos_group}.png")
            plt.close()