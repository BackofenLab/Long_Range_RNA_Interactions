from Codes.parameters import write_static_parameters, create_parameter_table
from Codes.evaluation import draw_lineplots, draw_energy_histo
from Codes.IntaRNA_code import main_intarna
import pandas as pd


database_path = "Data/Flavivirus_NCBI/Flavivirus_RefSeq_20220621"
static_param_file = "Data/static_params.csv"
parameter_table_file = "Data/parameter_table.csv"
raw_IntaRNA_output = "Results/IntaRNA_raw_output.txt"
IntaRNA_output = "Results/IntaRNA_output.csv"

## Outputs
energy_histo = "energy_histo.png"
line_plot = "interaction_lineplot.png"


static_d = {"energyVRNA": ["Data/rna_andronescu2007.par"],
            "intLenMax": [20],
            "extra_bases": [200],
            "extra_bases_roi": [100],
            "seedlength": [7]}

if __name__ == "__main__":
    write_static_parameters(static_d, static_param_file)
    create_parameter_table(database_path, static_d["extra_bases"][0], parameter_table_file)
    main_intarna(database_path, static_param_file, parameter_table_file,
                 IntaRNA_output, raw_IntaRNA_output)
    df = pd.read_csv("Results/IntaRNA_output.csv")
    draw_lineplots(df, static_d["extra_bases"][0], line_plot)
    draw_energy_histo(df, energy_histo)
