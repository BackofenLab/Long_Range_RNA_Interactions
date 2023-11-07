from Codes.parameters import write_static_parameters, create_parameter_table
from Codes.evaluation import draw_lineplots, draw_energy_histo
from Codes.IntaRNA_code import main_intarna
import pandas as pd


database_path = "Data/Flavivirus_NCBI/Flavivirus_RefSeq_20220621"
static_param_path = "Data/static_parameter.cfg"
parameter_table_file = "Data/parameter_table.csv"
raw_IntaRNA_output = "Results/IntaRNA_raw_output.txt"
IntaRNA_output = "Results/IntaRNA_output.csv"

## Outputs
energy_histo = "Results/energy_histo.png"
line_plot = "Results/interaction_lineplot.png"

## Static Parameters
extra_bases = 200
extra_bases_roi = 100
## IntaRNA specific:
static_d = {"energyVRNA": "Data/rna_andronescu2007.par",
            "intLenMax": 20,
            "seedBP": 7}


if __name__ == "__main__":
    write_static_parameters(static_d, static_param_path)
    create_parameter_table(database_path, extra_bases, parameter_table_file)
    main_intarna(database_path, static_param_path, extra_bases, extra_bases_roi,
                 parameter_table_file, IntaRNA_output, raw_IntaRNA_output)
    df = pd.read_csv("Results/IntaRNA_output.csv")
    draw_lineplots(df, extra_bases_roi, line_plot)
    draw_energy_histo(df, energy_histo)
