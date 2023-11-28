from Codes.parameters import write_static_parameters, create_parameter_table
from Codes.evaluation import draw_lineplots, draw_energy_histo, draw_energy_histo_subopt
from Codes.IntaRNA_code import main_intarna
from Codes.covariance import create_cms, cm_search
from Codes.alignment import get_alignment_sequences
from Codes.meme import get_meme_sequences
import pandas as pd

## Input IntaRNA Paths:
database_path = "Data/Flavivirus_NCBI/Flavivirus_RefSeq_20231111"

## Output IntaRNA Paths:
static_param_path = "Data/static_parameter.cfg"
parameter_table_file = "Data/parameter_table.csv"
raw_IntaRNA_output = "Results/IntaRNA_raw_output.txt"
IntaRNA_output = "Results/IntaRNA_output.csv"

## Input Covariance Model Paths:
stockholm_directory = "Data/Flavivirus_Stockholm"

## Output Covariance Model Paths:
covariance_dir = "Data/Flavivirus_Covariance"
cm_output = "Results/cm_search"

## Outputs
energy_histo = "Results/energy_histo.png"
line_plot = "Results/interaction_lineplot.png"

meme_output = "Results/MEME"

## Static Parameters
extra_bases = 200
extra_bases_roi = 100
outNumber = 4 # Allow n-1 subops, must be at least 1

## IntaRNA specific:
static_d = {"energyVRNA": "Data/rna_andronescu2007.par",
            "intLenMax": 20,
            "seedBP": 5,
            "accW": 50,
            "accL": 50,
            }


if __name__ == "__main__":
    write_static_parameters(static_d, static_param_path)
    create_parameter_table(database_path, extra_bases, parameter_table_file)
    main_intarna(static_param_path, extra_bases, extra_bases_roi,
                 parameter_table_file, IntaRNA_output, raw_IntaRNA_output,
                 outNumber)
    df = pd.read_csv(IntaRNA_output)
    ###create_cms(stockholm_directory, covariance_dir) ## Do not uncomment unless new data.
    df = cm_search(df, covariance_dir, database_path, cm_output)
    draw_lineplots(df, extra_bases_roi, line_plot)
    draw_energy_histo(df, energy_histo)
    draw_energy_histo_subopt(df, energy_histo)
    ##get_alignment_sequences(parameter_table_file, IntaRNA_output, extra_bases)
    get_meme_sequences(parameter_table_file, meme_output, IntaRNA_output, extra_bases)
