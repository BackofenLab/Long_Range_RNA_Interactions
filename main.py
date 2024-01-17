from Codes.parameters import write_static_parameters, create_parameter_table
from Codes.evaluation import draw_lineplots, draw_energy_histo, draw_energy_histo_subopt
from Codes.IntaRNA_code import main_intarna
from Codes.covariance import create_cms, cm_search
from Codes.meme import get_meme_sequences
from Codes.locarna import main_locarna
from Codes.meme_to_lineplot import meme_to_lineplot, glam2_to_lineplot
from Codes.locarna_with_mrri import main_mrri, main_loc_with_mrri
from Codes.cds_to_protein import cds_to_proteins
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

## MRRI+locARNA paths:
raw_MRRI_output = "Results/MRRI_raw_output.txt"
mrri_file_path = "Results/MRRI_output.csv"
output_loc_mmri_path = "Results/locARNA_with_MRRI"

## Outputs
energy_histo = "Results/energy_histo.png"
line_plot_output = "Results/interaction_lineplot.png"

meme_output = "Results/MEME"
locarna_output = "Results/locARNA"
amino_acids_output = "Results/AminoAcids.fa"

## Static Parameters
extra_bases = 200
extra_bases_roi = 100
outNumber = 4 # Allow n-1 subops, must be at least 1

CDS_left = 30+0
CDS_right = 70+0
CMHit_left = 30+0
CMHit_right = 30+0

## IntaRNA specific:
static_d = {"energyVRNA": "Data/rna_andronescu2007.par",
            "intLenMax": 20,
            "seedBP": 5,
            "accW": 50,
            "accL": 50,
            }


if __name__ == "__main__":
    #write_static_parameters(static_d, static_param_path)
    #create_parameter_table(database_path, extra_bases, parameter_table_file)
    #main_intarna(static_param_path, extra_bases, extra_bases_roi,
    #             parameter_table_file, IntaRNA_output, raw_IntaRNA_output,
    #             outNumber)
    #df = pd.read_csv(IntaRNA_output)
    ###create_cms(stockholm_directory, covariance_dir) ## Do not uncomment unless new data.
    #cm_search(df, covariance_dir, database_path, cm_output)
    #df = pd.read_csv(f"{cm_output}/Inta_plus_CM.csv")
    #draw_lineplots(df, extra_bases_roi, line_plot_output)
    #draw_energy_histo(df, energy_histo)
    #draw_energy_histo_subopt(df, energy_histo)
    #get_meme_sequences(parameter_table_file, f"{cm_output}/Inta_plus_CM.csv", meme_output)    
    #meme_to_lineplot(df, extra_bases_roi, meme_output)
    #glam2_to_lineplot(df, extra_bases_roi, meme_output) # HIGHLY improvised..
    #main_locarna(parameter_table_file, f"{cm_output}/Inta_plus_CM.csv", locarna_output, , CDS_left, CDS_right, CMHit_left, CMHit_right)
    main_mrri(parameter_table_file, static_param_path, extra_bases, extra_bases_roi, mrri_file_path, raw_MRRI_output)
    main_loc_with_mrri(mrri_file_path, f"{cm_output}/Inta_plus_CM.csv", parameter_table_file, output_loc_mmri_path, CDS_left, CDS_right, CMHit_left, CMHit_right)
    #cds_to_proteins(parameter_table_file, amino_acids_output, extra_bases_roi)