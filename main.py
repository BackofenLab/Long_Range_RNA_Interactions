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

tasks = { # Note: You cannot run later tasks without running the earlier ones at least once.
        "create_parameter_tables" : 1,
        "run_IntaRNA"             : 1,
        "CREATE_CMs"              : 0, ##!## Takes very long. Do not set to true unless new data.
        "run_CM_search"           : 1,
        "draw_IntaRNA_plots"      : 1,
        "MEME+GLAM2"              : 1,
        "run_locARNA"             : 1,
        "run_MRRI"                : 1,
        "locARNA+MRRI"            : 1,
        "draw_MRRI_plots"         : 1,
        "CDS_to_proteins"         : 1,
        }

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
cm_output_dir = "Results/cm_search"
cm_search_file = f"{cm_output_dir}/Inta_plus_CM.csv"

## MRRI+locARNA paths:
static_param_path_MRRI = "Data/static_parameter_MRRI.cfg"
raw_MRRI_output = "Results/MRRI_raw_output.txt"
mrri_file_path = "Results/MRRI_output.csv"
output_loc_mmri_path = "Results/locARNA_with_MRRI"
mrri_lineplot_path = "Results/interaction_lineplot_MRRI.png"

## Outputs
energy_histo = "Results/energy_histo.png"
lineplot_output = "Results/interaction_lineplot.png"

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
    if tasks["create_parameter_tables"]:
        write_static_parameters(static_d, static_param_path)
        create_parameter_table(database_path, extra_bases, parameter_table_file)
    if tasks["run_IntaRNA"]:
        main_intarna(static_param_path, extra_bases, extra_bases_roi,
                     parameter_table_file, IntaRNA_output, raw_IntaRNA_output,
                     outNumber)
    if tasks["CREATE_CMs"]:  ## Do not execute unless new data.
        create_cms(stockholm_directory, covariance_dir)
    if tasks["run_CM_search"]:
        df = pd.read_csv(IntaRNA_output)
        cm_search(df, covariance_dir, database_path, cm_output_dir, cm_search_file)
    if tasks["draw_IntaRNA_plots"]:
        df = pd.read_csv(cm_search_file)
        draw_lineplots(df, extra_bases_roi, lineplot_output)
        draw_energy_histo(df, energy_histo)
        draw_energy_histo_subopt(df, energy_histo)
    if tasks["MEME+GLAM2"]:
        get_meme_sequences(parameter_table_file, cm_search_file, meme_output)    
        meme_to_lineplot(df, extra_bases_roi, meme_output)
        glam2_to_lineplot(df, extra_bases_roi, meme_output) # HIGHLY improvised..
    if tasks["run_locARNA"]:
        main_locarna(parameter_table_file, cm_search_file, locarna_output, CDS_left, CDS_right, CMHit_left, CMHit_right)
    if tasks["run_MRRI"]:
        static_d_mrri = dict(static_d)
        if "intLenMax" in static_d_mrri:
            del static_d_mrri["intLenMax"]
        write_static_parameters(static_d_mrri, static_param_path_MRRI)
        main_mrri(parameter_table_file, static_param_path_MRRI, extra_bases, extra_bases_roi, mrri_file_path, raw_MRRI_output)
    if tasks["locARNA+MRRI"]:
        main_loc_with_mrri(mrri_file_path, cm_search_file, parameter_table_file, output_loc_mmri_path, CDS_left, CDS_right, CMHit_left, CMHit_right)
    if tasks["draw_MRRI_plots"]:
        mrri_df = pd.read_csv(mrri_file_path)
        cmdf = pd.read_csv(cm_search_file)
        mrri_df["cm_hit_f"] = pd.Series(cmdf["cm_hit_f"])
        mrri_df["cm_hit_t"] = pd.Series(cmdf["cm_hit_t"])
        mrri_df["cm_hit_src"] = pd.Series(cmdf["cm_hit_src"])
        draw_lineplots(mrri_df, extra_bases_roi, mrri_lineplot_path, nosubopts=True)
    if tasks["CDS_to_proteins"]:
        cds_to_proteins(parameter_table_file, amino_acids_output, extra_bases_roi)