from Codes.parameters import write_static_parameters, create_parameter_table
from Codes.evaluation import draw_lineplots, draw_energy_histo, draw_energy_histo_subopt
from Codes.IntaRNA_code import main_intarna
from Codes.covariance import create_cms, cm_search
from Codes.meme import get_meme_sequences
from Codes.locarna import main_locarna
from Codes.meme_to_lineplot import meme_to_lineplot, glam2_to_lineplot
from Codes.locarna_with_mrri import main_mrri, main_loc_with_mrri
from Codes.cds_to_protein import cds_to_proteins
from Codes.locarna_consensus import get_all_locarna_consensus
from Codes.energy_histos import plot_energy_histos
import pandas as pd
import os

tasks = { # Note: You cannot run later tasks without running the earlier ones at least once.
        "create_parameter_tables" : 0,
        "run_IntaRNA"             : 0,
        "CREATE_CMs"              : 0, ##!## Takes very long. Do not set to true unless new data.
        "run_CM_search"           : 0,
        "draw_IntaRNA_plots"      : 0,
        "MEME+GLAM2"              : 0,
        "run_locARNA"             : 0,
        "run_MRRI_1"              : 0, ## Restrictions like "run_IntaRNA"
        "run_MRRI_2"              : 0, ## Further restrictions
        "locARNA+MRRI"            : 0,
        "locARNA+MRRI+CARNA"      : 1, ## Old version of locARNA that also uses CARNA
        "draw_MRRI_plots"         : 0,
        "CDS_to_proteins"         : 0,
        }

## Input IntaRNA Paths:
database_path = "Data/Flavivirus_NCBI/Flavivirus_RefSeq_20231111"

## General results directory
results = "Results"
os.makedirs(results, exist_ok=True)

## Output IntaRNA Paths:
static_param_path = "Data/static_parameter.cfg"
parameter_table_file = "Data/parameter_table.csv"
raw_IntaRNA_output = f"{results}/IntaRNA_raw_output.txt"
IntaRNA_output = f"{results}/IntaRNA_output.csv"

## Input Covariance Model Paths:
stockholm_directory = "Data/Flavivirus_Stockholm"

## Output Covariance Model Paths:
covariance_dir = "Data/Flavivirus_Covariance"
cm_output_dir = f"{results}/cm_search"
cm_search_file = f"{cm_output_dir}/Inta_plus_CM.csv"

## MRRI+locARNA paths:
raw_MRRI_output_1 = f"{results}/MRRI_raw_output_1.txt"
raw_MRRI_output_2 = f"{results}/MRRI_raw_output_2.txt"
mrri_file_path_1 = f"{results}/MRRI_output_1.csv"
mrri_file_path_2 = f"{results}/MRRI_output_2.csv"
output_loc_mmri_path = f"{results}/locARNA_with_MRRI"
output_loc_mmri_carna_path = f"{results}/locARNA_with_MRRI_only_cm_pos"
output_loc_mmri_carna_path_mode_3 = f"{results}/locARNA_with_MRRI_only_inter"
mrri_lineplot_path_1 = f"{results}/interaction_lineplot_MRRI_1.png"
mrri_lineplot_path_2 = f"{results}/interaction_lineplot_MRRI_2.png"

## Outputs
energy_histo = f"{results}/old_energy_histo.png"
lineplot_output = f"{results}/interaction_lineplot.png"

meme_output = f"{results}/MEME"
locarna_output = f"{results}/locARNA"
amino_acids_output = f"{results}/AminoAcids.fa"

## Static Parameters
extra_bases = 200
extra_bases_roi = 100
outNumber = 4 # Allow n-1 subops, must be at least 1

## Region for locARNA
CDS_left = 40
CDS_right = 70
CMHit_left = 30
CMHit_right = 30

## IntaRNA specific:
static_d = {"energyVRNA": "Data/rna_andronescu2007.par",
            "intLenMax": 20,
            "seedBP": 5,
            "accW": 50,
            "accL": 50,
            "temperature": 18,
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
        df = pd.read_csv(cm_search_file)
        get_meme_sequences(parameter_table_file, cm_search_file, meme_output)    
        meme_to_lineplot(df, extra_bases_roi, meme_output)
        glam2_to_lineplot(df, extra_bases_roi, meme_output) # HIGHLY improvised..
    if tasks["run_locARNA"]:
        main_locarna(parameter_table_file, cm_search_file, locarna_output, CDS_left, CDS_right, CMHit_left, CMHit_right)
    if tasks["run_MRRI_1"]:
        param_mode = 1
        main_mrri(parameter_table_file, static_param_path, extra_bases, extra_bases_roi, mrri_file_path_1, raw_MRRI_output_1, param_mode)
    if tasks["run_MRRI_2"]:
        param_mode = 2
        main_mrri(parameter_table_file, static_param_path, extra_bases, extra_bases_roi, mrri_file_path_2, raw_MRRI_output_2, param_mode)
    if tasks["locARNA+MRRI"]:
        main_loc_with_mrri(mrri_file_path_2, cm_search_file, parameter_table_file, output_loc_mmri_path, CDS_left, CDS_right, CMHit_left, CMHit_right, cm_output_dir, use_carna=False)
        #get_all_locarna_consensus(output_loc_mmri_path)
    if tasks["locARNA+MRRI+CARNA"]:
        #main_loc_with_mrri(mrri_file_path_2, cm_search_file, parameter_table_file, output_loc_mmri_carna_path, CDS_left, CDS_right, CMHit_left, CMHit_right, cm_output_dir, use_carna=False, skip_FS=True, mode=2)
        main_loc_with_mrri(mrri_file_path_2, cm_search_file, parameter_table_file, output_loc_mmri_carna_path_mode_3, CDS_left, CDS_right, CMHit_left, CMHit_right, cm_output_dir, use_carna=False, skip_FS=True, mode=3)
        #get_all_locarna_consensus(output_loc_mmri_carna_path)
    if tasks["draw_MRRI_plots"]:
        cmdf = pd.read_csv(cm_search_file)
        if os.path.isfile(mrri_file_path_1):
            mrri_df_1 = pd.read_csv(mrri_file_path_1)
            mrri_df_1["cm_hit_f"] = pd.Series(cmdf["cm_hit_f"])
            mrri_df_1["cm_hit_t"] = pd.Series(cmdf["cm_hit_t"])
            mrri_df_1["cm_hit_src"] = pd.Series(cmdf["cm_hit_src"])
            draw_lineplots(mrri_df_1, extra_bases_roi, mrri_lineplot_path_1)
        if os.path.isfile(mrri_file_path_2):
            mrri_df_2 = pd.read_csv(mrri_file_path_2)
            mrri_df_2["cm_hit_f"] = pd.Series(cmdf["cm_hit_f"])
            mrri_df_2["cm_hit_t"] = pd.Series(cmdf["cm_hit_t"])
            mrri_df_2["cm_hit_src"] = pd.Series(cmdf["cm_hit_src"])
            draw_lineplots(mrri_df_2, extra_bases_roi, mrri_lineplot_path_2, draw_roi_box=True)
            plot_energy_histos(mrri_df_2, results)
    if tasks["CDS_to_proteins"]:
        cds_to_proteins(parameter_table_file, amino_acids_output, extra_bases_roi)
    