import subprocess
import os
import re


def run_cmbuild(input_stk, output_cm):
    """Create a covariance models from a stockholm file.
    
    input_stk (str): Filepath to a stockholm file
    output_cm (str): Filepath of the resulting CM file
    """
    print(f"cmbuild: {input_stk} `=> {output_cm}")
    cmd = ["cmbuild", output_cm, input_stk]
    p = subprocess.Popen(cmd, stdout=subprocess.PIPE)
    p.wait()

    
def calibrate_cm(cm):
    """Calibrate the created covariance model.
    Takes a long time and is expensive on CPU.
    
    cm (str): Filepath of the CM file to be calibrated
    """
    cmd = ["cmcalibrate", cm]
    p = subprocess.Popen(cmd)
    p.wait()

  
def create_cms(stk_dir, out_cm_dir, calibrate=True):
    """Creates a directory of covariance models out of a given directory of
    stockholm files and calibrates them afterwards.
    Warning: Calibrating is slow and expensive.
    
    stk_dir (str): Filepath to the directory of stockholm files.
    out_cm_dir (str): Filepath for the directory of finished cm files
    calibrate (bool): If True, calibrates cmfiles after creating them
    """
    os.makedirs(out_cm_dir, exist_ok=True)
    expression = re.compile(".*3SL.*")
    for root, dirs, files in os.walk(stk_dir):
        for file in files:
            if expression.match(file):
                cmfile_path = f"{out_cm_dir}/{file.split('.')[0]}_3SL.cm"
                run_cmbuild(f"{root}/{file}", cmfile_path)
                if calibrate:
                    calibrate_cm(cmfile_path)


def run_cm_search(cm_file, seq_file, output):
    """
    Runs cm_search with the given files.
    Should look like:
        cm_file = {cm_dir}/XXXX_3SL.cm
        seq_file = {seq_dir}/YYYY/XXXX/all_XXXX_3UTR.fa
    
    cm_file (str): Filepath for CM file
    seq_file (str): Filepath for Sequence file 
    output (str): Filepath of the resulting output table
    """
    cmd = ["cmsearch", "--tblout", output, "--toponly", cm_file, seq_file]
    p = subprocess.Popen(cmd, stdout=subprocess.PIPE)
    p.wait()
    return p
    
def cm_search(df, cm_dir, seq_dir, output_path):
    """Applies cm_search with a given cm directory on a given
    Sequence directory.

    df (df): Dataframe that will be updated with the locations of found matches.
    cm_dir (str): Path to the directory of .cm files. Files need to contain
                  a "_" after their main name. (eg.: "JEVG_3SL.cm")
                  This "main name" also should appear in the seq_dir directory.
                  This can be ignored if cm_dir was also built with this program.
    seq_dir (str): Path to the sequence directory/database.
    """
    os.makedirs(output_path, exist_ok=True)
    cm_dir_names = os.listdir(cm_dir)      ## This part is maybe unneccessary
    cm_end = cm_dir_names[0].split("_")[1] ## But I want to ensure modability
    cm_files = [i.split("_")[0] for i in os.listdir(cm_dir)]
    for file in os.listdir(seq_dir):
        if file.endswith("3UTR.fa"): ## Find the right file
        
        ####dir = root.split("/")[-1]
        ####if dir in cm_files:  # = Match in cm_dir
        ####for file in files:
        ####    if file.endswith("3UTR.fa"): 
            scores = {}
            for cm_file in cm_files:
                if not os.path.isfile(f"{cm_dir}/{cm_file}_{cm_end}"):
                    continue # Skip directories
                out_file = f"{output_path}/{cm_file}.cmout"
                p = run_cm_search(f"{cm_dir}/{cm_file}_{cm_end}", 
                                  f"{seq_dir}/{file}", out_file)
                f = open(out_file, "r")
                for line in f.readlines():
                    if line.startswith("#"):
                        continue
                    split_line = line.split()
                    id = split_line[17]
                    score = split_line[14]
                    if (not id in scores) or (scores[id] < score):
                        scores[id] = score # Update score if its better
                        f, t = split_line[7:9] # Save hit with best score
                        df.loc[df.id == id, "cm_hit_f"] = int(f)
                        df.loc[df.id == id, "cm_hit_t"] = int(t)
                        df.loc[df.id == id, "cm_hit_src"] = cm_file
    df.to_csv(f"{output_path}/Inta_plus_CM.csv", index=False)
    return df