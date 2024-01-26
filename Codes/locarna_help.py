import os
import sys
import subprocess
import ast
#from Codes.MRRI import MRRIHandler
import Codes.MRRI_main

def run_command(cmd):
    p = subprocess.run(cmd.split(" "))


def run_mlocarna(input_fasta, output_dir):
    """Run mlocarna on a given fasta file.
    
    input_fasta (str): Filepath to a fasta file to apply locarna on
    output_dir (str): Filepath of the resulting CM file
    """
    print(f"mlocarna: {input_fasta} `=> {output_dir}")
    cmd = ["conda", "run", "-n", "locarna",
           "mlocarna", input_fasta,
           #"--indel=-50", # Webserver parameter
           #"--indel-opening=-750", # Webserver parameter
           "--width=300",
           "--tgtdir", output_dir
           ]
    #print(" ".join(cmd))
    #raise
    p = subprocess.Popen(cmd, stdout=subprocess.PIPE)
    p.wait()


def run_rnaalifold(input_dir):
    print(f"RNAalifold : {input_dir}/result.aln => {input_dir}/alirna.ps+aln.ps")
    input_file = f"{input_dir}/result.aln"
    with open(input_file, "r") as f:
        for line in f.readlines():
            if line.startswith("#A1"):
                anchor_seq = line.split(" ")[-1]
                s1, s2 = anchor_seq.split("BBBBBBB")
                constraint = b"."+b"<"*(len(s1)-1) + b"xxxxxxx" + b">"*(len(s2)-1)
                break
    cmd = ["RNAalifold", input_file,
           "--aln", "--ribosum_scoring",
           "--cfactor", "0.6",
           "--nfactor", "0.5",
           "--color", "-C"
           ]
    p = subprocess.Popen(cmd, stdout=subprocess.PIPE, stdin=subprocess.PIPE, stderr=subprocess.PIPE)
    p.communicate(input=constraint)
    p.wait()
    os.rename("alirna.ps", f"{input_dir}/alirna.ps")
    os.rename("aln.ps", f"{input_dir}/aln.ps")


def run_ps_to_pdf(ps_file, output):
    print(f"ps2pdf: {ps_file} `=> {output}")
    cmd = ["ps2pdf", "-dEPSCrop", ps_file, output]
    p = subprocess.Popen(cmd, stdout=subprocess.PIPE)
    p.wait()

def hacked_MRRI_main(UTR5pCDS, UTR3pCDS, static_param_path):
    sys.argv = ["Codes/MRRI_main.py", "-q", UTR3pCDS, "-t", UTR5pCDS, "-p", static_param_path]
    MRRIHandler = Codes.MRRI_main.main()
    for qId in MRRIHandler.querySeq.keys():
        for tId in MRRIHandler.targetSeq.keys():
            BE = dict({ 'id1': tId, 'id2' : qId})
            B1 = MRRIHandler.runIntaRNA(BE)
            #print(B1)
            B2 = MRRIHandler.runIntaRNA(B1)
            #print(B2)
            B3 = MRRIHandler.runIntaRNA(B2)
            #print(B3)
            #raise
            return B1, [B2, B3]


def run_mrri(UTR5pCDS, UTR3pCDS, static_param_path):
    """Outdated. Currently not in use"""
    cmd = ["python", "Codes/MRRI_main.py", "-q", UTR3pCDS, "-t", UTR5pCDS, "-p", static_param_path]
    #print(" ".join(cmd))
    p = subprocess.Popen(cmd, stdout=subprocess.PIPE)
    p.wait()
    stdout = list(p.stdout)[0].decode("utf-8")
    d = ast.literal_eval(stdout)
    return stdout, d