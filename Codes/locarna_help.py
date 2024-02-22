import os
import sys
import subprocess
import ast
#from Codes.MRRI import MRRIHandler
import Codes.MRRI_main

def run_command(cmd):
    p = subprocess.run(cmd.split(" "))


def find_all(s, c):
    """Generator that finds and returns the locations 
    of all characters c in a string s.
    """
    idx = s.find(c)
    while idx != -1:
        yield idx
        idx = s.find(c, idx + 1)


def integrate_into(s1, s2, l):
    """
    s1 = "..AA..AA"
    s2 = "....CC.."
    l  = "C"
    res= "..AACCAA"
    """
    pos = find_all(s2, l)
    for i in pos:
        s1 = s1[:i] + l + s1[i+1:]
    return s1


def add_cm_hit_to_fs(fs, cmhit, cmstart):
    """
    fs    = "AAAAAA.NNN..a.aa.aaa......."
    cmhit = "<<-<->->>"
    cmstart= 8
    # Insert into our fs string
    result= "AAAAAA.NNN..a.aa.a<<-<->->>"
    # Replace replaced letters on the opposing site with "."
    result= "..AAAA.NNN..a.aa.a<<-<->->>"
    ..Should be good enough to replace the first n occurences of that letter?
    """
    from collection import Counter
    counts = Counter(fs[cmhit:]) # Count occurences of each letter thats being replaced
    for letter, amount in counts:
        if letter == ".":
            continue
        fs.replace(letter.upper(), ".", amount)
    return fs

def cm_compare(id, seq, cm_file):
    """- Open alignment file where the CMhit is from
    - find right sequence
    - read sequence from input and compare with cmhit
    - if match, insert respective part from cmhit
    - if not, shift cmhit sequence to the next base
    """
    #print(id)
    location = None
    result_string = ""
    name_found = False
    with open(cm_file, "r") as f:
        for line in f.readlines():
            if (not name_found) and line.startswith("#=GS") and line.split()[3] == id:
                name, location = line.split()[1].split("/")
                name_found = True
                #print(name, location)
            elif name_found and line.startswith(name):
                cm_seq = line.split()[1]
                cm_seq = cm_seq.upper().replace("U", "T")
            elif line.startswith("#=GC SS_cons"):
                align_cons = line.split()[-1]
            #elif id == "NC_018705.3":
            #    print(line)
    #print("AAAAAAAAAA")
    if not location:
        print(f"No CM alignment found for {id}")
        return None
    else:
        location_start, location_end = location.split("-")
        cut_seq = seq[int(location_start)+199:]
        cut_seq_pos = 0
        for i in range(len(cm_seq)):
            if cm_seq[i] == cut_seq[cut_seq_pos]:
                result_string += align_cons[i]
                cut_seq_pos += 1
        return result_string


def run_mlocarna(input_fasta, output_dir, use_carna=False):
    """Run mlocarna on a given fasta file.
    
    input_fasta (str): Filepath to a fasta file to apply locarna on
    output_dir (str): Filepath of the resulting CM file
    """
    
    if use_carna:
        print(f"mlocarna(+carna): {input_fasta} `=> {output_dir}")
        ## Doesnt work: :(
        conda_base = subprocess.Popen(["conda", "info", "--base"], stdout=subprocess.PIPE)
        conda_base.wait()
        carna_loc = list(conda_base.stdout)[0].decode("utf-8").rstrip("\n")
        carna_loc += "/envs/carna/bin/carna"
        #carna_loc = "/home/arkanini/miniconda3/envs/carna/bin/carna"
        cmd = ["conda", "run", "-n", "carna"]
    else:
        print(f"mlocarna: {input_fasta} `=> {output_dir}")
        cmd = ["conda", "run", "-n", "locarna"]
    cmd += ["mlocarna", input_fasta,
           #"--indel=-50", # Webserver parameter
           #"--indel-opening=-750", # Webserver parameter
           "--width=300",
           "--use-ribosum=true",
           "--tgtdir", output_dir
           ]
    if use_carna:
        cmd += [f"--pw-aligner={carna_loc}"]
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

def hacked_MRRI_main(UTR5pCDS, UTR3pCDS, static_param_path, param_mode):
    sys.argv = ["Codes/MRRI_main.py", "-q", UTR3pCDS, "-t", UTR5pCDS, "-p", static_param_path]
    MRRIHandler = Codes.MRRI_main.main()
    for qId in MRRIHandler.querySeq.keys():
        for tId in MRRIHandler.targetSeq.keys():
            BE = dict({ 'id1': tId, 'id2' : qId})
            B1 = MRRIHandler.runIntaRNA(BE, param_mode)
            #print(B1)
            if "id2" in B1:
                B2 = MRRIHandler.runIntaRNA(B1, param_mode)
            else:
                BErr = {"start1":0, "end1":0, "start2":0, "end2":0, "hybridDP":0}
                return [BErr] ## First IntaRNA round didnt find anything
            #print(B2)
            if "id2" in B2:
                B3 = MRRIHandler.runIntaRNA(B2, param_mode)
            else:
                return [B1]
            #print(B3)
            if "id2" in B3:
                return [B1, B2, B3]
            else:
                return [B1, B2]


def run_mrri(UTR5pCDS, UTR3pCDS, static_param_path):
    """Outdated. Currently not in use"""
    cmd = ["python", "Codes/MRRI_main.py", "-q", UTR3pCDS, "-t", UTR5pCDS, "-p", static_param_path]
    #print(" ".join(cmd))
    p = subprocess.Popen(cmd, stdout=subprocess.PIPE)
    p.wait()
    stdout = list(p.stdout)[0].decode("utf-8")
    d = ast.literal_eval(stdout)
    return stdout, d