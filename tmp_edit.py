
import subprocess
shift = 12
c = 0
boxshift = 68 # Oh god
boxshift_strength = 6 # Oh god
with open("aln_edited_2.ps", "w") as f:
    with open("aln_edited.ps", "r") as f2:
        for line in f2.readlines():
            c += 1
            if c >= 41 and c <= 1056: # Coloured boxes...
                sline = line.split(" ")
                sline[1] = str(float(sline[1]) + shift)
                sline[3] = str(float(sline[3]) + shift)
                if float(sline[0]) > 1000:
                    sline[0] = str(float(sline[0]) - boxshift*boxshift_strength)
                    sline[2] = str(float(sline[2]) - boxshift*boxshift_strength)
                f.write(" ".join(sline))
            elif c == 1058 or c == 1059: # CM only line
                sline = line.split(" ")
                sline[-2] = str(float(sline[-2]) + 1)
                f.write(" ".join(sline))
            elif c == 1060 or c == 1061: # Full alignment consensus line
                sline = line.split(" ")
                sline[-2] = str(float(sline[-2]) + 1)
                f.write(" ".join(sline))
            elif c >= 1062 and c <= 1137: # Strings
                sline = line.split(" ")
                sline[-2] = str(float(sline[-2]) + shift)
                if line.startswith("(258)"):
                    sline[1] = str(float(sline[1]) - boxshift*boxshift_strength)
                f.write(" ".join(sline))
            elif c >= 1139 and c <= 1245: # Grey stuff 
                sline = line.split(" ")
                sline[1] = str(float(sline[1]) + shift)
                sline[3] = str(float(sline[3]) + shift)
                f.write(" ".join(sline))
            elif c >= 1247 and c <= 1366: # Grey stuff part 2
                sline = line.split(" ")
                sline[1] = str(float(sline[1]) + shift)
                sline[3] = str(float(sline[3]) + shift)
                sline[0] = str(float(sline[0]) - 408)
                sline[2] = str(float(sline[2]) - 408)
                f.write(" ".join(sline))
            else: # Rest lines
                f.write(line)
                
p = subprocess.Popen(["ps2pdf", "-dEPSCrop", "aln_edited_2.ps"], stdout=subprocess.PIPE)
p.wait()