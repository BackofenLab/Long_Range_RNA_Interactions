from collections import defaultdict
import os

def get_consensus_base(position_d):
    best_k = "-"
    best_v = 0
    for k, v in position_d.items():
        if v > best_v:
            best_k = k
            best_v = v
    return best_k, best_v

def get_consensus(aln_path, output):
    d = defaultdict(lambda: defaultdict(int))
    consensus = ""
    consensus_prob = ""
    valid_sum = 0
    maxlen = 0
    with open(output, "w") as of:
        with open(aln_path, "r") as f:
            for line in f.readlines():
                of.write(line)
                if line.startswith(("CLUSTAL", "\n", "#")):
                    continue
                maxlen = len(line)
                line = line.split(" ")[-1].rstrip("\n")
                for i in range(0, len(line)):
                    d[i][line[i]] += 1
                valid_sum += 1
        for i in range(0, len(d)):
            position_d = d[i]
            cons_b, cons_p = get_consensus_base(position_d)
            consensus += cons_b
            cons_v = str(cons_p/valid_sum).split(".")[1]
            #cons_v = str(cons_v)[-1]
            consensus_prob += cons_v[0]
        #for key, value in d.items():
        #    d[key] = dict(d[key])
        #print(dict(d))
        of.write(f"#Consensus{' '*(maxlen - len(d) - 11)}{consensus}\n")
        of.write(f"#Consensus_Prob{' '*(maxlen - len(d) - 16)}{consensus_prob}\n")
    
def get_all_locarna_consensus(locarna_results_path):
    for root, dirs, files in os.walk(locarna_results_path):
        for file in files:
            if file == "result.aln":
                print(root)
                output = f"{locarna_results_path}/{root.split('/')[-2]}_cons.txt"
                get_consensus(f"{root}/{file}", output)