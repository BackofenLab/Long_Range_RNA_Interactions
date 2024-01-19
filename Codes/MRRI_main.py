import os
import sys
import argparse
import json
try:
    from Codes.MRRI import MRRI
except:
    from MRRI import MRRI
import re, csv
import numpy as np
from pathlib import Path

# ids of IntaRNA CSV output columns for which no :-separated value output is needed
uniqueValCols = ["id1","id2","RT"]



def main():
    '''
        Set defaults to arguments and get arguements from command line
    '''
    parser = argparse.ArgumentParser(
        description=" Multi Site RNA-RNA Interaction Prediction")
    parser.add_argument('-q', '--query', dest='query',
                        type=str, help='Query sequence', required=True)
    parser.add_argument('-t', '--target', dest='target',
                        type=str, help='Target sequence', required=True)
    parser.add_argument('-p', '--parameterFile', dest='parameterFile', default="",
                        help='Optional parameter file for IntaRNA provide further IntaRNA parameters and prediction constraints.')
    parser.add_argument( '--intarnaBin', type=Path, action="store", dest='intarnaBin', default=os.path.abspath(os.path.join(os.curdir, '..')),
                        help='Parameter that provides the IntaRNA binary (with path) to be called.')
    
    args = parser.parse_args()
    MRRIHandler = createMRRI(args)
    return MRRIHandler


def createMRRI(args):
    '''
        Create an Object of MRRI Class with the given Params
    '''
    args.intarnaBin = "IntaRNA" ####
    IntaRNA = {"exists": True}  #### findBinary(args.intarnaBin) ####
    if IntaRNA.get("exists") == True:
        MRRI_instance = MRRI(args)
        return MRRI_instance
    else:
        sys.exit(-1)

def findBinary(bin_name):
    '''
        Check if the Binary Exists. Return false if the Binary is absent
    '''
    exists = os.path.exists(bin_name) and os.path.isfile(bin_name)
       
    if not exists:
        sys.stderr.write("Cannot find IntaRNA via '"+str(bin_name)+"'.\n")
        sys.exit(-1)

    return {"exists": exists, "path":bin_name}

   
def printCSVRow(B1,B2, MRRIHandler):
    d1 = B1
    d2 = B2
    if (d2 == None):
        # use keys from d1 but set all values to empty string
        d2 = dict.fromkeys(d1,'')
    else: # correct B1 energies
        ED = MRRIHandler.getEDunconstraint(d1)
        d1['ED1'] = str(ED[0])
        d1['ED2'] = str(ED[1])
        d1['E'] = str(round(float(d1['E_hybrid']) + float(ED[0]) + float(ED[1]),2))

    line = []
    for key in sorted(d1.keys()):
        if key == 'E':
            E = round(float(d1['E']) + (0 if d2[key]=='' else float(d2['E'])),2)
            line.append(str(E))
        else:
            if uniqueValCols.__contains__(key) or str(d2[key]) == '':
                line.append(d1[key])
            else:
                line.append(str(d1[key])+':'+str(d2[key]))        
    ##print((';\t').join(line))
    return E ########################################

def csvHeader(B1):
    headerNames = []
    for element in sorted(B1.keys()):
        headerNames.append(element)
 #       if (element == 'E'):
  #          headerNames.append('Esite')
    print((';\t').join(headerNames))
       
def sameBlock( Block1, Block2 ):
    #print(B3)
    #print(B2)
    #B2['start1'] == B2['end1'] == B2['start2'] == B2['end2'] == B3['start1'] == B3['start2'] == B3['end1'] == B3['end2']
    result = all( Block1[k]==Block2[k] for k in ('start1', 'end1', 'start2', 'end2'))
    print(result)


if __name__ == '__main__':
    MRRIHandler = main()
    for qId in MRRIHandler.querySeq.keys():
        for tId in MRRIHandler.targetSeq.keys():
            B2 = dict({ 'id1': tId, 'id2' : qId})
            B1 = MRRIHandler.runIntaRNA(B2)
            B2 = None
            B3 = None
            ##print(B1) ### v1: Only 1 Round (like normal IntaRNA
            ##break
            #csvHeader(B1) ### These lines are only printing stuff and are useless
            iteration = 1
            E2 = 999 ###
            #E = 0    ###
            while True:
                #if iteration % 2 == 1:
                #    E = printCSVRow(B1, B2, MRRIHandler)
                #else:
                #    E = printCSVRow(B2, B1, MRRIHandler)
                B3 = MRRIHandler.runIntaRNA(B1)
                ####E = round(float(B1['E']) + float(B2['E']),2)
                #if(B2 != None and len(B1["hybridDP"]) >= len(B3["hybridDP"])):
                if(float(B1["E_hybrid"]) < float(B3["E_hybrid"])): 
                #if(float(B1["E"]) < float(B3["E"])):  ### v2: Directly compare energies
                #if E2 <= E: ### v3: Compare new energies # E2: Energy of B1, E: energy of B1
                #if( B2 == B3): ### v0: Currently doesnt work?
                    #sameBlock(B2, B3)
                    print(B1)
                    break
                else:
                    B2 = B1
                    B1 = B3
                    #E2 = E
                    #print(B1)
                iteration +=  1

                




    


   
