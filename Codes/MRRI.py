import os
import sys
import subprocess as s
from subprocess import Popen
import argparse
from typing import List, Tuple
import pandas as pd
import re
from itertools import groupby


class MRRI():

    def __init__(self, args):
        '''
            Constructor for the MRRI class with default value from args
        '''
        self.regex = "[ACGTUacgtu]+"
        self.args = args
        self.querySeq = self.parseFasta2dict( args.query, 'query' ) # dictionary of queryID 2 query sequence
        self.targetSeq = self.parseFasta2dict( args.target, 'target' ) # dictionary of targetID 2 target sequence
        self.b1 = None
        self.b2 = None
        self.b3 = None

    
    def get_param_mode(self, B1, param_mode=1):
        UTR5pCDS = self.targetSeq[B1['id1']]
        UTR3pCDS = self.querySeq[B1['id2']]
        tidxpos0 = -(len(UTR5pCDS)-200) # len of UTR5 alone
        qidxpos0 = -200 # = -200
        if param_mode == 1:
            tregion_s = tidxpos0 # = len(UTR5)
            tregion_e = 100 # = 100
            qidxpos0 = -200 # = -200
            qregion_s = -100 # = -100
            qregion_e = len(UTR3pCDS)-200 # len of UTR3 alone
        if param_mode == 2:
            tregion_s = -40 # = UTR5/CDS transition - 40
            tregion_e = +70 # = UTR5/CDS transition +70 40
            qregion_s = len(UTR3pCDS)-200-140
            qregion_e = len(UTR3pCDS)-200-50 # len of UTR3 alone
        tregion = f"{str(tregion_s)}-{str(tregion_e)}"
        qregion = f"{str(qregion_s)}-{str(qregion_e)}"
        return tidxpos0, qidxpos0, tregion, qregion

    def runIntaRNA(self, B1= None, param_mode=1):
        '''
            Run IntaRNA on given query and target
            params have value of start and end numbers
        
        '''
        complete = str(self.args.intarnaBin) + ' -q ' + self.querySeq[B1['id2']] + ' -t ' + self.targetSeq[B1['id1']] +' --outMode=C -n 1 '#--energyNoDangles'
        # add parameterFile to call if given
        if self.args.parameterFile :
            complete += " --parameterFile="+self.args.parameterFile
        # set require CSV columns
        complete += " --outCsvCols=id1,start1,end1,id2,start2,end2,subseqDP,hybridDP,E,E_hybrid,ED1,ED2"
        ########################################
        ############ Hacked in #################
        ########################################
        tidxpos0, qidxpos0, tregion, qregion = self.get_param_mode(B1, param_mode)
        complete += f" --tidxpos0 {str(tidxpos0)} " ## Transition UTR5-CDS
        complete += f"--qidxpos0 {str(qidxpos0)} "  ## Transition CDS-UTR3
        complete += f"--tregion {tregion} " ## UTR5start to end+100CDS
        complete += f"--qregion {qregion} " ## 100CDS to UTR3end ####
        ########################################
        ########################################
        ########################################
        # TODO (somewhen) parse parameterFile for outCsvCols and add user-requested csv-col ids not already within the list
        result = {"tAccConstr":"", "qAccConstr":""}
        if 'start1' in B1:
            result["tAccConstr"] += f"b:{B1['start1']}-{B1['end1']}"
            result["qAccConstr"] += f"b:{B1['start2']}-{B1['end2']}"
            if "tAccConstr" in B1 and B1["tAccConstr"]:
                result["tAccConstr"] += f",{B1['tAccConstr']}"
                result["qAccConstr"] += f",{B1['qAccConstr']}"
            complete += ' --tAccConstr="'+result["tAccConstr"]+'" --qAccConstr="'+result["qAccConstr"]+'" ' ## Add \xa0 before --tAccConst?
        #print("".join(complete))
        result.update(self.csv2dict(self.runCmdLine(complete)[0].replace("query",B1['id2']).replace("target",B1['id1'])))
        ### This might be neccessary? This actually works but something something utf-8.......
        #try:
        #    result.update(self.csv2dict(self.runCmdLine(complete)[0].replace("query",B1['id2']).replace("target",B1['id1'])))
        #except:
        #    pass
        #if "start1" in result:
        #    s1 = int(result["start1"])
        #    e1 = int(result["end1"])
        #    if s1 >= 0:
        #        #print(str(s1) + " ####################")
        #        result["start1"] = str(s1-1)
        #    if e1 >= 0:
        #        result["end1"] = str(e1-1)
        #print(result)
        return result


##################################
##### Currently unused stuff #####
##################################       


    def csv2dict( self, intarnaCsvOut ):
        if (intarnaCsvOut == None):
            return {}
        keys = intarnaCsvOut.split("\n")[0].split(";")
        values = intarnaCsvOut.split("\n")[1].split(";")
        op = zip(keys, values)
        return dict(op)

    def read_fasta_file(self, fastaname):
        """
        given a fasta file. yield tuples of id, sequence
        """
        fh = open(fastaname)
        fastaiter = (x[1] for x in groupby(fh, lambda line: line[0] == ">"))
        for idx in fastaiter:
            # drop the ">"
            idx = idx.next()[1:].strip()
            # join all sequence lines to one.
            seq = "".join(s.strip() for s in fastaiter.next())
            yield idx, seq

      

    def parseFasta2dict(self, sequenceInput, seqname ):
       
        if os.path.isfile(sequenceInput) == True:
            fastaReturn = self.read_fasta_file(sequenceInput)
            return fastaReturn
        else:
            matches = re.finditer(self.regex, sequenceInput, re.MULTILINE | re.IGNORECASE)
            tempDict = {}
            tempDict[seqname] = ''
            for match in matches:
                if match.group().strip("") != "":
                    tempDict[seqname] += str(match.group())
            return tempDict

    def getEDunconstraint(self, B1 ):
        complete = str(self.args.intarnaBin) + ' -q ' + self.querySeq[B1['id2']] + ' -t ' + self.targetSeq[B1['id1']] +' --out=/dev/null -n 0 '#--energyNoDangles'
        # add parameterFile to call if given
        if self.args.parameterFile :
            complete += " --parameterFile="+self.args.parameterFile
        # set require ED output
        complete += " --out=tAcc:STDOUT --out=qAcc:STDERR"
        # confine output to region of interest
        l1 = int(B1['end1'])-int(B1['start1'])+1
        l2 = int(B1['end2'])-int(B1['start2'])+1
        complete += " --tIntLenMax=" + str(l1) + " --qIntLenMax=" +str(l2)
        # run intarna
        outputED = self.runCmdLine( complete )
        # parse ED1 from stdout == outputED[0]
        # TODO CRASHES IF (tIdxPos0 != 1)
        ED1 = round(float(outputED[0].split('\n')[ int(B1['end1'])+1 ].split()[l1]),2)
        # parse ED2 from stdout == outputED[1]
        ED2 = round(float(outputED[1].split('\n')[ int(B1['end2'])+1 ].split()[l2]),2)
        # return combined constraints
        return [ED1,ED2]

    def runCmdLine(self,completeCall):
        #print(completeCall)
        ps = s.Popen(str(completeCall),  stdout = s.PIPE, stderr=s.PIPE, universal_newlines=True, shell=True)
        (stdout, stderr) = ps.communicate()
        if ps.returncode:  # If IntaRNA exits with a returncode != 0, skip this iteration
            sys.stderr.write("IntaRNA failed ({}) for call {}\n".format(stderr,completeCall))
            #exit(ps.returncode)
            return None
        return [stdout,stderr]
