from Codes.MRRI_main import main
import sys


def hacked_MRRI_main(UTR5pCDS, UTR3pCDS, static_param_path, param_mode):
    """
    param_mode (int): Decides region of interest for MRRI
                      1 - Whole sequence from 5' to 100 into CDS and last 100 to 3' end
                      2 - Limited to small area around 5'UTR/CDS transition and short area in 3' UTR
    """
    sys.argv = ["Codes/MRRI_main.py", "-q", UTR3pCDS, "-t", UTR5pCDS, "-p", static_param_path]
    MRRIHandler = main()
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
