#### this code adds the background text files
#### run: python3 stitchEBeamFiles.py <enter the options>
import os
import sys
import time
import pprint
import math
from ROOT import *
from collections import OrderedDict
import argparse


def main():
    
    
    parser = argparse.ArgumentParser(description='Code to get all BX track info')
    parser.add_argument('-out', action="store", dest="outFile", type=str, default="ePlusLaserBkgNewSamplesJan262021_AllBX_trackInfoClean.txt")
    parser.add_argument('-in', action="store", dest="inDir", type=str, default="NewSamplesEBeamOnlyFilesJan262021")
    parser.add_argument('-bx', action="store", dest="nbx", type=int, default=1.0)
    parser.add_argument('-ident', action="store", dest="identifier", type=str, default="ePlusLaserBkgNewSamplesJan262021")
    args = parser.parse_args()
    
    
    outFile       = open(args.outFile, "w")
    outFile.write("### bxNumber << pdg << track_id << det_id << xx << yy << eneg << ev_weight << vtx_x << vtx_y << vtx_z\n")
    inputDir      = args.inDir
    
    identifierFileName = args.identifier
    
    bxCounter = 0
    for bx in xrange(1, args.nbx+1):
        if(bx%20==0):print("bx processed: ", bx)
        bkgFileName   = open(inputDir+"/"+identifierFileName+"_DividedByBX"+str(bx)+"_trackInfoClean.txt")
        bxCounter    += 1
        ### write the bkg as it is
        for lines in bkgFileName.readlines():
            if '#' in lines:
                continue
            lines        = lines.rstrip()
            eachWord     = lines.split()
            modifiedLine = [str(bxCounter)] + eachWord[1:]
            writeLine    = ' '.join([elem for elem in modifiedLine])
            outFile.write(writeLine+"\n")
        
    outFile.close()
    
if __name__=="__main__":
    start = time.time()
    main()
    print("--- The time taken: ", time.time() - start, " s")
