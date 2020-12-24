### code to divide the background in separate BXs
### run: python getBXForEBeamOnly.py <list of background file names>
import os
import sys
import time
import pprint
import math
from ROOT import *
from collections import OrderedDict
import argparse


def main():
    inFile = open(sys.argv[1])
    outFile1 = open("EBeamOnlyWIS_DividedByBX1.txt", "w")
    outFile2 = open("EBeamOnlyWIS_DividedByBX2.txt", "w")
    outFile3 = open("EBeamOnlyWIS_DividedByBX3.txt", "w")
    outFile4 = open("EBeamOnlyWIS_DividedByBX4.txt", "w")
    outFile5 = open("EBeamOnlyWIS_DividedByBX5.txt", "w")
    totalNumber = 0
    for lines in inFile.readlines():
        lines = lines.rstrip()
        print("I am working on: ",lines)
        rootFile = TFile(lines, "READ")
        dirAddress = gDirectory.Get("hist")
        h0Value = dirAddress.Get("h0")
        eventNumber = h0Value.GetEntries()
        totalNumber += eventNumber
        print("totalNumber of electrons: ",totalNumber)
            
        if(totalNumber < 1500000000*0.826):outFile1.write(lines+"\n")
        elif(totalNumber > 1500000000*0.826 and totalNumber < 3000000000*0.826):outFile2.write(lines+"\n")
        elif(totalNumber > 3000000000*0.826 and totalNumber < 4500000000*0.826):outFile3.write(lines+"\n")
        elif(totalNumber > 4500000000*0.826 and totalNumber < 6000000000*0.826):outFile4.write(lines+"\n")
        else: outFile5.write(lines+"\n")
            
        
    outFile1.close()
    outFile2.close()
    outFile3.close()
    outFile4.close()
    outFile5.close()
    
            
            
    
if __name__=="__main__":
    start = time.time()
    main()
    print("--- The time taken: ", time.time() - start, " s")
