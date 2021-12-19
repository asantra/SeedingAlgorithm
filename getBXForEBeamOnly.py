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
    totalNumber = 0
    for lines in inFile.readlines():
        lines        = lines.rstrip()
        if('#' in lines): continue
        print("I am working on: ",lines)
        
        try:
            rootFile     = TFile(lines, "READ")
            dirAddress   = gDirectory.Get("hist")
            h0Value      = dirAddress.Get("h0")
            eventNumber  = h0Value.GetEntries()
            totalNumber += eventNumber
            #print("totalNumber of electrons: ",totalNumber)
            
            i = int(totalNumber) // 1500000000
            print("The division = ", i)
            outFile1 = open("gPlusLaser_list_root_55ae8938_DividedByBX"+str(i+1)+".txt", "a")
            outFile1.write(lines+"\n")
            
        except:
            print("something wrong here: ", lines)
        
        
    print("Total number of events: ",totalNumber)
        
            
    outFile1.close()
            
            
    
if __name__=="__main__":
    start = time.time()
    main()
    print("--- The time taken: ", time.time() - start, " s")
