import os
import sys
import time
import pprint
import math
from ROOT import *
import array
from makeTrackDiagrams import *
from collections import OrderedDict 


#### Z position of staves
z1inner = GetLayerZ(1000,0)
z2inner = GetLayerZ(1000,2)
z3inner = GetLayerZ(1000,4)
z4inner = GetLayerZ(1000,6)
z1outer = GetLayerZ(1000,1)
z2outer = GetLayerZ(1000,3)
z3outer = GetLayerZ(1000,5)
z4outer = GetLayerZ(1000,7)



# get the layer from the stave number and its z position in mm
def getStaveZ(stave):
    zPos = -99999.0
    if(stave==0):
        zPos = 3864.5125
    elif(stave==1):
        zPos = 3876.5125
    elif(stave==2):
        zPos = 3964.5125
    elif(stave==3):
        zPos = 3976.5125
    elif(stave==4):
        zPos = 4064.5125
    elif(stave==5):
        zPos = 4076.5125
    elif(stave==6):
        zPos = 4164.5125
    elif(stave==7):
        zPos = 4176.5125
    else:
        zPos = -99999.0
    return zPos


def main():
    # give the input text file containing all the track information
    inTextFile     = sys.argv[1]
    inputTrackInfo = open(inTextFile)
    ### open histogram to know the seed information
    plotSuffixName      = ""
    if (("_" in inTextFile) and ("WIS" in inTextFile)):
        eachName            = inTextFile.split('.')[0].split('_')
        suffixName          = "_".join(eachName[2:])
    else:
        suffixName          = inTextFile.split('.')[0]
    
    outFile           = TFile("seedingInformationSmallScript_"+suffixName+".root", "RECREATE")
    outFile.cd()
    hAllPossible      = TH1D("hAllPossible", "all possible track combination; bunch crossing; number of track combination", 9508, 0, 9508)
    hSeedPossible     = TH1D("hSeedPossible", "seed track combination; bunch crossing; number of seed track", 9508, 0, 9508)
    hSeedMultiplicity = TH1D("hSeedMultiplicity", "hSeedMultiplicity", 50, 0, 50)
    hSigEnergy        = TH1D("hSigEnergy", "hSigEnergy", 200, 0, 20)
    
    # all the track info is in the following list
    position = []
    # get the information from the text files
    for lines in inputTrackInfo.readlines():
        lines       = lines.rstrip()
        eachWord    = lines.split()
        bxNumber    = int(eachWord[0])
        trackId     = int(eachWord[2])
        pdgId       = int(eachWord[1])
        trackEnergy = float(eachWord[6])
        if(pdgId!=-11): continue
        ### select if only background or signal tracks wanted
        if(trackId!=1): continue
        position.append([bxNumber, trackId, int(eachWord[3])-1000, float(eachWord[4]), float(eachWord[5]), float(eachWord[6]), float(eachWord[7])])

    
    for bxCounter in range(1,9509):
        # separate each bx now
        eachBXValue = []
        for tracks in position:
            ### the below is needed for e+laser hics setup
            if tracks[0] == bxCounter:
                eachBXValue.append(tracks)

        
        ### fill up the x,y, z and E values from each of the tracker layers
        allR1Inner = []; 
        allR1Outer = [];
        
        for values in eachBXValue:
            zPosition = getStaveZ(values[2])
            ### x, y, z and E
            if (values[2] == 0):
                allR1Inner.append([values[3], values[4], zPosition, values[5], values[6]])
            elif (values[2] == 1):
                allR1Outer.append([values[3], values[4], zPosition, values[5], values[6]])
            else:
                print("stave not needed")
        
        ### removing the overlap region of inner and outer stave
        allR1Unique = allR1Inner
        
        for r1Out in allR1Outer:
            #### remove all points having an x overlap with inner stave: not 100% effective
            if r1Out[0] > (308.53 + 29.94176/2.):  ### the x position of last chip on inner stave layer 1+half of the x size 
                allR1Unique.append(r1Out)
            
        for r1 in allR1Unique:
            hSigEnergy.Fill(r1[3], r1[4])
        
    outFile.Write()
    outFile.Close()
    

if __name__ == "__main__":
    start = time.time()
    main()
    print("-------- The processing time: ",time.time() - start, " s")
