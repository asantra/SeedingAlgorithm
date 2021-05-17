#### this code mixes the background tracks and signal tracks per bunch crossing
#### run: python3 getSigBkgTextFile.py <bxNumberWanted>
import os
import sys
import time
import pprint
import math
from ROOT import *
from collections import OrderedDict
import argparse


def main():
    bxNumberWanted      = int(sys.argv[1])
    signalTracksNumber  = int(sys.argv[2])
    
    ### hics signal
    bkgFileName   = open("NewSamplesEBeamOnlyFilesMar62021/ePlusLaserBkgKaptonWindowNewSamplesMarch62021_DividedByBX"+str(bxNumberWanted)+"_trackInfoClean.txt")
    sigFileName   = open("../Outputfile/list_root_hics_165gev_w0_3000nm_jeti40_122020_9550dac4_trackInfoClean.txt")
    outFile       = open("BkgEBeam_SignalPositronhics3000nm_jeti40_122020_9550dac4_BX"+str(bxNumberWanted)+"_SignalTracks"+str(signalTracksNumber)+"_trackInfoClean.txt","w")
    
    
    ## outFile       = open("BkgEBeam_SignalHics5000nmProvisional_BX"+str(bxNumberWanted)+"_trackInfoClean.txt","w")
    ## outFile       = open("BkgEBeam_SignalHics8000nmProvisional_BX"+str(bxNumberWanted)+"_trackInfoClean.txt","w")
    ## outFile       = open("BkgEBeam_SignalHics3000nmOldReverse_BX"+str(bxNumberWanted)+"_trackInfoClean.txt","w")
    ## outFile       = open("BkgEBeam_SignalHics3000nmOldForSignalMultiplicityLessThan20Or5000nmForSignalMultiplicityMoreThan20_BX"+str(bxNumberWanted)+"_SignalTracks"+str(signalTracksNumber)+"_trackInfoClean.txt","w")
    
    
    #### bppp signal
    #bkgFileName   = open("NewSamplesGPlusLaserBkgFilesJan262021/gPlusLaserBkgNewSamplesJan262021_DividedByBX"+str(bxNumberWanted)+"_trackInfoClean.txt")
    #sigFileName   = open("list_root_hics_165gev_w0_3000nm_jeti40_122020_9550dac4_trackInfoClean.txt")
    #outFile       = open("BkgGBeam_SignalPositronbppp3000nmOr5000nm_BX"+str(bxNumberWanted)+"_SignalTracks"+str(signalTracksNumber)+"_trackInfoClean.txt","w")
    
    ## write the bkg as it is
    for lines in bkgFileName.readlines():
        outFile.write(lines)
    print("Printed bkg to mix")
    
    ### work on the signal, first collect all the signal lines in a dictionary
    position = []
    particleInEachBX = {}
    for lines2 in sigFileName.readlines():
        lines2   = lines2.rstrip()
        if "#" in lines2: continue
        eachWord = lines2.split()
        
        bxNumber  = int(eachWord[0])
        trackId   = int(eachWord[2])
        pdgId     = int(eachWord[1])
        staveId   = int(eachWord[3])
        
        position.append([bxNumber, pdgId, trackId, staveId])
        particleInEachBX.setdefault(bxNumber, []).append(lines2) ### this ensures that the new line for the same bxnumber does not erase the previous line, the new line will be appended as a list 
        
        
    ### count the number of signal positrons in each BX and then put BX: number of signals in a dictionary
    print("Now counting signals in each BX")
    outRootFile = TFile("signalTrackMultiplicityhics_w0_3000nm_jeti40_122020_9550dac4_BX"+str(bxNumberWanted)+".root", "RECREATE")
    outRootFile.cd()
    signalTrackMultiplicity = TH1F("signalTrackMultiplicity", "track multiplicity; number of tracks; BX",80, 0, 400)
    signalNumber     = {}
    for bx in range(0,1001+1):
        counter = 0
        if(bx%100==0): print("BX: ", bx)
        for tracks in position:
            if tracks[0] == bx:
                if(tracks[1] == -11 and tracks[2] == 1 and (tracks[3]==1000 or tracks[3]==1001)):
                    counter += 1
        signalNumber[bx] = counter ### store the number of signal for each BX
                
    ### sort the list according to the number of signals, descending order
    sortedSignalNumber = {k: v for k, v in sorted(signalNumber.items(), key=lambda item: item[1], reverse=True)}
    ### this is with ascending order
    #sortedSignalNumber = {k: v for k, v in sorted(signalNumber.items(), key=lambda item: item[1])}
    
    #print(sortedSignalNumber)
    ### write one bx of signal to the output file
    writeNumber = 0
    for key in sortedSignalNumber:
        if key==0:
            continue
            
        #print(key,sortedSignalNumber[key])
        
        #### select the number of tracks in each signal BX
        ##about 1 track, 5 tracks, 10 tracks, 20 tracks , 50 tracks, 100 tracks, 200 tracks, 500 tracks
        ### give a spread of 5 for each BX for low signal multiplicity
        if((signalTracksNumber - 200) < sortedSignalNumber[key] < (signalTracksNumber+400)):
            writeNumber += 1
        signalTrackMultiplicity.Fill(sortedSignalNumber[key])
        if(writeNumber!=bxNumberWanted): continue
        print("writing only the ",key," th signal BX, signalMultiplicity: ",sortedSignalNumber[key])
        #print(particleInEachBX[key])
        for eachLine in particleInEachBX[key]:
            #print("I am inside key loop: key is: ", key)
            findPDG     = int(eachLine.split()[1])
            findTrackId = int(eachLine.split()[2])
            #print("findPDG: ", findPDG, " findTrackId: ", findTrackId)
            if findPDG == -11 and findTrackId == 1:
                #print("I am writing")
                outFile.write(eachLine+"\n")
        break
            
    
    outFile.close()
    outRootFile.Write()
    outRootFile.Close()
                
            
    
if __name__=="__main__":
    start = time.time()
    main()
    print("--- The time taken: ", time.time() - start, " s")
