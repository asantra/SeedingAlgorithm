########################################################################################################################
##### code to make seeds from list of tracks                                                                       #####
##### The track list must have this order: <BX> <pdgID> <trackId> <stave> <x> <y> <E> <weight>                     ##### 
##### Run this code: python findSeed.py -l <trackList.txt> -s <bool, true for signal positron> -e <energyCutinKeV> #####
##### written by: arka.santra@cern.ch                                                                              #####
########################################################################################################################

import os
import sys
import time
import pprint
import math
from ROOT import *
from collections import OrderedDict
import argparse
### needed for plotting of tracks and dimension from the spreadsheet 
from makeTrackDiagrams import *

### the seed energy width
EseedMin = 1.0 # GeV
EseedMax = 16.5 # GeV

### length of the dipole in meters
LB       = 1.029
#//MeV mass of electron/positron
meMeV    = 0.5109989461 
meGeV    = meMeV/1000.
### me^2 in GeV
meGeV2   = meGeV*meGeV
### 1 Tesla magnetic field used for now
B        = 1.0 
### use only positive x, that's the positron side
side     = "Positron" 


### The seed lorentz vector, need for energy plotting
p       = TLorentzVector()
#### This is the cutflow dictionary
cutFlowDict = OrderedDict(
    [('noCut', 0), 
     ('x1Gtx4',0), 
     ('x1*x4Negative',0), 
     ('z1Eqz4', 0), 
     ('yDipoleExitGt10',0), 
     ('xDipoleExitLt5',0), 
     ('xDipoleExitWithinDipole', 0), 
     ('xDipoleExitLt0',0), 
     ('checkClusterTracksMiddleLayers', 0), 
     ('trackEnergy', 0)]
    )

# get the unit vector along one vector
def rUnit2(r1, r2):
    r = (r2-r1).Unit()
    return r
    

### check if there are hits along one track in layer 2 and layer 3
def check_clusters(r1, r4, allR2Inner, allR2Outer, allR3Inner, allR3Outer):
    yAbsMargins = 0.2 # mm (a "road" of 200 microns around the line between r4 and r1)
    xAbsMargins = 0.2 # mm (a "road" of 200 microns around the line between r4 and r1)
    
    r1min       = [r1[0]-xAbsMargins, r1[1]-yAbsMargins, r1[2]]
    r1max       = [r1[0]+xAbsMargins, r1[1]+yAbsMargins, r1[2]]
    r4min       = [r4[0]-xAbsMargins, r4[1]-yAbsMargins, r4[2]]
    r4max       = [r4[0]+xAbsMargins, r4[1]+yAbsMargins, r4[2]]


    # check possible clusters in layer 2, for both inner and outer stave
    y2minInner = yofz(r1min, r4min, z2inner)
    y2maxInner = yofz(r1max, r4max, z2inner)
    x2minInner = xofz(r1min, r4min, z2inner)
    x2maxInner = xofz(r1max, r4max, z2inner)
    
    y2minOuter = yofz(r1min, r4min, z2outer)
    y2maxOuter = yofz(r1max, r4max, z2outer)
    x2minOuter = xofz(r1min, r4min, z2outer)
    x2maxOuter = xofz(r1max, r4max, z2outer)
    
    ### separate out staves for allR2, work on the inner stave of layer 2
    accept2Inner = False
    accept2Outer = False
    for i2 in range(0, len(allR2Inner)):
        ### remember allR2 has x, y, z and E saved in the list
        accept2yzInner = ( (allR2Inner[i2][1] >= y2minInner and allR2Inner[i2][1] <= y2maxInner) )
        if(not accept2yzInner): continue
    
        accept2xzInner = ( ( allR2Inner[i2][0] >= x2minInner and allR2Inner[i2][0] <= x2maxInner ) )
        if(not accept2xzInner): continue
        accept2Inner = True
        break
    
    ### separate out staves for allR2
    ### only if there is no suitable track from layer 2 inner stave, go to layer 2 outer stave
    if(not accept2Inner):
        for i2 in range(0, len(allR2Outer)):
            ### remember allR2 has x, y, z and E saved in the list
            accept2yzOuter = ( (allR2Outer[i2][1] >= y2minOuter and allR2Outer[i2][1] <= y2maxOuter) )
            if(not accept2yzOuter): continue
        
            accept2xzOuter = ( ( allR2Outer[i2][0] >= x2minOuter and allR2Outer[i2][0] <= x2maxOuter ) )
            if(not accept2xzOuter): continue
            accept2Outer = True
            break
    
    #### return false if both of the inner and outer acceptances are false
    if(not (accept2Inner or accept2Outer)):
        return False
    
    
    
    #/// check possible clusters in layer 3, for both inner and outer stave
    y3minInner = yofz(r1min, r4min, z3inner) ### z for inner stave layer 3 is 4064.5125
    y3maxInner = yofz(r1max, r4max, z3inner)
    x3minInner = xofz(r1min, r4min, z3inner)
    x3maxInner = xofz(r1max, r4max, z3inner)
    
    y3minOuter = yofz(r1min, r4min, z3outer) ### z for Outer stave layer 3 is 4076.5125
    y3maxOuter = yofz(r1max, r4max, z3outer)
    x3minOuter = xofz(r1min, r4min, z3outer)
    x3maxOuter = xofz(r1max, r4max, z3outer)
    
    
    ### separate out staves for allR3, work on the inner stave of layer 3
    accept3Inner = False
    accept3Outer = False
    for i3 in range(0, len(allR3Inner)):
        ### remember allR3 has x, y, z and E saved in the list
        accept3yzInner = ( (allR3Inner[i3][1] >= y3minInner and allR3Inner[i3][1] <= y3maxInner) )
        if(not accept3yzInner): continue
    
        accept3xzInner = ( ( allR3Inner[i3][0] >= x3minInner and allR3Inner[i3][0] <= x3maxInner ) )
        if(not accept3xzInner): continue
        accept3Inner = True
        break
    
    ### separate out staves for allR3
    ### only if no track found in layer 3 inner stave, go to layer 3 outer stave
    if(not accept3Inner):
        for i3 in range(0, len(allR3Outer)):
            ### remember allR3 has x, y, z and E saved in the list
            accept3yzOuter = ( (allR3Outer[i3][1] >= y3minOuter and allR3Outer[i3][1] <= y3maxOuter) )
            if(not accept3yzOuter): continue
        
            accept3xzOuter = ( ( allR3Outer[i3][0] >= x3minOuter and allR3Outer[i3][0] <= x3maxOuter ) )
            if(not accept3xzOuter): continue
            accept3Outer = True
            break
    
    #### return false if both of the inner and outer acceptances are false
    if(not (accept3Inner or accept3Outer)):
        return False
    
    return True


### making the seeds from two tracks from innermost layer and outermost layer
def makeseed(r1, r4, allR2Inner, allR2Outer, allR3Inner, allR3Outer):
    cutFlowDict['noCut'] += 1
    
    if(abs(r1[0]) >= abs(r4[0])):
        return False  ### |x1| must be smaller than |x4|
    cutFlowDict['x1Gtx4'] += 1
    
    if(r1[0]*r4[0] < 0):
        return False  ### the x value should have the same sign,  i.e. on the same side of the beam
    cutFlowDict['x1*x4Negative'] += 1
    
    if(r1[2] == r4[2]):
        return False  ### if z1=z4..., this is also impossible
    cutFlowDict['z1Eqz4'] += 1
    
    yDipoleExit = yofz(r1, r4, zDipoleExit)
    xDipoleExit = xofz(r1, r4, zDipoleExit)
    ### The following cuts are coming because of the distribution of the signal
    ### the following cannot be 20mm/2 according to the signal tracks
    if(abs(yDipoleExit) > 10.0):
        return False
    cutFlowDict['yDipoleExitGt10'] += 1
    
    ### This is according to the signal x:y at the dipole exit
    if(abs(xDipoleExit) < 5.0):
        return False  # the track should point to |x|<~1.0 at the dipole exit
    cutFlowDict['xDipoleExitLt5'] += 1
    
    if(abs(xDipoleExit) > xDipoleWidth/2):
        return False
    cutFlowDict['xDipoleExitWithinDipole'] += 1

    ### select only positive x or negative x accroding to the particle
    if( (side=="Positron" and xDipoleExit < 0) or (side=="Electron" and xDipoleExit > 0)):
        return False
    cutFlowDict['xDipoleExitLt0'] += 1
    
    if(not check_clusters(r1, r4, allR2Inner, allR2Outer, allR3Inner, allR3Outer)):
        return False #// minimum one cluster at layer 2 and one at layer 3
    cutFlowDict['checkClusterTracksMiddleLayers'] += 1
    
    #### find the energy of the tracks
    x0        = 0
    z0        = zofx(r1, r4, x0)
    xExit     = abs(xofz(r1, r4, zDipoleExit))
    H         = abs((zDipoleExit-z0))/1000.0 ### converting H from mm to m
    xExitInM  = xExit/1000.0                 ### converting xExit from mm to m
    R         = H*(LB)/xExitInM + xExitInM   ### // This is the radius of curvature for a track
    P         = 0.3*B*R                      ### here B in Tesla and R in m
    v1        = TVector2(r1[2], r1[1])
    v4        = TVector2(r4[2], r4[1])
    u         = rUnit2(v1, v4)
    uz        = u.X()
    uy        = u.Y()
    px        = 0    
    py        = P*uy
    pz        = P*uz
    p.SetPxPyPzE(px, py, pz, math.sqrt(px*px + py*py + pz*pz + meGeV2))
    
    ### checking the track energy
    if(p.E() < EseedMin or p.E() > EseedMax): 
        return False
    cutFlowDict['trackEnergy'] += 1

    return True





def main():
    # give the input text file containing all the track information
    parser = argparse.ArgumentParser(description='Code to find seed tracks')
    parser.add_argument('-l', action="store", dest="inFile", default="list_root_hics_165gev_w0_3000nm_WIS_trackInfoClean.txt")
    parser.add_argument('-s', action="store_true", dest="needSignal", default=False)
    parser.add_argument('-e', action="store", dest="eCut", type=float, default=0.0)
    args = parser.parse_args()
    
    ### open the file containing track file
    inTextFile      = args.inFile
    inputTrackInfo  = open(inTextFile)
    energyCut       = args.eCut
    energyCutSuffix = str(energyCut)+"KeVCut"
    
    signalCutSuffix = ""
    if(args.needSignal):
        signalCutSuffix = "OnlySignal"
    else:
        signalCutSuffix = "SignalAndBackground"
    
    ### open histogram to know the seed information
    if (("_" in inTextFile) and ("hics" in inTextFile)):
        eachName            = inTextFile.split('.')[0].split('_')
        suffixName          = "_".join(eachName[2:])
    else:
        suffixName          = inTextFile.split('.')[0]
    
    
    outFile             = TFile("seedingInformation_"+suffixName+"_"+energyCutSuffix+"_"+signalCutSuffix+".root", "RECREATE")
    outFile.cd()
    
    hAllPossible        = TH1D("hAllPossible", "all possible track combination; bunch crossing; number of track combination", 9508, 0, 9508)
    hSeedPossible       = TH1D("hSeedPossible", "seed track combination; bunch crossing; number of seed track", 9508, 0, 9508)
    hSeedMultiplicity   = TH1D("hSeedMultiplicity", "seed multiplicity; number of seeds; BX", 50, 0, 50)
    hSignalMultiplicity = TH1D("hSignalMultiplicity", "number of signals; number of signals; BX", 50, 0, 50)
    hSigEnergy          = TH1D("hSigEnergy", "signal energy; Energy [GeV]; Entries", 200, 0, 20)
    hSeedEnergy         = TH1D("hSeedEnergy", "seed energy; Energy [GeV]; Entries", 200, 0, 20)
    
    

    #### select the number of bunch crossing for different files, BX is useful only for hics setup
    if 'hics' in inTextFile:
        nBX          = 9508
        checkBXMatch = True
    else:
        nBX          = 1
        checkBXMatch = False
       
    # all the track info is in the following list
    position = []
    # get the information from the text files
    for lines in inputTrackInfo.readlines():
        lines = lines.rstrip()
        eachWord = lines.split()
        bxNumber = int(eachWord[0])
        trackId  = int(eachWord[2])
        pdgId    = int(eachWord[1])
        ### here the preliminary selection of tracks can be done based on trackid and pdgid
        ### do not want photons, as they will not leave track
        if(pdgId==22): continue
        ### if needed, select only the signal  
        if(args.needSignal and not(pdgId==-11 and trackId==1)): continue
        position.append([bxNumber, trackId, int(eachWord[3])-1000, float(eachWord[4]), float(eachWord[5]), float(eachWord[6]), float(eachWord[7])])

    
    for bxCounter in range(1, nBX+1):
        # separate each bx now
        eachBXValue = []
        for tracks in position:
            ### the below is needed for e+laser hics setup
            if (not checkBXMatch):
                eachBXValue.append(tracks)
            elif(checkBXMatch and (bxCounter == tracks[0])):
                eachBXValue.append(tracks)
            else:
                continue

        #pprint.pprint(eachBXValue)
        ### fill up the x,y, z and E values from each of the tracker layers
        allR1Inner = []; allR2Inner = []; allR3Inner = []; allR4Inner = []
        allR1Outer = []; allR2Outer = []; allR3Outer = []; allR4Outer = []
        
        for values in eachBXValue:
            ### x, y, z and E
            if (values[2] == 0):
                allR1Inner.append([values[3], values[4], z1inner, values[5], values[6]])
            elif (values[2] == 1):
                allR1Outer.append([values[3], values[4], z1outer, values[5], values[6]])
            elif (values[2] == 2):
                allR2Inner.append([values[3], values[4], z2inner, values[5], values[6]])
            elif (values[2] == 3):
                allR2Outer.append([values[3], values[4], z2outer, values[5], values[6]])
            elif (values[2] == 4):
                allR3Inner.append([values[3], values[4], z3inner, values[5], values[6]])
            elif (values[2] == 5):
                allR3Outer.append([values[3], values[4], z3outer, values[5], values[6]])
            elif (values[2] == 6):
                allR4Inner.append([values[3], values[4], z4inner, values[5], values[6]])
            elif (values[2] == 7):
                allR4Outer.append([values[3], values[4], z4outer, values[5], values[6]])
            else:
                print("ERROR::Suitable staves not found!")
                quit()
        
        ### removing the overlap region of inner and outer stave
        allR1Unique = allR1Inner
        allR4Unique = allR4Outer
        
        for r1Out in allR1Outer:
            #### remove all points having an x overlap with inner stave: not 100% effective
            if r1Out[0] > (xInnerStaveLastChip + xSizeInnerStaveLastChip/2.):  ### the x position of last chip on inner stave layer 1+half of the x size 
                allR1Unique.append(r1Out)
            
        for r1 in allR1Unique:
            hSigEnergy.Fill(r1[3])
            
        hSignalMultiplicity.Fill(len(allR1Unique))
        
        for r4In in allR4Inner:
            #### remove all points having an x overlap with outer stave: not 100% effective
            if r4In[0] < (xOuterStaveFirstChip - xSizeOuterStaveFirstChip/2.):   ### the x position of last chip on outer stave layer 4-half of the x size 
                allR4Unique.append(r4In)
                
        #pprint.pprint(allR1Unique)
        #pprint.pprint(allR4Unique)
        #### loop over points in r4:
        counter    = 0
        allCounter = 0
        print("For the BX: ", bxCounter," the number of layer 1 tracks: ", len(allR1Unique)," and layer 4 tracks: ", len(allR4Unique))
        allSeedsTrackLines = []
        for r4 in allR4Unique:
            for r1 in allR1Unique:
                seed = makeseed(r1, r4, allR2Inner, allR2Outer, allR3Inner, allR3Outer)
                allCounter += 1
                if(seed):
                    #print("I found a seed here: ", r1, " and ", r4)
                    hSeedEnergy.Fill(p.E())
                    allSeedsTrackLines.append(GetExtendedTrackLine([r1, r4]))
                    counter += 1

        print("Number of seed: ", counter)
        print("Number of combination: ", allCounter)
        hAllPossible.SetBinContent(bxCounter, allCounter)
        hSeedPossible.SetBinContent(bxCounter, counter)
        hSeedMultiplicity.Fill(counter)
        
    outFile.Write()
    outFile.Close()
    inputTrackInfo.close()
    
    ### get the dipole
    dipole = GetDipole()
    ### get the sensors
    sensors = []
    for index, row in df.iterrows():
        detid   = row["detid"]
        layerid = row["layerid"]
        sensors.append( GetSensor(detid,layerid) )

    ### draw
    cnv = TCanvas("cnv","",500,500)
    view = TView.CreateView(1)
    view.SetRange(-600,-500,2500, +600,+500,4500)
    view.ShowAxis()
    dipole.Draw()
    #trackline.Draw()
    #trkpoints.Draw()
    #extendedline.Draw()
    for seedTrack in allSeedsTrackLines:
        seedTrack.Draw()
    for sensor in sensors: sensor.Draw()
    cnv.SaveAs("trackingDiagram_"+suffixName+"_"+energyCutSuffix+"_"+signalCutSuffix+".pdf")
    cnv.SaveAs("trackingDiagram_"+suffixName+"_"+energyCutSuffix+"_"+signalCutSuffix+".root")
    pprint.pprint(cutFlowDict)


if __name__ == "__main__":
    start = time.time()
    main()
    print("-------- The processing time: ",time.time() - start, " s")
