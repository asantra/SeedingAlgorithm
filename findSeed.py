import os
import sys
import time
import pprint
import math
from ROOT import *
import array
from makeTrackDiagrams import *

EseedMin = 1.0; #// GeV
EseedMax = 16.5; #// GeV
LB       = 1.029; ### what is the length of the dipole in meters?
meMeV    = 0.5109989461; #//MeV mass of electron/positron, what to do with photons?
meGeV    = meMeV/1000.;
meGeV2   = meGeV*meGeV;
B        = 1.0  ### 1 Tesla magnetic field used for now
side     = "Positron"
# get the unit vector along one vector
def rUnit2(r1, r2):
    r = (r2-r1).Unit()
    return r
    

def check_clusters(r1, r4, allR2Inner, allR2Outer, allR3Inner, allR3Outer):
    yAbsMargins = 0.2 #// mm (a "road" of 200 microns around the line between r4 and r1)
    xAbsMargins = 0.2 #// mm (a "road" of 200 microns around the line between r4 and r1)
    
    r1min       = [r1[0]-xAbsMargins, r1[1]-yAbsMargins, r1[2]]
    r1max       = [r1[0]+xAbsMargins, r1[1]+yAbsMargins, r1[2]]
    r4min       = [r4[0]-xAbsMargins, r4[1]-yAbsMargins, r4[2]]
    r4max       = [r4[0]+xAbsMargins, r4[1]+yAbsMargins, r4[2]]


    #/// check possible clusters in layer 2, for both inner and outer stave
    y2minInner = yofz(r1min, r4min, 3964.5125) ### zposition of layer 2 inner stave: 3964.5135 mm
    y2maxInner = yofz(r1max, r4max, 3964.5125)
    x2minInner = xofz(r1min, r4min, 3964.5125)
    x2maxInner = xofz(r1max, r4max, 3964.5125)
    
    y2minOuter = yofz(r1min, r4min, 3976.5125) ### zposition of layer 2 outer stave: 3976.5135 mm
    y2maxOuter = yofz(r1max, r4max, 3976.5125)
    x2minOuter = xofz(r1min, r4min, 3976.5125)
    x2maxOuter = xofz(r1max, r4max, 3976.5125)
    
    accept2Inner = False
    accept2Outer = False
    
    ### separate out staves for allR2
    for i2 in range(0, len(allR2Inner)):
        ### remember allR2 has x, y, z and E saved in the list
        accept2yzInner = ( (allR2Inner[i2][1] >= y2minInner and allR2Inner[i2][1] <= y2maxInner) )
        if(not accept2yzInner): continue
    
        accept2xzInner = ( ( allR2Inner[i2][0] >= x2minInner and allR2Inner[i2][0] <= x2maxInner ) )
        if(not accept2xzInner): continue
        accept2Inner = True
        break
    
    ### separate out staves for allR2
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
    y3minInner = yofz(r1min, r4min, 4064.5125) ### z for inner stave layer 3 is 4064.5125
    y3maxInner = yofz(r1max, r4max, 4064.5125)
    x3minInner = xofz(r1min, r4min, 4064.5125)
    x3maxInner = xofz(r1max, r4max, 4064.5125)
    
    y3minOuter = yofz(r1min, r4min, 4076.5125) ### z for Outer stave layer 3 is 4076.5125
    y3maxOuter = yofz(r1max, r4max, 4076.5125)
    x3minOuter = xofz(r1min, r4min, 4076.5125)
    x3maxOuter = xofz(r1max, r4max, 4076.5125)
    
    accept3Inner = False
    accept3Outer = False
    
    ### separate out staves for allR3
    for i3 in range(0, len(allR3Inner)):
        ### remember allR3 has x, y, z and E saved in the list
        accept3yzInner = ( (allR3Inner[i3][1] >= y3minInner and allR3Inner[i3][1] <= y3maxInner) )
        if(not accept3yzInner): continue
    
        accept3xzInner = ( ( allR3Inner[i3][0] >= x3minInner and allR3Inner[i3][0] <= x3maxInner ) )
        if(not accept3xzInner): continue
        accept3Inner = True
        break
    
    ### separate out staves for allR3
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


# making the seeds
def makeseed(r1, r4, allR2Inner, allR2Outer, allR3Inner, allR3Outer):
    
    if(abs(r1[0]) >= abs(r4[0])):
        return False  # // |x1| must be smaller than |x4|
    if(r1[0]*r4[0] < 0):
        return False  ### the x value should have the same sign,  i.e. on the same side of the beam
    if(r1[2] == r4[2]):
        return False  # // if z1=z4..., this is also impossible
    yDipoleExit = yofz(r1, r4, zDipoleExit)
    xDipoleExit = xofz(r1, r4, zDipoleExit)
    ### is this right?
    #if(abs(yDipoleExit) > yDipoleHeight/2):
        #return False
        
    ### The following cuts are coming because of the distribution of the signal
    if(abs(yDipoleExit) > 10.0):
        return False
    ##### why is this? This is accroding to the signal x:y at the dipole exit
    if(abs(xDipoleExit) < 5.0):
        return False  # // the track should point to |x|<~1.0 at the dipole exit
    if(abs(xDipoleExit) > xDipoleWidth/2):
        return False
    if( (side=="Positron" and xDipoleExit < 0) or (side=="Electron" and xDipoleExit > 0)):
        return False
    if(not check_clusters(r1, r4, allR2Inner, allR2Outer, allR3Inner, allR3Outer)):
        return False #// minimum one cluster at layer 2 and one at layer 3

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
    p         = TLorentzVector()
    p.SetPxPyPzE(px, py, pz, math.sqrt(px*px + py*py + pz*pz + meGeV2))
    # // if(i4==0 and side=="Eside") cout << "px=" << px << ", py=" << py << ", pz=" << pz << endl
    # // cout << "side=" << side << ", px=" << px << ", py=" << py << ", pz=" << pz << endl
    # EseedMax = (process=="bppp") ? EseedMaxBPPP : EseedMaxTRIDENT	// GeV
    
    
    ### checking the track energy
    if(p.E() < EseedMin or p.E() > EseedMax): return False

    return True



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
    
    outFile           = TFile("seedingInformation_"+suffixName+".root", "RECREATE")
    outFile.cd()
    hAllPossible      = TH1D("hAllPossible", "all possible track combination; bunch crossing; number of track combination", 9508, 0, 9508)
    hSeedPossible     = TH1D("hSeedPossible", "seed track combination; bunch crossing; number of seed track", 9508, 0, 9508)
    hSeedMultiplicity = TH1D("hSeedMultiplicity", "hSeedMultiplicity", 50, 0, 50)
    
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
        ### select if only background or signal tracks wanted
        if(pdgId==22): continue
        #if(trackId==1): continue                                                                                                                                                                                                                                                                                                                        
        position.append([bxNumber, trackId, int(eachWord[3])-1000, float(eachWord[4]), float(eachWord[5]), float(eachWord[6])])

    
    for bxCounter in range(1,9509):
        # separate each bx now
        eachBXValue = []
        for tracks in position:
            #if True: #
            ### the below is needed for e+laser hics setup
            if tracks[0] == bxCounter:
                eachBXValue.append(tracks)

        
        ### fill up the x,y, z and E values from each of the tracker layers
        allR1Inner = []; allR2Inner = []; allR3Inner = []; allR4Inner = []
        allR1Outer = []; allR2Outer = []; allR3Outer = []; allR4Outer = []
        
        for values in eachBXValue:
            zPosition = getStaveZ(values[2])
            ### x, y, z and E
            if (values[2] == 0):
                allR1Inner.append([values[3], values[4], zPosition])
            elif (values[2] == 1):
                allR1Outer.append([values[3], values[4], zPosition])
            elif (values[2] == 2):
                allR2Inner.append([values[3], values[4], zPosition])
            elif (values[2] == 3):
                allR2Outer.append([values[3], values[4], zPosition])
            elif (values[2] == 4):
                allR3Inner.append([values[3], values[4], zPosition])
            elif (values[2] == 5):
                allR3Outer.append([values[3], values[4], zPosition])
            elif (values[2] == 6):
                allR4Inner.append([values[3], values[4], zPosition])
            elif (values[2] == 7):
                allR4Outer.append([values[3], values[4], zPosition])
            else:
                print("ERROR::Suitable staves not found!")
                quit()
        
        ### removing the overlap region of inner and outer stave
        allR1Unique = allR1Inner
        allR4Unique = allR4Outer
        
        for r1Out in allR1Outer:
            #### remove all points having an x overlap with inner stave: not 100% effective
            if r1Out[0] > (308.53 + 29.94176/2.):  ### the x position of last chip on inner stave layer 1+half of the x size 
                allR1Unique.append(r1Out)
                
        for r4In in allR4Inner:
            #### remove all points having an x overlap with outer stave: not 100% effective
            if r4In[0] < (298.53 - 29.94176/2.):   ### the x position of last chip on outer stave layer 4-half of the x size 
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
                    print("I found a seed here: ", r1, " and ", r4)
                    allSeedsTrackLines.append(GetExtendedTrackLine([r1, r4]))
                    counter += 1

        print("Number of seed: ", counter)
        print("Number of combination: ", allCounter)
        hAllPossible.SetBinContent(bxCounter, allCounter)
        hSeedPossible.SetBinContent(bxCounter, counter)
        hSeedMultiplicity.Fill(counter)
        
    outFile.Write()
    outFile.Close()
    
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
    cnv.SaveAs("trackingDiagram_"+suffixName+".pdf")
    cnv.SaveAs("trackingDiagram_"+suffixName+".root")

if __name__ == "__main__":
    start = time.time()
    main()
    print("-------- The processing time: ",time.time() - start, " s")
