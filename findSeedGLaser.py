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
import csv, array
import numpy as np
### needed for plotting of tracks and dimension from the spreadsheet 
from makeTrackDiagrams import *

### the seed energy width
EseedMin       = 2.0  # GeV
EseedMax       = 15.0 # GeV
EseedMinPrelim = 0.2  # GeV
EseedMaxPrelim = 18.0 # GeV

### the seed pY width
PseedMin = -0.01 # GeV
PseedMax = 0.01 # GeV

###
xDipoleExitMax = 280.
xDipoleExitMin = 40.

### cut on nseeds
nseedsNoFitMax  = 5
nseedsNoFitMax2 = 30

#//MeV mass of electron/positron
meMeV    = 0.5109989461 
meGeV    = meMeV/1000.
### me^2 in GeV
meGeV2   = meGeV*meGeV

### the width of the road around the pivot tracks
yAbsMargins        = 0.13 # mm (a "road" of 200 microns around the line between r4 and r1)
xAbsMargins        = 0.13 # mm (a "road" of 200 microns around the line between r4 and r1)

zPositionListInner = [z2inner, z3inner]
zPositionListOuter = [z2outer, z3outer]
ddList = []



#### This is the cutflow dictionary
cutFlowDict = OrderedDict(
    [('noCut', 0), 
     ('x1Gtx4',0), 
     ('x1*x4Negative',0), 
     ('z1Eqz4', 0), 
     ('yDipoleExitGt5p4',0), 
     ('xDipoleExitLt25',0), 
     ('xDipoleExitGt165', 0), 
     ('xDipoleExitLt0',0), 
     ('seedEnergy', 0),
     ('checkClusterTracksMiddleLayers', 0),
     ('checkClusterFit', 0),
     ('trackEnergy', 0),
     ('checkClusterXDistance', 0),
     ('checkClusterTrackPy', 0),
     ('checkClusterTrackPyLoose', 0),
     ('checkClusterTrackPyTight', 0)]
    )

histos =  {'hSVDValues2Global': TH1D("hSVDValues2Global", "output of SVD[1]; fit quality [1]; Events", 200, 0, 0.3), 
           'hSVDValues3Global': TH1D("hSVDValues3Global", "output of SVD[2]; fit quality [2]; Events", 1000, 0, 0.1),
           'hSVDValues2LooseGlobal': TH1D("hSVDValues2LooseGlobal", "output of SVD[1]; fit quality [1]; Events", 200, 0, 0.3), 
           'hSVDValues3LooseGlobal': TH1D("hSVDValues3LooseGlobal", "output of SVD[2]; fit quality [2]; Events", 1000, 0, 0.1),
           'hSVDValues2TightGlobal': TH1D("hSVDValues2TightGlobal", "output of SVD[1]; fit quality [1]; Events", 200, 0, 0.3), 
           'hSVDValues3TightGlobal': TH1D("hSVDValues3TightGlobal", "output of SVD[2]; fit quality [2]; Events", 1000, 0, 0.1),
           'hSeedDistanceGlobal': TH1D("hSeedDistanceGlobal", "seed distance wrt analytical line; d [m]; Entries", 500, 0, 0.1),
           'hSeedDistanceLooseGlobal': TH1D("hSeedDistanceLooseGlobal", "seed distance wrt analytical line (loose); d [m]; Entries", 500, 0, 0.10),
           'hSeedDistanceTightGlobal': TH1D("hSeedDistanceTightGlobal", "seed distance wrt analytical line (tight); d [m]; Entries", 500, 0, 0.10),
           'hSVDValues2TightVsEnergyGlobal': TH2D("hSVDValues2TightVsEnergyGlobal", "output of SVD[1]; Energy [GeV]; fit quality [1]; Events", 35, 0, 17.5, 100, 0, 0.3),
           'hSeedDistanceTightVsEnergyGlobal': TH2D("hSeedDistanceTightVsEnergyGlobal", "seed distance wrt analytical line (tight) vs Energy; Energy [GeV]; d[m]; Entries", 35, 0, 17.5, 100, 0, 0.10)
           }



### read the dEdX vs E curve for the electron
results = csv.reader(open("dEdXSilicon_ForElectron.csv"), delimiter=",")
xdEdX   = array.array('f',[]); ydEdX = array.array('f', [])
for eachRow in results:
    xdEdX.append(float(eachRow[0]))
    ydEdX.append(float(eachRow[1]))
    

### This is to find the dEdx
def getdEdX(energy):
    ### the dEdX vs E curve, x axis energy is in MeV, dEdX in y is in keV/um, remember the unit difference
    mevEnergy        = energy*1e3
    interpolateddEdX = np.interp(mevEnergy, xdEdX, ydEdX)
    ### interpolateddEdX in keV/um
    return interpolateddEdX
    
    
### This is to find the absorbed energy in the FPC layer and cooling region
def energyAbsorbed(staveId, energy, vtxZ):
    if(staveId == 0 or staveId == 8):
        zPos = z1inner
    elif(staveId == 1 or staveId == 9):
        zPos = z1outer
    elif(staveId == 2 or staveId == 10):
        zPos = z2inner
    elif(staveId == 3 or staveId == 11):
        zPos = z2outer
    elif(staveId == 4 or staveId == 12):
        zPos = z3inner
    elif(staveId == 5 or staveId == 13):
        zPos = z3outer
    elif(staveId == 6 or staveId == 14):
        zPos = z4inner
    elif(staveId == 7 or staveId == 15):
        zPos = z4outer
    else:
        print("No suitable stave found!!! Exiting")
        quit()
        
    ### front side
    if((zPos - vtxZ) > 0):
        absorbingFactor = 150
    ### back side
    else:
        absorbingFactor = 250
    
    ### this dEdX is in keV/um, multiply by the width to get the energy loss
    dEdX            = getdEdX(energy)
    absorbed        = dEdX*absorbingFactor
    remainingEnergy = energy*1e6 - absorbed
    ### print("staveId: ", staveId, " zPos: ", zPos, " vtxZ: ", vtxZ, " energy (keV): ", energy*1e6, "dEdX (keV/um): ", dEdX, " factor: ", absorbingFactor, " absorbed (keV): ", absorbed, " remainingEnergy (keV): ", remainingEnergy)
    return remainingEnergy
    
    
        
### functions needed for distance cut from the analytical line in x1 vs x4 plane
## get R in m given momentum p
def getR(p):
    R = p/(0.3*B)
    return R

### get R^2 in m^2 given momentum p
def getR2(p):
    R = getR(p)
    return R*R

### get the xTangent in m given ZTangent in m and momentum p
def getXTangent(ZT,p):
    R2 = getR2(p)
    XT = (R2/LB-ZT)*(LB/math.sqrt(R2-LB2))
    return XT

### get x in Sasha's coordinate in m given zTangent and momentum p
def getx(ZT,p):
    XT = getXTangent(ZT,p)
    R = getR(p)
    x = R-XT
    return x

### get r1 and r2 (x,y --> 2 points) given the Z0 and Z4. This is in m
def getR1R2(p1, p2):
    r1    = [getx(Z0,p1), getx(Z4,p1)]
    r2    = [getx(Z0,p2), getx(Z4, p2)]
    #print("r1: ",r1, " r2: ",r2)
    #quit()
    return r1, r2

#### find the distance from the analytical line to the r1 and r2 point in m
def Distance(x0Test,x4Test,r1,r2):
    aOb = -(r1[1]-r2[1])/(r1[0]-r2[0])
    cOb = (r1[1]-r2[1])/(r1[0]-r2[0])*r1[0] - r1[1]
    ### now find the distance
    # d = |a*x0Test + b*x4Test + c|/sqrt(a^2 + b^2) --> d = |(a/b)*x0Test + x4Test + (c/b)|/sqrt((a/b)^2 + 1^2)
    d = abs(aOb*x0Test + x4Test + (cOb))/math.sqrt(aOb*aOb + 1)
    return d

    
    
# get the unit vector along one vector
def rUnit2(r1, r2):
    r = (r2-r1).Unit()
    return r

### get the 3d line
def line3d(z,m_xz,c_xz,m_yz,c_yz):
   # the intersection of those two planes and
   # the function for the line would be:
   # or:
   x = (z - c_xz)/m_xz
   y = (z - c_yz)/m_yz
   return x,y


### get the TLorentzVector of the seed
def getSeedMomentum(r1, r4):
    #### find the energy of the tracks
    x0        = 0
    z0        = zofx(r1, r4, x0)
    xExit     = abs(xofz(r1, r4, zDipoleActiveExit))
    yExit     = yofz(r1, r4, zDipoleActiveExit)
    H         = abs((zDipoleActiveExit-z0))*mm2m ### converting H from mm to m
    xExitInM  = xExit*mm2m                       ### converting xExit from mm to m
    R         = H*(LB)/xExitInM + xExitInM       ### // This is the radius of curvature for a track
    P         = 0.3*B*R                          ### here B in Tesla and R in m
    v1        = TVector2(r1[2], r1[1])
    v4        = TVector2(r4[2], r4[1])
    u         = rUnit2(v1, v4)
    uz        = u.X()
    uy        = u.Y()
    px        = 0    
    py        = P*uy
    pz        = P*uz
    ### The seed lorentz vector, need for energy plotting
    p         = TLorentzVector()
    p.SetPxPyPzE(px, py, pz, math.sqrt(px*px + py*py + pz*pz + meGeV2))
    
    return p, xExit, yExit



#### fit with chi2
def seed3dChi2fit(r1,r2Inner, r2Outer, r3Inner, r3Outer, r4):
   
   xList = [r1[0],r4[0]]
   yList = [r1[1],r4[1]]
   zList = [r1[2],r4[2]]
   
   if len(r2Inner)>0:
       xList.append(r2Inner[0])
       yList.append(r2Inner[1])
       zList.append(r2Inner[2])
       
   if len(r2Outer)>0:
       xList.append(r2Outer[0])
       yList.append(r2Outer[1])
       zList.append(r2Outer[2])
       
   if len(r3Inner)>0:
       xList.append(r3Inner[0])
       yList.append(r3Inner[1])
       zList.append(r3Inner[2])
       
   if len(r3Outer)>0:
       xList.append(r3Outer[0])
       yList.append(r3Outer[1])
       zList.append(r3Outer[2])
       
   x = np.array(xList)
   y = np.array(yList)
   z = np.array(zList)
   
   # this will find the slope and x-intercept of a plane
   # parallel to the y-axis that best fits the data
   A_xz         = np.vstack((x, np.ones(len(x)))).T
   result_xz    = np.linalg.lstsq(A_xz, z,rcond=None)
   m_xz, c_xz   = result_xz[0]
   residuals_xz = result_xz[1]

   # again for a plane parallel to the x-axis
   A_yz         = np.vstack((y, np.ones(len(y)))).T
   result_yz    = np.linalg.lstsq(A_yz, z,rcond=None)
   m_yz, c_yz   = result_yz[0]
   residuals_yz = result_yz[1]
   
   zz           = np.array([zDipoleActiveExit, r4[2]])
   xx,yy        = line3d(zz, m_xz, c_xz, m_yz, c_yz)
   lfit         = TPolyLine3D()
   for i in range(2):
       lfit.SetNextPoint(zz[i],yy[i],xx[i])
   
   return residuals_xz,residuals_yz, lfit
   
   
### use svd fit
def seed3dfitSVDWithList(r1,r2Inner, r2Outer, r3Inner, r3Outer, r4):
   
   xList = [r1[0],r4[0]]
   yList = [r1[1],r4[1]]
   zList = [r1[2],r4[2]]
   
   if len(r2Inner)>0:
       xList.append(r2Inner[0])
       yList.append(r2Inner[1])
       zList.append(r2Inner[2])
       
   if len(r2Outer)>0:
       xList.append(r2Outer[0])
       yList.append(r2Outer[1])
       zList.append(r2Outer[2])
       
   if len(r3Inner)>0:
       xList.append(r3Inner[0])
       yList.append(r3Inner[1])
       zList.append(r3Inner[2])
       
   if len(r3Outer)>0:
       xList.append(r3Outer[0])
       yList.append(r3Outer[1])
       zList.append(r3Outer[2])
       
   x = np.array(xList)
   y = np.array(yList)
   z = np.array(zList)


   data = np.concatenate((x[:, np.newaxis], 
                          y[:, np.newaxis], 
                          z[:, np.newaxis]), 
                         axis=1)

   # Calculate the mean of the points, i.e. the 'center' of the cloud
   datamean = data.mean(axis=0)

   # Do an SVD on the mean-centered data (Singular Value Decomposition)
   uu, dd, vv = np.linalg.svd(data - datamean) 

   # Now vv[0] contains the first principal component, i.e. the direction
   # vector of the 'best fit' line in the least squares sense.

   # Now generate some points along this best fit line, for plotting.

   # I use -7, 7 since the spread of the data is roughly 14
   # and we want it to have mean 0 (like the points we did
   # the svd on). Also, it's a straight line, so we only need 2 points.
   # linepts = vv[0] * np.mgrid[-7:7:2j][:, np.newaxis]
   linepts = vv[0] * np.mgrid[-170:170:2j][:, np.newaxis]

   # shift by the mean to get the line in the right place
   linepts += datamean
   #print("linepoints: ", linepts)
   return linepts, dd ## dd is a 1D array of the data singular values


### draw the SVD fit
def drawFit(name,linepts,hits):
    g = TGraph2D()
    g.SetMarkerSize(3)
    g.SetMarkerStyle(20)
    g.SetMarkerColor(ROOT.kRed)
    g.SetLineColor(ROOT.kRed)
    for i in range(len(hits)):
        g.SetPoint(i,hits[i][0],hits[i][1],hits[i][2])
    
    g.GetXaxis().SetRangeUser(0, 650)
    g.GetYaxis().SetRangeUser(-10,+10)
    g.GetZaxis().SetRangeUser(3900, 4400)
    
    
    lfit = TPolyLine3D()
    for point in linepts:
        lfit.SetNextPoint(point[0],point[1],point[2])
    lfit.SetLineColor(ROOT.kBlue)
    cnv = TCanvas("","",2000,2000)
    view = TView.CreateView(1)
    xviewmin = 3900 #if(isel(r1[0])) else xPsideL
    xviewmax = 4200 #xEsideR if(isel(r1[0])) else 0
    #view.SetAutoRange()
    view.SetRange(0,-0.8, xviewmin , 650,+0.8,xviewmax)
    view.ShowAxis()
    g.Draw("p0")
    lfit.Draw("same")
    cnv.SaveAs(name)
   
### get the expected number of hits in the layer 2 and layer 3 of the tracker given the position of the track
def getExpectedHits(r1, r4, energy, side):
    global xAbsMargins, yAbsMargins
    
    if(energy < 4.0):
        xAbsMargins = 0.22
        yAbsMargins = 0.22
    else:
        xAbsMargins = 0.13
        yAbsMargins = 0.13
    
    
    expectedHits2Inner = 0
    expectedHits2Outer = 0
    
    expectedHits3Inner = 0
    expectedHits3Outer = 0
    
    r1min              = [r1[0]-xAbsMargins, r1[1]-yAbsMargins, r1[2]]
    r1max              = [r1[0]+xAbsMargins, r1[1]+yAbsMargins, r1[2]]
    r4min              = [r4[0]-xAbsMargins, r4[1]-yAbsMargins, r4[2]]
    r4max              = [r4[0]+xAbsMargins, r4[1]+yAbsMargins, r4[2]]
    
    # check possible clusters in layer 2, for both inner and outer stave
    x2minInner         = xofz(r1min, r4min, z2inner)
    x2maxInner         = xofz(r1max, r4max, z2inner)
    
    x2minOuter         = xofz(r1min, r4min, z2outer)
    x2maxOuter         = xofz(r1max, r4max, z2outer)
    
    #/// check possible clusters in layer 3, for both inner and outer stave
    x3minInner         = xofz(r1min, r4min, z3inner)
    x3maxInner         = xofz(r1max, r4max, z3inner)

    x3minOuter         = xofz(r1min, r4min, z3outer)
    x3maxOuter         = xofz(r1max, r4max, z3outer)
    
    
    x2Inner            = xofz(r1, r4, z2inner)
    x2Outer            = xofz(r1, r4, z2outer)
    x3Inner            = xofz(r1, r4, z3inner)
    x3Outer            = xofz(r1, r4, z3outer)
    
    for chipName, boundaries in xBoundaries.items(): 
        #print(chipName)
        if(side=="Positron"):
            if "layerid2" in chipName and x2Inner > boundaries[0] and x2Inner < boundaries[1]:
                expectedHits2Inner += 1
            if "layerid3" in chipName and x2Outer > boundaries[0] and x2Outer < boundaries[1]:
                expectedHits2Outer += 1
            if "layerid4" in chipName and x3Inner > boundaries[0] and x3Inner < boundaries[1]:
                expectedHits3Inner += 1
            if "layerid5" in chipName and x3Outer > boundaries[0] and x3Outer < boundaries[1]:
                expectedHits3Outer += 1
        else:
            if "layerid10" in chipName and x2Inner > boundaries[0] and x2Inner < boundaries[1]:
                expectedHits2Inner += 1
            if "layerid11" in chipName and x2Outer > boundaries[0] and x2Outer < boundaries[1]:
                expectedHits2Outer += 1
            if "layerid12" in chipName and x3Inner > boundaries[0] and x3Inner < boundaries[1]:
                expectedHits3Inner += 1
            if "layerid13" in chipName and x3Outer > boundaries[0] and x3Outer < boundaries[1]:
                expectedHits3Outer += 1
    
    
    return expectedHits2Inner, expectedHits2Outer, expectedHits3Inner, expectedHits3Outer
        
    
    

### check if there are hits along one track in layer 2 and layer 3
def check_clusters(r1, r4, allR2Inner, allR2Outer, allR3Inner, allR3Outer, side, energy):
    global ddList, xAbsMargins, yAbsMargins
    expectedHits2Inner, expectedHits2Outer, expectedHits3Inner, expectedHits3Outer = getExpectedHits(r1, r4, energy, side)
    hitsOnRoad2Inner = 0   ### how many tracks accepted in the road along the r1 and r4
    hitsOnRoad2Outer = 0
    hitsOnRoad3Inner = 0   ### how many tracks accepted in the road along the r1 and r4
    hitsOnRoad3Outer = 0
    
    
    if(energy < 4.0):
        xAbsMargins = 0.22
        yAbsMargins = 0.22
    else:
        xAbsMargins = 0.13
        yAbsMargins = 0.13
    
    r1min            = [r1[0]-xAbsMargins, r1[1]-yAbsMargins, r1[2]]
    r1max            = [r1[0]+xAbsMargins, r1[1]+yAbsMargins, r1[2]]
    r4min            = [r4[0]-xAbsMargins, r4[1]-yAbsMargins, r4[2]]
    r4max            = [r4[0]+xAbsMargins, r4[1]+yAbsMargins, r4[2]]


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
    innerR2FromMatching = []
    outerR2FromMatching = []
    for i2 in range(0, len(allR2Inner)):
        ### remember allR2 has x, y, z and E saved in the list
        accept2yzInner = ( (allR2Inner[i2][1] >= y2minInner and allR2Inner[i2][1] <= y2maxInner) )
        if(not accept2yzInner): continue
        
    
        accept2xzInner = ( ( allR2Inner[i2][0] >= x2minInner and allR2Inner[i2][0] <= x2maxInner ) )
        if(not accept2xzInner): continue
        accept2Inner        = True
        hitsOnRoad2Inner   += 1
        innerR2FromMatching.append(allR2Inner[i2]) ### collect all matched tracks from layer 2, remove the break
    
    ### separate out staves for allR2
    for i2 in range(0, len(allR2Outer)):
        ### remember allR2 has x, y, z and E saved in the list
        accept2yzOuter = ( (allR2Outer[i2][1] >= y2minOuter and allR2Outer[i2][1] <= y2maxOuter) )
        if(not accept2yzOuter): continue
    
        accept2xzOuter = ( ( allR2Outer[i2][0] >= x2minOuter and allR2Outer[i2][0] <= x2maxOuter ) )
        if(not accept2xzOuter): continue
        accept2Outer      = True
        hitsOnRoad2Outer += 1
        outerR2FromMatching.append(allR2Outer[i2]) ### collect all matched tracks from layer 2, remove the break
    
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
    innerR3FromMatching = []
    outerR3FromMatching = []
    for i3 in range(0, len(allR3Inner)):
        ### remember allR3 has x, y, z and E saved in the list
        accept3yzInner = ( (allR3Inner[i3][1] >= y3minInner and allR3Inner[i3][1] <= y3maxInner) )
        if(not accept3yzInner): continue
    
        accept3xzInner = ( ( allR3Inner[i3][0] >= x3minInner and allR3Inner[i3][0] <= x3maxInner ) )
        if(not accept3xzInner): continue
        accept3Inner        = True
        hitsOnRoad3Inner   += 1
        innerR3FromMatching.append(allR3Inner[i3]) ### collect all matched tracks from layer 3, remove the break
        #break
    
    ### separate out staves for allR3
    for i3 in range(0, len(allR3Outer)):
        ### remember allR3 has x, y, z and E saved in the list
        accept3yzOuter = ( (allR3Outer[i3][1] >= y3minOuter and allR3Outer[i3][1] <= y3maxOuter) )
        if(not accept3yzOuter): continue
    
        accept3xzOuter = ( ( allR3Outer[i3][0] >= x3minOuter and allR3Outer[i3][0] <= x3maxOuter ) )
        if(not accept3xzOuter): continue
        accept3Outer        = True
        hitsOnRoad3Outer   += 1
        outerR3FromMatching.append(allR3Outer[i3]) ### collect all matched tracks from layer 3, remove the break
    
    ### find the number of actual matched tracks
    
    nMatched = (expectedHits2Inner<=hitsOnRoad2Inner)+(expectedHits2Outer<=hitsOnRoad2Outer)+(expectedHits3Inner<=hitsOnRoad3Inner)+(expectedHits3Outer<=hitsOnRoad3Outer)
    ### find the number of expected tracks
    nExpected = expectedHits2Inner + expectedHits2Outer + expectedHits3Inner + expectedHits3Outer
    #print("nMatched: ", nMatched, " nExpected: ", nExpected)
    ### only if we have hit in one of the two inner layer of the tracker, we select the seed track, otherwise we reject the track
    if(nMatched>=3):
        return nMatched, nExpected, innerR2FromMatching, outerR2FromMatching, innerR3FromMatching, outerR3FromMatching
    else:
        return nMatched, nExpected, [], [], [], []
        
        
#### prepare the SVD fit for the seed track
def makeSeedFit(r1, r4, nMatched, nExpected, innerR2FromMatching, outerR2FromMatching, innerR3FromMatching, outerR3FromMatching, nseedsNoFit, side="Positron"):
    ### now check if the SVD match works
    allDDList = {'dd':[], 'linepts':[], 'i2Inner':[], 'i2Outer':[], 'i3Inner':[], 'i3Outer':[]}
    #### number of matches tracks in each of the staves of the layer 2 and layer 3
    n2i = len(innerR2FromMatching)
    n2o = len(outerR2FromMatching)
    n3i = len(innerR3FromMatching)
    n3o = len(outerR3FromMatching)
    
    
    if(n2i> 0) :
        for i2Inner in range(n2i):
            r2Inner = innerR2FromMatching[i2Inner]
            if(n2o > 0):
                for i2Outer in range(n2o):
                    r2Outer = outerR2FromMatching[i2Outer]
                    if(n3i > 0):
                        for i3Inner in range(n3i):
                            r3Inner = innerR3FromMatching[i3Inner]
                            if(n3o > 0):
                                for i3Outer in range(n3o):
                                    r3Outer = outerR3FromMatching[i3Outer]
                                    lfitpts, dd = seed3dfitSVDWithList(r1,r2Inner, r2Outer, r3Inner, r3Outer, r4)
                                    allDDList['dd'].append(dd.tolist())
                                    allDDList['i2Inner'].append(i2Inner)
                                    allDDList['i2Outer'].append(i2Outer)
                                    allDDList['i3Inner'].append(i3Inner)
                                    allDDList['i3Outer'].append(i3Outer)
                                    allDDList['linepts'].append(lfitpts)
                            else:
                                    r3Outer = []
                                    lfitpts, dd = seed3dfitSVDWithList(r1,r2Inner, r2Outer, r3Inner, r3Outer, r4)
                                    allDDList['dd'].append(dd.tolist())
                                    allDDList['i2Inner'].append(i2Inner)
                                    allDDList['i2Outer'].append(i2Outer)
                                    allDDList['i3Inner'].append(i3Inner)
                                    allDDList['i3Outer'].append(-1)
                                    allDDList['linepts'].append(lfitpts)
                                    
                    else:
                        r3Inner = []
                        if(n3o > 0):
                            for i3Outer in range(n3o):
                                r3Outer = outerR3FromMatching[i3Outer]
                                lfitpts, dd = seed3dfitSVDWithList(r1,r2Inner, r2Outer, r3Inner, r3Outer, r4)
                                allDDList['dd'].append(dd.tolist())
                                allDDList['i2Inner'].append(i2Inner)
                                allDDList['i2Outer'].append(i2Outer)
                                allDDList['i3Inner'].append(-1)
                                allDDList['i3Outer'].append(i3Outer)
                                allDDList['linepts'].append(lfitpts)
                        else:
                                r3Outer = []
                                lfitpts, dd = seed3dfitSVDWithList(r1,r2Inner, r2Outer, r3Inner, r3Outer, r4)
                                allDDList['dd'].append(dd.tolist())
                                allDDList['i2Inner'].append(i2Inner)
                                allDDList['i2Outer'].append(i2Outer)
                                allDDList['i3Inner'].append(-1)
                                allDDList['i3Outer'].append(-1)
                                allDDList['linepts'].append(lfitpts)
                                
                                
            else:
                r2Outer = []
                if(n3i > 0):
                    for i3Inner in range(n3i):
                        r3Inner = innerR3FromMatching[i3Inner]
                        if(n3o > 0):
                            for i3Outer in range(n3o):
                                r3Outer = outerR3FromMatching[i3Outer]
                                lfitpts, dd = seed3dfitSVDWithList(r1,r2Inner, r2Outer, r3Inner, r3Outer, r4)
                                allDDList['dd'].append(dd.tolist())
                                allDDList['i2Inner'].append(i2Inner)
                                allDDList['i2Outer'].append(-1)
                                allDDList['i3Inner'].append(i3Inner)
                                allDDList['i3Outer'].append(i3Outer)
                                allDDList['linepts'].append(lfitpts)
                        else:
                                r3Outer = []
                                lfitpts, dd = seed3dfitSVDWithList(r1,r2Inner, r2Outer, r3Inner, r3Outer, r4)
                                allDDList['dd'].append(dd.tolist())
                                allDDList['i2Inner'].append(i2Inner)
                                allDDList['i2Outer'].append(-1)
                                allDDList['i3Inner'].append(i3Inner)
                                allDDList['i3Outer'].append(-1)
                                allDDList['linepts'].append(lfitpts)
                                
                else:
                    r3Inner = []
                    if(n3o > 0):
                        for i3Outer in range(n3o):
                            r3Outer = outerR3FromMatching[i3Outer]
                            lfitpts, dd = seed3dfitSVDWithList(r1,r2Inner, r2Outer, r3Inner, r3Outer, r4)
                            allDDList['dd'].append(dd.tolist())
                            allDDList['i2Inner'].append(i2Inner)
                            allDDList['i2Outer'].append(-1)
                            allDDList['i3Inner'].append(-1)
                            allDDList['i3Outer'].append(i3Outer)
                            allDDList['linepts'].append(lfitpts)
                    else:
                            r3Outer = []
                            lfitpts, dd = seed3dfitSVDWithList(r1,r2Inner, r2Outer, r3Inner, r3Outer, r4)
                            allDDList['dd'].append(dd.tolist())
                            allDDList['i2Inner'].append(i2Inner)
                            allDDList['i2Outer'].append(-1)
                            allDDList['i3Inner'].append(-1)
                            allDDList['i3Outer'].append(-1)
                            allDDList['linepts'].append(lfitpts) 
                            
                            
    else:
        r2Inner = []
        if(n2o > 0):
            for i2Outer in range(n2o):
                r2Outer = outerR2FromMatching[i2Outer]
                if(n3i > 0):
                    for i3Inner in range(n3i):
                        r3Inner = innerR3FromMatching[i3Inner]
                        if(n3o > 0):
                            for i3Outer in range(n3o):
                                r3Outer = outerR3FromMatching[i3Outer]
                                lfitpts, dd = seed3dfitSVDWithList(r1,r2Inner, r2Outer, r3Inner, r3Outer, r4)
                                allDDList['dd'].append(dd.tolist())
                                allDDList['i2Inner'].append(-1)
                                allDDList['i2Outer'].append(i2Outer)
                                allDDList['i3Inner'].append(i3Inner)
                                allDDList['i3Outer'].append(i3Outer)
                                allDDList['linepts'].append(lfitpts)
                        else:
                                r3Outer = []
                                lfitpts, dd = seed3dfitSVDWithList(r1,r2Inner, r2Outer, r3Inner, r3Outer, r4)
                                allDDList['dd'].append(dd.tolist())
                                allDDList['i2Inner'].append(-1)
                                allDDList['i2Outer'].append(i2Outer)
                                allDDList['i3Inner'].append(i3Inner)
                                allDDList['i3Outer'].append(-1)
                                allDDList['linepts'].append(lfitpts)
                                
                else:
                    r3Inner = []
                    if(n3o > 0):
                        for i3Outer in range(n3o):
                            r3Outer = outerR3FromMatching[i3Outer]
                            lfitpts, dd = seed3dfitSVDWithList(r1,r2Inner, r2Outer, r3Inner, r3Outer, r4)
                            allDDList['dd'].append(dd.tolist())
                            allDDList['i2Inner'].append(-1)
                            allDDList['i2Outer'].append(i2Outer)
                            allDDList['i3Inner'].append(-1)
                            allDDList['i3Outer'].append(i3Outer)
                            allDDList['linepts'].append(lfitpts)
                    else:
                            r3Outer = []
                            lfitpts, dd = seed3dfitSVDWithList(r1,r2Inner, r2Outer, r3Inner, r3Outer, r4)
                            allDDList['dd'].append(dd.tolist())
                            allDDList['i2Inner'].append(-1)
                            allDDList['i2Outer'].append(i2Outer)
                            allDDList['i3Inner'].append(-1)
                            allDDList['i3Outer'].append(-1)
                            allDDList['linepts'].append(lfitpts)
                            
                            
        else:
            r2Outer = []
            if(n3i > 0):
                for i3Inner in range(n3i):
                    r3Inner = innerR3FromMatching[i3Inner]
                    if(n3o > 0):
                        for i3Outer in range(n3o):
                            r3Outer = outerR3FromMatching[i3Outer]
                            lfitpts, dd = seed3dfitSVDWithList(r1,r2Inner, r2Outer, r3Inner, r3Outer, r4)
                            allDDList['dd'].append(dd.tolist())
                            allDDList['i2Inner'].append(-1)
                            allDDList['i2Outer'].append(-1)
                            allDDList['i3Inner'].append(i3Inner)
                            allDDList['i3Outer'].append(i3Outer)
                            allDDList['linepts'].append(lfitpts)
                    else:
                            r3Outer = []
                            lfitpts, dd = seed3dfitSVDWithList(r1,r2Inner, r2Outer, r3Inner, r3Outer, r4)
                            allDDList['dd'].append(dd.tolist())
                            allDDList['i2Inner'].append(-1)
                            allDDList['i2Outer'].append(-1)
                            allDDList['i3Inner'].append(i3Inner)
                            allDDList['i3Outer'].append(-1)
                            allDDList['linepts'].append(lfitpts)
                            
            else:
                r3Inner = []
                if(n3o > 0):
                    for i3Outer in range(n3o):
                        r3Outer = outerR3FromMatching[i3Outer]
                        lfitpts, dd = seed3dfitSVDWithList(r1,r2Inner, r2Outer, r3Inner, r3Outer, r4)
                        allDDList['dd'].append(dd.tolist())
                        allDDList['i2Inner'].append(-1)
                        allDDList['i2Outer'].append(-1)
                        allDDList['i3Inner'].append(-1)
                        allDDList['i3Outer'].append(i3Outer)
                        allDDList['linepts'].append(lfitpts)
                else:
                        return False, {}
    
    ### start with a very high value, in the end we will take the smallest one
    ddValue0 = 1e11
    ddValue1 = 1e11
    ddValue2 = 1e11
    ### the index of the best matched track
    iWinner = -1
    #print("ddList:", allDDList)
    for i in range(len(allDDList['dd'])):
        ### calculate the seed energy
        xDipoleFromFit = xofz(allDDList['linepts'][i][0], allDDList['linepts'][i][1], zDipoleExit)*mm2m
        xLayer4FromFit = xofz(allDDList['linepts'][i][0], allDDList['linepts'][i][1], r4[2])*mm2m
        
        yDipoleFromFit = yofz(allDDList['linepts'][i][0], allDDList['linepts'][i][1], zDipoleExit)*mm2m
        yLayer4FromFit = yofz(allDDList['linepts'][i][0], allDDList['linepts'][i][1], r4[2])*mm2m

        ### here r1 and r4 must be in mm
        R1 = [xDipoleFromFit*m2mm, yDipoleFromFit*m2mm, zDipoleExit]
        R4 = [xLayer4FromFit*m2mm, yLayer4FromFit*m2mm, r4[2]]
        
        #### re-evaluate the seed energy from the best fit track
        pSeed, xExit, yExit = getSeedMomentum(R1, R4)
        
        
        ddValue = allDDList['dd'][i]
        #print("ddValue1", ddValue)
        histos['hSVDValues2Global'].Fill(ddValue[1])
        histos['hSVDValues3Global'].Fill(ddValue[2])
        if(nMatched==4):
            histos['hSVDValues2TightGlobal'].Fill(ddValue[1])
            histos['hSVDValues3TightGlobal'].Fill(ddValue[2])
            histos['hSVDValues2TightVsEnergyGlobal'].Fill(pSeed.E(),ddValue[1])
        else:
            histos['hSVDValues2LooseGlobal'].Fill(ddValue[1])
            histos['hSVDValues3LooseGlobal'].Fill(ddValue[2])
            
        
        ### printing out py
        #print("-->before cut: ddValue1: ", ddValue[1], " ddValue2 :", ddValue[2], " and pSeed.E() ", pSeed.E(), " nseedsNoFit: ", nseedsNoFit, " nMatched: ",nMatched)
        ### high multiplicity
        if(nseedsNoFit > nseedsNoFitMax2):
            ### energy < 4 GeV
            if pSeed.E() < 4.0:
                if ((0.0 < ddValue[1] < 0.12) and (0.0 < ddValue[2] < 0.02)):
                    if ddValue[1] < ddValue1:
                        ddValue0 = ddValue[0]
                        ddValue1 = ddValue[1] 
                        ddValue2 = ddValue[2]
                        iWinner  = i
            ### energy > 4 GeV
            else:
                if ((0.0 < ddValue[1] < 0.08) and (0.0 < ddValue[2] < 0.05)):
                    if ddValue[1] < ddValue1:
                        ddValue0 = ddValue[0]
                        ddValue1 = ddValue[1] 
                        ddValue2 = ddValue[2]
                        iWinner  = i
                        
        #### medium multiplicity
        elif(nseedsNoFitMax < nseedsNoFit <= nseedsNoFitMax2):
            ### energy < 4 GeV
            if(pSeed.E()) < 4:
                if ((0.0 < ddValue[1] < 0.08) and (0.0 < ddValue[2] < 0.05)):
                    if ddValue[1] < ddValue1:
                        ddValue0 = ddValue[0]
                        ddValue1 = ddValue[1] 
                        ddValue2 = ddValue[2]
                        iWinner  = i
            ### energy > 4 GeV
            else:
                if ((0.0 < ddValue[1] < 0.08) and (0.0 < ddValue[2] < 0.06)):
                    if ddValue[1] < ddValue1:
                        ddValue0 = ddValue[0]
                        ddValue1 = ddValue[1] 
                        ddValue2 = ddValue[2]
                        iWinner  = i
                        
        ### low multiplicity
        else:
            ### tight tracks
            if(nMatched == 4):
                ### energy < 4 GeV
                if(pSeed.E() < 4):
                    if(side=="Electron"):
                        #if ((0.0 < ddValue[1] < 0.1) and (0.0 < ddValue[2] < 0.004)):
                        if ((0.0 < ddValue[1] < 0.1) and (0.0 < ddValue[2] < 0.1)):
                            if ddValue[1] < ddValue1:
                                ddValue0 = ddValue[0]
                                ddValue1 = ddValue[1] 
                                ddValue2 = ddValue[2]
                                iWinner  = i
                    else:
                        if ((0.0 < ddValue[1] < 0.1) and (0.0 < ddValue[2] < 0.1)):
                            if ddValue[1] < ddValue1:
                                ddValue0 = ddValue[0]
                                ddValue1 = ddValue[1] 
                                ddValue2 = ddValue[2]
                                iWinner  = i
                            
                ### energy > 4 GeV
                else:
                    if(side=="Electron"):
                        #if ((0.0 < ddValue[1] < 0.04) and (0.0 < ddValue[2] < 0.1)):
                        if ((0.0 < ddValue[1] < 0.1) and (0.0 < ddValue[2] < 0.1)):
                            if ddValue[1] < ddValue1:
                                ddValue0 = ddValue[0]
                                ddValue1 = ddValue[1] 
                                ddValue2 = ddValue[2]
                                iWinner  = i
                    else:
                        if ((0.0 < ddValue[1] < 0.1) and (0.0 < ddValue[2] < 0.1)):
                            if ddValue[1] < ddValue1:
                                ddValue0 = ddValue[0]
                                ddValue1 = ddValue[1] 
                                ddValue2 = ddValue[2]
                                iWinner  = i
                            
            ### loose tracks
            else:
                ### energy < 4 GeV
                if(pSeed.E() < 4.):
                    if ((0.0 < ddValue[1] < 0.2) and (0.0 < ddValue[2] < 0.05)):
                        if ddValue[1] < ddValue1:
                            ddValue0 = ddValue[0]
                            ddValue1 = ddValue[1] 
                            ddValue2 = ddValue[2]
                            iWinner  = i
                ### energy > 4 GeV
                else:
                    if ((0.0 < ddValue[1] < 0.065) and (0.0 < ddValue[2] < 0.065)):
                        if ddValue[1] < ddValue1:
                            ddValue0 = ddValue[0]
                            ddValue1 = ddValue[1] 
                            ddValue2 = ddValue[2]
                            iWinner  = i
                            
        #print("-->after cut: ddValue1: ", ddValue[1], " ddValue2 :", ddValue[2], " and pSeed.E() ", pSeed.E(), " nseedsNoFit: ", nseedsNoFit, " nMatched: ",nMatched)
    
    ### only if a good fit is available, return True
    if iWinner >= 0:
        #print("ddValue0: ", ddValue0, "ddValue1: ",ddValue1, "ddValue2: ", ddValue2)
        return True, {"ddValue0": ddValue0, "ddValue1":ddValue1, "ddValue2": ddValue2, "linepts":allDDList['linepts'][iWinner], "r2Inner":allDDList['i2Inner'][iWinner], "r2Outer":allDDList['i2Outer'][iWinner], "r3Inner":allDDList['i3Inner'][iWinner], "r3Outer":allDDList['i3Outer'][iWinner]}
    else:
        return False, {}



### making the seeds from two tracks from innermost layer and outermost layer, check the seeding cuts in the process
def makeseedNoFit(r1, r4, allR2Inner, allR2Outer, allR3Inner, allR3Outer, side):
    
    if(abs(r1[0]) >= abs(r4[0])):
        return False, {}  ### |x1| must be smaller than |x4|
    
    
    if(r1[0]*r4[0] < 0):
        return False, {}  ### the x value should have the same sign,  i.e. on the same side of the beam
    
    
    if(r1[2] == r4[2]):
        return False, {}  ### if z1=z4..., this is also impossible
    
    
    yDipoleExit = yofz(r1, r4, zDipoleExit)
    xDipoleExit = xofz(r1, r4, zDipoleExit)
    ### The following cuts are coming because of the distribution of the signal
    ### the following cannot be 10.8mm/2 according to the signal tracks
    if(abs(yDipoleExit) > 10.8/2):
        return False, {}
    
    
    ### This is according to the signal x:y at the dipole exit
    if(abs(xDipoleExit) < xDipoleExitMin/2):
        return False, {}  # the track should point to |x|<~1.0 at the dipole exit
    
    
    if(abs(xDipoleExit) > xDipoleExitMax/2):
        return False, {}
    

    ### select only positive x or negative x accroding to the particle
    if( (side=="Positron" and xDipoleExit < 0) or (side=="Electron" and xDipoleExit > 0)):
        return False, {}
    
    
    ### get the momentum of the seed from the r1 and r4
    p, xExit, yExit = getSeedMomentum(r1, r4)
    
    ### checking the track energy
    if(p.E() < EseedMinPrelim or p.E() > EseedMaxPrelim): 
        return False, {}
    
    ### add a 2 GeV cut
    if(p.E() < 2.0):
        return False, {}
    
    nMatched, nExpected, innerR2FromMatching, outerR2FromMatching, innerR3FromMatching, outerR3FromMatching = check_clusters(r1, r4, allR2Inner, allR2Outer, allR3Inner, allR3Outer, side, p.E())
    
    if nMatched < 3:
        return False, {}
    
    return True, {}
    
    

### making the seeds from two tracks from innermost layer and outermost layer, check the seeding cuts in the process
def makeseed(r1, r4, allR2Inner, allR2Outer, allR3Inner, allR3Outer, side, r1GeV, r10GeV, useFit=1, nseedsNoFit=0):
    cutFlowDict['noCut'] += 1
    
    if(abs(r1[0]) >= abs(r4[0])):
        return False, {}  ### |x1| must be smaller than |x4|
    cutFlowDict['x1Gtx4'] += 1
    
    if(r1[0]*r4[0] < 0):
        return False, {}  ### the x value should have the same sign,  i.e. on the same side of the beam
    cutFlowDict['x1*x4Negative'] += 1
    
    if(r1[2] == r4[2]):
        return False, {}  ### if z1=z4..., this is also impossible
    cutFlowDict['z1Eqz4'] += 1
    
    yDipoleExit = yofz(r1, r4, zDipoleExit)
    xDipoleExit = xofz(r1, r4, zDipoleExit)
    ### The following cuts are coming because of the distribution of the signal
    ### the following cannot be 10.8mm/2 according to the signal tracks
    if(abs(yDipoleExit) > 10.8/2):
        return False, {}
    cutFlowDict['yDipoleExitGt5p4'] += 1
    
    ### This is according to the signal x:y at the dipole exit
    if(abs(xDipoleExit) < xDipoleExitMin/2):
        return False, {}  # the track should point to |x|<~1.0 at the dipole exit
    cutFlowDict['xDipoleExitLt25'] += 1
    
    if(abs(xDipoleExit) > xDipoleExitMax/2):
        return False, {}
    cutFlowDict['xDipoleExitGt165'] += 1

    ### select only positive x or negative x accroding to the particle
    if( (side=="Positron" and xDipoleExit < 0) or (side=="Electron" and xDipoleExit > 0)):
        return False, {}
    cutFlowDict['xDipoleExitLt0'] += 1
    
    ### get the momentum of the seed from the r1 and r4
    p, xExit, yExit = getSeedMomentum(r1, r4)
    
    ### checking the track energy
    if(p.E() < EseedMinPrelim or p.E() > EseedMaxPrelim): 
        return False, {}
    cutFlowDict['seedEnergy'] += 1
    
    nMatched, nExpected, innerR2FromMatching, outerR2FromMatching, innerR3FromMatching, outerR3FromMatching = check_clusters(r1, r4, allR2Inner, allR2Outer, allR3Inner, allR3Outer, side, p.E())
    
    if nMatched < 3:
        return False, {}
    
    cutFlowDict['checkClusterTracksMiddleLayers'] += 1
    
    
    passFit, winnerFit = makeSeedFit(r1, r4, nMatched, nExpected, innerR2FromMatching, outerR2FromMatching, innerR3FromMatching, outerR3FromMatching, nseedsNoFit, side)
    
    ### add the nMatched and nExpected to the winnerFit dictionary
    winnerFit.update({'nMatched':nMatched, 'nExpected':nExpected, 'xExit':xExit, 'yExit':yExit, 'passFit':passFit, "pSeedPreFit":p, 'innerR2FromMatching':innerR2FromMatching, 'outerR2FromMatching':outerR2FromMatching, 'innerR3FromMatching':innerR3FromMatching, 'outerR3FromMatching':outerR3FromMatching})
    
    ### decide when we need fit or not
    if((useFit==1) and (not passFit)):
        return False, {}
    
    cutFlowDict['checkClusterFit'] += 1
    
    ### if fit is done, then take the x and z from the fit
    if(passFit):
        xDipoleFromFit = xofz(winnerFit['linepts'][0], winnerFit['linepts'][1], zDipoleExit)*mm2m
        xLayer4FromFit = xofz(winnerFit['linepts'][0], winnerFit['linepts'][1], r4[2])*mm2m
        
        yDipoleFromFit = yofz(winnerFit['linepts'][0], winnerFit['linepts'][1], zDipoleExit)*mm2m
        yLayer4FromFit = yofz(winnerFit['linepts'][0], winnerFit['linepts'][1], r4[2])*mm2m

        ### here r1 and r4 must be in mm
        R1 = [xDipoleFromFit*m2mm, yDipoleFromFit*m2mm, zDipoleExit]
        R4 = [xLayer4FromFit*m2mm, yLayer4FromFit*m2mm, r4[2]]
        
        #### re-evaluate the seed energy from the best fit track
        pSeed, xExit, yExit = getSeedMomentum(R1, R4)
        
        ### checking the track energy
        if(pSeed.E() < EseedMin or pSeed.E() > EseedMax): 
            return False, {}
        #print("Energy: ", pSeed.E())
        cutFlowDict['trackEnergy'] += 1
        
        
        ### update the dictionary
        winnerFit.update({"pSeedPostFit":pSeed})
        winnerFit["xExit"] = xExit
        winnerFit["yExit"] = yExit
        
        #### calculate the distance from x4:x_exit analytical line
        d              = Distance(xDipoleFromFit,xLayer4FromFit,r1GeV,r10GeV)
        #print("The cluster x distance: ",d)
        histos['hSeedDistanceGlobal'].Fill(d)
        if(nMatched==4):
            histos['hSeedDistanceTightGlobal'].Fill(d)
            histos['hSeedDistanceTightVsEnergyGlobal'].Fill(pSeed.E(),d)
        else:

            histos['hSeedDistanceLooseGlobal'].Fill(d)
        
        
        #print("nseedsNoFit: ", nseedsNoFit, " pSeed.E(): ", pSeed.E(), " d: ", d, " pSeed.Py(): ", pSeed.Py())
        
        ### medium and high multiplicity
        if(nseedsNoFit > nseedsNoFitMax2):
            if(pSeed.E() < 4.0):
                if(d > 10*mm2m):
                    return False, {}
            else:
                if(d > 5*mm2m):
                    return False, {}
            
        ### medium multiplicity
        elif(nseedsNoFitMax < nseedsNoFit <= nseedsNoFitMax2):
            if(nMatched == 4):
                if(pSeed.E() < 4.0):
                    if(d > 2*mm2m):
                        return False, {}
                else:
                    if(d > 6*mm2m):
                        return False, {}
            else:
                if(pSeed.E() < 4.0):
                    if(d > 7*mm2m):
                        return False, {}
                else:
                    if(d > 5*mm2m):
                        return False, {}
        ### low multiplicity
        else:
            ### tight tracks
            if(nMatched == 4):
                if(pSeed.E() < 4.0):
                    if(d > 6*mm2m):
                        return False, {}
                else:
                    if(d > 4*mm2m):
                        return False, {}
            ### loose tracks
            else:
                if(pSeed.E() < 4.0):
                    if(d > 5*mm2m):
                        return False, {}
                else:
                    if(d > 5*mm2m):
                        return False, {}
                
        cutFlowDict['checkClusterXDistance'] += 1
        winnerFit.update({"distance":d})
        
        
        
        
        
        ### checking the track pY
        if(nseedsNoFit > nseedsNoFitMax):
            if(abs(pSeed.Py()) > PseedMax*2): 
                return False, {}
        else:
            if(pSeed.E() < 4):
                if(nMatched == 4):
                    if(abs(pSeed.Py()) > PseedMax*2): 
                        return False, {}
                else:
                    if(abs(pSeed.Py()) > PseedMax*1.9): 
                        return False, {}
            else:
                if(nMatched == 4):
                    if(abs(pSeed.Py()) > PseedMax*7): 
                        return False, {}
                else:
                    if(abs(pSeed.Py()) > PseedMax*6.2): 
                        return False, {}
            
        cutFlowDict['checkClusterTrackPy'] += 1
        
        #print("$$$$: nseedsNoFit: ", nseedsNoFit, " pSeed.E(): ", pSeed.E(), " pSeed.Py(): ", pSeed.Py(), " d: ", d )
        #print("ddValue0: ", winnerFit["ddValue0"], " ddValue1: ", winnerFit["ddValue1"], " ddValue2: ", winnerFit["ddValue2"])
        
        
        if(winnerFit['nMatched'] == 4):
            cutFlowDict['checkClusterTrackPyTight'] += 1
        else:
            cutFlowDict['checkClusterTrackPyLoose'] += 1

    return True, winnerFit



#### add diagonal vertex cut
#def checkVtxCut(vtxx, vtxz):
    #if((-45.0 < vtxx <= -25.5 ) and (3800 < vtxz <= 4200)):
        #return True
    #elif((-55.0 < vtxx <= -35.0) and (4200 < vtxz <= 4500)):
        #return True
    #elif((-60.0 < vtxx <= -45.0) and (4500 < vtxz <= 5000)):
        #return True
    #elif((-70 < vtxx <= -50.0) and (5000 < vtxz <= 5500)):
        #return True
    #elif((-75 < vtxx <= -60.0) and (5500 < vtxz <= 6000)):
        #return True 
    #elif((-80 < vtxx <= -65.0) and (5500 < vtxz <= 6000)):
        #return True
    #elif((-90 < vtxx <= -70.0) and (6000 < vtxz <= 6500)):
        #return True
    #elif((-100.0 < vtxx <= -80.0) and (6500 < vtxz <= 7000)):
        #return True
    #else:
        #return False
    
### anything below x, remove it cut
def checkVtxCut(vtxx, vtxz):
    if((vtxx <= -25.5 ) and (3800 < vtxz <= 4200)):
        return True
    elif((vtxx <= -35.0) and (4200 < vtxz <= 4500)):
        return True
    elif((vtxx <= -45.0) and (4500 < vtxz <= 5000)):
        return True
    elif((vtxx <= -50.0) and (5000 < vtxz <= 5500)):
        return True
    elif((vtxx <= -60.0) and (5500 < vtxz <= 6000)):
        return True 
    elif((vtxx <= -65.0) and (5500 < vtxz <= 6000)):
        return True
    elif((vtxx <= -70.0) and (6000 < vtxz <= 6500)):
        return True
    elif((vtxx <= -80.0) and (6500 < vtxz <= 7000)):
        return True
    else:
        return False

### The main function
def main():
    # give the input text file containing all the track information
    parser = argparse.ArgumentParser(description='Code to find seed tracks')
    parser.add_argument('-l', action="store", dest="inFile", default="list_root_hics_165gev_w0_5000nm_provisional_10xi3_6cee466a_trackInfoClean.txt")
    parser.add_argument('-s', action="store_true", dest="needSignal", default=False)
    parser.add_argument('-f', action="store", dest="needFit", type=int, default=1)
    parser.add_argument('-e', action="store", dest="eCut", type=float, default=0.0)
    parser.add_argument('-p', action="store", dest="Particle", type=str, default="Positron")
    args = parser.parse_args()
    
    
    for index, row in df.iterrows():
        detid   = row["detid"]
        layerid = row["layerid"]
        xMin, xMax = GetSensorXBoundaries(detid, layerid)
        xBoundaries.update({"detid"+str(detid)+"_layerid"+str(layerid):[xMin, xMax]})
    
    
    ### open the file containing track file
    inTextFile      = args.inFile
    inputTrackInfo  = open(inTextFile)
    energyCutSuffix = "VariableEnergyCut"
    side            = args.Particle
    
    signalCutSuffix = ""
    if(args.needSignal):
        signalCutSuffix = "OnlySignal"
    else:
        signalCutSuffix = "SignalAndBackground"
        
    particleSuffix = ""
    if(side=="Positron"):
        particleSuffix="PositronSide"
    else:
        particleSuffix="ElectronSide"
        
        
    ### open histogram to know the seed information
    if (("_" in inTextFile) and ("hics" in inTextFile)):
        eachName            = inTextFile.split('.')[0].split('_')
        suffixName          = "_".join(eachName[2:])
    else:
        suffixName          = inTextFile.split('.')[0]
    
    
    outFile                   = TFile("seedingInformation_"+suffixName+"_"+energyCutSuffix+"_"+signalCutSuffix+"_"+particleSuffix+".root", "RECREATE")
    outFile.cd()
    
    hAllPossible              = TH1D("hAllPossible", "all possible track combination; bunch crossing; number of track combination", 9508, 0, 9508)
    hSeedPossible             = TH1D("hSeedPossible", "seed track combination; bunch crossing; number of seed track", 9508, 0, 9508)
    hSeedMultiplicity         = TH1D("hSeedMultiplicity", "seed multiplicity; number of seeds; BX", 2000, 0, 2000)
    hSeedMultiplicityPrelim   = TH1D("hSeedMultiplicityPrelim", "seed multiplicity prelim; number of seeds; BX", 2000, 0, 2000)
    hSignalMultiplicity       = TH1D("hSignalMultiplicity", "number of signals; number of signals; BX", 2000, 0, 2000)
    hXLayer4XDipole           = TH2D("hXLayer4XDipole", "number of signal distribution; x_{Dipole} [mm]; x_{Layer4}",330,0,330,650,0,650)
    hXExitYExit               = TH2D("hXExitYExit", "track distribution; x_{Exit} [mm]; y_{Exit} [mm]",650,0,650,40,-10,10)
    hSigEnergy                = TH1D("hSigEnergy", "signal energy; Energy [GeV]; Entries", 200, 0, 20)
    
    hSeedEnergy               = TH1D("hSeedEnergy", "seed energy; Energy [GeV]; Entries", 200, 0, 20)
    hSeedEnergyLoose          = TH1D("hSeedEnergyLoose", "seed energy (Loose); Energy [GeV]; Entries", 200, 0, 20)
    hSeedEnergyTight          = TH1D("hSeedEnergyTight", "seed energy (Tight); Energy [GeV]; Entries", 200, 0, 20)
    hSeedPy                   = TH1D("hSeedPy", "p'_{Y}; p'_{Y} [GeV]; Entries", 100, -0.02, 0.02)
    hSeedPyLoose              = TH1D("hSeedPyLoose", "p'_{Y} (Loose); p'_{Y} [GeV]; Entries", 100, -0.02, 0.02)
    hSeedPyTight              = TH1D("hSeedPyTight", "p'_{Y} (Tight); p'_{Y} [GeV]; Entries", 100, -0.02, 0.02)
    
    hSeedDistance             = TH1D("hSeedDistance", "seed distance wrt analytical line; d [m]; Entries", 200, 0, 0.005)
    hSeedDistanceLoose        = TH1D("hSeedDistanceLoose", "seed distance wrt analytical line (loose); d [m]; Entries", 200, 0, 0.005)
    hSeedDistanceTight        = TH1D("hSeedDistanceTight", "seed distance wrt analytical line (tight); d [m]; Entries", 200, 0, 0.005)
    hSVDValues1               = TH1D("hSVDValues1", "output of SVD[0]; fit quality [0]; Events", 50, 200, 250)
    hSVDValues2               = TH1D("hSVDValues2", "output of SVD[1]; fit quality [1]; Events", 50, 0, 0.1)
    hSVDValues3               = TH1D("hSVDValues3", "output of SVD[2]; fit quality [2]; Events", 25, 0, 0.01)
    hSVDValues1Tight          = TH1D("hSVDValues1Tight", "output of SVD[0]; fit quality [0]; Events", 50, 200, 250)
    hSVDValues2Tight          = TH1D("hSVDValues2Tight", "output of SVD[1], tight; fit quality [1]; Events", 50, 0, 0.1)
    hSVDValues3Tight          = TH1D("hSVDValues3Tight", "output of SVD[2], tight; fit quality [2]; Events", 25, 0, 0.01)
    hSVDValues1Loose          = TH1D("hSVDValues1Loose", "output of SVD[0]; fit quality [0]; Events", 50, 200, 250)
    hSVDValues2Loose          = TH1D("hSVDValues2Loose", "output of SVD[1], loose; fit quality [1]; Events", 50, 0, 0.1)
    hSVDValues3Loose          = TH1D("hSVDValues3Loose", "output of SVD[2], loose; fit quality [2]; Events", 25, 0, 0.01)
    
    

    #### select the number of bunch crossing for different files, BX is useful only for hics setup, 494 BX for 5000 nm new signal, 9508 BX for 3000 nm old signal
    if 'hics' in inTextFile:
        nBX          = 494
        checkBXMatch = True
    ### g-laser signal
    elif 'bppp_165gev_w0_3000nm' in inTextFile:
        nBX          = 1001
        checkBXMatch = True
    elif 'bppp_165gev_w0_5000nm' in inTextFile:
        nBX          = 1981
        checkBXMatch = True
    elif 'bppp_165gev_w0_8000nm' in inTextFile:
        nBX          = 1943
        checkBXMatch = True
    ### e-beam only background
    elif 'EBeamOnlyNewSamples' in inTextFile:
        nBX          = 160
        checkBXMatch = True
    elif 'ePlusLaserBkgNewSamples' in inTextFile:
        nBX          = 139
        checkBXMatch = True
    ### g+laser background
    elif 'gPlusLaserBkgNewSamples' in inTextFile:
        nBX          = 129
        checkBXMatch = True
    ### all other cases
    else:
        nBX          = 1
        checkBXMatch = False
        
    ### get the tracks given the momentum, we get 1 GeV and 10 GeV momentum tracks  
    r1GeV, r10GeV = getR1R2(1, 10)
       
    # all the track info is in the following list
    position = []
    # get the information from the text files
    for lines in inputTrackInfo.readlines():
        lines = lines.rstrip()
        if "#" in lines: continue
        eachWord  = lines.split()
        bxNumber  = int(eachWord[0])
        trackId   = int(eachWord[2])
        pdgId     = int(eachWord[1])
        staveId   = int(eachWord[3])-1000
        xPos      = float(eachWord[4])
        yPos      = float(eachWord[5])
        energyVal = float(eachWord[6])
        weight    = float(eachWord[7])
        vtx_x     = float(eachWord[8])
        vtx_y     = float(eachWord[9])
        vtx_z     = float(eachWord[10])
        ### here the preliminary selection of tracks can be done based on trackid and pdgid
        ### do not want photons, as they will not leave track
        if(pdgId==22): continue
        ### if needed, select only the signal  
        if( args.needSignal and ( (side=="Positron" and not(pdgId==-11 and trackId==1)) or (side=="Electron" and not(pdgId==11 and trackId==1)) ) ):
            continue
        failVtxCut = checkVtxCut(vtx_x, vtx_z)
        ### required to suppress unwanted background
        #if(vtx_x < -25.5 and (vtx_z > 3600 and vtx_z < 4600)): continue
        #if(failVtxCut): continue
        ### ask for variable energy cut, energyAfterCut is in keV
        energyAfterCut = energyAbsorbed(staveId, energyVal, vtx_z)
        #if(float(eachWord[6]) > args.eCut*1e-6):
        #### taking tracks with more than 2 GeV of energy
        if(energyAfterCut > 0):
            #### bxNumber, trackId, staveId, x, y, E, weight
            position.append([bxNumber, trackId, staveId, xPos, yPos, energyVal, weight, vtx_x, vtx_y, vtx_z])

    allSeedsTrackLines = []
    #### run the seeding algorithm per BX
    for bxCounter in range(1, nBX+1):
        # separate each bx now
        #if(bxCounter!=58):continue
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
            if(side=="Positron"):
                if (values[2] == 0):
                    allR1Inner.append([values[3], values[4], z1inner, values[5], values[6], values[7], values[8], values[9]])
                elif (values[2] == 1):
                    allR1Outer.append([values[3], values[4], z1outer, values[5], values[6], values[7], values[8], values[9]])
                elif (values[2] == 2):
                    allR2Inner.append([values[3], values[4], z2inner, values[5], values[6], values[7], values[8], values[9]])
                elif (values[2] == 3):
                    allR2Outer.append([values[3], values[4], z2outer, values[5], values[6], values[7], values[8], values[9]])
                elif (values[2] == 4):
                    allR3Inner.append([values[3], values[4], z3inner, values[5], values[6], values[7], values[8], values[9]])
                elif (values[2] == 5):
                    allR3Outer.append([values[3], values[4], z3outer, values[5], values[6], values[7], values[8], values[9]])
                elif (values[2] == 6):
                    allR4Inner.append([values[3], values[4], z4inner, values[5], values[6], values[7], values[8], values[9]])
                elif (values[2] == 7):
                    allR4Outer.append([values[3], values[4], z4outer, values[5], values[6], values[7], values[8], values[9]])
                else:
                    continue
            else:
                if (values[2] == 8):
                    allR1Inner.append([values[3], values[4], z1inner, values[5], values[6], values[7], values[8], values[9]])
                elif (values[2] == 9):
                    allR1Outer.append([values[3], values[4], z1outer, values[5], values[6], values[7], values[8], values[9]])
                elif (values[2] == 10):
                    allR2Inner.append([values[3], values[4], z2inner, values[5], values[6], values[7], values[8], values[9]])
                elif (values[2] == 11):
                    allR2Outer.append([values[3], values[4], z2outer, values[5], values[6], values[7], values[8], values[9]])
                elif (values[2] == 12):
                    allR3Inner.append([values[3], values[4], z3inner, values[5], values[6], values[7], values[8], values[9]])
                elif (values[2] == 13):
                    allR3Outer.append([values[3], values[4], z3outer, values[5], values[6], values[7], values[8], values[9]])
                elif (values[2] == 14):
                    allR4Inner.append([values[3], values[4], z4inner, values[5], values[6], values[7], values[8], values[9]])
                elif (values[2] == 15):
                    allR4Outer.append([values[3], values[4], z4outer, values[5], values[6], values[7], values[8], values[9]])
                else:
                    continue
                
        
                
        
        ### removing the overlap region of inner and outer stave
        allR1Unique = allR1Inner
        allR4Unique = allR4Outer
        
        for r1Out in allR1Outer:
            #### remove all points having an x overlap with inner stave: not 100% effective
            if(side=="Positron" and r1Out[0] > (xInnerStaveLastChip + xSizeInnerStaveLastChip/2.)):  ### the x position of last chip on inner stave layer 1+half of the x size 
                allR1Unique.append(r1Out)
            elif(side=="Electron" and r1Out[0] < -(xInnerStaveLastChip + xSizeInnerStaveLastChip/2.)):
                allR1Unique.append(r1Out)
            else:
                continue
                
            
        for r1 in allR1Unique:
            if(r1[3]>2.0): hSigEnergy.Fill(r1[3])
        
        ### count seeds more than 2 GeV 
        trueSignalMoreThan2GeV = 0
        for r1 in allR1Unique:
            if r1[3] >= 2.0:
                trueSignalMoreThan2GeV += 1
            
        hSignalMultiplicity.Fill(trueSignalMoreThan2GeV)
        
        for r4In in allR4Inner:
            #### remove all points having an x overlap with outer stave: not 100% effective
            if(side=="Positron" and r4In[0] < (xOuterStaveFirstChip - xSizeOuterStaveFirstChip/2.)):   ### the x position of last chip on outer stave layer 4-half of the x size 
                allR4Unique.append(r4In)
            elif(side=="Electron" and r4In[0] > -(xOuterStaveFirstChip - xSizeOuterStaveFirstChip/2.)):
                allR4Unique.append(r4In)
            else:
                continue
                
        #### loop over points in r4:
        counter    = 0
        allCounter = 0
        print("For the BX: ", bxCounter," the number of layer 1 tracks: ", len(allR1Unique)," and layer 4 tracks: ", len(allR4Unique))
        
        
        ### don't need to run the seeding algorithm if it is signal
        if(args.needSignal):
            print("+++++  actual signals                 : ", len(allR1Unique), " +++++") 
            print("+++++  actual signals after energy cut: ", trueSignalMoreThan2GeV, " +++++")
            
        seedCounterNoFit = 0
        for r4 in allR4Unique[:]:
            for r1 in allR1Unique[:]:
                seed, winnerDict = makeseedNoFit(r1, r4, allR2Inner, allR2Outer, allR3Inner, allR3Outer, side)
                if(seed):
                    seedCounterNoFit += 1
        
        
        
        for r4 in allR4Unique[:]:
            for r1 in allR1Unique[:]:
                ### This is the return of makeseed function
                #True, {"linepts":allDDList['linepts'][iWinner], "r2Inner":allDDList['i2Inner'][iWinner], "r2Outer":allDDList['i2Outer'][iWinner], "r3Inner":allDDList['i3Inner'][iWinner], "r3Outer":allDDList['i3Outer'][iWinner], "pSeed":p.E(), "nMatched": nMatched, "nExpected":nExpected}
                
                seed, winnerDict = makeseed(r1, r4, allR2Inner, allR2Outer, allR3Inner, allR3Outer, side, r1GeV, r10GeV, args.needFit, seedCounterNoFit)
                allCounter += 1
                
                if(seed):
                    ### values which can be obtained even before fitting, but with prefit results
                    hXExitYExit.Fill(winnerDict['xExit'],winnerDict['yExit'])
                    ### select the pSeed depending on the fitting
                    if winnerDict['passFit']:
                        pSeedName = "pSeedPostFit"
                    else:
                        pSeedName = "pSeedPreFit"
                    
                    hSeedEnergy.Fill(winnerDict[pSeedName].E())
                    hSeedPy.Fill(winnerDict[pSeedName].Py())
                    if(winnerDict['nMatched'] == 4):
                        hSeedEnergyTight.Fill(winnerDict[pSeedName].E())
                        hSeedPyTight.Fill(winnerDict[pSeedName].Py())
                    else:
                        hSeedEnergyLoose.Fill(winnerDict[pSeedName].E())
                        hSeedPyLoose.Fill(winnerDict[pSeedName].Py())
                    
                    ### values which can be obtained only after fitting
                    if winnerDict['passFit']:
                        hSeedDistance.Fill(winnerDict['distance'])
                        if(winnerDict['nMatched'] == 4 or winnerDict['nMatched'] == 3):
                            allSeedsTrackLines.append(GetExtendedTrackLine([r1, r4]))
                        xDipoleFromFit = xofz(winnerDict['linepts'][0], winnerDict['linepts'][1], zDipoleExit)
                        xLayer4FromFit = xofz(winnerDict['linepts'][0], winnerDict['linepts'][1], r4[2])
                        hXLayer4XDipole.Fill(xDipoleFromFit, xLayer4FromFit)
                        if(False and winnerDict['nMatched'] >= 3):
                            
                            print("nMatched: ", winnerDict['nMatched'], " nExpected: ", winnerDict['nExpected'])
                            print("vertex positions first layer: xpos ", r1[0], " ypos ", r1[1], " energy ", r1[3], " vtx_x ", r1[5], " vtx_y ", r1[6], " vtx_z ", r1[7])
                            if(len(winnerDict['innerR2FromMatching']) >=1 and len(winnerDict['innerR2FromMatching'][0])>=8):
                                print("vertex positions second inner layer: xpos ",winnerDict['innerR2FromMatching'][0][0], " ypos ", winnerDict['innerR2FromMatching'][0][1]," energy ", winnerDict['innerR2FromMatching'][0][3]," vtx_x ", winnerDict['innerR2FromMatching'][0][5], " vtx_y ", winnerDict['innerR2FromMatching'][0][6], " vtx_z ", winnerDict['innerR2FromMatching'][0][7])
                                
                            if(len(winnerDict['outerR2FromMatching']) >=1 and len(winnerDict['outerR2FromMatching'][0])>=8):
                                print("vertex positions second outer layer: xpos ",winnerDict['outerR2FromMatching'][0][0], " ypos ", winnerDict['outerR2FromMatching'][0][1], " energy ", winnerDict['outerR2FromMatching'][0][3]," vtx_x ", winnerDict['outerR2FromMatching'][0][5], " vtx_y ", winnerDict['outerR2FromMatching'][0][6], " vtx_z ", winnerDict['outerR2FromMatching'][0][7])
                             
                            if(len(winnerDict['innerR3FromMatching']) >=1 and len(winnerDict['innerR3FromMatching'][0])>=8):
                                print("vertex positions third inner layer: xpos ", winnerDict['innerR3FromMatching'][0][0], " ypos ", winnerDict['innerR3FromMatching'][0][1], " energy ", winnerDict['innerR3FromMatching'][0][3]," vtx_x ", winnerDict['innerR3FromMatching'][0][5], " vtx_y ", winnerDict['innerR3FromMatching'][0][6], " vtx_z ", winnerDict['innerR3FromMatching'][0][7])
                            if(len(winnerDict['outerR3FromMatching']) >=1 and len(winnerDict['outerR3FromMatching'][0])>=8):
                                print("vertex positions third outer layer: xpos ",winnerDict['outerR3FromMatching'][0][0], " ypos ", winnerDict['outerR3FromMatching'][0][1], " energy ", winnerDict['outerR3FromMatching'][0][3]," vtx_x ", winnerDict['outerR3FromMatching'][0][5], " vtx_y ", winnerDict['outerR3FromMatching'][0][6], " vtx_z ", winnerDict['outerR3FromMatching'][0][7])
                             
                             
                            print("vertex positions last layer: xpos ", r4[0], " ypos ", r4[1], " energy ", r4[3], " vtx_x ", r4[5], " vtx_y ", r4[6], " vtx_z ", r4[7])
                            
                            
                        hSVDValues1.Fill(winnerDict["ddValue0"])
                        hSVDValues2.Fill(winnerDict["ddValue1"])
                        hSVDValues3.Fill(winnerDict["ddValue2"])
                        if(winnerDict['nMatched'] == 4):
                            hSeedDistanceTight.Fill(winnerDict['distance'])
                            hSVDValues1Tight.Fill(winnerDict["ddValue0"])
                            hSVDValues2Tight.Fill(winnerDict["ddValue1"])
                            hSVDValues3Tight.Fill(winnerDict["ddValue2"])
                        else:
                            hSeedDistanceLoose.Fill(winnerDict['distance'])
                            hSVDValues1Loose.Fill(winnerDict["ddValue0"])
                            hSVDValues2Loose.Fill(winnerDict["ddValue1"])
                            hSVDValues3Loose.Fill(winnerDict["ddValue2"])
                        counter += 1

        print("Number of seed: ", counter)
        print("Number of combination: ", allCounter)
        print("Preliminary multiplicity: ",seedCounterNoFit)
        hSeedMultiplicityPrelim.Fill(seedCounterNoFit)
        hAllPossible.SetBinContent(bxCounter, allCounter)
        hSeedPossible.SetBinContent(bxCounter, counter)
        hSeedMultiplicity.Fill(counter)
     
    outFile.cd()
    for hName, hist in histos.items():
        hist.Write()
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
    cnv.SaveAs("trackingDiagram_"+suffixName+"_"+energyCutSuffix+"_"+signalCutSuffix+"_"+particleSuffix+".pdf")
    cnv.SaveAs("trackingDiagram_"+suffixName+"_"+energyCutSuffix+"_"+signalCutSuffix+"_"+particleSuffix+".root")
    pprint.pprint(cutFlowDict)


#### call the main function
if __name__ == "__main__":
    start = time.time()
    main()
    print("-------- The processing time: ",time.time() - start, " s")
