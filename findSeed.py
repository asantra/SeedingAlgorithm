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
EseedMin = 0.5 # GeV
EseedMax = 17.5 # GeV




#//MeV mass of electron/positron
meMeV    = 0.5109989461 
meGeV    = meMeV/1000.
### me^2 in GeV
meGeV2   = meGeV*meGeV

zPositionListInner = [z2inner, z3inner]
zPositionListOuter = [z2outer, z3outer]
ddList = []



#### This is the cutflow dictionary
cutFlowDict = OrderedDict(
    [('noCut', 0), 
     ('x1Gtx4',0), 
     ('x1*x4Negative',0), 
     ('z1Eqz4', 0), 
     ('yDipoleExitGt15',0), 
     ('xDipoleExitLt25',0), 
     ('xDipoleExitGt165', 0), 
     ('xDipoleExitLt0',0), 
     ('trackEnergy', 0),
     ('checkClusterTracksMiddleLayers', 0),
     ('checkClusterFit', 0),
     ('checkClusterXDistance', 0)]
    )

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
    
    
        
  
def getR(p):
   R = p/(0.3*B)
   return R

def getR2(p):
   R = getR(p)
   return R*R

def getXTangent(ZT,p):
   R2 = getR2(p)
   XT = (R2/LB-ZT)*(LB/math.sqrt(R2-LB2))
   return XT

def getx(ZT,p):
   XT = getXTangent(ZT,p)
   R = getR(p)
   x = R-XT
   return x

def getR1R2(p1, p2):
    r1    = [getx(Z0,p1), getx(Z4,p1)]
    r2    = [getx(Z0,p2), getx(Z4, p2)]
    #print("r1: ",r1, " r2: ",r2)
    #quit()
    return r1, r2

def Distance(x0Test,x4Test,r1,r2):
   ### general straight line equation
   # a*x0 + b*x4 + c = 0 --> (a/b)*x0 + x4 + (c/b) = 0
   # (a/b)*r1[0] + r1[1] + (c/b) = 0
   # (a/b)*r2[0] + r2[1] + (c/b) = 0  -->  (a/b)*(r1[0]-r2[0]) + (r1[1]-r2[1]) = 0 --> a/b = -(r1[1]-r2[1])/(r1[0]-r2[0])
   # -(r1[1]-r2[1])/(r1[0]-r2[0])*r1[0] + r1[1] + (c/b) = 0 --> c/b = (r1[1]-r2[1])/(r1[0]-r2[0])*r1[0] - r1[1]
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


### use svd fit
def seed3dfitSVD3Points(r1,r2,r3, name="", dodraw=False):
   g = TGraph2D()
   g.SetMarkerSize(3)
   g.SetMarkerStyle(20)
   g.SetMarkerColor(ROOT.kRed)
   g.SetLineColor(ROOT.kRed)
   g.SetPoint(0,r1[0],r1[1],r1[2])
   g.SetPoint(1,r2[0],r2[1],r2[2])
   g.SetPoint(2,r3[0],r3[1],r3[2])
   g.GetXaxis().SetRangeUser(0, 650)
   g.GetYaxis().SetRangeUser(-10,+10)
   g.GetZaxis().SetRangeUser(3900, 4400) 
   
    
   x = np.array([r1[0],r2[0],r3[0]])
   y = np.array([r1[1],r2[1],r3[1]])
   z = np.array([r1[2],r2[2],r3[2]])

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
   
   if(dodraw):
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
      cnv.SaveAs(name+".pdf")
   
   return linepts, dd ## dd is a 1D array of the data singular values
   
    



### use svd fit
def seed3dfitSVD(r1,r2,r3,r4, name="", dodraw=False):
   g = TGraph2D()
   g.SetMarkerSize(3)
   g.SetMarkerStyle(20)
   g.SetMarkerColor(ROOT.kRed)
   g.SetLineColor(ROOT.kRed)
   g.SetPoint(0,r1[0],r1[1],r1[2])
   g.SetPoint(1,r2[0],r2[1],r2[2])
   g.SetPoint(2,r3[0],r3[1],r3[2])
   g.SetPoint(3,r4[0],r4[1],r4[2])
   g.GetXaxis().SetRangeUser(0, 650)
   g.GetYaxis().SetRangeUser(-10,+10)
   g.GetZaxis().SetRangeUser(3900, 4400) 
   
    
   x = np.array([r1[0],r2[0],r3[0],r4[0]])
   y = np.array([r1[1],r2[1],r3[1],r4[1]])
   z = np.array([r1[2],r2[2],r3[2],r4[2]])

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
   
   if(dodraw):
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
      cnv.SaveAs(name+".pdf")
   
   return linepts, dd ## dd is a 1D array of the data singular values
   
   
   
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
   
def getExpectedHits(r1, r4):
    expectedHits2Inner = 0
    expectedHits2Outer = 0
    
    expectedHits3Inner = 0
    expectedHits3Outer = 0
    
    yAbsMargins    = 0.13 # mm (a "road" of 200 microns around the line between r4 and r1)
    xAbsMargins    = 0.13 # mm (a "road" of 200 microns around the line between r4 and r1)
    
    r1min          = [r1[0]-xAbsMargins, r1[1]-yAbsMargins, r1[2]]
    r1max          = [r1[0]+xAbsMargins, r1[1]+yAbsMargins, r1[2]]
    r4min          = [r4[0]-xAbsMargins, r4[1]-yAbsMargins, r4[2]]
    r4max          = [r4[0]+xAbsMargins, r4[1]+yAbsMargins, r4[2]]
    
    # check possible clusters in layer 2, for both inner and outer stave
    x2minInner = xofz(r1min, r4min, z2inner)
    x2maxInner = xofz(r1max, r4max, z2inner)
    
    x2minOuter = xofz(r1min, r4min, z2outer)
    x2maxOuter = xofz(r1max, r4max, z2outer)
    
    #/// check possible clusters in layer 3, for both inner and outer stave
    x3minInner = xofz(r1min, r4min, z3inner)
    x3maxInner = xofz(r1max, r4max, z3inner)

    x3minOuter = xofz(r1min, r4min, z3outer)
    x3maxOuter = xofz(r1max, r4max, z3outer)
    
    
    x2Inner = xofz(r1, r4, z2inner)
    x2Outer = xofz(r1, r4, z2outer)
    x3Inner = xofz(r1, r4, z3inner)
    x3Outer = xofz(r1, r4, z3outer)
    
    for chipName, boundaries in xBoundaries.items(): 
        if "layerid2" in chipName and x2Inner > boundaries[0] and x2Inner < boundaries[1]:
            expectedHits2Inner += 1
        if "layerid3" in chipName and x2Outer > boundaries[0] and x2Outer < boundaries[1]:
            expectedHits2Outer += 1
        if "layerid4" in chipName and x3Inner > boundaries[0] and x3Inner < boundaries[1]:
            expectedHits3Inner += 1
        if "layerid5" in chipName and x3Outer > boundaries[0] and x3Outer < boundaries[1]:
            expectedHits3Outer += 1
    
    
    return expectedHits2Inner, expectedHits2Outer, expectedHits3Inner, expectedHits3Outer
        
    
    

### check if there are hits along one track in layer 2 and layer 3
def check_clusters(r1, r4, allR2Inner, allR2Outer, allR3Inner, allR3Outer, side):
    global ddList
    expectedHits2Inner, expectedHits2Outer, expectedHits3Inner, expectedHits3Outer = getExpectedHits(r1, r4)
    #print(expectedHits2Inner, expectedHits2Outer,  expectedHits3Inner, expectedHits3Outer)
    hitsOnRoad2Inner = 0   ### how many tracks accepted in the road along the r1 and r4
    hitsOnRoad2Outer = 0
    hitsOnRoad3Inner = 0   ### how many tracks accepted in the road along the r1 and r4
    hitsOnRoad3Outer = 0
    yAbsMargins    = 0.13 # mm (a "road" of 200 microns around the line between r4 and r1)
    xAbsMargins    = 0.13 # mm (a "road" of 200 microns around the line between r4 and r1)
    
    r1min          = [r1[0]-xAbsMargins, r1[1]-yAbsMargins, r1[2]]
    r1max          = [r1[0]+xAbsMargins, r1[1]+yAbsMargins, r1[2]]
    r4min          = [r4[0]-xAbsMargins, r4[1]-yAbsMargins, r4[2]]
    r4max          = [r4[0]+xAbsMargins, r4[1]+yAbsMargins, r4[2]]


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
        #break
    
    ### separate out staves for allR2
    ### only if there is no suitable track from layer 2 inner stave, go to layer 2 outer stave
    #if(not accept2Inner):
    if True:
        for i2 in range(0, len(allR2Outer)):
            ### remember allR2 has x, y, z and E saved in the list
            accept2yzOuter = ( (allR2Outer[i2][1] >= y2minOuter and allR2Outer[i2][1] <= y2maxOuter) )
            if(not accept2yzOuter): continue
        
            accept2xzOuter = ( ( allR2Outer[i2][0] >= x2minOuter and allR2Outer[i2][0] <= x2maxOuter ) )
            if(not accept2xzOuter): continue
            accept2Outer      = True
            hitsOnRoad2Outer += 1
            outerR2FromMatching.append(allR2Outer[i2]) ### collect all matched tracks from layer 2, remove the break
            #break
    
    ##### return false if both of the inner and outer acceptances are false
    #if(not (accept2Inner or accept2Outer)):
        #return False
    
    
    
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
    ### only if no track found in layer 3 inner stave, go to layer 3 outer stave
    #if(not accept3Inner):
    ### go to layer 3
    if True:
        for i3 in range(0, len(allR3Outer)):
            ### remember allR3 has x, y, z and E saved in the list
            accept3yzOuter = ( (allR3Outer[i3][1] >= y3minOuter and allR3Outer[i3][1] <= y3maxOuter) )
            if(not accept3yzOuter): continue
        
            accept3xzOuter = ( ( allR3Outer[i3][0] >= x3minOuter and allR3Outer[i3][0] <= x3maxOuter ) )
            if(not accept3xzOuter): continue
            accept3Outer        = True
            hitsOnRoad3Outer   += 1
            outerR3FromMatching.append(allR3Outer[i3]) ### collect all matched tracks from layer 3, remove the break
            #break
    
    ##### return false if both of the inner and outer acceptances are false
    #if(not (accept3Inner or accept3Outer)):
        #return False
    
    
    #### return false if both of the inner and outer acceptances are false
    #if(not (accept2Inner or accept2Outer or accept3Inner or accept3Outer)):
        #return False
    
    nMatched = (expectedHits2Inner<=hitsOnRoad2Inner)+(expectedHits2Outer<=hitsOnRoad2Outer)+(expectedHits3Inner<=hitsOnRoad3Inner)+(expectedHits3Outer<=hitsOnRoad3Outer)
    nExpected = expectedHits2Inner + expectedHits2Outer + expectedHits3Inner + expectedHits3Outer
    
    if(nMatched>=3):
        return nMatched, nExpected, innerR2FromMatching, outerR2FromMatching, innerR3FromMatching, outerR3FromMatching
    else:
        return nMatched, nExpected, [], [], [], []
        
        
        
def makeSeedFit(r1, r4, nMatched, nExpected, innerR2FromMatching, outerR2FromMatching, innerR3FromMatching, outerR3FromMatching):
    ### now check if the SVD match works
    allDDList = {'dd':[], 'linepts':[], 'i2Inner':[], 'i2Outer':[], 'i3Inner':[], 'i3Outer':[]}
    
    #### if we have 4 hits in 4 layers
    #name = "SVDFit_nMatched"+str(nMatched)+"_nExpected"+str(nExpected)+".pdf"
    #cnv = TCanvas("","",2000,2000)
    #cnv.SaveAs(name+"(")
    #print("inner r2 length: ", len(innerR2FromMatching))
    #print("outer r2 length: ", len(outerR2FromMatching))
    #print("inner r3 length: ", len(innerR3FromMatching))
    #print("outer r3 length: ", len(outerR3FromMatching))
    
    n2i = len(innerR2FromMatching)
    n2o = len(outerR2FromMatching)
    n3i = len(innerR3FromMatching)
    n3o = len(outerR3FromMatching)
    
    allHits = {}
    for i2Inner in range(len(innerR2FromMatching)):
        allHits.update({"Layer2Inner_"+str(i2Inner):innerR2FromMatching[i2Inner]})
    for i2Outer in range(len(outerR2FromMatching)):
        allHits.update({"Layer2Outer_"+str(i2Outer):outerR2FromMatching[i2Outer]})
    for i3Inner in range(len(innerR3FromMatching)):
        allHits.update({"Layer3Inner_"+str(i3Inner):innerR3FromMatching[i3Inner]})
    for i3Outer in range(len(outerR3FromMatching)):
        allHits.update({"Layer3Outer_"+str(i3Outer):outerR3FromMatching[i3Outer]})
    
    
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
                    
                        #r3Outer = []
                        #lfitpts, dd = seed3dfitSVDWithList(r1,r2Inner, r2Outer, r3Inner, r3Outer, r4)
                        #allDDList['dd'].append(dd.tolist())
                        #allDDList['i2Inner'].append(-1)
                        #allDDList['i2Outer'].append(-1)
                        #allDDList['i3Inner'].append(-1)
                        #allDDList['i3Outer'].append(-1)
                        #allDDList['linepts'].append(lfitpts)
                        
    
    ddValue1 = 10e10
    ddValue2 = 10e10
    iWinner = -1
    #print("ddList:", allDDList)
    for i in range(len(allDDList['dd'])):
        ddValue = allDDList['dd'][i]
        #print("ddValue1", ddValue)
        if ((0.0 < ddValue[1] < 0.1) and (0 < ddValue[2] < 0.05)):
            if ddValue[1] < ddValue1:
                ddValue1 = ddValue[1] 
                ddValue2 = ddValue[2]
                iWinner  = i
    
    if iWinner >= 0:
        return True, {"ddValue1":ddValue1, "ddValue2": ddValue2, "linepts":allDDList['linepts'][iWinner], "r2Inner":allDDList['i2Inner'][iWinner], "r2Outer":allDDList['i2Outer'][iWinner], "r3Inner":allDDList['i3Inner'][iWinner], "r3Outer":allDDList['i3Outer'][iWinner], 'nMatched':nMatched, 'nExpected':nExpected}
    else:
        return False, {}


### making the seeds from two tracks from innermost layer and outermost layer
def makeseed(r1, r4, allR2Inner, allR2Outer, allR3Inner, allR3Outer, side, r1GeV, r10GeV):
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
    ### the following cannot be 20mm/2 according to the signal tracks
    if(abs(yDipoleExit) > 10.8/2):
        return False, {}
    cutFlowDict['yDipoleExitGt15'] += 1
    
    ### This is according to the signal x:y at the dipole exit
    if(abs(xDipoleExit) < 20.0):
        return False, {}  # the track should point to |x|<~1.0 at the dipole exit
    cutFlowDict['xDipoleExitLt25'] += 1
    
    if(abs(xDipoleExit) > 330.0/2):
        return False, {}
    cutFlowDict['xDipoleExitGt165'] += 1

    ### select only positive x or negative x accroding to the particle
    if( (side=="Positron" and xDipoleExit < 0) or (side=="Electron" and xDipoleExit > 0)):
        return False, {}
    cutFlowDict['xDipoleExitLt0'] += 1
    
    #### find the energy of the tracks
    x0        = 0
    z0        = zofx(r1, r4, x0)
    #xExit     = abs(xofz(r1, r4, zDipoleActiveExit))
    #H         = abs((zDipoleActiveExit-z0))/1000.0 ### converting H from mm to m
    xExit     = abs(xofz(r1, r4, zDipoleActiveExit))
    H         = abs((zDipoleActiveExit-z0))/1000.0 ### converting H from mm to m
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
    ### The seed lorentz vector, need for energy plotting
    p         = TLorentzVector()
    p.SetPxPyPzE(px, py, pz, math.sqrt(px*px + py*py + pz*pz + meGeV2))
    
    ### checking the track energy
    if(p.E() < EseedMin or p.E() > EseedMax): 
        return False, {}
    cutFlowDict['trackEnergy'] += 1
    
    
    nMatched, nExpected, innerR2FromMatching, outerR2FromMatching, innerR3FromMatching, outerR3FromMatching = check_clusters(r1, r4, allR2Inner, allR2Outer, allR3Inner, allR3Outer, side)
    
    if nMatched < 3:
        return False, {}
    
    cutFlowDict['checkClusterTracksMiddleLayers'] += 1
    
    #expectedHits2Inner, expectedHits2Outer, expectedHits3Inner, expectedHits3Outer = getExpectedHits(r1, r4)
    #print("2time:", expectedHits2Inner, expectedHits2Outer,  expectedHits3Inner, expectedHits3Outer)
    
    passFit, winnerFit = makeSeedFit(r1, r4, nMatched, nExpected, innerR2FromMatching, outerR2FromMatching, innerR3FromMatching, outerR3FromMatching)
    
    if(not passFit):
        return False, {}
    
    cutFlowDict['checkClusterFit'] += 1
    winnerFit.update({"pSeed":p})
    
    xDipoleFromFit = xofz(winnerFit['linepts'][0], winnerFit['linepts'][1], zDipoleExit)*mm2m
    xLayer4FromFit = xofz(winnerFit['linepts'][0], winnerFit['linepts'][1], r4[2])*mm2m
    d = Distance(xDipoleFromFit,xLayer4FromFit,r1GeV,r10GeV)
    #print("The cluster x distance: ",d)
    if(d > 5*mm2m):
        return False, {}
    
    cutFlowDict['checkClusterXDistance'] += 1
    winnerFit.update({"distance":d})
    
    
    
    #if()

    return True, winnerFit





def main():
    # give the input text file containing all the track information
    parser = argparse.ArgumentParser(description='Code to find seed tracks')
    parser.add_argument('-l', action="store", dest="inFile", default="list_root_hics_165gev_w0_5000nm_provisional_10xi3_6cee466a_trackInfoClean.txt")
    parser.add_argument('-s', action="store_true", dest="needSignal", default=False)
    parser.add_argument('-e', action="store", dest="eCut", type=float, default=0.0)
    parser.add_argument('-p', action="store", dest="Particle", type=str, default="Positron")
    args = parser.parse_args()
    
    
    for index, row in df.iterrows():
        detid   = row["detid"]
        layerid = row["layerid"]
        if(layerid>7):
            continue
        xMin, xMax = GetSensorXBoundaries(detid, layerid)
        xBoundaries.update({"detid"+str(detid)+"_layerid"+str(layerid):[xMin, xMax]})
    
    #print(xBoundaries)
    
    ### open the file containing track file
    inTextFile      = args.inFile
    inputTrackInfo  = open(inTextFile)
    #energyCutSuffix = str(args.eCut)+"KeVCut"
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
    
    
    outFile             = TFile("seedingInformation_"+suffixName+"_"+energyCutSuffix+"_"+signalCutSuffix+"_"+particleSuffix+".root", "RECREATE")
    outFile.cd()
    
    hAllPossible        = TH1D("hAllPossible", "all possible track combination; bunch crossing; number of track combination", 9508, 0, 9508)
    hSeedPossible       = TH1D("hSeedPossible", "seed track combination; bunch crossing; number of seed track", 9508, 0, 9508)
    hSeedMultiplicity   = TH1D("hSeedMultiplicity", "seed multiplicity; number of seeds; BX", 2000, 0, 2000)
    hSignalMultiplicity = TH1D("hSignalMultiplicity", "number of signals; number of signals; BX", 2000, 0, 2000)
    hXLayer4XDipole     = TH2D("hXLayer4XDipole", "number of signal distribution; x_{Dipole} [mm]; x_{Layer4}",330,0,330,650,0,650)
    hSigEnergy          = TH1D("hSigEnergy", "signal energy; Energy [GeV]; Entries", 200, 0, 20)
    hSeedEnergy         = TH1D("hSeedEnergy", "seed energy; Energy [GeV]; Entries", 200, 0, 20)
    hSeedEnergyLoose    = TH1D("hSeedEnergyLoose", "seed energy (Loose); Energy [GeV]; Entries", 200, 0, 20)
    hSeedEnergyTight    = TH1D("hSeedEnergyTight", "seed energy (Tight); Energy [GeV]; Entries", 200, 0, 20)
    hSeedDistance       = TH1D("hSeedDistance", "seed distance wrt analytical line; d [m]; Entries", 200, 0, 0.005)
    hSeedDistanceLoose  = TH1D("hSeedDistanceLoose", "seed distance wrt analytical line (loose); d [m]; Entries", 200, 0, 0.005)
    hSeedDistanceTight  = TH1D("hSeedDistanceTight", "seed distance wrt analytical line (tight); d [m]; Entries", 200, 0, 0.005)
    hSVDValues1         = TH1D("hSVDValues1", "output of SVD[0]; fit quality [0]; Events", 2000, 100, 300)
    hSVDValues2         = TH1D("hSVDValues2", "output of SVD[1]; fit quality [1]; Events", 240, 0, 0.8)
    hSVDValues3         = TH1D("hSVDValues3", "output of SVD[2]; fit quality [2]; Events", 240, 0, 0.8)
    
    

    #### select the number of bunch crossing for different files, BX is useful only for hics setup
    if 'hics' in inTextFile:
        nBX          = 100
        checkBXMatch = True
    else:
        nBX          = 1
        checkBXMatch = False
        
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
        if(args.needSignal and not(pdgId==-11 and trackId==1)): continue
        ### if there is an energy cut
        
        ### ask for variable energy cut, energyAfterCut is in keV
        energyAfterCut = energyAbsorbed(staveId, energyVal, vtx_z)
        #if(float(eachWord[6]) > args.eCut*1e-6):
        if(energyAfterCut > 0):
        #if True:
            #### bxNumber, trackId, staveId, x, y, E, weight
            position.append([bxNumber, trackId, staveId, xPos, yPos, energyVal, weight])

    #pprint.pprint(position)
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
            if(side=="Positron"):
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
                    continue
            else:
                if (values[2] == 8):
                    allR1Inner.append([values[3], values[4], z1inner, values[5], values[6]])
                elif (values[2] == 9):
                    allR1Outer.append([values[3], values[4], z1outer, values[5], values[6]])
                elif (values[2] == 10):
                    allR2Inner.append([values[3], values[4], z2inner, values[5], values[6]])
                elif (values[2] == 11):
                    allR2Outer.append([values[3], values[4], z2outer, values[5], values[6]])
                elif (values[2] == 12):
                    allR3Inner.append([values[3], values[4], z3inner, values[5], values[6]])
                elif (values[2] == 13):
                    allR3Outer.append([values[3], values[4], z3outer, values[5], values[6]])
                elif (values[2] == 14):
                    allR4Inner.append([values[3], values[4], z4inner, values[5], values[6]])
                elif (values[2] == 15):
                    allR4Outer.append([values[3], values[4], z4outer, values[5], values[6]])
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
            hSigEnergy.Fill(r1[3])
            
        hSignalMultiplicity.Fill(len(allR1Unique))
        
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
        allSeedsTrackLines = []
        for r4 in allR4Unique[:]:
            for r1 in allR1Unique[:]:
                #True, {"linepts":allDDList['linepts'][iWinner], "r2Inner":allDDList['i2Inner'][iWinner], "r2Outer":allDDList['i2Outer'][iWinner], "r3Inner":allDDList['i3Inner'][iWinner], "r3Outer":allDDList['i3Outer'][iWinner], "pSeed":p.E()}
                seed, winnerDict = makeseed(r1, r4, allR2Inner, allR2Outer, allR3Inner, allR3Outer, side, r1GeV, r10GeV)
                allCounter += 1
                if(seed):
                    #print("I found a seed here: ", r1, " and ", r4)
                    hSeedEnergy.Fill(winnerDict["pSeed"].E())
                    
                    allSeedsTrackLines.append(GetExtendedTrackLine([r1, r4]))
                    xDipoleFromFit = xofz(winnerDict['linepts'][0], winnerDict['linepts'][1], zDipoleExit)
                    xLayer4FromFit = xofz(winnerDict['linepts'][0], winnerDict['linepts'][1], r4[2])
                    hXLayer4XDipole.Fill(xDipoleFromFit, xLayer4FromFit)
                    hSeedDistance.Fill(winnerDict['distance'])
                    
                    if(winnerDict['nMatched'] == 4):
                        hSeedEnergyTight.Fill(winnerDict["pSeed"].E())
                        hSeedDistanceTight.Fill(winnerDict['distance'])
                    else:
                        hSeedEnergyLoose.Fill(winnerDict["pSeed"].E())
                        hSeedDistanceLoose.Fill(winnerDict['distance'])
                    
                    counter += 1
                    #allR1Unique.remove(r1)
                    #print(ddList)
                    for ddValue in ddList:
                        #hSVDValues1.Fill(ddValue[0])
                        hSVDValues2.Fill(winnerDict["ddValue1"])
                        hSVDValues3.Fill(winnerDict["ddValue2"])
                    #print("svdOutput: ", svdOutput)
                    ##print(len(svdOutput), " and ", svdOutput[0])
                    #hSVDValues.Fill(svdOutput)

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
    cnv.SaveAs("trackingDiagram_"+suffixName+"_"+energyCutSuffix+"_"+signalCutSuffix+"_"+particleSuffix+".pdf")
    cnv.SaveAs("trackingDiagram_"+suffixName+"_"+energyCutSuffix+"_"+signalCutSuffix+"_"+particleSuffix+".root")
    pprint.pprint(cutFlowDict)


if __name__ == "__main__":
    start = time.time()
    main()
    print("-------- The processing time: ",time.time() - start, " s")
