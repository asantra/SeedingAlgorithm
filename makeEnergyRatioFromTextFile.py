#### this code mixes the background tracks and signal tracks per bunch crossing
#### run: python3 makeEnergyRatioFromTextFile.py <bxNumberWanted>
import os
import sys
import time
import pprint
import math
from ROOT import *
from collections import OrderedDict
import argparse

### global variables

m2mm          = 1000
mm2m          = 1./1000

### these are in mm
xDipoleWidth  = 330.0
yDipoleHeight = 108.0
zDipoleExit   = 2748.0

### length of the dipole in meters, needed in meters because pf the p = 0.3*B*R formula
LB            = 1.029
LB2           = LB*LB
LD            = 1.396
zDipoleActiveExit = zDipoleExit - (LD - LB)*m2mm/2  ### in mm, to be consistent with Sasha's values

### 1 Tesla magnetic field used for now
B        = 1.0

def xofz(r1,r2,z):
   dz = r2[2]-r1[2]
   dx = r2[0]-r1[0]
   if(dz==0):
   	print("ERROR in xofz: dz=0")
   	quit()
   a = dx/dz
   b = r1[0]-a*r1[2]
   x = a*z+b
   return x


### get radius of curvature in m
def getRadius(p):
    B = 1.0
    R = p/(0.3*B)
    return R


### get the xTangent in m given ZTangent in m and momentum p
def getXTangent(ZT,p):
    R  = getRadius(p)
    R2 = R*R
    XT = (R2/LB-ZT)*(LB/math.sqrt(R2-LB2))
    return XT

### get x in Sasha's coordinate in m given zTangent and momentum p
def getx(ZT,p):
    XT = getXTangent(ZT,p)
    R  = getRadius(p)
    x  = R-XT
    return x



def main():
    
    # give the signal sample you want to use, old with less signal tracks or new with more signal tracks
    parser = argparse.ArgumentParser(description='Code to find energy ratio at IP and at tracks')
    parser.add_argument('-in', action="store", dest="inFile", default="list_root_hics_165gev_w0_3000nm_jeti40_122020_9550dac4SignalPositrons_trackInfoClean.txt")
    args = parser.parse_args()
    
    
    ### hics signal
    sigFileName      = open("../Outputfile/"+args.inFile)
    
    suffixName       = args.inFile.split('.')[0]
    
    rootFile         = TFile("../Outputfile/energyRatio_"+suffixName+".root", "RECREATE")
    rootFile.cd()
    hEnergyIP        = TH1D("hEnergyIP", "Energy at IP; E_{IP} [GeV]; Signal Positrons/BX", 200, 0, 20)
    hEnergyTracker   = TH1D("hEnergyTracker", "Energy at Tracker; E_{Tracker} [GeV]; Signal Positrons/BX", 200, 0, 20)
    hEnergyRatio     = TH1D("hEnergyRatio", "Ratio of energy; E_{Tracker}/E_{IP}; Signal Positrons/BX", 120, 0.0, 1.2)
    hEnergyDiff      = TH1D("hEnergyDiff", "energy diff; (E_{IP} - E_{Tracker}}/E_{IP}; Signal Positrons/BX", 120, 0.0, 1.2)
    hEnergyDiffVsEip = TH2D("hEnergyDiffVsEip", "energy diff vs e_ip; E_{IP} [GeV]; (E_{IP} - E_{Tracker}}/E_{IP}",140,0,14,105,0.0,1.05)
    
    hDeltaX          = TH1D("hDeltaX","#Delta(x) of tracker and analytical calculation from fields; (x_{Field} - x_{Tracker}) [mm]; Particles/BX", 180, -30, 30)
    hDeltaY          = TH1D("hDeltaY","#Delta(y) of tracker and analytical calculation from fields; (y_{Field} - y_{Tracker}) [mm]; Particles/BX", 180, -10, 10)
    hDeltaXFraction  = TH1D("hDeltaXFraction","(x_{Field} - x_{TrackerFace})/x_{Field}; (x_{Field} - x_{Tracker})/x_{Field} ; Particles/BX", 1000, -0.1, 0.1)
    hDeltaYFraction  = TH1D("hDeltaYFraction","(y_{Field} - y_{TrackerFace})/y_{Field}; (y_{Field} - y_{Tracker})/y_{Field} ; Particles/BX", 400, -0.5, 0.5)
    
    hDeltaXFractionVsXField = TH2D("hDeltaXFractionVsXField","(x_{Field} - x_{TrackerFace})/x_{Field} vs x_{Field}; x_{Field} [mm]; (x_{Field} - x_{Tracker})/x_{Field}; Particles/BX", 600, 0, 600, 1000, -0.1, 0.1)
    
    hXSimulationVsXField = TH2D("hXSimulationVsXField","x_{TrackerFace} vs x_{Field}; x_{Field} [mm];  x_{Tracker} [mm]; Particles/BX", 600, 0, 600, 600, 0, 600)
    
    hStaveZVsXField = TH2D("hStaveZVsXField","Z position vs x_{Field}; x_{Field} [mm];  z_Stave [mm]; Particles/BX", 600, 0, 600, 50, 3850, 3900)
    
    hDeltaXFractionVsXSimulation = TH2D("hDeltaXFractionVsXSimulation","(x_{Field} - x_{TrackerFace})/x_{Field} vs x_{TrackerFace}; x_{TrackerFace} [mm]; (x_{Field} - x_{Tracker})/x_{Field}; Particles/BX", 600, 0, 600, 1000, -0.1, 0.1)
    hDeltaYFractionVsYField = TH2D("hDeltaYFractionVsYField","(y_{Field} - y_{TrackerFace})/y_{Field} vs y_{Field}; y_{Field} [mm]; (y_{Field} - y_{Tracker})/y_{Field}; Particles/BX", 15, -7.5, 7.5, 400, -0.5, 0.5)
    
    
    hDeltaXMagnet    = TH1D("hDeltaXMagnet","#Delta(x) at magnet exit; x_{Extrapolate} - x_{Magnet} [mm]; Particles/BX", 100, -200, 200)
    
    ### work on the signal, first collect all the signal lines in a dictionary
    ### the text file contains these values
    ### bxNumber << pdg << track_id << det_id << xx << yy << eneg << ev_weight << vtx_x << vtx_y << vtx_z << parentid
    
    position = []
    particleInEachBX = {}
    for lines2 in sigFileName.readlines():
        lines2   = lines2.rstrip()
        if "#" in lines2: continue
        eachWord = lines2.split()
        
        bxNumber  = int(eachWord[0])
        pdgId     = int(eachWord[1])
        trackId   = int(eachWord[2])
        staveId   = int(eachWord[3])
        x0        = float(eachWord[4])
        y0        = float(eachWord[5])
        energy    = float(eachWord[6])
        weight    = float(eachWord[7])
        parentid  = int(eachWord[11])
        pxx       = float(eachWord[12])
        pyy       = float(eachWord[13])
        pzz       = float(eachWord[14])
        
        
        
        if not ((trackId==1) and ((staveId == -1) or (staveId >=1000 and staveId <=1001) or (staveId >=1006 and staveId<=1007))): continue
        
        keyVal    = str(bxNumber)+'_'+str(parentid)
        p         = math.sqrt(pxx*pxx+pyy*pyy+pzz*pzz)
        
        ### this is in the coordinate where the center of curved track is at (0,0), in m unit
        
        
        particleInEachBX.setdefault(keyVal, []).append({"pdgId": pdgId, "trackId": trackId, "staveId":staveId, "energy":energy, "weight":weight, "p":p, "xFile":x0, "yFile":y0, "pyy":pyy, "pzz":pzz}) 
        ### this ensures that the new line for the same bxnumber does not erase the previous line, the new line will be appended as a list 
        
     
    print("Total bx:", bxNumber)
    print("Scale: ", 1.0/float(bxNumber))
    #pprint.pprint(particleInEachBX)
    
    counter = 0
    for keys in particleInEachBX:
        if(len(particleInEachBX[keys]) < 2): continue
        #print(keys, ":", particleInEachBX[keys][0])
        if(counter%100000==0):print("processed: ",counter)
        E_IP           = 0
        weight_IP      = 1
        E_Tracker      = 0
        weight_Tracker = 1
        xField         = -999999
        xTracker       = -999999
        p_IP           = -999999
        pYY            = -999999
        pZZ            = -999999
        xInitial       = -999999
        yInitial       = -999999
        stave          = -999999
        xSimulation    = -999999
        ySimulation    = -999999
        xMagnet        = -999999
        xExtrapolate   = -999999
        R1             = []
        R4             = []
        
        ### get values for each stave
        lengthList = len(particleInEachBX[keys])
        for i in  range(lengthList):
            if(particleInEachBX[keys][i]["staveId"]==-1):
                E_IP      = particleInEachBX[keys][i]["energy"]
                weight_IP = particleInEachBX[keys][i]["weight"]
                p_IP      = particleInEachBX[keys][i]["p"]
                xInitial  = particleInEachBX[keys][i]["xFile"]
                yInitial  = particleInEachBX[keys][i]["yFile"]
                pYY       = particleInEachBX[keys][i]["pyy"]
                pZZ       = particleInEachBX[keys][i]["pzz"]
                
                
            if(particleInEachBX[keys][i]["staveId"]==1000):  ## and particleInEachBX[keys][i]["xFile"] < (308.53+29.94176/2.0)
                E_Tracker      = particleInEachBX[keys][i]["energy"]
                weight_Tracker = particleInEachBX[keys][i]["weight"]
                stave          = particleInEachBX[keys][i]["staveId"]
                xSimulation    = particleInEachBX[keys][i]["xFile"]
                ySimulation    = particleInEachBX[keys][i]["yFile"]
            elif(particleInEachBX[keys][i]["staveId"]==1001): ##  and particleInEachBX[keys][i]["xFile"] > (308.53+29.94176/2.0)
                E_Tracker      = particleInEachBX[keys][i]["energy"]
                weight_Tracker = particleInEachBX[keys][i]["weight"]
                stave          = particleInEachBX[keys][i]["staveId"]
                xSimulation    = particleInEachBX[keys][i]["xFile"]
                ySimulation    = particleInEachBX[keys][i]["yFile"]
                
            if(particleInEachBX[keys][i]["staveId"]==1000):
                R1 = [particleInEachBX[keys][i]["xFile"], particleInEachBX[keys][i]["yFile"], 3864.5125]
            if(particleInEachBX[keys][i]["staveId"]==1001):
                R1 = [particleInEachBX[keys][i]["xFile"], particleInEachBX[keys][i]["yFile"], 3876.5125]
            if(particleInEachBX[keys][i]["staveId"]==1006):
                R4 = [particleInEachBX[keys][i]["xFile"], particleInEachBX[keys][i]["yFile"], 4164.5125]
            if(particleInEachBX[keys][i]["staveId"]==1007):
                R4 = [particleInEachBX[keys][i]["xFile"], particleInEachBX[keys][i]["yFile"], 4176.5125]
                
        if(E_IP==0):
            print("Something is wrong in this key: ", keys)
        else:
            ratio = (E_Tracker*weight_Tracker)/(E_IP*weight_IP)
            diff  = 1.0 - (ratio)
            
            if(diff < 0.0):
                print("E_Tracker: ", E_Tracker, " E_IP: ", E_IP, " diff: ", E_IP-E_Tracker)
                
            hEnergyIP.Fill(E_IP*weight_IP)
            hEnergyTracker.Fill(E_Tracker*weight_Tracker)
            hEnergyRatio.Fill(ratio)
            hEnergyDiff.Fill(diff)
            hEnergyDiffVsEip.Fill(E_IP, diff)
            
            if(stave==1000): ### and xSimulation < (308.53+29.94176/2.0)
                zStave       = (3864.5125-1533.5)*mm2m
                yValue       = (pYY/pZZ)*(3864.5125)
            if(stave==1001): ###  and xSimulation > (308.53+29.94176/2.0)
                zStave       = (3876.5125-1533.5)*mm2m
                yValue       = (pYY/pZZ)*(3876.5125)
                
            #### get the x value from the extrapolation of tangent tracks, in m
            if(p_IP > 0 and xSimulation != -999999):
                xValue = getx(zStave,p_IP)
                xValue = xValue*m2mm
                
                
                
                ### adding the initial position of signal particles
                xCalculation   = xInitial+xValue ###xInitial
                yCalculation   = yInitial+yValue
                deltaX         = (xCalculation - xSimulation)
                deltaXFraction = deltaX/xCalculation
                
                deltaY         = (yCalculation - ySimulation)
                deltaYFraction = deltaY/yCalculation
                
                hDeltaX.Fill(deltaX)
                hDeltaXFraction.Fill(deltaXFraction)
                hDeltaXFractionVsXField.Fill(xCalculation, deltaXFraction)
                hDeltaXFractionVsXSimulation.Fill(xSimulation, deltaXFraction)
                
                hXSimulationVsXField.Fill(xCalculation, xSimulation)
                
                
                #print(xSimulation, xCalculation)
                zValue = zStave*m2mm + 1533.5
                hStaveZVsXField.Fill(xCalculation, zValue)
                hDeltaY.Fill(deltaY)
                hDeltaYFraction.Fill(deltaYFraction)
                hDeltaYFractionVsYField.Fill(yCalculation, deltaYFraction)
                
            ### get values at magnet exit window
            if(p_IP > 0 and R1 and R4):
                xMagnet      = getx(2748.0-1533.5, p_IP)
                xExtrapolate = xofz(R1, R4, 2748.0)
                hDeltaXMagnet.Fill(xExtrapolate - xMagnet)
            
        counter+=1
        
    
    hEnergyIP.Scale(1.0/bxNumber)
    hEnergyTracker.Scale(1.0/bxNumber)
    hEnergyRatio.Scale(1.0/bxNumber)
    hEnergyDiff.Scale(1.0/bxNumber)
    hEnergyDiffVsEip.Scale(1.0/bxNumber)
    
    hDeltaX.Scale(1.0/bxNumber)
    hDeltaXFraction.Scale(1.0/bxNumber)
    hDeltaXFractionVsXField.Scale(1.0/bxNumber)
    hDeltaXFractionVsXSimulation.Scale(1.0/bxNumber)
    hDeltaY.Scale(1.0/bxNumber)
    hDeltaYFraction.Scale(1.0/bxNumber)
    hDeltaYFractionVsYField.Scale(1.0/bxNumber)
    hXSimulationVsXField.Scale(1.0/bxNumber)
    hStaveZVsXField.Scale(1.0/bxNumber)
    
    
    hDeltaXMagnet.Scale(1.0/bxNumber)
    
    
    
    rootFile.Write()
    rootFile.Close()
                
            
    
if __name__=="__main__":
    start = time.time()
    main()
    print("--- The time taken: ", time.time() - start, " s")
