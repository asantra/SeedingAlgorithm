#### this code mixes the background tracks and signal tracks per bunch crossing
#### run: python3 getSigBkgTextFile.py <bxNumberWanted>
import os
import sys
import time
import pprint
import math, array
from ROOT import *
from collections import OrderedDict
import argparse

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
    
def energyBins():
  #//// variable binned in X axis histograms
  logmin           = -6;
  logmax           = 1;
  nbins            = 450;
  logbinwidth      = (logmax-logmin)/float(nbins);
  xpoints          = []
  
  for i in xrange(0, nbins+1):
      #print((logmin + i*logbinwidth), pow( 10,(logmin + i*logbinwidth) ))
      xpoints.append(pow( 10,(logmin + i*logbinwidth) ))
                     
  xpoints.append(2*pow( 10,1))
  return xpoints

def main():
    
    #outFile  = TFile("AllBackground2DDistribution_WithVtxCut.root", "RECREATE")
    outFile  = TFile("AllGplusLaserBackground2DDistribution_WithoutVtxCut.root", "RECREATE")
    outFile.cd()
    xbins  = energyBins()
    xarray = array.array('d',xbins)
    #print(xarray, " ", len(xarray))
    
    allHistoDict  = {}
    for i in xrange(0, 16):
        allHistoDict.update({"h2DVtxxVtxz_Electron_"+str(i):TH2D("tracking_planes_background_vtxz_vtxx_electrons_"+str(i),"Electron vertex for stave"+str(i)+"; vtxz [mm]; vtxx [mm]", 6000, 1000, 7000, 1000, -500, 500)})
        allHistoDict.update({"h2DVtxxVtxz_Positron_"+str(i):TH2D("tracking_planes_background_vtxz_vtxx_positrons_"+str(i),"Positron vertex for stave"+str(i)+"; vtxz [mm]; vtxx [mm]", 6000, 1000, 7000, 1000, -500, 500)})
        allHistoDict.update({"h2DVtxyVtxz_Electron_"+str(i):TH2D("tracking_planes_background_vtxz_vtxy_electrons_"+str(i),"Electron vertex for stave"+str(i)+"; vtxz [mm]; vtxy [mm]", 6000, 1000, 7000, 100, -50, 50)})
        allHistoDict.update({"h2DVtxyVtxz_Positron_"+str(i):TH2D("tracking_planes_background_vtxz_vtxy_positrons_"+str(i),"Positron vertex for stave"+str(i)+"; vtxz [mm]; vtxy [mm]", 6000, 1000, 7000, 100, -50, 50)})
        allHistoDict.update({"h2DVtxyVtxx_Electron_"+str(i):TH2D("tracking_planes_background_vtxx_vtxy_electrons_"+str(i),"Electron vertex for stave"+str(i)+"; vtxx [mm]; vtxy [mm]", 1000, -500, 500, 100, -50, 50)})
        allHistoDict.update({"h2DVtxyVtxx_Positron_"+str(i):TH2D("tracking_planes_background_vtxx_vtxy_positrons_"+str(i),"Positron vertex for stave"+str(i)+"; vtxx [mm]; vtxy [mm]", 1000, -500, 500, 100, -50, 50)})
        allHistoDict.update({"h2DEVtxx_Positron_"+str(i):TH2D("tracking_planes_background_vtx_x_track_e_positrons_log_"+str(i),"Positron E vs vtxx for stave"+str(i)+"; vtxx [mm]; E [GeV]", 1000, -500, 500, len(xbins)-1, xarray)})
        allHistoDict.update({"h2DEVtxz_Positron_"+str(i):TH2D("tracking_planes_background_vtx_z_track_e_positrons_log_"+str(i),"Positron E vs vtxz for stave"+str(i)+"; vtxz [mm]; E [GeV]", 6000, 1000, 7000, len(xbins)-1, xarray)})
        allHistoDict.update({"h2DEVtxx_Electron_"+str(i):TH2D("tracking_planes_background_vtx_x_track_e_electrons_log_"+str(i),"Electron E vs vtxx for stave"+str(i)+"; vtxx [mm]; E [GeV]", 1000, -500, 500, len(xbins)-1, xarray)})
        allHistoDict.update({"h2DEVtxz_Electron_"+str(i):TH2D("tracking_planes_background_vtx_z_track_e_electrons_log_"+str(i),"Electron E vs vtxz for stave"+str(i)+"; vtxz [mm]; E [GeV]", 6000, 1000, 7000, len(xbins)-1, xarray)})
        
        #allHistoDict.update({"h2DEVtxx_Positron_"+str(i):TH2D("h2DEVtxx_Positron_"+str(i),"Positron E vs vtxx for stave"+str(i)+"; vtxx [mm]; E [GeV]", 1000, 500, 500, 1000, 0, 10)})
        #allHistoDict.update({"h2DEVtxz_Positron_"+str(i):TH2D("h2DEVtxz_Positron_"+str(i),"Positron E vs vtxz for stave"+str(i)+"; vtxz [mm]; E [GeV]", 6000, 1000, 7000, 1000, 0, 10)})
        #allHistoDict.update({"h2DEVtxx_Electron_"+str(i):TH2D("h2DEVtxx_Electron_"+str(i),"Electron E vs vtxx for stave"+str(i)+"; vtxx [mm]; E [GeV]", 1000, 500, 500, 1000, 0, 10)})
        #allHistoDict.update({"h2DEVtxz_Electron_"+str(i):TH2D("h2DEVtxz_Electron_"+str(i),"Electron E vs vtxz for stave"+str(i)+"; vtxz [mm]; E [GeV]", 6000, 1000, 7000, 1000, 0, 10)})
    
    
    #bkgFileName   = open("EBeamOnlyNewSamples_AllBX_trackInfoClean.txt")
    bkgFileName   = open("gPlusLaserBkgNewSamples_AllBX_trackInfoClean.txt")
    lineCounter   = 0
    ### write the bkg as it is
    for lines in bkgFileName.readlines():
        if '#' in lines:
            continue
        lineCounter += 1
        if(lineCounter%10000==0): print("processed: ", lineCounter)
        lines        = lines.rstrip()
        eachWord     = lines.split()
        bxNumber     = int(eachWord[0])
        trackId      = int(eachWord[2])
        pdgId        = int(eachWord[1])
        staveId      = int(eachWord[3])-1000
        xPos         = float(eachWord[4])
        yPos         = float(eachWord[5])
        energyVal    = float(eachWord[6])
        weight       = float(eachWord[7])
        vtx_x        = float(eachWord[8])
        vtx_y        = float(eachWord[9])
        vtx_z        = float(eachWord[10])
        
        failVtxCut = checkVtxCut(vtx_x, vtx_z)
        #if(vtx_x < -25.5 and (vtx_z > 3600 and vtx_z < 4600)): continue
        #if(failVtxCut): continue
    
        if(pdgId == 11 and trackId!=1):
            if(energyVal > 10): print("Electron: ", energyVal, " ", vtx_x, " ", vtx_z)
            allHistoDict["h2DVtxxVtxz_Electron_"+str(staveId)].Fill(vtx_z, vtx_x, weight)
            allHistoDict["h2DVtxyVtxz_Electron_"+str(staveId)].Fill(vtx_z, vtx_y, weight)
            allHistoDict["h2DVtxyVtxx_Electron_"+str(staveId)].Fill(vtx_x, vtx_y, weight)
            allHistoDict["h2DEVtxx_Electron_"+str(staveId)].Fill(vtx_x, energyVal, weight)
            allHistoDict["h2DEVtxz_Electron_"+str(staveId)].Fill(vtx_z, energyVal, weight)
        if(pdgId == -11 and trackId!=1):
            if(energyVal > 10): print("Positron: ", energyVal, " ", vtx_x, " ", vtx_z)
            allHistoDict["h2DVtxxVtxz_Positron_"+str(staveId)].Fill(vtx_z, vtx_x, weight)
            allHistoDict["h2DVtxyVtxz_Positron_"+str(staveId)].Fill(vtx_z, vtx_y, weight)
            allHistoDict["h2DVtxyVtxx_Positron_"+str(staveId)].Fill(vtx_x, vtx_y, weight)
            allHistoDict["h2DEVtxx_Positron_"+str(staveId)].Fill(vtx_x, energyVal, weight)
            allHistoDict["h2DEVtxz_Positron_"+str(staveId)].Fill(vtx_z, energyVal, weight)
        
    for keys in allHistoDict:
        #allHistoDict[keys].Scale(1./160.)
        allHistoDict[keys].Scale(1./30.)
        allHistoDict[keys].Write()
    outFile.Close()
    
if __name__=="__main__":
    start = time.time()
    main()
    print("--- The time taken: ", time.time() - start, " s")
