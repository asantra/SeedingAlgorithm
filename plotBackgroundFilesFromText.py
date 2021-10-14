#### this code mixes the background tracks and signal tracks per bunch crossing
#### run: python3 plotBackgroundFilesFromText.py -l  <event file in text mode>
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
    
def energyBins(nbins):
  #//// variable binned in X axis histograms
  logmin           = -6;
  logmax           = 1;
  logbinwidth      = (logmax-logmin)/float(nbins);
  xpoints          = []
  
  for i in range(0, nbins+1):
      #print((logmin + i*logbinwidth), pow( 10,(logmin + i*logbinwidth) ))
      xpoints.append(pow( 10,(logmin + i*logbinwidth) ))
                     
  xpoints.append(2*pow( 10,1))
  return xpoints

def main():
    
    parser = argparse.ArgumentParser(description='Code to get 2D plots')
    parser.add_argument('-l', action="store", dest="inFile", type=str, default="ePlusLaserBackground_list_root_7671ee4c_First3BX_trackInfo.txt")
    args = parser.parse_args()
    
    inDir = '/Users/arkasantra/arka/Sasha_Work/OutputFile'
    
    bkgFileName   = open(inDir+'/'+args.inFile)
    withoutText   = args.inFile.split('.txt')[0]
    rootFile      = withoutText+".root"    
    
    if 'ePlusLaserBkgNewSamples' in args.inFile:
        nbx = 298.0
    elif 'gPlusLaserBkgNewSamples' in args.inFile:
        nbx = 80.0
    elif 'fast_9f6b6590_al_window' in args.inFile:
        nbx = 23.4
    elif 'fast_9f6b6590_34_35' in args.inFile:
        nbx = 16.56
    elif 'ePlusLaserKaptonWindow_hics_background_fast_9f6b6590_0_16' in args.inFile:
        nbx = 120.0
    elif 'hics_background_fast_9f6b6590_17_33' in args.inFile:
        nbx = 141.38
    elif 'ePlusLaserKaptonWindow' in args.inFile:
        nbx = 278.0
    elif 'ePlusLaserBackground' in args.inFile:
        nbx = 2.13
    else:
        nbx = 1.0
    
    print('BX selected: ',nbx)
    
    outFile       = TFile(inDir+'/'+rootFile, "RECREATE")
    outFile.cd()
    nbins  = 450
    xbins  = energyBins(nbins)
    xarray = array.array('d',xbins)
    
    #print(xarray, " ", len(xarray))
    
    allHistoDict  = {}
    for i in range(0, 16):
        #### gamma
        allHistoDict.update({"tracking_planes_background_track_x_gamma_"+str(i):TH1D("tracking_planes_background_track_x_gamma_"+str(i),"tracking_planes_background_track_x_gamma_"+str(i),1300, -650.0, 650.0)})
        allHistoDict.update({"tracking_planes_background_track_x_gamma_sumE_"+str(i):TH1D("tracking_planes_background_track_x_gamma_sumE_"+str(i),"tracking_planes_background_track_x_gamma_sumE_"+str(i),1300, -650.0, 650.0)})
        allHistoDict.update({"tracking_planes_background_track_x_gamma_1GeVCut_"+str(i):TH1D("tracking_planes_background_track_x_gamma_1GeVCut_"+str(i),"tracking_planes_background_track_x_gamma_1GeVCut_"+str(i),1300, -650.0, 650.0)})
        allHistoDict.update({"tracking_planes_background_track_e_gamma_log_"+str(i):TH1D("tracking_planes_background_track_e_gamma_log_"+str(i),"tracking_planes_background_track_e_gamma_log_"+str(i),nbins+1, xarray)})
        allHistoDict.update({"tracking_planes_background_vtx_z_track_e_gamma_log_"+str(i):TH2D("tracking_planes_background_vtx_z_track_e_gamma_log_"+str(i),"tracking_planes_background_vtx_z_track_e_gamma_log_"+str(i),4000, -5000, 15000, nbins+1, xarray)})
        allHistoDict.update({"tracking_planes_background_vtx_x_track_e_gamma_log_"+str(i):TH2D("tracking_planes_background_vtx_x_track_e_gamma_log_"+str(i),"tracking_planes_background_vtx_x_track_e_gamma_log_"+str(i),6000, -3000, 3000, nbins+1, xarray)})
        allHistoDict.update({"tracking_planes_background_vtxz_vtxx_gamma_"+str(i):TH2D("tracking_planes_background_vtxz_vtxx_gamma_"+str(i), "tracking_planes_background_vtxz_vtxx_gamma_"+str(i), 4000, -5000, 15000, 6000, -3000, 3000)})
        
        
                             
        #### positrons
        allHistoDict.update({"tracking_planes_background_track_x_positrons_"+str(i):TH1D("tracking_planes_background_track_x_positrons_"+str(i),"tracking_planes_background_track_x_positrons_"+str(i),1300, -650.0, 650.0)})
        allHistoDict.update({"tracking_planes_background_track_x_positrons_sumE_"+str(i):TH1D("tracking_planes_background_track_x_positrons_sumE_"+str(i),"tracking_planes_background_track_x_positrons_sumE_"+str(i),1300, -650.0, 650.0)})
        allHistoDict.update({"tracking_planes_background_track_x_positrons_1GeVCut_"+str(i):TH1D("tracking_planes_background_track_x_positrons_1GeVCut_"+str(i),"tracking_planes_background_track_x_positrons_1GeVCut_"+str(i),1300, -650.0, 650.0)})
        allHistoDict.update({"tracking_planes_background_track_e_positrons_log_"+str(i):TH1D("tracking_planes_background_track_e_positrons_log_"+str(i),"tracking_planes_background_track_e_positrons_log_"+str(i),nbins+1, xarray)})
        allHistoDict.update({"tracking_planes_background_vtx_z_track_e_positrons_log_"+str(i):TH2D("tracking_planes_background_vtx_z_track_e_positrons_log_"+str(i),"tracking_planes_background_vtx_z_track_e_positrons_log_"+str(i),4000, -5000, 15000, nbins+1, xarray)})
        allHistoDict.update({"tracking_planes_background_vtx_x_track_e_positrons_log_"+str(i):TH2D("tracking_planes_background_vtx_x_track_e_positrons_log_"+str(i),"tracking_planes_background_vtx_x_track_e_positrons_log_"+str(i),6000, -3000, 3000, nbins+1, xarray)})
        allHistoDict.update({"tracking_planes_background_vtxz_vtxx_positrons_"+str(i):TH2D("tracking_planes_background_vtxz_vtxx_positrons_"+str(i), "tracking_planes_background_vtxz_vtxx_positrons_"+str(i), 4000, -5000, 15000, 6000, -3000, 3000)})
        
        
        #### electrons
        allHistoDict.update({"tracking_planes_background_track_x_electrons_"+str(i):TH1D("tracking_planes_background_track_x_electrons_"+str(i),"tracking_planes_background_track_x_electrons_"+str(i),1300, -650.0, 650.0)})
        allHistoDict.update({"tracking_planes_background_track_x_electrons_sumE_"+str(i):TH1D("tracking_planes_background_track_x_electrons_sumE_"+str(i),"tracking_planes_background_track_x_electrons_sumE_"+str(i),1300, -650.0, 650.0)})
        allHistoDict.update({"tracking_planes_background_track_x_electrons_1GeVCut_"+str(i):TH1D("tracking_planes_background_track_x_electrons_1GeVCut_"+str(i),"tracking_planes_background_track_x_electrons_1GeVCut_"+str(i),1300, -650.0, 650.0)})
        allHistoDict.update({"tracking_planes_background_track_e_electrons_log_"+str(i):TH1D("tracking_planes_background_track_e_electrons_log_"+str(i),"tracking_planes_background_track_e_electrons_log_"+str(i),nbins+1, xarray)})
        allHistoDict.update({"tracking_planes_background_vtx_z_track_e_electrons_log_"+str(i):TH2D("tracking_planes_background_vtx_z_track_e_electrons_log_"+str(i),"tracking_planes_background_vtx_z_track_e_electrons_log_"+str(i),4000, -5000, 15000, nbins+1, xarray)})
        allHistoDict.update({"tracking_planes_background_vtx_x_track_e_electrons_log_"+str(i):TH2D("tracking_planes_background_vtx_x_track_e_electrons_log_"+str(i),"tracking_planes_background_vtx_x_track_e_electrons_log_"+str(i),6000, -3000, 3000, nbins+1, xarray)})
        allHistoDict.update({"tracking_planes_background_vtxz_vtxx_electrons_"+str(i):TH2D("tracking_planes_background_vtxz_vtxx_electrons_"+str(i), "tracking_planes_background_vtxz_vtxx_electrons_"+str(i), 4000, -5000, 15000, 6000, -3000, 3000)})
        
        
        
        #### protons
        allHistoDict.update({"tracking_planes_background_track_x_protons_"+str(i):TH1D("tracking_planes_background_track_x_protons_"+str(i),"tracking_planes_background_track_x_protons_"+str(i),1300, -650.0, 650.0)})
        allHistoDict.update({"tracking_planes_background_track_x_protons_sumE_"+str(i):TH1D("tracking_planes_background_track_x_protons_sumE_"+str(i),"tracking_planes_background_track_x_protons_sumE_"+str(i),1300, -650.0, 650.0)})
        allHistoDict.update({"tracking_planes_background_track_x_protons_1GeVCut_"+str(i):TH1D("tracking_planes_background_track_x_protons_1GeVCut_"+str(i),"tracking_planes_background_track_x_protons_1GeVCut_"+str(i),1300, -650.0, 650.0)})
        allHistoDict.update({"tracking_planes_background_track_e_protons_log_"+str(i):TH1D("tracking_planes_background_track_e_protons_log_"+str(i),"tracking_planes_background_track_e_protons_log_"+str(i),nbins+1, xarray)})
        allHistoDict.update({"tracking_planes_background_vtx_z_track_e_protons_log_"+str(i):TH2D("tracking_planes_background_vtx_z_track_e_protons_log_"+str(i),"tracking_planes_background_vtx_z_track_e_protons_log_"+str(i),4000, -5000, 15000, nbins+1, xarray)})
        allHistoDict.update({"tracking_planes_background_vtx_x_track_e_protons_log_"+str(i):TH2D("tracking_planes_background_vtx_x_track_e_protons_log_"+str(i),"tracking_planes_background_vtx_x_track_e_protons_log_"+str(i),6000, -3000, 3000, nbins+1, xarray)})
        allHistoDict.update({"tracking_planes_background_vtxz_vtxx_protons_"+str(i):TH2D("tracking_planes_background_vtxz_vtxx_protons_"+str(i), "tracking_planes_background_vtxz_vtxx_protons_"+str(i), 4000, -5000, 15000, 6000, -3000, 3000)})
        
        
        
        #### neutrons
        allHistoDict.update({"tracking_planes_background_track_x_neutrons_"+str(i):TH1D("tracking_planes_background_track_x_neutrons_"+str(i),"tracking_planes_background_track_x_neutrons_"+str(i),1300, -650.0, 650.0)})
        allHistoDict.update({"tracking_planes_background_track_x_neutrons_sumE_"+str(i):TH1D("tracking_planes_background_track_x_neutrons_sumE_"+str(i),"tracking_planes_background_track_x_neutrons_sumE_"+str(i),1300, -650.0, 650.0)})
        allHistoDict.update({"tracking_planes_background_track_x_neutrons_1GeVCut_"+str(i):TH1D("tracking_planes_background_track_x_neutrons_1GeVCut_"+str(i),"tracking_planes_background_track_x_neutrons_1GeVCut_"+str(i),1300, -650.0, 650.0)})
        allHistoDict.update({"tracking_planes_background_track_e_neutrons_log_"+str(i):TH1D("tracking_planes_background_track_e_neutrons_log_"+str(i),"tracking_planes_background_track_e_neutrons_log_"+str(i),nbins+1, xarray)})
        allHistoDict.update({"tracking_planes_background_vtx_z_track_e_neutrons_log_"+str(i):TH2D("tracking_planes_background_vtx_z_track_e_neutrons_log_"+str(i),"tracking_planes_background_vtx_z_track_e_neutrons_log_"+str(i),4000, -5000, 15000, nbins+1, xarray)})
        allHistoDict.update({"tracking_planes_background_vtx_x_track_e_neutrons_log_"+str(i):TH2D("tracking_planes_background_vtx_x_track_e_neutrons_log_"+str(i),"tracking_planes_background_vtx_x_track_e_neutrons_log_"+str(i),6000, -3000, 3000, nbins+1, xarray)})
        allHistoDict.update({"tracking_planes_background_vtxz_vtxx_neutrons_"+str(i):TH2D("tracking_planes_background_vtxz_vtxx_neutrons_"+str(i), "tracking_planes_background_vtxz_vtxx_neutrons_"+str(i), 4000, -5000, 15000, 6000, -3000, 3000)})
        
        
        
        
        #### muons
        allHistoDict.update({"tracking_planes_background_track_x_muons_"+str(i):TH1D("tracking_planes_background_track_x_muons_"+str(i),"tracking_planes_background_track_x_muons_"+str(i),1300, -650.0, 650.0)})
        allHistoDict.update({"tracking_planes_background_track_x_muons_sumE_"+str(i):TH1D("tracking_planes_background_track_x_muons_sumE_"+str(i),"tracking_planes_background_track_x_muons_sumE_"+str(i),1300, -650.0, 650.0)})
        allHistoDict.update({"tracking_planes_background_track_x_muons_1GeVCut_"+str(i):TH1D("tracking_planes_background_track_x_muons_1GeVCut_"+str(i),"tracking_planes_background_track_x_muons_1GeVCut_"+str(i),1300, -650.0, 650.0)})
        allHistoDict.update({"tracking_planes_background_track_e_muons_log_"+str(i):TH1D("tracking_planes_background_track_e_muons_log_"+str(i),"tracking_planes_background_track_e_muons_log_"+str(i),nbins+1, xarray)})
        allHistoDict.update({"tracking_planes_background_vtx_z_track_e_muons_log_"+str(i):TH2D("tracking_planes_background_vtx_z_track_e_muons_log_"+str(i),"tracking_planes_background_vtx_z_track_e_muons_log_"+str(i),4000, -5000, 15000, nbins+1, xarray)})
        allHistoDict.update({"tracking_planes_background_vtx_x_track_e_muons_log_"+str(i):TH2D("tracking_planes_background_vtx_x_track_e_muons_log_"+str(i),"tracking_planes_background_vtx_x_track_e_muons_log_"+str(i),6000, -3000, 3000, nbins+1, xarray)})
        allHistoDict.update({"tracking_planes_background_vtxz_vtxx_muons_"+str(i):TH2D("tracking_planes_background_vtxz_vtxx_muons_"+str(i), "tracking_planes_background_vtxz_vtxx_muons_"+str(i), 4000, -5000, 15000, 6000, -3000, 3000)})
        
        
        
        
        #### pions
        allHistoDict.update({"tracking_planes_background_track_x_pions_"+str(i):TH1D("tracking_planes_background_track_x_pions_"+str(i),"tracking_planes_background_track_x_pions_"+str(i),1300, -650.0, 650.0)})
        allHistoDict.update({"tracking_planes_background_track_x_pions_sumE_"+str(i):TH1D("tracking_planes_background_track_x_pions_sumE_"+str(i),"tracking_planes_background_track_x_pions_sumE_"+str(i),1300, -650.0, 650.0)})
        allHistoDict.update({"tracking_planes_background_track_x_pions_1GeVCut_"+str(i):TH1D("tracking_planes_background_track_x_pions_1GeVCut_"+str(i),"tracking_planes_background_track_x_pions_1GeVCut_"+str(i),1300, -650.0, 650.0)})
        allHistoDict.update({"tracking_planes_background_track_e_pions_log_"+str(i):TH1D("tracking_planes_background_track_e_pions_log_"+str(i),"tracking_planes_background_track_e_pions_log_"+str(i),nbins+1, xarray)})
        allHistoDict.update({"tracking_planes_background_vtx_z_track_e_pions_log_"+str(i):TH2D("tracking_planes_background_vtx_z_track_e_pions_log_"+str(i),"tracking_planes_background_vtx_z_track_e_pions_log_"+str(i),4000, -5000, 15000, nbins+1, xarray)})
        allHistoDict.update({"tracking_planes_background_vtx_x_track_e_pions_log_"+str(i):TH2D("tracking_planes_background_vtx_x_track_e_pions_log_"+str(i),"tracking_planes_background_vtx_x_track_e_pions_log_"+str(i),6000, -3000, 3000, nbins+1, xarray)})
        allHistoDict.update({"tracking_planes_background_vtxz_vtxx_pions_"+str(i):TH2D("tracking_planes_background_vtxz_vtxx_pions_"+str(i), "tracking_planes_background_vtxz_vtxx_pions_"+str(i), 4000, -5000, 15000, 6000, -3000, 3000)})
        
        
        
        #### pions 0
        allHistoDict.update({"tracking_planes_background_track_x_pi0_"+str(i):TH1D("tracking_planes_background_track_x_pi0_"+str(i),"tracking_planes_background_track_x_pi0_"+str(i),1300, -650.0, 650.0)})
        allHistoDict.update({"tracking_planes_background_track_x_pi0_sumE_"+str(i):TH1D("tracking_planes_background_track_x_pi0_sumE_"+str(i),"tracking_planes_background_track_x_pi0_sumE_"+str(i),1300, -650.0, 650.0)})
        allHistoDict.update({"tracking_planes_background_track_x_pi0_1GeVCut_"+str(i):TH1D("tracking_planes_background_track_x_pi0_1GeVCut_"+str(i),"tracking_planes_background_track_x_pi0_1GeVCut_"+str(i),1300, -650.0, 650.0)})
        allHistoDict.update({"tracking_planes_background_track_e_pi0_log_"+str(i):TH1D("tracking_planes_background_track_e_pi0_log_"+str(i),"tracking_planes_background_track_e_pi0_log_"+str(i),nbins+1, xarray)})
        allHistoDict.update({"tracking_planes_background_vtx_z_track_e_pi0_log_"+str(i):TH2D("tracking_planes_background_vtx_z_track_e_pi0_log_"+str(i),"tracking_planes_background_vtx_z_track_e_pi0_log_"+str(i),4000, -5000, 15000, nbins+1, xarray)})
        allHistoDict.update({"tracking_planes_background_vtx_x_track_e_pi0_log_"+str(i):TH2D("tracking_planes_background_vtx_x_track_e_pi0_log_"+str(i),"tracking_planes_background_vtx_x_track_e_pi0_log_"+str(i),6000, -3000, 3000, nbins+1, xarray)})
        allHistoDict.update({"tracking_planes_background_vtxz_vtxx_pi0_"+str(i):TH2D("tracking_planes_background_vtxz_vtxx_pi0_"+str(i), "tracking_planes_background_vtxz_vtxx_pi0_"+str(i), 4000, -5000, 15000, 6000, -3000, 3000)})
    
    
    
    lineCounter   = 0
    ### write the bkg as it is
    for lines in bkgFileName.readlines():
        if '#' in lines:
            continue
        lineCounter += 1
        if(lineCounter%100000==0): print("processed: ", lineCounter)
        #if(lineCounter > 200000):break
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
        ###electrons
        if(pdgId == 11 and trackId!=1):
            allHistoDict["tracking_planes_background_track_x_electrons_"+str(staveId)].Fill(xPos, weight)
            if(energyVal>1.0):
                allHistoDict["tracking_planes_background_track_x_electrons_1GeVCut_"+str(staveId)].Fill(xPos, weight)
            allHistoDict["tracking_planes_background_track_x_electrons_sumE_"+str(staveId)].Fill(xPos, energyVal*weight)
            
            allHistoDict["tracking_planes_background_track_e_electrons_log_"+str(staveId)].Fill(energyVal, weight)
            allHistoDict["tracking_planes_background_vtx_z_track_e_electrons_log_"+str(staveId)].Fill(vtx_z, energyVal, weight)
            allHistoDict["tracking_planes_background_vtx_x_track_e_electrons_log_"+str(staveId)].Fill(vtx_x, energyVal, weight)
            allHistoDict["tracking_planes_background_vtxz_vtxx_electrons_"+str(staveId)].Fill(vtx_z, vtx_x, weight)
        
        ###protons and anti-protons
        if((pdgId == 2212 or pdgId == -2212) and trackId!=1):
            allHistoDict["tracking_planes_background_track_x_protons_"+str(staveId)].Fill(xPos, weight)
            if(energyVal>1.0):
                allHistoDict["tracking_planes_background_track_x_protons_1GeVCut_"+str(staveId)].Fill(xPos, weight)
            allHistoDict["tracking_planes_background_track_x_protons_sumE_"+str(staveId)].Fill(xPos, energyVal*weight)
            
            allHistoDict["tracking_planes_background_track_e_protons_log_"+str(staveId)].Fill(energyVal, weight)
            allHistoDict["tracking_planes_background_vtx_z_track_e_protons_log_"+str(staveId)].Fill(vtx_z, energyVal, weight)
            allHistoDict["tracking_planes_background_vtx_x_track_e_protons_log_"+str(staveId)].Fill(vtx_x, energyVal, weight)
            allHistoDict["tracking_planes_background_vtxz_vtxx_protons_"+str(staveId)].Fill(vtx_z, vtx_x, weight)
            
        ### neutrons
        if(pdgId == 2112 and trackId!=1):
            allHistoDict["tracking_planes_background_track_x_neutrons_"+str(staveId)].Fill(xPos, weight)
            if(energyVal>1.0):
                allHistoDict["tracking_planes_background_track_x_neutrons_1GeVCut_"+str(staveId)].Fill(xPos, weight)
            allHistoDict["tracking_planes_background_track_x_neutrons_sumE_"+str(staveId)].Fill(xPos, energyVal*weight)
            
            allHistoDict["tracking_planes_background_track_e_neutrons_log_"+str(staveId)].Fill(energyVal, weight)
            allHistoDict["tracking_planes_background_vtx_z_track_e_neutrons_log_"+str(staveId)].Fill(vtx_z, energyVal, weight)
            allHistoDict["tracking_planes_background_vtx_x_track_e_neutrons_log_"+str(staveId)].Fill(vtx_x, energyVal, weight)
            allHistoDict["tracking_planes_background_vtxz_vtxx_neutrons_"+str(staveId)].Fill(vtx_z, vtx_x, weight)
            
        ### muons and anti-muons
        if((pdgId == 13 or pdgId == -13) and trackId!=1):
            allHistoDict["tracking_planes_background_track_x_muons_"+str(staveId)].Fill(xPos, weight)
            if(energyVal>1.0):
                allHistoDict["tracking_planes_background_track_x_muons_1GeVCut_"+str(staveId)].Fill(xPos, weight)
            allHistoDict["tracking_planes_background_track_x_muons_sumE_"+str(staveId)].Fill(xPos, energyVal*weight)
            
            allHistoDict["tracking_planes_background_track_e_muons_log_"+str(staveId)].Fill(energyVal, weight)
            allHistoDict["tracking_planes_background_vtx_z_track_e_muons_log_"+str(staveId)].Fill(vtx_z, energyVal, weight)
            allHistoDict["tracking_planes_background_vtx_x_track_e_muons_log_"+str(staveId)].Fill(vtx_x, energyVal, weight)
            allHistoDict["tracking_planes_background_vtxz_vtxx_muons_"+str(staveId)].Fill(vtx_z, vtx_x, weight)
            
        ### pi+ or pi-
        if((pdgId == 211 or pdgId == -211) and trackId!=1):
            allHistoDict["tracking_planes_background_track_x_pions_"+str(staveId)].Fill(xPos, weight)
            if(energyVal>1.0):
                allHistoDict["tracking_planes_background_track_x_pions_1GeVCut_"+str(staveId)].Fill(xPos, weight)
            allHistoDict["tracking_planes_background_track_x_pions_sumE_"+str(staveId)].Fill(xPos, energyVal*weight)
            
            allHistoDict["tracking_planes_background_track_e_pions_log_"+str(staveId)].Fill(energyVal, weight)
            allHistoDict["tracking_planes_background_vtx_z_track_e_pions_log_"+str(staveId)].Fill(vtx_z, energyVal, weight)
            allHistoDict["tracking_planes_background_vtx_x_track_e_pions_log_"+str(staveId)].Fill(vtx_x, energyVal, weight)
            allHistoDict["tracking_planes_background_vtxz_vtxx_pions_"+str(staveId)].Fill(vtx_z, vtx_x, weight)
            
        ### pi0
        if((pdgId == 111 or pdgId == -111) and trackId!=1):
            allHistoDict["tracking_planes_background_track_x_pi0_"+str(staveId)].Fill(xPos, weight)
            if(energyVal>1.0):
                allHistoDict["tracking_planes_background_track_x_pi0_1GeVCut_"+str(staveId)].Fill(xPos, weight)
            allHistoDict["tracking_planes_background_track_x_pi0_sumE_"+str(staveId)].Fill(xPos, energyVal*weight)
            
            allHistoDict["tracking_planes_background_track_e_pi0_log_"+str(staveId)].Fill(energyVal, weight)
            allHistoDict["tracking_planes_background_vtx_z_track_e_pi0_log_"+str(staveId)].Fill(vtx_z, energyVal, weight)
            allHistoDict["tracking_planes_background_vtx_x_track_e_pi0_log_"+str(staveId)].Fill(vtx_x, energyVal, weight)
            allHistoDict["tracking_planes_background_vtxz_vtxx_pi0_"+str(staveId)].Fill(vtx_z, vtx_x, weight)
            
        ### positrons
        if(pdgId == -11 and trackId!=1):
            allHistoDict["tracking_planes_background_track_x_positrons_"+str(staveId)].Fill(xPos, weight)
            if(energyVal>1.0):
                allHistoDict["tracking_planes_background_track_x_positrons_1GeVCut_"+str(staveId)].Fill(xPos, weight)
            allHistoDict["tracking_planes_background_track_x_positrons_sumE_"+str(staveId)].Fill(xPos, energyVal*weight)
            allHistoDict["tracking_planes_background_track_e_positrons_log_"+str(staveId)].Fill(energyVal, weight)
            allHistoDict["tracking_planes_background_vtx_z_track_e_positrons_log_"+str(staveId)].Fill(vtx_z, energyVal, weight)
            allHistoDict["tracking_planes_background_vtx_x_track_e_positrons_log_"+str(staveId)].Fill(vtx_x, energyVal, weight)
            allHistoDict["tracking_planes_background_vtxz_vtxx_positrons_"+str(staveId)].Fill(vtx_z, vtx_x, weight)
            
            
        ### photon
        if(pdgId == 22 and trackId!=1):
            allHistoDict["tracking_planes_background_track_x_gamma_"+str(staveId)].Fill(xPos, weight)
            if(energyVal>1.0):
                allHistoDict["tracking_planes_background_track_x_gamma_1GeVCut_"+str(staveId)].Fill(xPos, weight)
            allHistoDict["tracking_planes_background_track_x_gamma_sumE_"+str(staveId)].Fill(xPos, energyVal*weight)
            allHistoDict["tracking_planes_background_track_e_gamma_log_"+str(staveId)].Fill(energyVal, weight)
            allHistoDict["tracking_planes_background_vtx_z_track_e_gamma_log_"+str(staveId)].Fill(vtx_z, energyVal, weight)
            allHistoDict["tracking_planes_background_vtx_x_track_e_gamma_log_"+str(staveId)].Fill(vtx_x, energyVal, weight)
            allHistoDict["tracking_planes_background_vtxz_vtxx_gamma_"+str(staveId)].Fill(vtx_z, vtx_x, weight)
        
        
    for keys in allHistoDict:
        allHistoDict[keys].Scale(1./nbx)
        allHistoDict[keys].Write()
    outFile.Close()
    
if __name__=="__main__":
    start = time.time()
    main()
    print("--- The time taken: ", time.time() - start, " s")
