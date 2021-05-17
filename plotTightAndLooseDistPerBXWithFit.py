import os, sys, glob, time
from ROOT import *
import argparse
from copy import copy, deepcopy
sys.path.insert(0, '/Users/arkasantra/arka/include')
from Functions import *
import pprint



def DrawHistsRatio(FirstTH1, LegendName, PlotColor, xrange1down, xrange1up, yrange1down, yrange1up, xaxisTitle, CanvasName, h2, yline1low, yline1up, drawline=False, logy=False, LatexName='', LatexName2='', TeVTag=False, doSumw2=False, doAtlas=False, doLumi=False, noRatio=False, do80=False, do59=False):
   Tex = MakeLatex(0.44,0.79,LatexName)
   Tex2 = MakeLatex(0.48,0.70,LatexName2)
   c = TCanvas("c","c",900, 900)
   #c.SetGrid(0)
   pad1 = TPad('pad1','pad1',0,0.25,1,1)
   pad1.SetBottomMargin(0.0)
   pad1.SetFillStyle(0)
   #pad1.SetGrid(0)
   pad1.Draw()
   pad1.cd()
   gStyle.SetOptStat(0)
   if(logy):
     pad1.SetLogy()
     
   line = MakeLine(xrange1down,yline1low,xrange1up,yline1up)
   legend1 = LegendMaker()
   tex1 = TLatex(); tex2 = TLatex(); tex3 = TLatex()
   L = [tex1, tex2, tex3]
   TexMaker(L, doAtlas, doLumi, noRatio, do80, do59)
   FirstTH1[0].GetYaxis().SetRangeUser(yrange1down,yrange1up)
   FirstTH1[0].GetXaxis().SetRangeUser(xrange1down,xrange1up)
   WaterMark = TexWaterMark('Preliminary')
   for i in range(0, len(FirstTH1)):
     FirstTH1[i].SetTitle("")
     FirstTH1[i].GetYaxis().SetTitle(FirstTH1[i].GetYaxis().GetTitle())#(FirstTH1[i].GetYaxis().GetTitle())
     FirstTH1[i].GetYaxis().SetTitleOffset(0.95)
     FirstTH1[i].GetYaxis().SetTitleSize(0.05)
     FirstTH1[i].GetYaxis().SetLabelSize(0.045)
     FirstTH1[i].GetXaxis().SetTitle("")
     FirstTH1[i].GetXaxis().SetLabelSize(0.15)
     
     #w = FirstTH1[i].Integral()
     ### comment in if you want only shape distribution
     #FirstTH1[i].Scale(1.0/w)
     FirstTH1[i] = SetHistColorEtc(FirstTH1[i], PlotColor[i])
     FirstTH1[i].SetMarkerSize(1)
     FirstTH1[i].SetMarkerStyle(20)
     if(i==0):legend1.AddEntry(FirstTH1[i],LegendName[i], "f")
     else: legend1.AddEntry(FirstTH1[i],LegendName[i], "lp")
   
   FirstTH1[0].GetXaxis().SetLabelSize(0.0);
   FirstTH1[0].GetYaxis().SetTitle("Tracks/BX");
   FirstTH1[0].SetFillColor(PlotColor[0])
   FirstTH1[0].Draw("hist")
   #WaterMark.Draw("sames")
   gPad.SetTickx()
   gPad.SetTicky()
   gPad.Modified(); gPad.Update();
   for i in xrange(1, len(FirstTH1)):
     FirstTH1[i].Draw("lp sames")  
     gPad.Modified(); gPad.Update();
   #FirstTH1[0].Draw("ep sames")
   gPad.RedrawAxis()
   L[1].Draw("sames")
   L[2].Draw("sames")
   if(TeVTag):
       TexTeV = TLatex(0.892,0.914,"#sqrt{s}=13 TeV")#36.075 fb^{-1} (13 TeV)
       TexTeV.SetNDC()
       TexTeV.SetTextAlign(31)
       TexTeV.SetTextFont(42)
       TexTeV.SetTextSize(0.037)
       TexTeV.SetLineWidth(2) 
       TexTeV.Draw()
   else:
       L[0].Draw()
   LUXELabel(0.25,0.85,"CDR")
   legend1.Draw("sames")
   Tex.Draw("sames")
   Tex2.Draw("sames")
   
   c.cd()
   pad2 = TPad("pad2", "pad2",0.0,0.0,1.0,0.245)
   pad2.SetTopMargin(0.0)
   pad2.SetBottomMargin(0.3)
   pad2.SetFillStyle(0)
   pad2.SetGrid(0)
   pad2.Draw()
   pad2.cd()
   gStyle.SetOptStat(0)
   gPad.SetTickx()
   gPad.SetTicky()
   
   if(doSumw2):
       for i in xrange(0, len(FirstTH1)):
          FirstTH1[i].Sumw2() 
   
   h2.Divide(FirstTH1[1],FirstTH1[0])
   h2.GetXaxis().SetRangeUser(xrange1down,xrange1up)
   h2.GetXaxis().SetTitle(xaxisTitle)
   h2.GetYaxis().SetTitle(h2.GetYaxis().GetTitle())
   h2.GetYaxis().CenterTitle()
   h2.SetTitle('')
   h2 = OneRatioAxisSize(h2)
   h2.GetYaxis().SetRangeUser(0., 2.4)
   h2.GetYaxis().SetNdivisions(6)
   h2.SetLineColor(kBlack)
   h2.SetMarkerColor(kBlack)
   h2.SetMarkerStyle(20)
   h2.SetMarkerSize(1)
   h2.SetLineStyle(1)
   h2.Draw("ep")
   if(drawline):
     line.Draw()
     
   SaveFile(c, CanvasName)
   c.Clear()
     
     

def main():
    gROOT.LoadMacro("LuxeStyle.C")
    gROOT.LoadMacro("LuxeLabels.C")
    gROOT.SetBatch()
    SetLuxeStyle()
    gROOT.SetBatch()
    directory = "RecoAndTruthEnergy_Jan142021"  ## _LooseSVDanddCuts_December29_2020
    if not os.path.exists(directory):
        os.makedirs(directory)
    # give the signal sample you want to use, old with less signal tracks or new with more signal tracks
    parser = argparse.ArgumentParser(description='Code to find seed tracks')
    parser.add_argument('-n', action="store", dest="nTracks", type=str, default="1")
    parser.add_argument('-p', action="store", dest="postFix", type=str, default="WithFit3or4HitsTracksAndDistanceCut")
    parser.add_argument('-f', action="store", dest="postFixNoFit", type=str, default="WithoutFit3or4HitsTracksAndDistanceCut")
    
    args = parser.parse_args()
    inDirectory = "."
    
    
    #### signal and background file, with fit
    signalAndBackgroundBX1 = TFile(inDirectory+"/seedingInformation_BkgEBeam_SignalHics3000nmOldForSignalMultiplicityLessThan20Or5000nmForSignalMultiplicityMoreThan20_BX1_SignalTracks"+args.nTracks+"_trackInfoClean_VariableEnergyCut_SignalAndBackground_PositronSide_"+args.postFix+".root", "READ")
    
    signalAndBackgroundBX2 = TFile(inDirectory+"/seedingInformation_BkgEBeam_SignalHics3000nmOldForSignalMultiplicityLessThan20Or5000nmForSignalMultiplicityMoreThan20_BX2_SignalTracks"+args.nTracks+"_trackInfoClean_VariableEnergyCut_SignalAndBackground_PositronSide_"+args.postFix+".root", "READ")
    
    signalAndBackgroundBX3 = TFile(inDirectory+"/seedingInformation_BkgEBeam_SignalHics3000nmOldForSignalMultiplicityLessThan20Or5000nmForSignalMultiplicityMoreThan20_BX3_SignalTracks"+args.nTracks+"_trackInfoClean_VariableEnergyCut_SignalAndBackground_PositronSide_"+args.postFix+".root", "READ")
    
    signalAndBackgroundBX4 = TFile(inDirectory+"/seedingInformation_BkgEBeam_SignalHics3000nmOldForSignalMultiplicityLessThan20Or5000nmForSignalMultiplicityMoreThan20_BX4_SignalTracks"+args.nTracks+"_trackInfoClean_VariableEnergyCut_SignalAndBackground_PositronSide_"+args.postFix+".root", "READ")
    
    
    
    
    #### signal and background file, without fit
    signalAndBackgroundWithoutFitBX1 = TFile(inDirectory+"/seedingInformation_BkgEBeam_SignalHics3000nmOldForSignalMultiplicityLessThan20Or5000nmForSignalMultiplicityMoreThan20_BX1_SignalTracks"+args.nTracks+"_trackInfoClean_VariableEnergyCut_signalAndBackground_PositronSide_"+args.postFixNoFit+".root", "READ")
    
    signalAndBackgroundWithoutFitBX2 = TFile(inDirectory+"/seedingInformation_BkgEBeam_SignalHics3000nmOldForSignalMultiplicityLessThan20Or5000nmForSignalMultiplicityMoreThan20_BX2_SignalTracks"+args.nTracks+"_trackInfoClean_VariableEnergyCut_signalAndBackground_PositronSide_"+args.postFixNoFit+".root", "READ")
    
    signalAndBackgroundWithoutFitBX3 = TFile(inDirectory+"/seedingInformation_BkgEBeam_SignalHics3000nmOldForSignalMultiplicityLessThan20Or5000nmForSignalMultiplicityMoreThan20_BX3_SignalTracks"+args.nTracks+"_trackInfoClean_VariableEnergyCut_signalAndBackground_PositronSide_"+args.postFixNoFit+".root", "READ")
    
    signalAndBackgroundWithoutFitBX4 = TFile(inDirectory+"/seedingInformation_BkgEBeam_SignalHics3000nmOldForSignalMultiplicityLessThan20Or5000nmForSignalMultiplicityMoreThan20_BX4_SignalTracks"+args.nTracks+"_trackInfoClean_VariableEnergyCut_signalAndBackground_PositronSide_"+args.postFixNoFit+".root", "READ")
    





    #### only signal file
    signalFromBkgFileBX1 = TFile(inDirectory+"/seedingInformation_BkgEBeam_SignalHics3000nmOldForSignalMultiplicityLessThan20Or5000nmForSignalMultiplicityMoreThan20_BX1_SignalTracks"+args.nTracks+"_trackInfoClean_VariableEnergyCut_OnlySignal_PositronSide_"+args.postFix+".root", "READ")
    signalFromBkgFileBX2 = TFile(inDirectory+"/seedingInformation_BkgEBeam_SignalHics3000nmOldForSignalMultiplicityLessThan20Or5000nmForSignalMultiplicityMoreThan20_BX2_SignalTracks"+args.nTracks+"_trackInfoClean_VariableEnergyCut_OnlySignal_PositronSide_"+args.postFix+".root", "READ")
    signalFromBkgFileBX3 = TFile(inDirectory+"/seedingInformation_BkgEBeam_SignalHics3000nmOldForSignalMultiplicityLessThan20Or5000nmForSignalMultiplicityMoreThan20_BX3_SignalTracks"+args.nTracks+"_trackInfoClean_VariableEnergyCut_OnlySignal_PositronSide_"+args.postFix+".root", "READ")
    signalFromBkgFileBX4 = TFile(inDirectory+"/seedingInformation_BkgEBeam_SignalHics3000nmOldForSignalMultiplicityLessThan20Or5000nmForSignalMultiplicityMoreThan20_BX4_SignalTracks"+args.nTracks+"_trackInfoClean_VariableEnergyCut_OnlySignal_PositronSide_"+args.postFix+".root", "READ")
    
    
    plotSuffix = "BkgEBeam_SignalHics3000nmOldForSignalMultiplicityLessThan20Or5000nmForSignalMultiplicityMoreThan20_SVDFit"
    latexName2  = 'signal tracks ~ '+args.nTracks




    ### with fit, energy from signal+background, inclusive
    hSeedEnergy_signalAndBackgroundBX1  = signalAndBackgroundBX1.Get("hSeedEnergy")
    hSeedEnergy_signalAndBackgroundBX2  = signalAndBackgroundBX2.Get("hSeedEnergy")
    hSeedEnergy_signalAndBackgroundBX3  = signalAndBackgroundBX3.Get("hSeedEnergy")
    hSeedEnergy_signalAndBackgroundBX4  = signalAndBackgroundBX4.Get("hSeedEnergy")

    hSeedEnergy_signalAndBackgroundBX1.Add(hSeedEnergy_signalAndBackgroundBX2)
    hSeedEnergy_signalAndBackgroundBX1.Add(hSeedEnergy_signalAndBackgroundBX3)
    hSeedEnergy_signalAndBackgroundBX1.Add(hSeedEnergy_signalAndBackgroundBX4)
    
    hSeedEnergy_signalAndBackgroundBX1.Rebin(4)
    
    
    ### with fit, energy from signal+background, tight
    hSeedEnergyTight_signalAndBackgroundBX1  = signalAndBackgroundBX1.Get("hSeedEnergyTight")
    hSeedEnergyTight_signalAndBackgroundBX2  = signalAndBackgroundBX2.Get("hSeedEnergyTight")
    hSeedEnergyTight_signalAndBackgroundBX3  = signalAndBackgroundBX3.Get("hSeedEnergyTight")
    hSeedEnergyTight_signalAndBackgroundBX4  = signalAndBackgroundBX4.Get("hSeedEnergyTight")

    hSeedEnergyTight_signalAndBackgroundBX1.Add(hSeedEnergyTight_signalAndBackgroundBX2)
    hSeedEnergyTight_signalAndBackgroundBX1.Add(hSeedEnergyTight_signalAndBackgroundBX3)
    hSeedEnergyTight_signalAndBackgroundBX1.Add(hSeedEnergyTight_signalAndBackgroundBX4)
    
    hSeedEnergyTight_signalAndBackgroundBX1.Rebin(4)
    
    
    
    ### with fit, energy from signal+background, Loose
    hSeedEnergyLoose_signalAndBackgroundBX1  = signalAndBackgroundBX1.Get("hSeedEnergyLoose")
    hSeedEnergyLoose_signalAndBackgroundBX2  = signalAndBackgroundBX2.Get("hSeedEnergyLoose")
    hSeedEnergyLoose_signalAndBackgroundBX3  = signalAndBackgroundBX3.Get("hSeedEnergyLoose")
    hSeedEnergyLoose_signalAndBackgroundBX4  = signalAndBackgroundBX4.Get("hSeedEnergyLoose")

    hSeedEnergyLoose_signalAndBackgroundBX1.Add(hSeedEnergyLoose_signalAndBackgroundBX2)
    hSeedEnergyLoose_signalAndBackgroundBX1.Add(hSeedEnergyLoose_signalAndBackgroundBX3)
    hSeedEnergyLoose_signalAndBackgroundBX1.Add(hSeedEnergyLoose_signalAndBackgroundBX4)
    
    hSeedEnergyLoose_signalAndBackgroundBX1.Rebin(4)
    
    
    
    
    
    
    
    
    
    ### with fit, energy from signal+background, inclusive
    hSeedEnergy_signalAndBackgroundWithoutFitBX1  = signalAndBackgroundWithoutFitBX1.Get("hSeedEnergy")
    hSeedEnergy_signalAndBackgroundWithoutFitBX2  = signalAndBackgroundWithoutFitBX2.Get("hSeedEnergy")
    hSeedEnergy_signalAndBackgroundWithoutFitBX3  = signalAndBackgroundWithoutFitBX3.Get("hSeedEnergy")
    hSeedEnergy_signalAndBackgroundWithoutFitBX4  = signalAndBackgroundWithoutFitBX4.Get("hSeedEnergy")

    hSeedEnergy_signalAndBackgroundWithoutFitBX1.Add(hSeedEnergy_signalAndBackgroundWithoutFitBX2)
    hSeedEnergy_signalAndBackgroundWithoutFitBX1.Add(hSeedEnergy_signalAndBackgroundWithoutFitBX3)
    hSeedEnergy_signalAndBackgroundWithoutFitBX1.Add(hSeedEnergy_signalAndBackgroundWithoutFitBX4)
    
    hSeedEnergy_signalAndBackgroundWithoutFitBX1.Rebin(4)
    
    
    ### with fit, energy from signal+background, tight
    hSeedEnergyTight_signalAndBackgroundWithoutFitBX1  = signalAndBackgroundWithoutFitBX1.Get("hSeedEnergyTight")
    hSeedEnergyTight_signalAndBackgroundWithoutFitBX2  = signalAndBackgroundWithoutFitBX2.Get("hSeedEnergyTight")
    hSeedEnergyTight_signalAndBackgroundWithoutFitBX3  = signalAndBackgroundWithoutFitBX3.Get("hSeedEnergyTight")
    hSeedEnergyTight_signalAndBackgroundWithoutFitBX4  = signalAndBackgroundWithoutFitBX4.Get("hSeedEnergyTight")

    hSeedEnergyTight_signalAndBackgroundWithoutFitBX1.Add(hSeedEnergyTight_signalAndBackgroundWithoutFitBX2)
    hSeedEnergyTight_signalAndBackgroundWithoutFitBX1.Add(hSeedEnergyTight_signalAndBackgroundWithoutFitBX3)
    hSeedEnergyTight_signalAndBackgroundWithoutFitBX1.Add(hSeedEnergyTight_signalAndBackgroundWithoutFitBX4)
    
    hSeedEnergyTight_signalAndBackgroundWithoutFitBX1.Rebin(4)
    
    
    
    ### with fit, energy from signal+background, Loose
    hSeedEnergyLoose_signalAndBackgroundWithoutFitBX1  = signalAndBackgroundWithoutFitBX1.Get("hSeedEnergyLoose")
    hSeedEnergyLoose_signalAndBackgroundWithoutFitBX2  = signalAndBackgroundWithoutFitBX2.Get("hSeedEnergyLoose")
    hSeedEnergyLoose_signalAndBackgroundWithoutFitBX3  = signalAndBackgroundWithoutFitBX3.Get("hSeedEnergyLoose")
    hSeedEnergyLoose_signalAndBackgroundWithoutFitBX4  = signalAndBackgroundWithoutFitBX4.Get("hSeedEnergyLoose")

    hSeedEnergyLoose_signalAndBackgroundWithoutFitBX1.Add(hSeedEnergyLoose_signalAndBackgroundWithoutFitBX2)
    hSeedEnergyLoose_signalAndBackgroundWithoutFitBX1.Add(hSeedEnergyLoose_signalAndBackgroundWithoutFitBX3)
    hSeedEnergyLoose_signalAndBackgroundWithoutFitBX1.Add(hSeedEnergyLoose_signalAndBackgroundWithoutFitBX4)
    
    hSeedEnergyLoose_signalAndBackgroundWithoutFitBX1.Rebin(4)
    
    


    #### actual signal
    hSignalEnergyBX1                     = signalFromBkgFileBX1.Get("hSigEnergy")
    hSignalEnergyBX2                     = signalFromBkgFileBX2.Get("hSigEnergy")
    hSignalEnergyBX3                     = signalFromBkgFileBX3.Get("hSigEnergy")
    hSignalEnergyBX4                     = signalFromBkgFileBX4.Get("hSigEnergy")

    hSignalEnergyBX1.Add(hSignalEnergyBX2)
    hSignalEnergyBX1.Add(hSignalEnergyBX3)
    hSignalEnergyBX1.Add(hSignalEnergyBX4)
    
    hSignalEnergyBX1.Rebin(4)
    
    
    #### seed energy, matched with signal
    hSeedEnergyMatchedBX1                     = signalFromBkgFileBX1.Get("hSeedEnergy")
    hSeedEnergyMatchedBX2                     = signalFromBkgFileBX2.Get("hSeedEnergy")
    hSeedEnergyMatchedBX3                     = signalFromBkgFileBX3.Get("hSeedEnergy")
    hSeedEnergyMatchedBX4                     = signalFromBkgFileBX4.Get("hSeedEnergy")

    hSeedEnergyMatchedBX1.Add(hSeedEnergyMatchedBX2)
    hSeedEnergyMatchedBX1.Add(hSeedEnergyMatchedBX3)
    hSeedEnergyMatchedBX1.Add(hSeedEnergyMatchedBX4)
    
    hSeedEnergyMatchedBX1.Rebin(4)

    
    #if(args.nTracks=="1"):
        #yrange1up   = 3
    #elif(args.nTracks=="5"):
        #yrange1up   = 3
    #elif(args.nTracks=="10"):
        #yrange1up   = 6
    #elif(args.nTracks=="20"):
        #yrange1up   = 8
    #elif(args.nTracks=="50"):
        #yrange1up   = 15
    #elif(args.nTracks=="100"):
        #yrange1up   = 20
    #elif(args.nTracks=="200"):
        #yrange1up   = 50
    #else:
        #yrange1up   = 500

    yrange1up   = 1e4
    xAxisName   = "E [GeV]"
    yAxisName   = "Tracks/BX"
    xrange1down = 0
    xrange1up   = 17
    yrange1down = 0.001
    
    yline1low   = 1
    yline1up    = 1
    
    drawline    = False
    logy        = True
    
    latexName   = "Data equivalent to 4 BXs"
    latexName3  = "sig from e+laser"
    latexName4  = "bkg from beam-only"
    latexName5  = ""
    
    leftLegend  = False
    doAtlas     = False
    doLumi      = False
    noRatio     = False
    do80        = False
    do59        = False
    drawPattern = ""
    logz        = False
    logx        = False
    drawline    = True
    energyPlots = True
    
    FirstTH1 = [hSignalEnergyBX1, hSeedEnergy_signalAndBackgroundBX1, hSeedEnergyTight_signalAndBackgroundBX1, hSeedEnergyLoose_signalAndBackgroundBX1, hSeedEnergy_signalAndBackgroundWithoutFitBX1, hSeedEnergyTight_signalAndBackgroundWithoutFitBX1, hSeedEnergyLoose_signalAndBackgroundWithoutFitBX1]
    
    #### put the x axis label and y axis label
    for i in xrange(0,len(FirstTH1)):
        FirstTH1[i].Scale(1./4.0)
        FirstTH1[i].GetYaxis().SetTitle(yAxisName)
        FirstTH1[i].GetXaxis().SetTitle(xAxisName)

    

    h2 = hSeedEnergy_signalAndBackgroundBX1.Clone("h2")
    h2.Reset()
    h2.GetYaxis().SetTitle("#frac{seed energy (reco, sig+bkg)}{energy (sig, from Geant4)}")

    h3 = hSeedEnergy_signalAndBackgroundBX1.Clone("h3")
    h3.Reset()
    h3.GetYaxis().SetTitle("#frac{seed energy (reco, sig matched)}{energy (sig, from Geant4)}")
    
    h4 = hSeedEnergy_signalAndBackgroundBX1.Clone("h4")
    h4.Reset()
    h4.GetYaxis().SetTitle("#frac{seed energy (reco, sig+bkg)}{energy (sig, from Geant4)}")

    h5 = hSeedEnergy_signalAndBackgroundBX1.Clone("h5")
    h5.Reset()
    h5.GetYaxis().SetTitle("#frac{seed energy (reco, sig matched)}{energy (sig, from Geant4)}")



    SecondTH1  = [hSignalEnergyBX1, hSeedEnergy_signalAndBackgroundWithoutFitBX1, hSeedEnergy_signalAndBackgroundBX1]

    PlotColor  = [kGray, kBlack, kBlack]
    LegendName = ['True signal', 'Inclusive', 'Inclusive+Fit']
    
    DrawHistsRatioFour(SecondTH1, LegendName, PlotColor, xrange1down, xrange1up, yrange1down, yrange1up, directory+"/energyDistributionSeedAndSignalInclusive_WithRatio_nTracks"+args.nTracks+"_"+plotSuffix, h2, h3, yline1low, yline1up, drawline, logy, False, False, latexName, latexName2, latexName3, latexName4, latexName5, h4, h5, energyPlots)
    
    
    
    SecondTH1 = [hSignalEnergyBX1, hSeedEnergyTight_signalAndBackgroundWithoutFitBX1, hSeedEnergyTight_signalAndBackgroundBX1]

    PlotColor = [kGray, 2, 2]
    LegendName = ['True signal','Tight', 'Tight+Fit']
    
    DrawHistsRatioFour(SecondTH1, LegendName, PlotColor, xrange1down, xrange1up, yrange1down, yrange1up, directory+"/energyDistributionSeedAndSignalTight_WithRatio_nTracks"+args.nTracks+"_"+plotSuffix, h2, h3, yline1low, yline1up, drawline, logy, False, False, latexName, latexName2, latexName3, latexName4, latexName5, h4, h5, energyPlots)
    
    
    SecondTH1 = [hSignalEnergyBX1, hSeedEnergyLoose_signalAndBackgroundWithoutFitBX1, hSeedEnergyLoose_signalAndBackgroundBX1]

    PlotColor = [kGray, 4, 4]
    LegendName = ['True signal', 'Loose', 'Loose+Fit']
    
    DrawHistsRatioFour(SecondTH1, LegendName, PlotColor, xrange1down, xrange1up, yrange1down, yrange1up, directory+"/energyDistributionSeedAndSignalLoose_WithRatio_nTracks"+args.nTracks+"_"+plotSuffix, h2, h3, yline1low, yline1up, drawline, logy, False, False, latexName, latexName2, latexName3, latexName4, latexName5, h4, h5, energyPlots)
    
    
    
    SecondTH1 = [hSignalEnergyBX1, hSeedEnergyTight_signalAndBackgroundBX1]

    PlotColor = [kGray, kBlack]
    LegendName = ['Signal positrons','Reconstructed tracks']
    
    #DrawHistsRatioFour(SecondTH1, LegendName, PlotColor, xrange1down, xrange1up, yrange1down, yrange1up, directory+"/energyDistributionSeedAndSignalTightReco_WithRatio_nTracks"+args.nTracks+"_"+plotSuffix, h2, h3, yline1low, yline1up, drawline, logy, False, False, latexName, latexName2, latexName3, latexName4, latexName5, h4, h5, energyPlots)
    
    
    
    DrawHistsRatio(SecondTH1, LegendName, PlotColor, xrange1down, xrange1up, 2e-1, 2e2, "E [GeV]", directory+"/energyDistributionSeedAndSignalTightReco_WithRatio_nTracks"+args.nTracks+"_"+plotSuffix, h2, yline1low, yline1up, drawline, logy, "e+laser (5 #mum)")
    
    
    
    
    
    
    ##### make distribution of the distance parameter
    
    ### with fit, distance from signal+background, inclusive
    hSeedDistance_signalAndBackgroundBX1  = signalAndBackgroundBX1.Get("hSeedDistance")
    hSeedDistance_signalAndBackgroundBX2  = signalAndBackgroundBX2.Get("hSeedDistance")
    hSeedDistance_signalAndBackgroundBX3  = signalAndBackgroundBX3.Get("hSeedDistance")
    hSeedDistance_signalAndBackgroundBX4  = signalAndBackgroundBX4.Get("hSeedDistance")

    hSeedDistance_signalAndBackgroundBX1.Add(hSeedDistance_signalAndBackgroundBX2)
    hSeedDistance_signalAndBackgroundBX1.Add(hSeedDistance_signalAndBackgroundBX3)
    hSeedDistance_signalAndBackgroundBX1.Add(hSeedDistance_signalAndBackgroundBX4)
    
    hSeedDistance_signalAndBackgroundBX1.Rebin(4)
    
    
    ### with fit, distance from signal+background, tight
    hSeedDistanceTight_signalAndBackgroundBX1  = signalAndBackgroundBX1.Get("hSeedDistanceTight")
    hSeedDistanceTight_signalAndBackgroundBX2  = signalAndBackgroundBX2.Get("hSeedDistanceTight")
    hSeedDistanceTight_signalAndBackgroundBX3  = signalAndBackgroundBX3.Get("hSeedDistanceTight")
    hSeedDistanceTight_signalAndBackgroundBX4  = signalAndBackgroundBX4.Get("hSeedDistanceTight")

    hSeedDistanceTight_signalAndBackgroundBX1.Add(hSeedDistanceTight_signalAndBackgroundBX2)
    hSeedDistanceTight_signalAndBackgroundBX1.Add(hSeedDistanceTight_signalAndBackgroundBX3)
    hSeedDistanceTight_signalAndBackgroundBX1.Add(hSeedDistanceTight_signalAndBackgroundBX4)
    
    hSeedDistanceTight_signalAndBackgroundBX1.Rebin(4)
    
    
    ### with fit, distance from signal+background, loose
    hSeedDistanceLoose_signalAndBackgroundBX1  = signalAndBackgroundBX1.Get("hSeedDistanceLoose")
    hSeedDistanceLoose_signalAndBackgroundBX2  = signalAndBackgroundBX2.Get("hSeedDistanceLoose")
    hSeedDistanceLoose_signalAndBackgroundBX3  = signalAndBackgroundBX3.Get("hSeedDistanceLoose")
    hSeedDistanceLoose_signalAndBackgroundBX4  = signalAndBackgroundBX4.Get("hSeedDistanceLoose")

    hSeedDistanceLoose_signalAndBackgroundBX1.Add(hSeedDistanceLoose_signalAndBackgroundBX2)
    hSeedDistanceLoose_signalAndBackgroundBX1.Add(hSeedDistanceLoose_signalAndBackgroundBX3)
    hSeedDistanceLoose_signalAndBackgroundBX1.Add(hSeedDistanceLoose_signalAndBackgroundBX4)
    
    hSeedDistanceLoose_signalAndBackgroundBX1.Rebin(4)
    
    
    
    #### seed energy, matched with signal
    hSeedDistanceMatchedBX1                     = signalFromBkgFileBX1.Get("hSeedDistance")
    hSeedDistanceMatchedBX2                     = signalFromBkgFileBX2.Get("hSeedDistance")
    hSeedDistanceMatchedBX3                     = signalFromBkgFileBX3.Get("hSeedDistance")
    hSeedDistanceMatchedBX4                     = signalFromBkgFileBX4.Get("hSeedDistance")

    hSeedDistanceMatchedBX1.Add(hSeedDistanceMatchedBX2)
    hSeedDistanceMatchedBX1.Add(hSeedDistanceMatchedBX3)
    hSeedDistanceMatchedBX1.Add(hSeedDistanceMatchedBX4)
    hSeedDistanceMatchedBX1.Rebin(4)
    
    FirstTH1   = [hSeedDistanceMatchedBX1, hSeedDistanceTight_signalAndBackgroundBX1, hSeedDistanceLoose_signalAndBackgroundBX1]
    
    LegendName = ["sig-only seeds", "sig+bkg, tight", "sig+bkg, loose"]
    PlotColor  = [kGray, 2, 4]
    xAxisName  = "Distance from analytical line in the x_{4}:x_{exit} plane [m]"
    yAxisName  = "Entries/BX"
    xrange1down = 0.0
    xrange1up   = 0.005
    yrange1down = 1e-1
    yrange1up   = 1e3
    logy        = True
    
    #### put the x axis label and y axis label
    for i in xrange(0,len(FirstTH1)):
        FirstTH1[i].Scale(1./4.0)
        #FirstTH1[i].GetYaxis().SetTitle(yAxisName)
        #FirstTH1[i].GetXaxis().SetTitle(xAxisName)
    
    DrawHists(FirstTH1, LegendName, PlotColor,xAxisName, yAxisName, xrange1down, xrange1up, yrange1down, yrange1up, directory+"/distanceDistributionSeedAndSignal_nTracks"+args.nTracks+"_"+plotSuffix, yline1low, yline1up, drawline, logy, latexName, latexName2, latexName3, leftLegend, doAtlas, doLumi, noRatio, do80, do59, drawPattern, logz, logx, latexName4)
    
    
    
    
    ### with fit, distance from signal+background, inclusive
    hSVDValues1_signalAndBackgroundBX1  = signalAndBackgroundBX1.Get("hSVDValues1")
    hSVDValues1_signalAndBackgroundBX2  = signalAndBackgroundBX2.Get("hSVDValues1")
    hSVDValues1_signalAndBackgroundBX3  = signalAndBackgroundBX3.Get("hSVDValues1")
    hSVDValues1_signalAndBackgroundBX4  = signalAndBackgroundBX4.Get("hSVDValues1")

    hSVDValues1_signalAndBackgroundBX1.Add(hSVDValues1_signalAndBackgroundBX2)
    hSVDValues1_signalAndBackgroundBX1.Add(hSVDValues1_signalAndBackgroundBX3)
    hSVDValues1_signalAndBackgroundBX1.Add(hSVDValues1_signalAndBackgroundBX4)
    
    #hSVDValues1_signalAndBackgroundBX1.Rebin(4)
    
    
    ### with fit, distance from signal+background, tight
    hSVDValues1Tight_signalAndBackgroundBX1  = signalAndBackgroundBX1.Get("hSVDValues1Tight")
    hSVDValues1Tight_signalAndBackgroundBX2  = signalAndBackgroundBX2.Get("hSVDValues1Tight")
    hSVDValues1Tight_signalAndBackgroundBX3  = signalAndBackgroundBX3.Get("hSVDValues1Tight")
    hSVDValues1Tight_signalAndBackgroundBX4  = signalAndBackgroundBX4.Get("hSVDValues1Tight")

    hSVDValues1Tight_signalAndBackgroundBX1.Add(hSVDValues1Tight_signalAndBackgroundBX2)
    hSVDValues1Tight_signalAndBackgroundBX1.Add(hSVDValues1Tight_signalAndBackgroundBX3)
    hSVDValues1Tight_signalAndBackgroundBX1.Add(hSVDValues1Tight_signalAndBackgroundBX4)
    
    #hSVDValues1Tight_signalAndBackgroundBX1.Rebin(4)
    
    
    ### with fit, distance from signal+background, loose
    hSVDValues1Loose_signalAndBackgroundBX1  = signalAndBackgroundBX1.Get("hSVDValues1Loose")
    hSVDValues1Loose_signalAndBackgroundBX2  = signalAndBackgroundBX2.Get("hSVDValues1Loose")
    hSVDValues1Loose_signalAndBackgroundBX3  = signalAndBackgroundBX3.Get("hSVDValues1Loose")
    hSVDValues1Loose_signalAndBackgroundBX4  = signalAndBackgroundBX4.Get("hSVDValues1Loose")

    hSVDValues1Loose_signalAndBackgroundBX1.Add(hSVDValues1Loose_signalAndBackgroundBX2)
    hSVDValues1Loose_signalAndBackgroundBX1.Add(hSVDValues1Loose_signalAndBackgroundBX3)
    hSVDValues1Loose_signalAndBackgroundBX1.Add(hSVDValues1Loose_signalAndBackgroundBX4)
    
    #hSVDValues1Loose_signalAndBackgroundBX1.Rebin(4)
    
    
    
    #### seed energy, matched with signal
    hSVDValues1MatchedBX1                     = signalFromBkgFileBX1.Get("hSVDValues1")
    hSVDValues1MatchedBX2                     = signalFromBkgFileBX2.Get("hSVDValues1")
    hSVDValues1MatchedBX3                     = signalFromBkgFileBX3.Get("hSVDValues1")
    hSVDValues1MatchedBX4                     = signalFromBkgFileBX4.Get("hSVDValues1")

    hSVDValues1MatchedBX1.Add(hSVDValues1MatchedBX2)
    hSVDValues1MatchedBX1.Add(hSVDValues1MatchedBX3)
    hSVDValues1MatchedBX1.Add(hSVDValues1MatchedBX4)
    #hSVDValues1MatchedBX1.Rebin(4)
    
    FirstTH1   = [hSVDValues1MatchedBX1, hSVDValues1Tight_signalAndBackgroundBX1, hSVDValues1Loose_signalAndBackgroundBX1]
    
    LegendName = ["sig-only seeds", "sig+bkg, tight", "sig+bkg, loose"]
    PlotColor  = [kGray, 2, 4]
    xAxisName  = "First SVD value"
    yAxisName  = "Entries/BX"
    xrange1down = 200
    xrange1up   = 250
    yrange1down = 1e-1
    yrange1up   = 1e5
    logy        = True
    #### put the x axis label and y axis label
    for i in xrange(0,len(FirstTH1)):
        FirstTH1[i].Scale(1./4.0)
    
    DrawHists(FirstTH1, LegendName, PlotColor,xAxisName, yAxisName, xrange1down, xrange1up, yrange1down, yrange1up, directory+"/svdValues1DistributionSeedAndSignal_nTracks"+args.nTracks+"_"+plotSuffix, yline1low, yline1up, drawline, logy, latexName, latexName2, latexName3, leftLegend, doAtlas, doLumi, noRatio, do80, do59, drawPattern, logz, logx, latexName4)
    
    
    
    
    
    
    ##### make distribution of the distance parameter
    
    ### with fit, distance from signal+background, inclusive
    hSVDValues2_signalAndBackgroundBX1  = signalAndBackgroundBX1.Get("hSVDValues2")
    hSVDValues2_signalAndBackgroundBX2  = signalAndBackgroundBX2.Get("hSVDValues2")
    hSVDValues2_signalAndBackgroundBX3  = signalAndBackgroundBX3.Get("hSVDValues2")
    hSVDValues2_signalAndBackgroundBX4  = signalAndBackgroundBX4.Get("hSVDValues2")

    hSVDValues2_signalAndBackgroundBX1.Add(hSVDValues2_signalAndBackgroundBX2)
    hSVDValues2_signalAndBackgroundBX1.Add(hSVDValues2_signalAndBackgroundBX3)
    hSVDValues2_signalAndBackgroundBX1.Add(hSVDValues2_signalAndBackgroundBX4)
    
    #hSVDValues2_signalAndBackgroundBX1.Rebin(4)
    
    
    ### with fit, distance from signal+background, tight
    hSVDValues2Tight_signalAndBackgroundBX1  = signalAndBackgroundBX1.Get("hSVDValues2Tight")
    hSVDValues2Tight_signalAndBackgroundBX2  = signalAndBackgroundBX2.Get("hSVDValues2Tight")
    hSVDValues2Tight_signalAndBackgroundBX3  = signalAndBackgroundBX3.Get("hSVDValues2Tight")
    hSVDValues2Tight_signalAndBackgroundBX4  = signalAndBackgroundBX4.Get("hSVDValues2Tight")

    hSVDValues2Tight_signalAndBackgroundBX1.Add(hSVDValues2Tight_signalAndBackgroundBX2)
    hSVDValues2Tight_signalAndBackgroundBX1.Add(hSVDValues2Tight_signalAndBackgroundBX3)
    hSVDValues2Tight_signalAndBackgroundBX1.Add(hSVDValues2Tight_signalAndBackgroundBX4)
    
    #hSVDValues2Tight_signalAndBackgroundBX1.Rebin(4)
    
    
    ### with fit, distance from signal+background, loose
    hSVDValues2Loose_signalAndBackgroundBX1  = signalAndBackgroundBX1.Get("hSVDValues2Loose")
    hSVDValues2Loose_signalAndBackgroundBX2  = signalAndBackgroundBX2.Get("hSVDValues2Loose")
    hSVDValues2Loose_signalAndBackgroundBX3  = signalAndBackgroundBX3.Get("hSVDValues2Loose")
    hSVDValues2Loose_signalAndBackgroundBX4  = signalAndBackgroundBX4.Get("hSVDValues2Loose")

    hSVDValues2Loose_signalAndBackgroundBX1.Add(hSVDValues2Loose_signalAndBackgroundBX2)
    hSVDValues2Loose_signalAndBackgroundBX1.Add(hSVDValues2Loose_signalAndBackgroundBX3)
    hSVDValues2Loose_signalAndBackgroundBX1.Add(hSVDValues2Loose_signalAndBackgroundBX4)
    
    #hSVDValues2Loose_signalAndBackgroundBX1.Rebin(4)
    
    
    
    #### seed energy, matched with signal
    hSVDValues2MatchedBX1                     = signalFromBkgFileBX1.Get("hSVDValues2")
    hSVDValues2MatchedBX2                     = signalFromBkgFileBX2.Get("hSVDValues2")
    hSVDValues2MatchedBX3                     = signalFromBkgFileBX3.Get("hSVDValues2")
    hSVDValues2MatchedBX4                     = signalFromBkgFileBX4.Get("hSVDValues2")

    hSVDValues2MatchedBX1.Add(hSVDValues2MatchedBX2)
    hSVDValues2MatchedBX1.Add(hSVDValues2MatchedBX3)
    hSVDValues2MatchedBX1.Add(hSVDValues2MatchedBX4)
    #hSVDValues2MatchedBX1.Rebin(4)
    
    FirstTH1   = [hSVDValues2MatchedBX1, hSVDValues2Tight_signalAndBackgroundBX1, hSVDValues2Loose_signalAndBackgroundBX1]
    
    LegendName = ["sig-only seeds", "sig+bkg, tight", "sig+bkg, loose"]
    PlotColor  = [kGray, 2, 4]
    xAxisName  = "Second SVD value"
    yAxisName  = "Entries/BX"
    xrange1down = 0.0
    xrange1up   = 0.1
    yrange1down = 1e-1
    yrange1up   = 1e3
    logy        = True
    #### put the x axis label and y axis label
    for i in xrange(0,len(FirstTH1)):
        FirstTH1[i].Scale(1./4.0)
    
    DrawHists(FirstTH1, LegendName, PlotColor,xAxisName, yAxisName, xrange1down, xrange1up, yrange1down, yrange1up, directory+"/svdValues2DistributionSeedAndSignal_nTracks"+args.nTracks+"_"+plotSuffix, yline1low, yline1up, drawline, logy, latexName, latexName2, latexName3, leftLegend, doAtlas, doLumi, noRatio, do80, do59, drawPattern, logz, logx, latexName4)
    
    
    
    
    
    
    
    
    ##### make distribution of the distance parameter
    
    ### with fit, distance from signal+background, inclusive
    hSVDValues3_signalAndBackgroundBX1  = signalAndBackgroundBX1.Get("hSVDValues3")
    hSVDValues3_signalAndBackgroundBX2  = signalAndBackgroundBX2.Get("hSVDValues3")
    hSVDValues3_signalAndBackgroundBX3  = signalAndBackgroundBX3.Get("hSVDValues3")
    hSVDValues3_signalAndBackgroundBX4  = signalAndBackgroundBX4.Get("hSVDValues3")

    hSVDValues3_signalAndBackgroundBX1.Add(hSVDValues3_signalAndBackgroundBX2)
    hSVDValues3_signalAndBackgroundBX1.Add(hSVDValues3_signalAndBackgroundBX3)
    hSVDValues3_signalAndBackgroundBX1.Add(hSVDValues3_signalAndBackgroundBX4)
    
    #hSVDValues3_signalAndBackgroundBX1.Rebin(4)
    
    
    ### with fit, distance from signal+background, tight
    hSVDValues3Tight_signalAndBackgroundBX1  = signalAndBackgroundBX1.Get("hSVDValues3Tight")
    hSVDValues3Tight_signalAndBackgroundBX2  = signalAndBackgroundBX2.Get("hSVDValues3Tight")
    hSVDValues3Tight_signalAndBackgroundBX3  = signalAndBackgroundBX3.Get("hSVDValues3Tight")
    hSVDValues3Tight_signalAndBackgroundBX4  = signalAndBackgroundBX4.Get("hSVDValues3Tight")

    hSVDValues3Tight_signalAndBackgroundBX1.Add(hSVDValues3Tight_signalAndBackgroundBX2)
    hSVDValues3Tight_signalAndBackgroundBX1.Add(hSVDValues3Tight_signalAndBackgroundBX3)
    hSVDValues3Tight_signalAndBackgroundBX1.Add(hSVDValues3Tight_signalAndBackgroundBX4)
    
    #hSVDValues3Tight_signalAndBackgroundBX1.Rebin(4)
    
    
    ### with fit, distance from signal+background, loose
    hSVDValues3Loose_signalAndBackgroundBX1  = signalAndBackgroundBX1.Get("hSVDValues3Loose")
    hSVDValues3Loose_signalAndBackgroundBX2  = signalAndBackgroundBX2.Get("hSVDValues3Loose")
    hSVDValues3Loose_signalAndBackgroundBX3  = signalAndBackgroundBX3.Get("hSVDValues3Loose")
    hSVDValues3Loose_signalAndBackgroundBX4  = signalAndBackgroundBX4.Get("hSVDValues3Loose")

    hSVDValues3Loose_signalAndBackgroundBX1.Add(hSVDValues3Loose_signalAndBackgroundBX2)
    hSVDValues3Loose_signalAndBackgroundBX1.Add(hSVDValues3Loose_signalAndBackgroundBX3)
    hSVDValues3Loose_signalAndBackgroundBX1.Add(hSVDValues3Loose_signalAndBackgroundBX4)
    
    #hSVDValues3Loose_signalAndBackgroundBX1.Rebin(4)
    
    
    
    #### seed energy, matched with signal
    hSVDValues3MatchedBX1                     = signalFromBkgFileBX1.Get("hSVDValues3")
    hSVDValues3MatchedBX2                     = signalFromBkgFileBX2.Get("hSVDValues3")
    hSVDValues3MatchedBX3                     = signalFromBkgFileBX3.Get("hSVDValues3")
    hSVDValues3MatchedBX4                     = signalFromBkgFileBX4.Get("hSVDValues3")

    hSVDValues3MatchedBX1.Add(hSVDValues3MatchedBX2)
    hSVDValues3MatchedBX1.Add(hSVDValues3MatchedBX3)
    hSVDValues3MatchedBX1.Add(hSVDValues3MatchedBX4)
    #hSVDValues3MatchedBX1.Rebin(4)
    #c1 = TCanvas()
    #hSVDValues3_signalAndBackgroundBX1.Draw()
    #c1.SaveAs("test.pdf")
    FirstTH1   = [hSVDValues3MatchedBX1, hSVDValues3Tight_signalAndBackgroundBX1, hSVDValues3Loose_signalAndBackgroundBX1]
    
    LegendName = ["sig-only seeds", "sig+bkg, tight", "sig+bkg, loose"]
    PlotColor  = [kGray, 2, 4]
    xAxisName  = "Third SVD value"
    yAxisName  = "Entries/BX"
    xrange1down = 0.0
    xrange1up   = 0.01
    yrange1down = 1e-1
    yrange1up   = 1e3
    logy        = True
    #### put the x axis label and y axis label
    for i in xrange(0,len(FirstTH1)):
        FirstTH1[i].Scale(1./4.0)
    
    #c1 = TCanvas()
    #hSVDValues3Loose_signalAndBackgroundBX1.Draw()
    #c1.SaveAs("test.pdf")
    
    DrawHists(FirstTH1, LegendName, PlotColor,xAxisName, yAxisName, xrange1down, xrange1up, yrange1down, yrange1up, directory+"/svdValues3DistributionSeedAndSignal_nTracks"+args.nTracks+"_"+plotSuffix, yline1low, yline1up, drawline, logy, latexName, latexName2, latexName3, leftLegend, doAtlas, doLumi, noRatio, do80, do59, drawPattern, logz, logx, latexName4)
    
    
    
    
    ##### py
    
    ### with fit, distance from signal+background, inclusive
    hSeedPy_signalAndBackgroundBX1  = signalAndBackgroundBX1.Get("hSeedPy")
    hSeedPy_signalAndBackgroundBX2  = signalAndBackgroundBX2.Get("hSeedPy")
    hSeedPy_signalAndBackgroundBX3  = signalAndBackgroundBX3.Get("hSeedPy")
    hSeedPy_signalAndBackgroundBX4  = signalAndBackgroundBX4.Get("hSeedPy")

    hSeedPy_signalAndBackgroundBX1.Add(hSeedPy_signalAndBackgroundBX2)
    hSeedPy_signalAndBackgroundBX1.Add(hSeedPy_signalAndBackgroundBX3)
    hSeedPy_signalAndBackgroundBX1.Add(hSeedPy_signalAndBackgroundBX4)
    
    #hSeedPy_signalAndBackgroundBX1.Rebin(4)
    
    
    ### with fit, distance from signal+background, tight
    hSeedPyTight_signalAndBackgroundBX1  = signalAndBackgroundBX1.Get("hSeedPyTight")
    hSeedPyTight_signalAndBackgroundBX2  = signalAndBackgroundBX2.Get("hSeedPyTight")
    hSeedPyTight_signalAndBackgroundBX3  = signalAndBackgroundBX3.Get("hSeedPyTight")
    hSeedPyTight_signalAndBackgroundBX4  = signalAndBackgroundBX4.Get("hSeedPyTight")

    hSeedPyTight_signalAndBackgroundBX1.Add(hSeedPyTight_signalAndBackgroundBX2)
    hSeedPyTight_signalAndBackgroundBX1.Add(hSeedPyTight_signalAndBackgroundBX3)
    hSeedPyTight_signalAndBackgroundBX1.Add(hSeedPyTight_signalAndBackgroundBX4)
    
    #hSeedPyTight_signalAndBackgroundBX1.Rebin(4)
    
    
    ### with fit, distance from signal+background, loose
    hSeedPyLoose_signalAndBackgroundBX1  = signalAndBackgroundBX1.Get("hSeedPyLoose")
    hSeedPyLoose_signalAndBackgroundBX2  = signalAndBackgroundBX2.Get("hSeedPyLoose")
    hSeedPyLoose_signalAndBackgroundBX3  = signalAndBackgroundBX3.Get("hSeedPyLoose")
    hSeedPyLoose_signalAndBackgroundBX4  = signalAndBackgroundBX4.Get("hSeedPyLoose")

    hSeedPyLoose_signalAndBackgroundBX1.Add(hSeedPyLoose_signalAndBackgroundBX2)
    hSeedPyLoose_signalAndBackgroundBX1.Add(hSeedPyLoose_signalAndBackgroundBX3)
    hSeedPyLoose_signalAndBackgroundBX1.Add(hSeedPyLoose_signalAndBackgroundBX4)
    
    #hSeedPyLoose_signalAndBackgroundBX1.Rebin(4)
    
    
    
    #### seed energy, matched with signal
    hSeedPyMatchedBX1                     = signalFromBkgFileBX1.Get("hSeedPy")
    hSeedPyMatchedBX2                     = signalFromBkgFileBX2.Get("hSeedPy")
    hSeedPyMatchedBX3                     = signalFromBkgFileBX3.Get("hSeedPy")
    hSeedPyMatchedBX4                     = signalFromBkgFileBX4.Get("hSeedPy")

    hSeedPyMatchedBX1.Add(hSeedPyMatchedBX2)
    hSeedPyMatchedBX1.Add(hSeedPyMatchedBX3)
    hSeedPyMatchedBX1.Add(hSeedPyMatchedBX4)
    #hSeedPyMatchedBX1.Rebin(4)
    
    FirstTH1   = [hSeedPyMatchedBX1, hSeedPyTight_signalAndBackgroundBX1, hSeedPyLoose_signalAndBackgroundBX1]
    
    LegendName = ["sig-only seeds", "sig+bkg, tight", "sig+bkg, loose"]
    PlotColor  = [kGray, 2, 4]
    xAxisName  = "p_{Y} [GeV]"
    yAxisName  = "Entries/BX"
    xrange1down = -0.02
    xrange1up   = 0.02
    yrange1down = 1e-1
    yrange1up   = 1e3
    logy        = True
    #### put the x axis label and y axis label
    for i in xrange(0,len(FirstTH1)):
        FirstTH1[i].Scale(1./4.0)
    
    DrawHists(FirstTH1, LegendName, PlotColor,xAxisName, yAxisName, xrange1down, xrange1up, yrange1down, yrange1up, directory+"/SeedPyDistributionSeedAndSignal_nTracks"+args.nTracks+"_"+plotSuffix, yline1low, yline1up, drawline, logy, latexName, latexName2, latexName3, leftLegend, doAtlas, doLumi, noRatio, do80, do59, drawPattern, logz, logx, latexName4)
    
    
    
    
    
    
    ##### py
    
    ### with fit, distance from signal+background, inclusive
    hXExitYExit_signalAndBackgroundBX1  = signalAndBackgroundBX1.Get("hXExitYExit")
    hXExitYExit_signalAndBackgroundBX2  = signalAndBackgroundBX2.Get("hXExitYExit")
    hXExitYExit_signalAndBackgroundBX3  = signalAndBackgroundBX3.Get("hXExitYExit")
    hXExitYExit_signalAndBackgroundBX4  = signalAndBackgroundBX4.Get("hXExitYExit")

    hXExitYExit_signalAndBackgroundBX1.Add(hXExitYExit_signalAndBackgroundBX2)
    hXExitYExit_signalAndBackgroundBX1.Add(hXExitYExit_signalAndBackgroundBX3)
    hXExitYExit_signalAndBackgroundBX1.Add(hXExitYExit_signalAndBackgroundBX4)
    
    #hXExitYExit_signalAndBackgroundBX1.Rebin(4)
    
    
    
    
    
    #### seed energy, matched with signal
    hXExitYExitMatchedBX1                     = signalFromBkgFileBX1.Get("hXExitYExit")
    hXExitYExitMatchedBX2                     = signalFromBkgFileBX2.Get("hXExitYExit")
    hXExitYExitMatchedBX3                     = signalFromBkgFileBX3.Get("hXExitYExit")
    hXExitYExitMatchedBX4                     = signalFromBkgFileBX4.Get("hXExitYExit")

    hXExitYExitMatchedBX1.Add(hXExitYExitMatchedBX2)
    hXExitYExitMatchedBX1.Add(hXExitYExitMatchedBX3)
    hXExitYExitMatchedBX1.Add(hXExitYExitMatchedBX4)
    #hXExitYExitMatchedBX1.Rebin(4)
    
    FirstTH1   = [hXExitYExit_signalAndBackgroundBX1]
    
    LegendName = ["sig-only seeds"]
    PlotColor  = [kGray]
    xAxisName  = "x_{Exit} [mm]"
    yAxisName  = "y_{Exit} [mm]"
    xrange1down = 0.0
    xrange1up   = 100.0
    yrange1down = -5
    yrange1up   = 5
    logy        = False
    #### put the x axis label and y axis label
    for i in xrange(0,len(FirstTH1)):
        FirstTH1[i].Scale(1./4.0)
    
    DrawHists(FirstTH1, LegendName, PlotColor,xAxisName, yAxisName, xrange1down, xrange1up, yrange1down, yrange1up, directory+"/XExitYExitDistributionSignalAndBackground_nTracks"+args.nTracks+"_"+plotSuffix, yline1low, yline1up, drawline, logy, latexName, latexName2, latexName3, leftLegend, doAtlas, doLumi, noRatio, do80, do59, "COLZ", logz, logx, latexName4)
    


if __name__=="__main__":
    start = time.time()
    main()
    print("The time taken: ", time.time() - start, " s")
