#!/usr/bin/python
import os
import math
import subprocess
import array
import numpy as np
import ast
import ROOT
from ROOT import TFile, TH1D, TH2D, TCanvas, TPolyLine, TPolyMarker, TRandom, TPad, TF1, TGraph

ROOT.gROOT.SetBatch(1)
ROOT.gStyle.SetOptFit(0);
ROOT.gStyle.SetOptStat(0);
ROOT.gStyle.SetPadBottomMargin(0.15)
ROOT.gStyle.SetPadLeftMargin(0.16)
ROOT.gStyle.SetLineWidth(1)
# storage =  ROOT.gSystem.ExpandPathName("$STORAGEDIR")
storage =  os.path.expandvars("$STORAGEDIR")


### globals
m2mm        = 1000
z4          = 4176.5125 ## [mm]
z0          = 2748 ## [mm] ## 2050
Ldipole     = 1396 ## [mm] --> check with Sasha
B           = 1.0 ## [T]
LB          = 1029 ## [mm]
zDipoleExit = z0-(Ldipole-LB)/2
D4          = (z4-zDipoleExit) ## [mm]
LB2         = LB*LB

ZE = LB
Z0 = (LB+Ldipole)/2
Z4 = LB+D4

histos = {}



def getR(p):
   R = p/(0.3*B)*m2mm
   return R

def getR2(p):
   R = getR(p)
   return R*R

def getXTangent(ZT,p):
   R2 = getR2(p)
   XT = (R2/LB-ZT)*(LB/math.sqrt(R2-LB2))
   return XT

def getx(ZT,p,side="Pside"):
   XT = getXTangent(ZT,p)
   R = getR(p)
   x = R-XT
   # return x if(side=="Eside") else -x
   return x

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

def Book():
   tf = TFile("/Users/hod/Downloads/seedingInformation_hics_165gev_w0_5000nm_provisional_10xi3_6cee466a_trackInfoClean_VariableEnergyCut_OnlySignal_PositronSide.root","READ")
   histos.update( { "h_x4_vs_xDipole" : tf.Get("hXLayer4XDipole") } )
   histos["h_x4_vs_xDipole"].SetDirectory(0)
   
   histos.update( { "h_distance" : TH1D("h_distance",";Distance to analytical line [mm];Number of seeds",200,0,100) } )
   histos["h_distance"].SetDirectory(0)
   print(histos)

def Graph():
   xmin = histos["h_x4_vs_xDipole"].GetXaxis().GetXmin()
   xmax = histos["h_x4_vs_xDipole"].GetXaxis().GetXmax()
   xlist = []
   ylist = []
   for p in np.arange(1.3, 11, 0.5):
      R     = getR(p)
      xExit = getx(ZE,p)
      x0    = getx(Z0,p)
      x4    = getx(Z4,p)
      # print("p:",p,"[GeV], R:",R,"[mm], xExit:",xExit,"[mm], x0:",x0,"[mm], x4:",x4,"[mm]")
      xlist.append(x0)
      ylist.append(x4)
   x  = array.array('d', xlist)
   y  = array.array('d', ylist)
   g = TGraph(len(x),x,y)
   g.SetLineColor(ROOT.kRed)
   g.SetMarkerColor(ROOT.kRed)
   g.SetMarkerStyle(20)
   g.SetMarkerSize(0.3)
   g.GetXaxis().SetTitle(histos["h_x4_vs_xDipole"].GetXaxis().GetTitle())
   g.GetYaxis().SetTitle(histos["h_x4_vs_xDipole"].GetYaxis().GetTitle())
   n = len(xlist)
   r1 = [xlist[0],ylist[0]]
   r2 = [xlist[n-1],ylist[n-1]]
   return g,r1,r2





###################################################
###################################################

## run
Book()
g,r1,r2 = Graph()

## draw
cnv = TCanvas("cnv","",500,500)
ROOT.gPad.SetTicks(1,1)
histos["h_x4_vs_xDipole"].Draw("col")
g.Draw("PC same")
cnv.SaveAs("h_x4_vs_xDipole.pdf(")

## check
for x0 in range(1,histos["h_x4_vs_xDipole"].GetNbinsX()+1):
   for x4 in range(1,histos["h_x4_vs_xDipole"].GetNbinsY()+1):
      if(histos["h_x4_vs_xDipole"].GetBinContent(x0,x4)>0):
         x0Test = histos["h_x4_vs_xDipole"].GetXaxis().GetBinCenter(x0)
         x4Test = histos["h_x4_vs_xDipole"].GetYaxis().GetBinCenter(x4)
         dTest = Distance(x0Test,x4Test,r1,r2)
         print("dTest=",dTest)
         histos["h_distance"].Fill(dTest)
cnv = TCanvas("cnv","",500,500)
ROOT.gPad.SetTicks(1,1)
ROOT.gPad.SetLogy()
histos["h_distance"].Draw("hist")
cnv.SaveAs("h_x4_vs_xDipole.pdf)")

