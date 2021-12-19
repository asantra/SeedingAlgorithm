#!/usr/bin/python3
import os
import math
import subprocess
import array
import numpy as np
import ROOT
from ROOT import TPolyMarker3D, TPolyLine3D, TPolyLine, TCanvas, TPad, TView, TLatex, TLegend, TGaxis
import pandas as pd


ROOT.gROOT.SetBatch(1)
ROOT.gStyle.SetOptFit(0);
ROOT.gStyle.SetOptStat(0);
ROOT.gStyle.SetPadBottomMargin(0.15)
ROOT.gStyle.SetPadLeftMargin(0.16)

### globals
m2mm        = 1000
mm2m        = 1./1000

### these are in mm
xDipoleWidth  = 330.0
yDipoleHeight = 108.0
zDipoleExit   = 2748.0
### length of the dipole in meters, needed in meters because pf the p = 0.3*B*R formula
LB       = 1.029
LD       = 1.396
zDipoleActiveExit = zDipoleExit - (LD - LB)*m2mm/2  ### in mm, to be consistent with Sasha's values

### 1 Tesla magnetic field used for now
B        = 1.0


z4       = 4176.5125 ## [mm]

D4       = (z4-zDipoleActiveExit)*mm2m ## [mm]
LB2      = LB*LB

### Z exit 
ZE       = LB

Z0       = (LB+LD)/2
Z4       = (LB+D4)

xBoundaries = {}


### read the data from Sasha's spreadhseet (need to convert to xlsx!)
xl = pd.ExcelFile("luxe_tracker_sensor_position.xlsx")
print(xl.sheet_names)
df = xl.parse("Sheet1")
print("Head of dataframe is:", df.head())




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
	
def yofz(r1,r2,z):
	dz = r2[2]-r1[2]
	dy = r2[1]-r1[1]
	if(dz==0):
		print("ERROR in yofz: dz=0")
		quit()
	a = dy/dz
	b = r1[1]-a*r1[2]
	y = a*z+b
	return y

def zofx(r1,r2,x):
	dz = r2[2]-r1[2]
	dx = r2[0]-r1[0]
	if(dx==0):
		print("ERROR in zofx: dx=0")
		quit()
	a = dz/dx
	b = r1[2]-a*r1[0]
	z = a*x+b
	return z



def GetLayerZ(detid,layerid):
   row = df[ (df['detid']==detid) & (df['layerid']==layerid) ]
   # print("row=",row)
   z = row["translation_z"].tolist()[0]
   return z


### get the z of the layers, global variables, can be used by other codes
z1inner = GetLayerZ(1000,0)
z2inner = GetLayerZ(1000,2)
z3inner = GetLayerZ(1000,4)
z4inner = GetLayerZ(1000,6)
z1outer = GetLayerZ(1000,1)
z2outer = GetLayerZ(1000,3)
z3outer = GetLayerZ(1000,5)
z4outer = GetLayerZ(1000,7)



   
def GetSensorXY(detid,layerid):
   row = df[ (df['detid']==detid) & (df['layerid']==layerid) ]
   # print("row=",row)
   x = row["translation_x"].tolist()[0]
   y = row["translation_y"].tolist()[0]
   return x,y
   
#def GetSensorSize(detid,layerid,axis):
   #row = df[ (df['detid']==detid) & (df['layerid']==layerid) ]
   ## print("row=",row)
   #size = row["size_x"].tolist()[0] if(axis=="x") else row["size_y"].tolist()[0]
   #return size

def GetDipole(color=ROOT.kBlack):
   xarr = array.array('d', [-xDipoleWidth/2,-xDipoleWidth/2,+xDipoleWidth/2,+xDipoleWidth/2,-xDipoleWidth/2])
   yarr = array.array('d', [-yDipoleHeight/2,+yDipoleHeight/2,+yDipoleHeight/2,-yDipoleHeight/2,-yDipoleHeight/2])
   zarr = array.array('d', [zDipoleExit,zDipoleExit,zDipoleExit,zDipoleExit,zDipoleExit])
   n = len(xarr)
   dipole = TPolyLine3D(n,xarr,yarr,zarr)
   dipole.SetLineColor(color)
   return dipole

### get the x,y position and x,y size of staves
def GetSensorSize(detid, layerid):
   row = df[ (df['detid']==detid) & (df['layerid']==layerid) ]
   # print("row=",row)
   x = row["translation_x"].tolist()[0]
   y = row["translation_y"].tolist()[0]
   z = row["translation_z"].tolist()[0]
   xsize = row["size_x"].tolist()[0]
   ysize = row["size_y"].tolist()[0]
   return x, y, xsize, ysize


### global variable, can be called by other functions
xInnerStaveLastChip,  yInnerStaveLastChip,  xSizeInnerStaveLastChip,  ySizeInnerStaveLastChip  = GetSensorSize(1008, 0)
xOuterStaveFirstChip, yOuterStaveFirstChip, xSizeOuterStaveFirstChip, ySizeOuterStaveFirstChip = GetSensorSize(1000, 7)


    
def GetSensor(detid,layerid,color=ROOT.kGreen+2):
   # print("looking for detid="+str(detid)+" and layerid="+str(layerid))
   row = df[ (df['detid']==detid) & (df['layerid']==layerid) ]
   # print("row=",row)
   x = row["translation_x"].tolist()[0]
   y = row["translation_y"].tolist()[0]
   z = row["translation_z"].tolist()[0]
   xsize = row["size_x"].tolist()[0]
   ysize = row["size_y"].tolist()[0]
   xhalf = xsize/2
   yhalf = ysize/2   
   # print("x=",x)
   # print("y=",y)
   # print("z=",z)
   # print("xsize=",xsize)
   # print("ysize=",ysize)
   xarr = array.array('d', [x-xhalf,x-xhalf,x+xhalf,x+xhalf,x-xhalf])
   yarr = array.array('d', [y-yhalf,y+yhalf,y+yhalf,y-yhalf,y-yhalf])
   zarr = array.array('d', [z,      z,      z,      z,      z])
   n = len(xarr)
   sensor = TPolyLine3D(n,xarr,yarr,zarr)
   sensor.SetLineColor(color)
   return sensor


def GetSensorXBoundaries(detid,layerid):
   # print("looking for detid="+str(detid)+" and layerid="+str(layerid))
   row = df[ (df['detid']==detid) & (df['layerid']==layerid) ]
   # print("row=",row)
   x = row["translation_x"].tolist()[0]
   y = row["translation_y"].tolist()[0]
   z = row["translation_z"].tolist()[0]
   xsize = row["size_x"].tolist()[0]
   ysize = row["size_y"].tolist()[0]
   xhalf = xsize/2
   return x-xhalf, x+xhalf

def GetTrackLine(rlist,color=ROOT.kRed):
   n = len(rlist)
   track = TPolyLine3D(n+1)
   track.SetLineColor(color)
   for i in range(n): track.SetPoint(i,rlist[i][0],rlist[i][1],rlist[i][2])
   return track
   
def GetTrackPoints(rlist,color=ROOT.kOrange):
   n = len(rlist)
   track = TPolyMarker3D(n+1)
   track.SetMarkerColor(color)
   for i in range(n): track.SetPoint(i,rlist[i][0],rlist[i][1],rlist[i][2])
   return track

def GetExtendedTrackLine(rlist,color=ROOT.kOrange):
   n = len(rlist)
   track = TPolyLine3D(n+1)
   track.SetLineColor(color)
   x = xofz(rlist[1],rlist[0],zDipoleExit)
   y = yofz(rlist[1],rlist[0],zDipoleExit)
   rDipoleExit = [x,y,zDipoleExit]
   rlistextended = [rlist[1], rDipoleExit]
   for i in range(n): track.SetPoint(i,rlistextended[i][0],rlistextended[i][1],rlistextended[i][2])
   return track


##############################################################################
##############################################################################
##############################################################################

def main():
    

    ### get the dipole
    dipole = GetDipole()

    
    rlist  = [[464.685, -2.89186,z1inner],[549.413, -0.553007,z4inner]]
    trkpoints = GetTrackPoints(rlist)
    trackline = GetTrackLine(rlist)
    extendedline = GetExtendedTrackLine(rlist)

    ### get the sensors
    sensors = []
    
    for index, row in df.iterrows():
        detid   = row["detid"]
        layerid = row["layerid"]
        sensors.append( GetSensor(detid,layerid) )

    print(xBoundaries)
    ### draw
    cnv = TCanvas("cnv","",500,500)
    view = TView.CreateView(1)
    view.SetRange(-600,-500,2500, +600,+500,4500)
    view.ShowAxis()
    for sensor in sensors: sensor.Draw()
    dipole.Draw()
    extendedline.Draw()
    trackline.Draw()
    trkpoints.Draw()
    
    cnv.SaveAs("sensors.pdf")
    cnv.SaveAs("sensors.root")

if __name__=="__main__":
    main()







