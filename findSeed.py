import os
import sys
import time
import pprint
import math
from ROOT import *
import array
from makeTrackDiagrams import *

EseedMin = 1.0; #// GeV
EseedMax = 16.0; #// GeV
LB = 1.029;
# get the unit vector along one vector
def rUnit2(r1, r2):
    r = (r2-r1).Unit()
    return r


# def check_clusters(i1, i4, side):
    # yAbsMargins = 0.02 #// cm (a "road" of 200 microns around the line between r4 and r1)
    # xAbsMargins = 0.02 #// cm (a "road" of 200 microns around the line between r4 and r1)
    #r1min[3] = { cached_clusters_xyz["x_L1_"+side][i1]-xAbsMargins, cached_clusters_xyz["y_L1_"+side][i1]-yAbsMargins, cached_clusters_xyz["z_L1_"+side][i1] }
    #r1max[3] = { cached_clusters_xyz["x_L1_"+side][i1]+xAbsMargins, cached_clusters_xyz["y_L1_"+side][i1]+yAbsMargins, cached_clusters_xyz["z_L1_"+side][i1] }
    #r4min[3] = { cached_clusters_xyz["x_L4_"+side][i4]-xAbsMargins, cached_clusters_xyz["y_L4_"+side][i4]-yAbsMargins, cached_clusters_xyz["z_L4_"+side][i4] }
    #r4max[3] = { cached_clusters_xyz["x_L4_"+side][i4]+xAbsMargins, cached_clusters_xyz["y_L4_"+side][i4]+yAbsMargins, cached_clusters_xyz["z_L4_"+side][i4] }

    # /// check possible clusters in layer 2
    #y2min = yofz(r1min,r4min,float(szlayers["L2"]))
    #y2max = yofz(r1max,r4max,float(szlayers["L2"]))
    #x2min = xofz(r1min,r4min,float(szlayers["L2"]))
    #x2max = xofz(r1max,r4max,float(szlayers["L2"]))
    #accept2 = False
    # for i2 in xrange(0, len(cached_clusters_xyz["x_L2_"+side])):
    #acceptyz = ( cached_clusters_xyz["y_L2_"+side][i2]>=y2min and cached_clusters_xyz["y_L2_"+side][i2]<=y2max )
    # if(not acceptyz): continue
    #acceptxz = ( cached_clusters_xyz["x_L2_"+side][i2]>=x2min and cached_clusters_xyz["x_L2_"+side][i2]<=x2max )
    # if(not acceptxz): continue
    #accept2 = True
    # break

    # if(not accept2) return False

    # y3min = yofz(r1min,r4min,(float)szlayers["L3"])
    # y3max = yofz(r1max,r4max,(float)szlayers["L3"])
    # x3min = xofz(r1min,r4min,(float)szlayers["L3"])
    # x3max = xofz(r1max,r4max,(float)szlayers["L3"])
    #accept3 = False
    # for i3 in xrange(0, len(cached_clusters_xyz["x_L3_"+side])):
    #acceptyz = ( cached_clusters_xyz["y_L3_"+side][i3]>=y3min and cached_clusters_xyz["y_L3_"+side][i3]<=y3max )
    # if(not acceptyz) continue
    #acceptxz = ( cached_clusters_xyz["x_L3_"+side][i3]>=x3min and cached_clusters_xyz["x_L3_"+side][i3]<=x3max )
    # if(not acceptxz) continue
    #accept3 = True
    # break

    # if(not accept3) return False

    # return True

# making the seeds
def makeseed(r1, r4, p):

    if(abs(r1[0]) >= abs(r4[0])):
        return False  # // |x1| must be smaller than |x4|
    if(r1[0] > 0 and r4[0] < 0):
        return False
    if(r1[0] < 0 and r4[0] > 0):
        return False
    if(r1[2] == r4[2]):
        return False  # // if z1=z4...
    yDipoleExitAbsMax = 0.7  # // cm
    xDipoleExitAbsMin = 4.  # // cm
    xDipoleExitAbsMax = 30.  # // cm
    yDipoleExit = yofz(r1, r4, zDipoleExit)
    xDipoleExit = xofz(r1, r4, zDipoleExit)
    if(abs(yDipoleExit) > yDipoleExitAbsMax):
        return False  # // the track should point to |y|<~0.2 at the dipole exit
    if(abs(xDipoleExit) < xDipoleExitAbsMin):
        return False  # // the track should point to |x|<~1.0 at the dipole exit
    if(abs(xDipoleExit) > xDipoleExitAbsMax):
        return False
    # if(not check_clusters(i1,i4,side)):     return False #// minimum one cluster at layer 2 and one at layer 3

    rnd = TRandom()
    rnd.SetSeed()
    posneg = rnd.Uniform(-1, +1)
    pxgaus = rnd.Gaus(7.2e-4, 5.0e-4)

    x0 = 0
    z0 = zofx(r1, r4, x0)
    xExit = abs(xofz(r1, r4, zDipoleExit))
    H = abs((zDipoleExit-z0))
    R = H*(LB)/xExit + xExit  # // look this up in my slides
    P = 0.3*B*R
    # // if(i4==0 and side=="Eside") cout << "z0=" << z0 << ", xExit=" << xExit << ", H=" << H << ", R=" << R << ", P=" << P << endl
    v1 = TVector2()
    v2 = TVector2()
    v1(r1[2], r1[1])
    v2(r4[2], r4[1])
    u = rUnit2(v1, v4)
    uz = u.X()
    uy = u.Y()
    px = -pxgaus
    py = P*uy
    pz = P*uz
    p.SetPxPyPzE(px, py, pz, math.sqrt(px*px + py*py + pz*pz + meGeV2))
    # // if(i4==0 and side=="Eside") cout << "px=" << px << ", py=" << py << ", pz=" << pz << endl
    # // cout << "side=" << side << ", px=" << px << ", py=" << py << ", pz=" << pz << endl
    # EseedMax = (process=="bppp") ? EseedMaxBPPP : EseedMaxTRIDENT	// GeV
    if(p.E() < EseedMin or p.E() > EseedMax): return False

    return True



# get the layer from the stave number
def getLayerId(stave):
    layerId = 0
    if ((stave == 1000) or (stave == 1001) or (stave == 1008) or (stave == 1009)):
        layerId = 1
    elif ((stave == 1002) or (stave == 1003) or (stave == 1010) or (stave == 1011)):
        layerId = 2
    elif ((stave == 1004) or (stave == 1005) or (stave == 1012) or (stave == 1013)):
        layerId = 3
    else:
        layerId = 4

    return layerId


def main():
    # give the input text file containing all the track information
    inputTrackInfo = open(sys.argv[1])
    # all the track info is in the following list
    position = []
    # get the information from the text files
    for lines in inputTrackInfo.readlines():
        lines = lines.rstrip()
        eachWord = lines.split()
        bxNumber = int(eachWord[0])
        layerNumber = getLayerId(int(eachWord[3]))
        trackId = int(eachWord[2])
        position.append([bxNumber, trackId, layerNumber, float(
            eachWord[4]), float(eachWord[5]), float(eachWord[6])])

    # separate each bx now
    eachBXValue = []
    bxCounter = 1
    for tracks in position:
        if tracks[0] == bxCounter:
            eachBXValue.append(tracks)

    pprint.pprint(eachBXValue)


if __name__ == "__main__":
    main()
