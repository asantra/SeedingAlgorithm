import csv, numpy, array

### read the dEdX vs E curve for the electron
results = csv.reader(open("dEdXSilicon_ForElectron.csv"), delimiter=",")
xdEdX = array.array('f',[]); ydEdX = array.array('f', [])
for eachRow in results:
    #print(eachRow)
    xdEdX.append(float(eachRow[0]))
    ydEdX.append(float(eachRow[1]))
mevEnergy = 0.01
interpolateddEdX = numpy.interp(mevEnergy, xdEdX, ydEdX)
print(interpolateddEdX)
