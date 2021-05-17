#### separate out bkg track info according to the bx number
#### usage: python separateBkgBX.py <inTextFile> <BXNumberWanted>

import os, sys, glob

def main():
    bkgFileName = sys.argv[1]
    bxWanted    = int(sys.argv[2])
    inFile      = open(bkgFileName)
    
    
    outDir      = "NewSamplesEBeamOnlyFilesMar62021/"
    
    if not os.path.exists(outDir):
        os.makedirs(outDir)
        
    outFileName = "ePlusLaserBkgKaptonWindowNewSamplesMarch62021_DividedByBX"+str(bxWanted)+"_trackInfoClean.txt" 
    
    outFile     = open(outDir+"/"+outFileName,"w")
    outFile.write('### bxNumber << pdg << track_id << det_id << xx << yy << eneg << ev_weight << vtx_x << vtx_y << vtx_z\n')
    
    for line in inFile.readlines():
        line     = line.rstrip()
        if "#" in line: continue
        eachWord = line.split()
        bxNumber = int(eachWord[0])
        if(bxNumber == bxWanted):
            #print(bxNumber)
            #print(bxWanted)
            outFile.write(line+'\n')
            
    outFile.close()
        
    
    
if __name__=="__main__":
    main()
