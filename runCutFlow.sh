#! /bin/bash

python3 findSeedCutFlowTrackPy.py -l EBeamOnlyNewSamples_AllBX_trackInfoClean.txt
# python3 findSeedCutFlowDistanceD.py -l EBeamOnlyNewSamples_AllBX_trackInfoClean.txt
# python3 findSeedCutFlowFitParameter.py -l EBeamOnlyNewSamples_AllBX_trackInfoClean.txt

python3 findSeedCutFlowTrackPy.py -l list_root_hics_165gev_w0_5000nm_provisional_10xi3_6cee466a_trackInfoClean.txt -s
# python3 findSeedCutFlowDistanceD.py -l list_root_hics_165gev_w0_5000nm_provisional_10xi3_6cee466a_trackInfoClean.txt -s
# python3 findSeedCutFlowFitParameter.py -l list_root_hics_165gev_w0_5000nm_provisional_10xi3_6cee466a_trackInfoClean.txt -s

# python3 findSeedCutFlowTrackPy.py -l gPlusLaserBkgNewSamples_AllBX_trackInfoClean.txt
# python3 findSeedCutFlowDistanceD.py -l gPlusLaserBkgNewSamples_AllBX_trackInfoClean.txt
# python3 findSeedCutFlowFitParameter.py -l gPlusLaserBkgNewSamples_AllBX_trackInfoClean.txt

# python3 findSeedCutFlowTrackPy.py -l gPlusLaserBkgNewSamples_AllBX_trackInfoClean.txt -p Electron
# python3 findSeedCutFlowDistanceD.py -l gPlusLaserBkgNewSamples_AllBX_trackInfoClean.txt -p Electron
# python3 findSeedCutFlowFitParameter.py -l gPlusLaserBkgNewSamples_AllBX_trackInfoClean.txt -p Electron
