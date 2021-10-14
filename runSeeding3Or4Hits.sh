#! /bin/zsh
postfix=""
needFitList="F"
needELaser=${1}
alias python3=/usr/local/bin/python3.9
for needFitIter in ${needFitList}; do
    if [[ ${needFitIter} == "F" ]]
    then
        postfix="WithFit3or4HitsTracksAndDistanceCut"
        needFit=1
    else
        postfix="WithoutFit3or4HitsTracksAndDistanceCut"
        needFit=0
    fi

    echo "!!!!!!#### The postfix to the root file: "${postfix}" and needFit: "${needFit}
    #signalTracksList="1 5 10 20 30 50 80 100 130 150 170 185 200 220"
    #250 300 500 750 
    signalTracksList="1000 1500 2000 2500 3000 3500"
    #signalTracksList="1 4 7 10 13 16"
    bxList="1 2 3 4"
    particleList="Positron"
     
    
    ### working with the signal+background case, 250 300  3000 3500
    ### 500 750 1000 1500 2000 
    for signalTracks in 250 300 500 750 1000 1500 2000 2500 3000; do
        echo "###############################################################"
        echo "########## working for nTracks: ${signalTracks} case ##########"
        echo "###############################################################"
        
        for bx in 1; do
            echo "######!!! BX: "${bx}"  !!!#########"
            for particle in ${particleList}; do
                echo "#####!!! Particle : "${particle}" !!!####"
                
                if [[ $1 == "E" ]]
                then
                
                    #### e+laser setup
                    
                    python3 findSeed.py -l BkgEBeam_Signal${particle}hics3000nm_jeti40_122020_9550dac4_BX${bx}_SignalTracks${signalTracks}_trackInfoClean.txt -s -f ${needFit} -p ${particle}
                    python3 findSeed.py -l BkgEBeam_Signal${particle}hics3000nm_jeti40_122020_9550dac4_BX${bx}_SignalTracks${signalTracks}_trackInfoClean.txt -f ${needFit} -p ${particle}
                
                
                    #### change the root file name prefix

                    mv seedingInformationFiles/seedingInformation_BkgEBeam_Signal${particle}hics3000nm_jeti40_122020_9550dac4_BX${bx}_SignalTracks${signalTracks}_trackInfoClean_VariableEnergyCut_OnlySignal_${particle}Side.root seedingInformationFiles/seedingInformation_BkgEBeam_Signal${particle}hics3000nm_jeti40_122020_9550dac4_BX${bx}_SignalTracks${signalTracks}_trackInfoClean_VariableEnergyCut_OnlySignal_${particle}Side_${postfix}.root

                

                    mv seedingInformationFiles/seedingInformation_BkgEBeam_Signal${particle}hics3000nm_jeti40_122020_9550dac4_BX${bx}_SignalTracks${signalTracks}_trackInfoClean_VariableEnergyCut_SignalAndBackground_${particle}Side.root seedingInformationFiles/seedingInformation_BkgEBeam_Signal${particle}hics3000nm_jeti40_122020_9550dac4_BX${bx}_SignalTracks${signalTracks}_trackInfoClean_VariableEnergyCut_SignalAndBackground_${particle}Side_${postfix}.root
                    
                    
                else
                    #### g+laser setup
    
                    python3 findSeedGLaser.py -l BkgGBeam_Signal${particle}bppp3000nmOr5000nm_BX${bx}_SignalTracks${signalTracks}_trackInfoClean.txt -s -f ${needFit} -p ${particle}
                    python3 findSeedGLaser.py -l BkgGBeam_Signal${particle}bppp3000nmOr5000nm_BX${bx}_SignalTracks${signalTracks}_trackInfoClean.txt -f ${needFit} -p ${particle}
                
                
                    #### change the root file name prefix

                    mv seedingInformationFiles/seedingInformation_BkgGBeam_Signal${particle}bppp3000nmOr5000nm_BX${bx}_SignalTracks${signalTracks}_trackInfoClean_VariableEnergyCut_OnlySignal_${particle}Side.root seedingInformationFiles/seedingInformation_BkgGBeam_Signal${particle}bppp3000nmOr5000nm_BX${bx}_SignalTracks${signalTracks}_trackInfoClean_VariableEnergyCut_OnlySignal_${particle}Side_${postfix}.root

                

                    mv seedingInformationFiles/seedingInformation_BkgGBeam_Signal${particle}bppp3000nmOr5000nm_BX${bx}_SignalTracks${signalTracks}_trackInfoClean_VariableEnergyCut_SignalAndBackground_${particle}Side.root seedingInformationFiles/seedingInformation_BkgGBeam_Signal${particle}bppp3000nmOr5000nm_BX${bx}_SignalTracks${signalTracks}_trackInfoClean_VariableEnergyCut_SignalAndBackground_${particle}Side_${postfix}.root
                fi
            done
        done
    done
    
    
# #     ### work with the background only file


#     ### work with the new background only file
#     if [[ $1 == "E" ]]
#     then
#         for bx in 1 2 3 4; do
#             for particle in ${particleList}; do
#                 python3 findSeed.py -l ePlusLaserBkgKaptonWindowNewSamplesMarch62021_DividedByBX${bx}_trackInfoClean.txt -f ${needFit} -p ${particle}
#                 mv seedingInformationFiles/seedingInformation_ePlusLaserBkgKaptonWindowNewSamplesMarch62021_DividedByBX${bx}_trackInfoClean_VariableEnergyCut_SignalAndBackground_${particle}Side.root seedingInformationFiles/seedingInformation_ePlusLaserBkgKaptonWindowNewSamplesMarch62021_DividedByBX${bx}_trackInfoClean_VariableEnergyCut_SignalAndBackground_${particle}Side_${postfix}.root
#             done
#         done
#     else
#         for bx in 1 2 3 4; do
#             python3 findSeedGLaser.py -l EBeamOnlyWIS_DividedByBX${bx}_trackInfoClean.txt -f ${needFit}
#             mv seedingInformationFiles/seedingInformation_EBeamOnlyWIS_DividedByBX${bx}_trackInfoClean_VariableEnergyCut_SignalAndBackground_PositronSide.root seedingInformationFiles/seedingInformation_EBeamOnlyWIS_DividedByBX${bx}_trackInfoClean_VariableEnergyCut_SignalAndBackground_PositronSide_${postfix}.root
#         done
#     fi
        

done
