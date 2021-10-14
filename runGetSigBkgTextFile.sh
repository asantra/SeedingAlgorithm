#! /bin/bash
signalTracksList="250 300 500 750 1000 1500 2000 2500 3000 3500"
bxList="1 2 3 4"
alias python3=/usr/local/bin/python3.9
python3 --version
#250 300 500 750 1000 1500 2000 2500  3000 3500 
for signalTracks in 3500; do
    for bx in 1 2 3 4; do
        echo "############# working with nTracks: "${signalTracks}" and BX: "${bx}" #########"
        python3 getSigBkgTextFile.py ${bx} ${signalTracks}
    done
done
