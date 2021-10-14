#! /bin/bash
signalTracksList="1 5 10 20 30 50 80 100 130 150 170 185 200 220"
for signalTracks in ${signalTracksList}; do
    echo "############# working with nTracks: "${signalTracks}" #########"
    python plotTightAndLooseDistPerBXWithFit.py -n ${signalTracks}
done
