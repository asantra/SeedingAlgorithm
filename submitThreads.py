#!/usr/bin/env python
### This is the script to run on multithreads
import os, sys
import subprocess
from subprocess import call
import Queue
import threading
import multiprocessing
import argparse
from random import seed, randint

def main():
    #parser = argparse.ArgumentParser(description='Code to find seed tracks')
    #parser.add_argument('-l', action="store", dest="inFile", default="list_root_hics_165gev_w0_5000nm_provisional_10xi3_6cee466a_trackInfoClean.txt")
    #parser.add_argument('-s', action="store_true", dest="needSignal", default=False)
    #parser.add_argument('-f', action="store", dest="needFit", type=int, default=1)
    #parser.add_argument('-e', action="store", dest="eCut", type=float, default=0.0)
    #parser.add_argument('-p', action="store", dest="Particle", type=str, default="Positron")
    #args = parser.parse_args()
    
    ## run!
    q = Queue.Queue()
   
    #command = 'bash runSeeding3Or4Hits.sh'
    command = 'python3 findSeed.py -l EBeamOnlyNewSamples_DividedByBX15_trackInfoClean.txt -f 1'
    print(command)
    q.put(command)

    def worker():
        while True:
            item = q.get()
            # execute a task: call a shell program and wait until it completes
            command  = str(item)
            p        = subprocess.Popen(command, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
            out, err = p.communicate()
            path     = "/Users/arkasantra/arka/Sasha_Work/seedingAlgorithm/multiCoreJobOutput"
            
            if not os.path.exists(path):
                os.makedirs(path)
                
            f = open(path+"/log_output.out", 'w')
            f.write(out)
            f.close()
            f = open(path+"/log_error.err", 'w')
            f.write(err)
            f.close()
            q.task_done()

    cpus=multiprocessing.cpu_count() #detect number of cores
    print("Creating %d threads" % cpus)
    for i in range(cpus):
        t = threading.Thread(target=worker)
        t.daemon = True
        t.start()

    q.join() # block until all tasks are done

if __name__=="__main__":
    main()
