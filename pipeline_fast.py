#! /usr/bin/env python3


# Protein Block Design, Master's degree (2nd year) - Structural Bioinformatics project
# Author : Damien GARCIA, M2BB
#
# This script is separated in 3 sections :
#   -   Initilization process :
#       The PDB file on command line argument is treated in order to generate the necessary
#       files (Amino acids sequence, PB sequence, forsa initial zScore).
#       Be aware that this method only works if you aim to generate potential homologues
#       for the chain A of the structure. If you wan't to work on any other chain, modify
#       the PDB file accordingly.
#   -   Generation cycle :
#       Sequences are generated following the length of the reference structure. The
#       generation follows a uniform law for the amino acid selection at each position to
#       avoid any a priori. Each sequence is scored using FORSA algorithm and is selected
#       only if the threshold (half the reference sequence zScore) is exceeded. This process
#       is repeated until a batch of 100 sequence is obtained.
#   -   Mutation cycle(s) :
#       Each sequence of the batch is cloned a 1000 times and mutated randomly with a
#       mutation rate set to 20%. All newly mutated sequence is scored and the top100 zScoring
#       sequences are selected for another cycle.
#       Cycles are repeated until either the number of repetition reach 5, or the reference
#       zScore is reached or exceeded by at least one sequence of the batch.
#
# Results are stored in csv files for the zScore and the corresponding sequences. A histogram
# graph of the zScore repartition of each batch of sequence is also provided.


# ------------------------------ Library imports ----------------------------- #


import os                               # File handling
import sys                              # Command line arguments
import time                             # Script execution elapsed time

import subprocess                       # Forsa and dssp algorithms execution
import random                           # Sequence generation and mutation rate
import pandas as pd                     # Dataframes
from matplotlib import pyplot as plt    # Plot generation

import uuid                             # Unique file name for parallelization process
from joblib import Parallel, delayed    # Parallelization process
from tqdm import tqdm                   # Progress Bar
from tqdm_joblib import tqdm_joblib     # Progress Bar for parallel execution

#import forsa                            #SWIG module from c
import pipeline                         #dont have to redefine functions



# ---------------------------------------------------------------------------- #
#                                     Main                                     #
# ---------------------------------------------------------------------------- #

if __name__ == "__main__":
    # Start timer
    start = time.time()

    # Print header
    pipeline.print_header()

    # Checking command line arguments
    # Structure of interest should be one chain only or the first chain of the pdb file.
    basename = ""
    if len(sys.argv) != 2 or not sys.argv[1].endswith(".pdb"):
        pipeline.print_help()
        exit(1)
    else:
        basename = os.path.basename(sys.argv[1])[:-4]

    # Dictionary structures containing all necessary paths
    path = {
        "forsa": "forsa_parameter/forsa_global",
        "dssp": "dssp/dssp-2.0.4-linux-i386",
        "dssp_sep": "dssp/dssp_separator.pl",
        "dssp2pb": "dssp/dssp_to_pb_tor_rmsda.pl",
        "pdb2aa": "dssp/pdb_to_aa.pl",

        "results": f"results_{basename}",
        "data": f"results_{basename}/data",
        "temp": f"results_{basename}/temp",
    }

    # Checking for files and directories existence
    if not os.path.exists(path['forsa']):
        print("[Err1] Missing required files to execute generation.")
        exit(1)
    if not os.path.exists(f"results_{basename}"):
        os.makedirs(path['results'])
        os.makedirs(path['data'])
        os.makedirs(path['temp'])
    elif not os.path.exists(path['data']):
        os.makedirs(path['data'])
    if not os.path.exists(path['temp']):
        os.makedirs(path['temp'])

    # Data types storing generated sequences and corresponding zScores
    initRun = pd.DataFrame({"zScore": pd.Series(dtype='float'),
                            "sequence": pd.Series(dtype='str')})
    newRun  = pd.DataFrame({"oriSeq": pd.Series(dtype='int'),
                            "zScore": pd.Series(dtype='float'),
                            "sequence": pd.Series(dtype='str')})

    # Generation constants
    pool = 100
    length, objective = pipeline.initialization(path, basename)
    threshold = objective/2

    # Generation and mutation cycles
    runIter, bestScore = 0, 0
    while bestScore < objective and runIter < 5:
        if runIter == 0:
            initRun = pipeline.init_process(path, basename, initRun, pool, threshold, length, bestScore)
            bestScore = max(initRun['zScore'].to_list())
        elif runIter == 1:
            sequences = pipeline.initRun['sequence'].to_list()
            newRun = pipeline.following_processes(path, basename, newRun, runIter, sequences, bestScore)
            bestScore = max(initRun['zScore'].to_list())
            print(newRun)
        else:
            sequences = newRun['sequence'].to_list()
            newRun = pipeline.following_processes(path, basename, newRun, runIter, sequences, bestScore)
            bestScore = max(newRun['zScore'].to_list())
            print(newRun)
        runIter += 1
    
    if len(os.listdir(path['temp'])) == 0:
        os.rmdir(path['temp'])
    
    # End timer
    end = time.time()
    elapsed = end - start
    
    print_results(runIter, elapsed, bestScore)
