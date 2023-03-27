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


# ---------------------------------------------------------------------------- #
#                                   Functions                                  #
# ---------------------------------------------------------------------------- #

# ------------------------------- Miscellaneous ------------------------------ #


def print_header():
    print('  ╔══════════════════════════════════════════════════════════════════════╗')
    print('  ║  ██████╗ ██████╗     ██████╗ ███████╗███████╗██╗ ██████╗ ███╗   ██╗  ║')
    print('  ║  ██╔══██╗██╔══██╗    ██╔══██╗██╔════╝██╔════╝██║██╔════╝ ████╗  ██║  ║')
    print('  ║  ██████╔╝██████╔╝    ██║  ██║█████╗  ███████╗██║██║  ███╗██╔██╗ ██║  ║')
    print('  ║  ██╔═══╝ ██╔══██╗    ██║  ██║██╔══╝  ╚════██║██║██║   ██║██║╚██╗██║  ║')
    print('  ║  ██║     ██████╔╝    ██████╔╝███████╗███████║██║╚██████╔╝██║ ╚████║  ║')
    print('  ║  ╚═╝     ╚═════╝     ╚═════╝ ╚══════╝╚══════╝╚═╝ ╚═════╝ ╚═╝  ╚═══╝  ║')
    print('  ╚══════════════════════════════════════════════════════════════════════╝\n')
                                                                  

def print_help():
    print('  ╔══════════════════════════════════════════════════════════════════════╗')
    print('  ║  Use : "./pipeline [pathToFile]"                                     ║')
    print('  ║        with [pathToFile] a PDB file containing one chain you wish    ║')
    print('  ║        to generate potential homologues of.                          ║')
    print('  ║                                                                      ║')
    print('  ║  For more informations, please visit:                                ║')
    print('  ║  https://github.com/Damien-Garcia-Bioinformatics/protein-design/     ║')
    print('  ╚══════════════════════════════════════════════════════════════════════╝\n')


def print_parameters(length, objective):
    print(f'  ╔══════════════════════════════════════════════════════════════════════╗')
    print(f'  ║  Selected parameters:                                                ║')
    print(f'  ║    - Batch size          = 100                                       ║')
    print(f'  ║    - sequences length    = {str(length):<42}║')
    print(f'  ║    - Selection threshold = {str(objective/2):<42}║')
    print(f'  ║    - Z-score objective   = {str(objective):<42}║')
    print(f'  ╚══════════════════════════════════════════════════════════════════════╝\n')


def print_results(runIter, elapsed, bestScore):
    print(f'\n  ╔══════════════════════════════════════════════════════════════════════╗')
    print(f'  ║  Script execution results:                                           ║')
    print(f'  ║    - Total number of cycles    = {str(runIter):<36}║')
    print(f'  ║    - Script execution time (s) = {str(round(elapsed, 2)):<36}║')
    print(f'  ║    - Highest scoring sequence  = {str(bestScore):<36}║')
    print(f'  ╚══════════════════════════════════════════════════════════════════════╝\n')

# Flatens list of lists to list
# Parallel processing appends batches for optimization purposes? --> Needs investigating (°_°)'
def flatten_list(list):
    return [item for sublist in list for item in sublist]


# ---------------------------- Sequence generation --------------------------- #

# Random amino acids sequence generator following uniform law.
def sequence_generator(length):
    aa = "ACDEFGHIKLMNPQRSTVWY"
    seq = ""
    for i in range(length):
        seq += aa[random.randint(0, 19)]
    return seq


def mutate(seq):
    rate = 0.20         # Mutation rate
    aa = "ACDEFGHIKLMNPQRSTVWY"
    newSeq = ""
    for i in range(len(seq)):
        if random.uniform(0, 1) < rate:
            newSeq += aa[random.randint(1, 19)]
        else:
            newSeq += seq[i]
    return newSeq


# ------------------------------- File Handling ------------------------------ #


# Write sequence in temporary file to feed forsa algorithm.
def write_sequence(seq, path):
    with open(path, "w+") as file:
        file.write(">tempSeq\n")
        file.write(seq)


# zScore extraction from forsa algorithm output.
def extract_zscore(path):
    with open(path) as file:
        for line in file:
            if line.startswith("raw score"):
                return float(line[line.rfind(':')+1:])


# Writes zScore and sequence from selected pool of sequence to csv file.
def write_csv(path, batch):
    print(f"\nWriting selected pool of sequences to '{path}'")
    batch.to_csv(path, sep=';', encoding='utf-8', index=False)


# --------------------------- Histogram generation --------------------------- #


def hist_generation(runIter, scores, bestScore):
    print(f"\nHistogram generated to 'gen{runIter}_hist.png'")
    if isinstance(scores, list):
        x = scores
    else:
        x = scores['zScore']
    plt.hist(x, density=True, bins=100)
    plt.xlabel('zScore')
    plt.ylabel('Percentage')
    if bestScore != 0:
        plt.axvline(bestScore, color='r', linestyle='--', label='Highest scoring sequence of previous cycle')
        plt.legend()
    plt.title(f'Generation{runIter} zScore distribution')
    plt.savefig(f'results_{basename}/gen{runIter}.png', dpi='figure', format='png')
    plt.clf()


# ------------------------------ Initialization ------------------------------ #


def initialization(path, basename):
    # Running PDB to AA program
    subprocess.run(f"perl {path['pdb2aa']} {basename}.pdb ./{path['data']}",
                   shell=True)
    length = 0
    with open(f"{path['data']}/{basename}.pdb.aa", "r+") as file:
        seq = file.read()
        file.seek(0)
        file.write(f">{basename}\n{seq}")
        length = len(seq)
    
    # Running DSSP program
    subprocess.run(f"./{path['dssp']} -i {basename}.pdb -o ./{path['data']}/{basename}.dssp",
                   shell=True)
    
    # Running DSSP_SEPARATOR program
    subprocess.run(f"perl {path['dssp_sep']} {path['data']}/{basename}.dssp ./{path['data']}",
                   shell=True)
    
    # Running DSSP to PB program
    subprocess.run(f"perl {path['dssp2pb']} {path['data']}/{basename}_A.dssp ./{path['data']}",
                   shell=True)
    with open(f"{path['data']}/{basename}_A.dssp.pb", "r+") as file:
        seq = file.read()
        file.seek(0)
        file.write(f">{basename}\n{seq}")

    # Running FORSA
    subprocess.run(f"./{path['forsa']} {path['data']}/{basename}.pdb.aa {path['data']}/{basename}_A.dssp.pb -5 > {path['data']}/{basename}.forsa",
                   shell=True)
    objective = extract_zscore(f"{path['data']}/{basename}.forsa")

    print_parameters(length, objective)

    return length, objective


# ---------------------------- Initial generation ---------------------------- #


# Parallelized process of generating sequence until initial validation of minimal score.
# Function called by init_process()
def init_generation(path, basename, threshold, length):
    seq = ""
    zScore = 0.0
    id = uuid.uuid4()
    while zScore < threshold:      # Initial sequence validation threshold is set to 3
        # Sequence generation
        seq = sequence_generator(length)
        write_sequence(seq, f"{path['temp']}/{id}.aaseq")
        
        # Forsa algorithm execution
        subprocess.run(
            f"./{path['forsa']} {path['temp']}/{id}.aaseq {path['data']}/{basename}_A.dssp.pb -5 > {path['temp']}/{id}.forsa",
            shell=True
        )
        zScore = extract_zscore(f"{path['temp']}/{id}.forsa")
    
    # Cleaning temporary files
    os.remove(f"{path['temp']}/{id}.aaseq")
    os.remove(f"{path['temp']}/{id}.forsa")

    return zScore, seq


# First process of sequences generation
def init_process(path, basename, initRun, pool, threshold, length, bestScore):
    # Parallelized generation and scoring of sequences until a valid pool of 100 sequences is obtained
    print("Generation of pool of sequences:")
    data = []
    with tqdm_joblib(desc="", total=pool) as progress_bar:
        data += Parallel(n_jobs=-1)(delayed(init_generation)(path, basename, threshold, length) for i in range(pool))
    
    for result in data:
        initRun.loc[len(initRun.index)] = [result[0], result[1]]
    
    print("\nInitial sequences generated:")
    print(initRun.sort_values(by=["zScore"], ascending=False))
    
    hist_generation(runIter, initRun, bestScore)
    write_csv(f"{path['results']}/gen0.csv", initRun)

    return initRun


# ----------------------------- Second generation ---------------------------- #


# Parallelized process of scoring sequences in all but first mutation/scoring cycles.
# Function called by following_processes().
def following_generation(path, basename, seq, seqN, nbCopies):
    id = uuid.uuid4()
    write_sequence(seq, f"{path['temp']}/{id}.aaseq")
    subprocess.run(
        f"./{path['forsa']} {path['temp']}/{id}.aaseq {path['data']}/{basename}_A.dssp.pb -5 > {path['temp']}/{id}.forsa",
        shell=True
    )
    oriSeq = seqN//nbCopies
    zScore = extract_zscore(f"{path['temp']}/{id}.forsa")

    os.remove(f"{path['temp']}/{id}.aaseq")
    os.remove(f"{path['temp']}/{id}.forsa")
    
    return oriSeq, zScore, seq


# All but first mutation/scoring cycles (Cycles 1 and over).
def following_processes(path, basename, newRun, runIter, sequences, bestScore):
    batch = []
    nbCopies = 100

    # Mutating sequences
    print("\nMutating sequences:")
    for i in tqdm(range(len(sequences))):
        # Generation of nbCopies mutated copies for each sequence
        batch.append(Parallel(n_jobs=-1)(delayed(mutate)(sequences[i]) for j in range(nbCopies)))
    batch = flatten_list(batch)

    # Forsa algorithm execution
    print("\nForsa zScore calculation:")
    data = []
    with tqdm_joblib(desc="", total=len(batch)) as progress_bar:
        data.append(Parallel(n_jobs=-1)(delayed(following_generation)(path, basename, batch[i], i, nbCopies) for i in range(len(batch))))
    data = flatten_list(data)

    for result in data:
        newRun.loc[len(newRun.index)] = [result[0], result[1], result[2]]
    
    hist_generation(runIter, newRun, bestScore)
    newRun = newRun.nlargest(100, 'zScore')
    write_csv(f"{path['results']}/gen{runIter}.csv", newRun)

    return newRun
    


# ---------------------------------------------------------------------------- #
#                                     Main                                     #
# ---------------------------------------------------------------------------- #

if __name__ == "__main__":
    # Start timer
    start = time.time()

    # Print header
    print_header()

    # Checking command line arguments
    # Structure of interest should be one chain only or the first chain of the pdb file.
    basename = ""
    if len(sys.argv) != 2 or not sys.argv[1].endswith(".pdb"):
        print_help()
        exit(1)
    else:
        basename = os.path.basename(sys.argv[1])[:-4]

    # Dictionary structures containing all necessary paths
    path = {
        "forsa": "forsa/forsa_global",
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
    length, objective = initialization(path, basename)
    threshold = objective/2

    # Generation and mutation cycles
    runIter, bestScore = 0, 0
    while bestScore < objective and runIter < 5:
        if runIter == 0:
            initRun = init_process(path, basename, initRun, pool, threshold, length, bestScore)
            bestScore = max(initRun['zScore'].to_list())
        elif runIter == 1:
            sequences = initRun['sequence'].to_list()
            newRun = following_processes(path, basename, newRun, runIter, sequences, bestScore)
            bestScore = max(initRun['zScore'].to_list())
            print(newRun)
        else:
            sequences = newRun['sequence'].to_list()
            newRun = following_processes(path, basename, newRun, runIter, sequences, bestScore)
            bestScore = max(newRun['zScore'].to_list())
            print(newRun)
        runIter += 1
    
    if len(os.listdir(path['temp'])) == 0:
        os.rmdir(path['temp'])
    
    # End timer
    end = time.time()
    elapsed = end - start
    
    print_results(runIter, elapsed, bestScore)
