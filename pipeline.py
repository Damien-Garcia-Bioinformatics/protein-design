#! /usr/bin/env python3


# Protein Block Design, Master's degree (2nd year) - Structural Bioinformatics project
# Author : Damien GARCIA, M2BB


# To do to dooo :
# Dans le rapport, ajouter les 5 meilleures séquences produites.
# Modélisation moléculaire d'au moins un séquence.


# ------------------------------ Library imports ----------------------------- #


import os                               # File handling
import sys                              # Command line arguments

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


# Prints formated header to console output 
def header():
    print("     ____   ____    ____              _               ")
    print("    |  _ \ | __ )  |  _ \   ___  ___ (_)  __ _  _ __  ")
    print("    | |_) ||  _ \  | | | | / _ \/ __|| | / _` || '_ \ ")
    print("    |  __/ | |_) | | |_| ||  __/\__ \| || (_| || | | |")
    print("    |_|    |____/  |____/  \___||___/|_| \__, ||_| |_|")
    print("                                         |___/        \n")


def help():
    print("Help function")


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
    print(f"Writing selected pool of sequences to '{path}'")
    batch.to_csv(path, sep=';', encoding='utf-8', index=False)


# --------------------------- Histogram generation --------------------------- #


def hist_init():
    return 0


def hist_generation(runIter, scores):
    print(f"Histogram generated to 'gen{runIter}_hist.png'")
    if isinstance(scores, list):
        x = scores
    else:
        x = scores['zScore']
    plt.hist(x, density=True, bins=100)
    plt.xlabel('zScore')
    plt.ylabel('Percentage')
    plt.axvline(3.0, color='r', linestyle='--', label='...')
    plt.title(f'Generation{runIter} zScore distribution')
    plt.savefig(f'results_{basename}/gen{runIter}.png', dpi='figure', format='png')
    plt.clf()


### WIP ###
def hist_diff(pool, scores1, scores2):
    nbCopies = 100
    s1 = scores1['zScore'].to_list()
    s2 = scores2['zScore'].to_list()

    result = []
    for i in range(len(s1)):
        count = 0
        for j in range(100+i*100):
            if s2[j] > s1[i]:
                count += 1
        result.append((count/100)*100)
        count
    # plt.bar(result)
    # plt.xlabel("Parent sequence ID")
    # plt.ylabel("Percentage")
    # plt.title("Percentage of sequence with higher zScore than their respective parent")
    # plt.show()
    # plt.savefig
    # plt.clf()
### WIP ###


# ------------------------------ Initialization ------------------------------ #


def initialization(path_progs, basename):
    # Running PDB to AA program
    subprocess.run(f"perl {path_progs['pdb2aa']} {basename}.pdb ./results_{basename}/data/",
                   shell=True)
    length = 0
    with open(f"results_{basename}/data/{basename}.pdb.aa", "r+") as file:
        seq = file.read()
        file.seek(0)
        file.write(f">{basename}\n{seq}")
        length = len(seq)
    
    # Running DSSP program
    subprocess.run(f"./{path_progs['dssp']} -i {basename}.pdb -o ./results_{basename}/data/{basename}.dssp",
                   shell=True)
    
    # Running DSSP_SEPARATOR program
    subprocess.run(f"perl {path_progs['dssp_sep']} results_{basename}/data/{basename}.dssp ./results_{basename}/data/",
                   shell=True)
    
    # Running DSSP to PB program
    subprocess.run(f"perl {path_progs['dssp2pb']} results_{basename}/data/{basename}_A.dssp ./results_{basename}/data/",
                   shell=True)
    with open(f"results_{basename}/data/{basename}_A.dssp.pb", "r+") as file:
        seq = file.read()
        file.seek(0)
        file.write(f">{basename}\n{seq}")

    # Running FORSA
    subprocess.run(f"./{path_progs['forsa']} results_{basename}/data/{basename}.pdb.aa results_{basename}/data/{basename}_A.dssp.pb -5 > results_{basename}/data/{basename}.forsa",
                   shell=True)
    objective = extract_zscore(f"results_{basename}/data/{basename}.forsa")

    print("Selected parameters are:")
    print(f"  - Batch size = 100")
    print(f"  - Seqeuences length = {length}")
    print(f"  - Selection threshold = {objective/2}")
    print(f"  - zScore objective = {objective}\n")

    return length, objective


# ---------------------------- Initial generation ---------------------------- #


# Parallelized process of generating sequence until initial validation of minimal score.
# Function called by init_process()
def init_generation(threshold, length):
    seq = ""
    zScore = 0.0
    id = uuid.uuid4()
    while zScore < threshold:      # Initial sequence validation threshold is set to 3
        # Sequence generation
        seq = sequence_generator(length)
        write_sequence(seq, f"results_{basename}/temp/{id}.aaseq")
        
        # Forsa algorithm execution
        subprocess.run(
            f"./{path_progs['forsa']} results_{basename}/temp/{id}.aaseq results_{basename}/data/{basename}_A.dssp.pb -5 > results_{basename}/temp/{id}.forsa",
            shell=True
        )
        zScore = extract_zscore(f"results_{basename}/temp/{id}.forsa")
    
    # Cleaning temporary files
    os.remove(f"results_{basename}/temp/{id}.aaseq")
    os.remove(f"results_{basename}/temp/{id}.forsa")

    return zScore, seq


# First process of sequences generation
def init_process(pool, initRun, threshold, length):
    # Parallelized generation and scoring of sequences until a valid pool of 100 sequences is obtained
    print("Generation of pool of sequences:")
    data = []
    with tqdm_joblib(desc="", total=pool) as progress_bar:
        data += Parallel(n_jobs=-1)(delayed(init_generation)(threshold, length) for i in range(pool))
    
    for result in data:
        initRun.loc[len(initRun.index)] = [result[0], result[1]]
    
    print("\nTop10 scoring sequences from generation:")
    print(initRun.sort_values(by=["zScore"], ascending=False))
    
    # hist_generation(runIter, allScores)
    write_csv(f"results_{basename}/gen0_data.csv", initRun)


# ----------------------------- Second generation ---------------------------- #


# Parallelized process of scoring sequences in all but first mutation/scoring cycles.
# Function called by following_processes().
def following_generation(seq, seqN, nbCopies):
    id = uuid.uuid4()
    write_sequence(seq, f"results_{basename}/temp/{id}.aaseq")
    subprocess.run(
        f"./{path_progs['forsa']} results_{basename}/temp/{id}.aaseq results_{basename}/data/{basename}_A.dssp.pb -5 > results_{basename}/temp/{id}.forsa",
        shell=True
    )
    oriSeq = seqN//nbCopies
    zScore = extract_zscore(f"results_{basename}/temp/{id}.forsa")

    os.remove(f"results_{basename}/temp/{id}.aaseq")
    os.remove(f"results_{basename}/temp/{id}.forsa")
    
    return oriSeq, zScore, seq


# All but first mutation/scoring cycles (Cycles 1 and over).
def following_processes(runIter, sequences, newRun):
    batch = []
    nbCopies = 100

    # Mutating sequences
    print("\nMutating sequences:")
    for i in tqdm(range(len(sequences))):
        # Generation of 1000 mutated copies for each sequence
        batch.append(Parallel(n_jobs=-1)(delayed(mutate)(sequences[i]) for j in range(nbCopies)))
    batch = flatten_list(batch)

    # Forsa algorithm execution
    print("\nForsa zScore calculation:")
    data = []
    with tqdm_joblib(desc="", total=len(batch)) as progress_bar:
        data.append(Parallel(n_jobs=-1)(delayed(following_generation)(batch[i], i, nbCopies) for i in range(len(batch))))
    data = flatten_list(data)

    for result in data:
        newRun.loc[len(newRun.index)] = [result[0], result[1], result[2]]
    
    hist_generation(runIter, newRun)
    newRun = newRun.nlargest(100, 'zScore')
    # print(newRun)
    write_csv(f"results_{basename}/gen{runIter}_data.csv", newRun)

    return newRun
    


# ---------------------------------------------------------------------------- #
#                                     Main                                     #
# ---------------------------------------------------------------------------- #

if __name__ == "__main__":
    # Checking command line arguments
    # Structure of interest should be one chain only or the first chain of the pdb file.
    basename = ""
    if len(sys.argv) != 2 or not sys.argv[1].endswith(".pdb"):
        help()
        exit(1)
    else:
        basename = os.path.basename(sys.argv[1])[:-4]

    # Dictionary structures containing all necessary paths
    path_progs = {
        "forsa":    "forsa/forsa_global",
        "dssp":     "dssp/dssp-2.0.4-linux-i386",
        "dssp_sep": "dssp/dssp_separator.pl",
        "dssp2pb":  "dssp/dssp_to_pb_tor_rmsda.pl",
        "pdb2aa":   "dssp/pdb_to_aa.pl",
    }
    # dir_path = {
    #     "data": f"results_{basename}/data/",
    #     "temp": f"results_{basename}/temp/",
    # }

    # path_files = {
    #     "pdb":   f"{basename}.pdb",
    #     "aaseq": f"results_{basename}/data/{basename}.pdb.aa",
    #     "dssp":  f"results_{basename}/data/{basename}.dssp",
    #     "pbseq": f"results_{basename}/data/{basename}.pbseq",
    #     "forsa": f"results_{basename}/data/{basename}_A.forsa",
    #     "temp":  f"results_{basename}/temp",
    # }

    # Checking for files and directories existence
    if not os.path.exists(path_progs['forsa']):
        print("[Err1] Missing required files to execute generation.")
        exit(1)
    if not os.path.exists(f"results_{basename}"):
        os.makedirs(f"results_{basename}")
        os.makedirs(f"results_{basename}/data")
        os.makedirs(f"results_{basename}/temp")
    elif not os.path.exists(f"results_{basename}/data"):
        os.makedirs(f"results_{basename}/data")
    if not os.path.exists(f"results_{basename}/temp"):
        os.makedirs(f"results_{basename}/temp")

    # Data types storing generated sequences and corresponding zScores
    initRun = pd.DataFrame({"zScore": pd.Series(dtype='float'),
                            "sequence": pd.Series(dtype='str')})
    newRun  = pd.DataFrame({"oriSeq": pd.Series(dtype='int'),
                            "zScore": pd.Series(dtype='float'),
                            "sequence": pd.Series(dtype='str')})
    
    # Print header
    header()

    # Generation constants
    pool = 100
    length, objective = initialization(path_progs, basename)
    threshold = objective/2

    # Generation and mutation cycles
    runIter, bestScore = 0, 0
    while runIter < 3 or bestScore < objective :
        if runIter == 0:
            init_process(pool, initRun, threshold, length)
        elif runIter == 1:
            sequences = initRun['sequence'].to_list()
            newRun = following_processes(runIter, sequences, newRun)
            bestScore = max(initRun['zScore'].to_list())
            print(newRun)
        else:
            sequences = newRun['sequence'].to_list()
            newRun = following_processes(runIter, sequences, newRun)
            bestScore = max(newRun['zScore'].to_list())
            print(newRun)
        print(bestScore)
        runIter += 1