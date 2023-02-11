#! /usr/bin/env python3


# Protein Block Design, Master's degree (2nd year) - Structural Bioinformatics project
# Author : Damien GARCIA, M2BB


# To do to dooo :
# Dans le rapport, ajouter les 5 meilleures séquences produites.
# Modélisation moléculaire d'au moins un séquence.
# Utiliser un générateur de nom aléatoire pour paralléliser sur des fichiers différents (forsa)

# Références à l'article de Forsa
# Article : A substitution matrix for strucral alphabet based on structural alignment of homolgous proteins ans tis applications.



import os                               # File handling
import subprocess                       # Forsa algorithm exection
import random                           # Sequence generation and mutation rate
import pandas as pd                     # Dataframes
from matplotlib import pyplot as plt    # Plot generation
from joblib import Parallel, delayed    # Parallelization process
import uuid                             # Unique file name for parallelization process


# ------------------------------- Miscellaneous ------------------------------ #


# Prints formated header to console output 
def header():
    print("     ____   ____    ____              _               ")
    print("    |  _ \ | __ )  |  _ \   ___  ___ (_)  __ _  _ __  ")
    print("    | |_) ||  _ \  | | | | / _ \/ __|| | / _` || '_ \ ")
    print("    |  __/ | |_) | | |_| ||  __/\__ \| || (_| || | | |")
    print("    |_|    |____/  |____/  \___||___/|_| \__, ||_| |_|")
    print("                                         |___/        \n")


# Code by 'Greenstick':
# https://stackoverflow.com/questions/3173320/text-progress-bar-in-terminal-with-block-characters
# Init =     printProgressBar(0, total, prefix='Progress:', suffix='Complete', length=50)
# In loop =  printProgressBar(i+1, pool, prefix='Progress:', suffix='Complete', length=50)
def printProgressBar(iteration, total, prefix='', suffix='', decimals=1,
                     length=100, fill='█', printEnd="\r"):
    """
    Call in a loop to create terminal progress bar
    @params:
        iteration   - Required  : current iteration (Int)
        total       - Required  : total iterations (Int)
        prefix      - Optional  : prefix string (Str)
        suffix      - Optional  : suffix string (Str)
        decimals    - Optional  : positive number of decimals in percent complete (Int)
        length      - Optional  : character length of bar (Int)
        fill        - Optional  : bar fill character (Str)
        printEnd    - Optional  : end character (e.g. "\r", "\r\n") (Str)
    """
    percent = ("{0:." + str(decimals) + "f}").format(100 * (iteration / float(total)))
    filledLength = int(length * iteration // total)
    bar = fill * filledLength + ' ' * (length - filledLength)
    print(f'\r{prefix} |{bar}| {percent}% {suffix}', end = printEnd)
    # Print New Line on Complete
    if iteration == total: 
        print()


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


def hist_generation(runIter, scores):
    print(f"Histogram generated to 'gen{runIter}_hist.png'")
    if isinstance(scores, list):
        x = scores
    else:
        x = scores['zScore']
    plt.hist(x, density=True, bins=100)
    plt.xlabel('zScore')
    plt.ylabel('Percentage')

    plt.title(f'Generation{runIter} zScore distribution')
    plt.savefig(f'gen{runIter}.png', dpi='figure', format='png')
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


# ---------------------------- Initial generation ---------------------------- #


def initial_generation(path, runIter, pool, initRun):
    allScores = []
    # Generation and scoring of sequences until a valid pool of 100 sequences is obtained
    print("Generation of pool of sequences:")
    printProgressBar(0, pool, prefix='Progress:', suffix='Complete', length=50)
    for i in range(pool):
        seq = ""
        zScore = 0.0
        while zScore < 3.0:      # Initial sequence validation threshold is set to 3 
            # Sequence generation
            seq = sequence_generator(length)
            write_sequence(seq, path['tempSeq'])
            
            # Forsa algorithm execution
            subprocess.run(
                f"./{path['forsa']} {path['tempSeq']} {path['pbseq']} -5 > {path['tempRes']}",
                shell=True
            )

            # Data Handling
            zScore = extract_zscore(path['tempRes'])
            allScores.append(zScore)
        
        initRun.loc[len(initRun.index)] = [zScore, seq]
        printProgressBar(i+1, pool, prefix='Progress:', suffix='Complete', length=50)

    print("\nPool selected from generation:\n")
    print(initRun.sort_values(by=["zScore"], ascending=False))
    
    print(f"\nTotal number of sequence generated = {len(allScores)}\n")

    hist_generation(runIter, allScores)
    write_csv("gen0_data.csv", initRun)


# ----------------------------- Second generation ---------------------------- #


def following_generation(path, runIter, sequences, newRun):
    batch = []
    nbCopies = 100
    # Mutating sequences
    print("\nMutating sequences:")
    printProgressBar(0, len(sequences), prefix='Progress:', suffix='Complete', length=50)
    for i in range(len(sequences)):
        # Generation of 1000 mutated copies for each sequence
        for j in range(nbCopies):
            batch.append(mutate(sequences[i]))
        printProgressBar(i+1, len(sequences), prefix='Progress:', suffix='Complete', length=50)

    # Forsa algorithm execution
    print("\nForsa zScore calculation:")
    printProgressBar(0, len(batch), prefix='Progress:', suffix='Complete', length=50)
    for i in range(len(batch)):
        write_sequence(batch[i], path['tempSeq'])
        subprocess.run(
            f"./{path['forsa']} {path['tempSeq']} {path['pbseq']} -5 > {path['tempRes']}",
            shell=True
        )
        newRun.loc[len(newRun.index)] = [i//nbCopies, extract_zscore(path['tempRes']), batch[i]]
        printProgressBar(i+1, len(batch), prefix='Progress:', suffix='Complete', length=50)
    
    hist_generation(runIter, newRun)
    newRun = newRun.nlargest(100, 'zScore')
    write_csv(f"gen{runIter}_data.csv", newRun)
    



# ----------------------------------- Main ----------------------------------- #


if __name__ == "__main__":
    # Dictionary structure containing all necessary paths
    path = {
        "forsa": "forsa/forsa_global",
        "pbseq": "2xiw.pb",
        "tempSeq": "temp_seq.txt",
        "tempRes": "temp_res.forsa",
        "temp": "temp",
        "result": "result"
    }

    # Checking for files and directories existence
    if not os.path.exists(path['forsa'] or not os.path.exists(path['pbseq'])):
        print("[Err1] Missing required files to execute generation.")
        exit(1)
    if not os.path.isdir(path['temp']):
        os.makedirs('temp')
    if not os.path.isdir(path['result']):
        os.makedirs('path')

    # Data types storing generated sequences and corresponding zScores
    initRun = pd.DataFrame({"zScore": pd.Series(dtype='float'),
                            "sequence": pd.Series(dtype='str')})
    newRun  = pd.DataFrame({"oriSeq": pd.Series(dtype='int'),
                            "zScore": pd.Series(dtype='float'),
                            "sequence": pd.Series(dtype='str')})

    # Generation constants
    length = 64         # Len of Sac7d protein used as reference for sequence generation.
    pool = 5            # Set to 100 but can be lowered for testing purposes.
    objective = 6.44    # zScore of reference for Sac7d 2XIW protein.

    # Print header
    header()

    # Generation and mutation cycles
    runIter, bestScore = 0, 0
    while runIter < 3 or bestScore < objective :
        if runIter == 0:
            initial_generation(path, runIter, pool, initRun)
        elif runIter == 1:
            sequences = initRun['sequence'].to_list()
            following_generation(path, runIter, sequences, newRun)
            bestScore = max(initRun['zScore'].to_list())
        else:
            sequences = newRun['sequence'].to_list()
            following_generation(path, runIter, sequences, newRun)
            bestScore = max(newRun['zScore'].to_list())
        print(bestScore)
        runIter += 1


    # hist_diff(pool, initRun, newRun)

    # Removing temporary files
    os.remove(path['tempRes'])
    os.remove(path['tempSeq'])