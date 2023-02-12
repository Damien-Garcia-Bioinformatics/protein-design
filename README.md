# Protein Block Design
Master's degree project aiming to design protein sequences folds using the "Protein Block" approach.

## About the project



## Getting started

### Installation
1. Clone the repo:
```
git clone https://github.com/Damien-Garcia-Bioinformatics/protein-design.git
```
2. Install requirements:
```
pip install -r requirements.txt
```

#### Full requirements list
```
contourpy==1.0.7
cycler==0.11.0
fonttools==4.38.0
joblib==1.2.0
kiwisolver==1.4.4
matplotlib==3.6.3
numpy==1.24.2
packaging==23.0
pandas==1.5.3
Pillow==9.4.0
pyparsing==3.0.9
python-dateutil==2.8.2
pytz==2022.7.1
six==1.16.0
tqdm==4.64.1
tqdm-joblib==0.0.2
```

### Usage

### Roadmap
- [x] Pipeline working with multiple cycles of mutation + scoring.
- [x] Parallelize initial generation of sequences, mutation and scoring cycles.
- [ ] Generate bar graphs for each generation cycle.
- [ ] First section of generation to get automaticaly the .pbseq used for forsa calculations.
- [ ] Command line arguments of fasta and pdb formats.
- [ ] Improve Command Line Output to easily assess the pipeline state.

### License
Distributed under the GNU GENERAL PUBLIC LICENSE (Version 3, 29 June 2007).
See ```LICENSE``` for more information.

### Contact
- Author: Damien Garcia
- eMail: damien.garcia@etu.univ-nantes.fr
- Project link: https://github.com/Damien-Garcia-Bioinformatics/protein-design


## Acknowledgments
- **FORSA:**
    - Use of a structural alphabet to find compatible folds for amino acid sequences.
    
    *Swapnil Mahajan, Alexandre G. de Brevern, Yves-Henri Sanejouand, Narayanaswamy Srinivasan, Bernard Offmann*

- **Protein Blocks:**
    - A substitution matrix for structural alphabet based on structural alignment of homologous proteins and its applications
    
    *Manoj Tyagi, Venkataraman S. Gowri, Narayanaswamy Srinivasan, Alexandre G. de Brevern, Bernard Offmann*.
    
    - Extension of a local backbone description using a structural alphabet: A new approach to the sequence-structure relationship
    
    *Alexandre G. de Brevern, Hélène Valadié, Serge Hazout, Catherine Etchebest*

- **DSSP:**
    - Knowledge-based protein secondary structure assignment
    *Dmitrij Frishman, Patrick Argos*
    ```
        SPDX-License-Identifier: BSD-2-Clause

        Copyright (c) 2020 NKI/AVL, Netherlands Cancer Institute
        Redistribution and use in source and binary forms, with or without
        modification, are permitted provided that the following conditions are met:
        1. Redistributions of source code must retain the above copyright notice, this
        list of conditions and the following disclaimer
        2. Redistributions in binary form must reproduce the above copyright notice,
        this list of conditions and the following disclaimer in the documentation
        and/or other materials provided with the distribution.

        THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
        ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
        WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
        DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE FOR
        ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
        (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
        LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
        ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
        (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
        SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
    ```