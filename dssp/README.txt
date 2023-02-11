Download DSSP source file from :
https://github.com/PDB-REDO/dssp

Install DSSP

Running the DSSP program is a "pre-processing" step required for PB Assignment.

Step1: DSSP Assignment

To generate a DSSP file, you could use the stand-alone version of the DSSP program or directly download the file from their webserver.
eg-
./dssp-2.0.4-linux-i386 1AQC.pdb >1AQC.dssp 
--------------------------------------------------------------------------------------------------------------------------------------

Step2: Chain separation

The script- dssp_separator.pl separates the DSSP file into individual chains.
eg-
perl dssp_separator.pl ./1AQC.dssp ./
produces- 1AQC_A.dssp 1AQC_A.dssp 1AQC_C.dssp 1AQC_D.dssp
--------------------------------------------------------------------------------------------------------------------------------------

Step3: PB assignment

For each of the chain separated DSSP files, PB assignment can be performed using the script dssp_to_pb_tor_rmsda.pl
eg-
perl dssp_to_pb_tor_rmsda.pl ./1AQC_A.dssp ./
produces-1AQC_A.dssp.pb 1AQC_A.dssp.rmsda 1AQC_A.dssp.torsion
The number of 'X' in the 1AQC_A.dssp.pb correspond to (n+4) where n is the number of missing residues as detected from the DSSP file.
