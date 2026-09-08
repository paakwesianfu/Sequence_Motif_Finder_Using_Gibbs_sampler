# DNA Motif Finder using Gibbs Sampling
A Gibbs sampling loop in written with python that identifies conserved motifs (patterns) from a set of DNA sequences.

Motif finding is an important problem in bioinformatics, often used to detect regulatory elements such as transcription factor binding sites.


# How It Works
The algorithim follows the ff. major steps 
1.Input parsing: 
    .Reads DNA sequences from FASTA format using standard input

2.Motif Insertion (for testing): 
    .Randomly inserts a known motif into the sequences (for testing)

3.Initializing: 
    .Selects initial motifs from each sequence at random

4.Gibbs Sampling Loop:
    .Removes one sequence
    .Builds a probability matrix (q-matrix) from the remaining motifs (using pseudo counts to avoid zero probabilities)
    .Computes overall nucleotide frequencies across all sequences creating a background matrix
    .Slides a window across the removed sequence and computes the likelihood ratio comparing the q-matrix and the background matrix and calculates a weight for each sequence
    .The weights are normalized into a probability distribution and a new motif is randomly selected based on this distribution

4.Tracks the best motif using a mismatch based scoring function

5.Outputs a consensus sequence from the best motifs 
NB: The script contains a block of code to which reruns the loop for as many times as specified by the user

# Usage 
Run the program using the command below:
python gibbs_sampler_v1.py a b s rerun < input.fasta

Arguments
Parameter       Description

a                : motif_length

b               : Number of Gibbs sampling iterations

s               : Motif sequence (for testing)

rerun           : Number of independent reruns (should have a value of at least 1)

input.fasta     : The fasta file to be used for the search 

                (note to use a correctly formatted fasta file with name starting with ">")
