import sys
import random
import numpy as np

############################## Defining Functions ###############################


def read_fasta():
    sequences = []
    current_seq = ''
    for line in sys.stdin:
        line = line.rstrip()
        if line.startswith(">"):
            if current_seq:
                sequences.append(current_seq)
                current_seq = ''
        else:
            current_seq += line
    if current_seq:
        sequences.append(current_seq)
    return sequences


def insert_sequence_randomly(sequences, site_sequence):
    inserted_sequences = []
    for seq in sequences:
        pos = random.randint(0, len(seq)-len(site_sequence))
        new_seq = seq[:pos] + site_sequence + seq[pos:]
        inserted_sequences.append(new_seq)
    return inserted_sequences


def calculate_q_matrix(motifs, pseudocount, alphabet):
    n_seq = len(motifs)
    n_sym = len(alphabet)
    width = len(motifs[0])
    count_matrix = [[0]*width for _ in range(n_sym)]
    for i, sym in enumerate(alphabet):
        for pos in range(width):
            count_matrix[i][pos] = sum(1 for m in motifs if m[pos] == sym)
    q_matrix = [[(count_matrix[i][pos]+pseudocount)/(n_seq+n_sym*pseudocount)
                 for pos in range(width)] for i in range(n_sym)]
    return count_matrix, q_matrix


def calculate_background_matrix(sequences, alphabet):
    total = sum(len(seq) for seq in sequences)
    counts = [sum(seq.count(sym) for seq in sequences) for sym in alphabet]
    return [c/total for c in counts]


def calculate_Ax(motif, q_matrix, background, alphabet):
    numerator = 1
    denominator = 1
    for i, char in enumerate(motif):
        idx = alphabet.index(char)
        numerator *= q_matrix[idx][i]
        denominator *= background[idx]
    return numerator/denominator


def Ax_to_probability(Ax):
    total = sum(Ax)
    return [x/total for x in Ax]


def score_motifs(motifs, alphabet):
    width = len(motifs[0])
    n_seq = len(motifs)
    n_sym = len(alphabet)
    counts = [[0]*width for _ in range(n_sym)]
    for pos in range(width):
        for i, sym in enumerate(alphabet):
            counts[i][pos] = sum(1 for m in motifs if m[pos] == sym)
    score = 0
    for pos in range(width):
        score += n_seq - max(counts[i][pos] for i in range(n_sym))
    return score


def consensus_sequence(motifs, alphabet):
    width = len(motifs[0])
    consensus = ''
    for pos in range(width):
        column = [m[pos] for m in motifs]
        consensus += max(alphabet, key=lambda x: column.count(x))
    return consensus

############################ Main Sampler Program ##################################


# Parsing command-line arguments
if len(sys.argv) < 5:
    print("Usage: python gibbs_sampler_v1.py a b s rerun < fastA_file")
    sys.exit(1)

a = int(sys.argv[1])
b = int(sys.argv[2])
s = sys.argv[3]
rerun = int(sys.argv[4])

# Read sequences from FASTA
sequences = read_fasta()

# Insert the sequence of width b randomly into each sequence
sequences_with_site = insert_sequence_randomly(sequences, s)

# Initialize random motifs from each sequence
alphabet = ['A', 'C', 'G', 'T']
motif_length = a
current_motifs = []
concencus_reruns = []

for u in range(rerun):
    for seq in sequences_with_site:
        start = random.randint(0, len(seq)-motif_length)
        current_motifs.append(seq[start:start+motif_length])

    best_motifs = current_motifs[:]
    best_score = score_motifs(current_motifs, alphabet)

    # Gibbs sampling loop
    for iteration in range(b):

        for i in range(len(sequences_with_site)):
            # Remove i-th sequence
            seq = sequences_with_site[i]
            other_motifs = current_motifs[:i] + current_motifs[i+1:]

            # Build q_matrix from other motifs
            _, q_matrix = calculate_q_matrix(other_motifs, 0.25, alphabet)
            b_matrix = calculate_background_matrix(
                sequences_with_site, alphabet)

            # Compute Ax for all positions in the removed sequence
            candidate_motifs = [seq[j:j+motif_length]
                                for j in range(len(seq)-motif_length+1)]
            Ax_values = [calculate_Ax(m, q_matrix, b_matrix, alphabet)

                         for m in candidate_motifs]
            Ax_prob = Ax_to_probability(Ax_values)

            # Sample new motif
            new_motif = str(np.random.choice(candidate_motifs, p=Ax_prob))
            current_motifs[i] = new_motif

        current_score = score_motifs(current_motifs, alphabet)
        if current_score < best_score:
            best_score = current_score
            best_motifs = current_motifs[:]

        print(f"\nIteration {iteration+1}:"
              f"Sampled motifs: {current_motifs}")

    # Output consensus
    final_consensus = consensus_sequence(best_motifs, alphabet)
    # print("\nFinal consensus sequence:", final_consensus)
    concencus_reruns.append(final_consensus)
consensus_from_reruns = consensus_sequence(concencus_reruns, alphabet)
print("\nConcences after", u, " reruns :", consensus_from_reruns)
print("\n", concencus_reruns)
