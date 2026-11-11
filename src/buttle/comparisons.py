lazy from collections import defaultdict
lazy from itertools import combinations
lazy from . import constants as btc
lazy from .models import molecules as btm
lazy from . import parsers as btp
lazy from . import validation as btv
lazy from . import buttle_rust as btrust


def motifs(protein1:btm.Protein, protein2:btm.Protein) -> dict:
    '''
    Given two protein sequences, this method looks for similar stretches fo amino acids, ranging from length 4 to length 40.

    Matches count as +1, while conservative substitutions count for +0.5. Any two peptides where total_score / length > threshold is considered a match.

    Threshold is more stringent for shorter peptides, and more lenient for longer peptides. Defined in constant THRESHOLD.
    '''

    # Determine shortest sequence and set that as seq1
    length1 = len(protein1.sequence)
    length2 = len(protein2.sequence)

    if length1 > length2:
        seq1 = protein2.sequence
        seq2 = protein1.sequence
    else:
        seq1 = protein1.sequence
        seq2 = protein2.sequence

    # Initialize the dictionary that will hold matches
    similars = defaultdict(list)

    # Loop that looks at peptides of length 4 to 40 in the two sequences, compares them, and looks for regions of similarity.
    for size in range(4, min(40, length1, length2)+1):

        for min_size, max_size, cutoff in btc.THRESHOLD:
            if min_size <= size <= max_size:
                threshold = size * cutoff
                break
        else:
            raise ValueError(f"No threshold defined for motif size {size}.")

        for n in range(0, len(seq1) - size + 1):
            peptide1 = seq1[n:n+size]
            for x in range(0, len(seq2) - size + 1):
                count = 0
                peptide2 = seq2[x:x+size]
                for z in range(size):
                    if peptide1[z] == peptide2[z]:
                        count += 1
                        continue
                    if peptide2[z] in btc.CONSERVATIVES[peptide1[z]]:
                        count += 0.5
                if count >= threshold:
                    similars[(size, threshold)].append(((peptide1, n+1, n + size), (peptide2, x+1, x + size), count))

    return similars


def identity(
        file:str
    ) -> dict[tuple[str,str], float]:

    sequences = btp.parse_fasta(file)

    if len(sequences) < 2:
            raise ValueError("Identity requires at least two sequences.")
    
    molecule_type = type(sequences[0]).__name__.lower()
    identities = {}

    for seq1, seq2 in combinations(sequences, 2):
        btv.check_alignment_size(len(seq1.sequence), len(seq2.sequence), backend="rust")
        align1, align2 = btrust._align(seq1.sequence, seq2.sequence, molecule_type)
        length = len(align1)
        identical = sum(1 for aa1, aa2 in zip(align1, align2, strict=True) if aa1==aa2 and (molecule_type != "dna" or aa1 != "N"))
        percent_identity = 100*(identical/length)
        identities[(seq1.name, seq2.name)] = percent_identity

    return identities