from . import constants as btc
from .models import molecules as btmm
from . import validation as btv
from math import inf
from operator import itemgetter


def pairwise_align(
        molecule1:btmm.Protein | btmm.DNA, 
        molecule2:btmm.Protein | btmm.DNA
    ) -> tuple[list[str], list[str]]:
    """
    Perform a global pairwise alignment using affine gap penalties.

    Both inputs must be of the same type (Protein or DNA).
    Protein alignments use BLOSUM62 substitution scores, while DNA
    alignments use fixed match and mismatch scores. Gap penalties
    are determined by the molecule type.

    Args:
        molecule1: First Protein or DNA object to align.
        molecule2: Second Protein or DNA object to align.

    Returns:
        A tuple containing two lists of aligned residues or bases,
        including '-' characters representing gaps.

    Raises:
        TypeError: If the inputs are of different or unsupported types.
    """

    if isinstance(molecule1, btmm.Protein) and isinstance(molecule2, btmm.Protein):
        gap_open = btc.PROTEIN_GAP_OPEN
        gap_extend = btc.PROTEIN_GAP_EXTEND
        sigma = _sigma_protein

    elif isinstance(molecule1, btmm.DNA) and isinstance(molecule2, btmm.DNA):
        gap_open = btc.DNA_GAP_OPEN
        gap_extend = btc.DNA_GAP_EXTEND
        sigma = _sigma_dna

    else:
        raise TypeError("Both sequences must be of the same type: Protein or DNA.")

    btv.check_alignment_size(len(molecule1.sequence), len(molecule2.sequence), backend="python")

    #Add a space at the beginning to represent the 0th position
    seq1 = " " + molecule1.sequence
    seq2 = " " + molecule2.sequence
    length1 = len(seq1)
    length2 = len(seq2)

    #Initialize the affine gap scoring matrices
    score_diag = {
        (0, 0): (0, "start")
    }
    score_vert = {
        (0, 0): (-inf, "impossible") 
    }
    score_hori = {
        (0, 0): (-inf, "impossible")
    }

    for j in range(1, length2):
        score_diag[(0, j)] = (-inf, "impossible")
        score_vert[(0, j)] = (-inf, "impossible")
        score_hori[(0, j)] = (gap_open + (j-1) * gap_extend, "diag" if j == 1 else "hori")

    for i in range(1, length1):
        score_diag[(i, 0)] = (-inf, "impossible")
        score_vert[(i, 0)] = (gap_open + (i-1) * gap_extend, "diag" if i == 1 else "vert")
        score_hori[(i, 0)] = (-inf, "impossible")

    #Iterate over the entire length of sequence 1
    for i in range(1, length1):

        #Iterate over the entire length of sequence 2
        for j in range(1, length2):

            #Choose best score for each of the three scoring matrices and add them
            #Current pair is a base-base alignment
            score = sigma(seq1[i], seq2[j])
            score_diag[(i, j)] = max(
                (score_diag[(i - 1, j - 1)][0] + score, "diag"),
                (score_hori[(i - 1, j - 1)][0] + score, "hori"),
                (score_vert[(i - 1, j - 1)][0] + score, "vert"),
                key = itemgetter(0)
            )
            
            #Current pair is a gap in sequence 2
            score_vert[(i, j)] = max(
                (score_diag[(i - 1, j)][0] + gap_open, "diag"),
                (score_vert[(i - 1, j)][0] + gap_extend, "vert"),
                key = itemgetter(0)
            )

            #Current pair is a gap in sequence 1
            score_hori[(i, j)] = max(
                (score_diag[(i, j - 1)][0] + gap_open, "diag"),
                (score_hori[(i, j - 1)][0] + gap_extend, "hori"),
                key = itemgetter(0)
            )

    return _global_backtrack(seq1, seq2, score_diag, score_hori, score_vert, length1-1, length2-1)
        

 
def _global_backtrack(
        sequence1:str, 
        sequence2:str, 
        score_diag:dict[tuple[int, int], tuple[float, str]], 
        score_hori:dict[tuple[int, int], tuple[float, str]], 
        score_vert:dict[tuple[int, int], tuple[float, str]], 
        i:int, 
        j:int
    ) -> tuple[list[str], list[str]]:
    '''
    Private method for taking pairwise score matrices and backtracking through them to assemble the optimal alignment.
    '''

    align1 = []
    align2 = []

    _, predecessor = max(
        (score_diag[(i, j)][0], "diag"),
        (score_hori[(i, j)][0], "hori"),
        (score_vert[(i, j)][0], "vert"),
        key = itemgetter(0)
    )

    while i > 0 or j > 0:

        if predecessor == "diag":
            predecessor = score_diag[(i, j)][1]
            align1.append(sequence1[i])
            align2.append(sequence2[j])
            i -= 1
            j -= 1

        elif predecessor == "hori":
            predecessor = score_hori[(i, j)][1]
            align1.append("-")
            align2.append(sequence2[j])
            j -= 1

        else:
            predecessor = score_vert[(i, j)][1]
            align1.append(sequence1[i])
            align2.append("-")
            i -= 1

    align1.reverse()
    align2.reverse()
    
    return align1, align2


def _sigma_dna(str1:str, str2:str) -> int:
    '''
    Returns the sigma score between two DNA bases.
    '''
    if str1 == str2 and str1 != "N":
        return btc.DNA_MATCH
    else:
        return btc.DNA_MISMATCH

def _sigma_protein(str1:str, str2:str) -> int:
    '''
    Returns the BLOSUM62 comparison score between two amino acids.
    '''
    return btc.BUTTLE_COMPARE[str1][str2]