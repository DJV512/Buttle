def _align(
        seq1: str,
        seq2: str,
        molecule_type: str
    ) -> tuple[str, str]:
    """
    Perform a global pairwise alignment using affine gap penalties.

    Both inputs must be of the same type (Protein sequence or DNA sequence).
    Protein alignments use BLOSUM62 substitution scores, while DNA
    alignments use fixed match and mismatch scores. Gap penalties
    are determined by the molecule type.

    Args:
        seq1: First Protein or DNA sequence to align.
        seq2: Second Protein or DNA sequence to align.

    Returns:
        A tuple containing two lists of aligned residues or bases,
        including '-' characters representing gaps.

    Warning:
        This function does not validate that the input sequences match
        the specified molecule type. Callers must perform this validation
        before invoking the function.
    """
    ...

def _merge(
        left_names:list[str], 
        left_sequences: list[str],
        right_names:list[str],
        right_sequences: list[str],
        molecule_type: str,
        gap_open:float,
        gap_extend:float,
    ) -> tuple[list[str], list[str]]:
    '''
    Globally align two profiles, made up of already aligned sequences, using affine gap penalties.

    Inputs are two Profiles to merge (separated out into their constituent lists of names and sequences to pass to rust),
    as well as the gap open and gap extend penalties to use. This internal function does not perform input validation.

    Returns a tuple of a list of sequence names, and a list of sequence alignments.
    '''
    ...