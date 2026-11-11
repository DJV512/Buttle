from .. import constants as btc
lazy from collections import Counter
from dataclasses import dataclass


@dataclass
class Sequence:
    name: str
    sequence: str


class Protein:

    def __init__(self, name:str, sequence:str, allow_x=False):
        self.name = name
        self.sequence = sequence.upper()
        self._validate_protein(allow_x = allow_x)


    def __str__(self):
        return self.name


    def __len__(self):
        return len(self.sequence)


    def __eq__(self, other):
        if isinstance(other, Protein):
            return self.name == other.name and self.sequence == other.sequence 


    def _validate_protein(self, allow_x):
        '''
        Only the normal 20 amino acids are allowed
        
        If allow_x is True, 'X' is also allowed.
        
        'X' is used as an alignment anchor in Multiple Sequence Alignment.
        '''
        if not self.sequence:
            raise ValueError("Empty sequence provided.")
        allowed_aas = set(btc.AMINO_ACIDS)
        if allow_x:
            allowed_aas.add("X")
        residues = set(self.sequence)
        if not residues.issubset(allowed_aas):
            raise ValueError(f"{self.sequence} is not a valid Protein sequence.")


    def charge(self) -> int:
        '''
        Calculates and returns its overall charge (at pH 7).
        '''

        charge = 0
        for aa in self.sequence:
            if aa in ["K", "R"]:
                charge += 1
            elif aa in ["D", "E"]:
                charge -= 1
        return charge


    def mw(self) -> float:
        '''
        Given a protein, calculates and returns its molecular weight.
        '''
        self._validate_protein(allow_x=False)
        weight = 18
        for aa in self.sequence:
            weight += btc.MW[aa]
        return weight


    def ec(self) -> int:
        '''
        Given a protein sequence, calculates and returns the extinction coefficient, in cm-1/M.
        '''

        ec = 0
        for aa in self.sequence:
            if aa == "W":
                ec += 5690
            elif aa == "Y":
                ec += 1280
            elif aa == "C":
                ec += 120
        return ec


    def glycosylation(self) -> list:
        '''
        Given a protein sequence, finds all potential N-linked glycosylations sites, denoted by the motif N-X-S/T, where X is any amino acid except proline.
        '''

        sites = []
        length = len(self.sequence)
        for i, aa in enumerate(self.sequence):
            if aa == "N" and i < length - 2:
                aa2 = self.sequence[i+1]
                if aa2 != "P":
                    aa3 = self.sequence[i+2]
                    if aa3 in ["T", "S"]:
                        sites.append(i)
        return sites



class DNA:

    def __init__(self, name: str, sequence: str, allow_n=False):
        self.name = name
        self.sequence = sequence.upper()
        self._validate_dna(allow_n=allow_n)


    def __str__(self):
        return self.name

    
    def __len__(self):
        return len(self.sequence)


    def __eq__(self, other):
        if isinstance(other, DNA):
            return self.name == other.name and self.sequence == other.sequence 


    def _validate_dna(self, allow_n):
        '''
        Only the normal 4 bases are allowed
        
        If allow_n is True, 'N' is also allowed.
        '''
        if not self.sequence:
            raise ValueError("Empty sequence provided.")
        allowed_bases = set("ACGT")
        if allow_n:
            allowed_bases.add("N")
        bases = set(self.sequence)
        if not bases.issubset(allowed_bases):
            raise ValueError(f"{self.sequence} is not a valid DNA sequence.")


    def transcribe(self) -> str:
        '''
        Returns the RNA transcription.
        '''

        rna = ""
        for base in self.sequence:
            if base == "A":
                rna += "A"
            elif base == "C":
                rna += "C"
            elif base == "G":
                rna += "G"
            elif base == "T":
                rna += "U"
            else:
                rna += "N"
        return rna


    def translate(self) -> tuple[str, str, str]:
        '''
        Returns the protein translations in all three forward frames.
        '''

        # Initialize the three translations with empty strings.
        proteinf1 = ""
        proteinf2 = ""
        proteinf3 = ""

        # Frame 1
        codons = [self.sequence[i:i+3] for i in range(0, len(self.sequence), 3)]
        for codon in codons:
            if len(codon) == 3:
                try:
                    proteinf1 += btc.CODON_TABLE[codon]
                except KeyError:
                    proteinf1 += "X"

        # Frame 2
        codons = [self.sequence[i:i+3] for i in range(1, len(self.sequence), 3)]
        for codon in codons:
            if len(codon) == 3:
                try:
                    proteinf2 += btc.CODON_TABLE[codon]
                except KeyError:
                    proteinf2 += "X"

        # Frame 3
        codons = [self.sequence[i:i+3] for i in range(2, len(self.sequence), 3)]
        for codon in codons:
            if len(codon) == 3:
                try:
                    proteinf3 += btc.CODON_TABLE[codon]
                except KeyError:
                    proteinf3 += "X"

        return proteinf1, proteinf2, proteinf3


    def reverse_comp(self) -> str:
        '''
        Returns the reverse complement DNA sequence.
        '''

        rev_comp_dna = ""
        rev_dna = self.sequence[::-1]
        for base in rev_dna:
            if base == "A":
                rev_comp_dna += "T"
            elif base == "C":
                rev_comp_dna += "G"
            elif base == "G":
                rev_comp_dna += "C"
            elif base == "T":
                rev_comp_dna += "A"
            else:
                rev_comp_dna += "N"

        return rev_comp_dna


    def gc(self) -> tuple[int, Counter, float, float]:
        '''
        Returns information about the base contents.
        '''

        length = len(self)
        contents = Counter(self.sequence)
        at_content = 100 * (contents["A"] + contents["T"]) / length
        gc_content = 100 * (contents["G"] + contents["C"]) / length
        return length, contents, at_content, gc_content
        