import buttle.pairwise as btpa
import buttle.models.molecules as btmm
import buttle.parsers as btparse
import buttle.constants as btc
from buttle import buttle_rust as btrust
import pytest
from hypothesis import given, settings, strategies as st
from itertools import combinations
from pathlib import Path


def test_sigma_dna_match():
    base1 = "A"
    base2 = "A"
    score = btpa._sigma_dna(base1, base2)
    assert score == 8


def test_sigma_dna_mismatch():
    base1 = "A"
    base2 = "G"
    score = btpa._sigma_dna(base1, base2)
    assert score == -5


def test_pairwise_dna_case_insensitive():
    dna1 = btmm.DNA("","acagtgca")
    dna2 = btmm.DNA("","actgtata")
    align1, align2 = btpa.pairwise_align(dna1, dna2)
    assert (align1, align2) == (["A", "C", "A", "G", "T", "G", "C", "A"],["A", "C", "T", "G", "T", "A", "T", "A"])


def test_pairwise_dna_leading_gap_first():
    dna1 = btmm.DNA("","AGTGTC")
    dna2 = btmm.DNA("","GTGATC")
    align1, align2 = btpa.pairwise_align(dna1, dna2)
    assert (align1, align2) == (["A", "G", "T", "G", "-", "T", "C"],["-", "G", "T", "G", "A", "T", "C"])


def test_pairwise_dna_leading_gap_second():
    dna2 = btmm.DNA("","AGTGTC")
    dna1 = btmm.DNA("","GTGATC")
    align1, align2 = btpa.pairwise_align(dna1, dna2)
    assert (align1, align2) == (["-", "G", "T", "G", "A", "T", "C"], ["A", "G", "T", "G", "-", "T", "C"])


def test_pairwise_dna_trailing_gap():
    dna1 = btmm.DNA("","AGTGTCA")
    dna2 = btmm.DNA("","AGTGATC")
    align1, align2 = btpa.pairwise_align(dna1, dna2)
    assert (align1, align2) == (["A", "G", "T", "G", "-", "T", "C", "A"],["A", "G", "T", "G", "A", "T", "C", "-"])


def test_pairwise_dna_identical():
    dna1 = btmm.DNA("","AGTGTC")
    dna2 = btmm.DNA("","AGTGTC")
    align1, align2 = btpa.pairwise_align(dna1, dna2)
    assert (align1, align2) == (["A", "G", "T", "G", "T", "C"],["A", "G", "T", "G", "T", "C"])


def test_pairwise_dna_single_mismatch():
    dna1 = btmm.DNA("","AGTGTC")
    dna2 = btmm.DNA("","AGAGTC")
    align1, align2 = btpa.pairwise_align(dna1, dna2)
    assert (align1, align2) == (["A", "G", "T", "G", "T", "C"],["A", "G", "A", "G", "T", "C"])


def test_pairwise_dna_keyword_arguments():
    dna1 = btmm.DNA("","AGTGTC")
    dna2 = btmm.DNA("","AGTGTC")
    align1, align2 = btpa.pairwise_align(molecule1=dna1, molecule2=dna2)
    assert (align1, align2) == (["A", "G", "T", "G", "T", "C"],["A", "G", "T", "G", "T", "C"])


def test_pairwise_protein_case_insesitive():
    protein1 = btmm.Protein("","avmyerf")
    protein2 = btmm.Protein("","llmfdkf")
    align1, align2 = btpa.pairwise_align(protein1, protein2)
    assert (align1, align2) == (["A", "V", "M", "Y", "E", "R", "F"],["L", "L", "M", "F", "D", "K", "F"])


def test_pairwise_protein_leading_gap_first():
    protein1 = btmm.Protein("","AVMYERF")
    protein2 = btmm.Protein("","QAVMFDKF")
    align1, align2 = btpa.pairwise_align(protein1, protein2)
    assert (align1, align2) == (["-","A", "V", "M", "Y", "E", "R", "F"],["Q", "A", "V", "M", "F", "D", "K", "F"])


def test_pairwise_protein_leading_gap_second():
    protein1 = btmm.Protein("","QAVMYERF")
    protein2 = btmm.Protein("","AVMFDKF")
    align1, align2 = btpa.pairwise_align(protein1, protein2)
    assert (align1, align2) == (["Q", "A", "V", "M", "Y", "E", "R", "F"],["-", "A", "V", "M", "F", "D", "K", "F"])


def test_pairwise_protein_single_point_mutation():
    protein1 = btmm.Protein("","AVMYERF")
    protein2 = btmm.Protein("","AVMYGRF")
    align1, align2 = btpa.pairwise_align(protein1, protein2)
    assert (align1, align2) == (["A", "V", "M", "Y", "E", "R", "F"],["A", "V", "M", "Y", "G", "R", "F"])


def test_pairwise_protein_both_gap_types():
    protein1 = btmm.Protein("","AQVMYERF")
    protein2 = btmm.Protein("","AVMYERPF")
    align1, align2 = btpa.pairwise_align(protein1, protein2)
    assert (align1, align2) == (["A", "Q", "V", "M", "Y", "E", "R", "-", "F"],["A", "-", "V", "M", "Y", "E", "R", "P", "F"])



protein_strategy = st.text(alphabet="ACDEFGHIKLMNPQRSTVWY",
                min_size = 1,
                max_size = 100)

dna_strategy = st.text(alphabet="ACGT",
                min_size = 1,
                max_size = 100)

@given(protein_strategy, protein_strategy)
@settings(max_examples=100, deadline = None)
def test_pairwise_align_protein_hypothesis(protein1, protein2):
    test_protein1 = btmm.Protein("", protein1)
    test_protein2 = btmm.Protein("", protein2)
    align1, align2 = btpa.pairwise_align(test_protein1, test_protein2)
    assert len(align1) == len(align2)
    assert "".join(align1).replace("-", "") == test_protein1.sequence
    assert "".join(align2).replace("-", "") == test_protein2.sequence


@given(dna_strategy, dna_strategy)
@settings(max_examples=100, deadline = None)
def test_pairwise_align_dna_hypothesis(dna1, dna2):
    test_dna1 = btmm.DNA("", dna1)
    test_dna2 = btmm.DNA("", dna2)
    align1, align2 = btpa.pairwise_align(test_dna1, test_dna2)
    assert len(align1) == len(align2)
    assert "".join(align1).replace("-", "") == test_dna1.sequence
    assert "".join(align2).replace("-", "") == test_dna2.sequence


@given(protein_strategy, protein_strategy)
@settings(max_examples=100, deadline=None)
def test_rust_matches_python_protein_hypothesis(seq1, seq2):
    protein1 = btmm.Protein("protein1", seq1)
    protein2 = btmm.Protein("protein2", seq2)

    python1, python2 = btpa.pairwise_align(protein1, protein2)
    rust1, rust2 = btrust._align(seq1, seq2, "protein")

    assert rust1 == "".join(python1)
    assert rust2 == "".join(python2)
    assert len(rust1) == len(rust2)
    assert rust1.replace("-", "") == seq1
    assert rust2.replace("-", "") == seq2


@given(dna_strategy, dna_strategy)
@settings(max_examples=100, deadline=None)
def test_rust_matches_python_dna_hypothesis(seq1, seq2):
    dna1 = btmm.DNA("dna1", seq1)
    dna2 = btmm.DNA("dna2", seq2)

    python1, python2 = btpa.pairwise_align(dna1, dna2)
    rust1, rust2 = btrust._align(seq1, seq2, "dna")

    assert rust1 == "".join(python1)
    assert rust2 == "".join(python2)
    assert len(rust1) == len(rust2)
    assert rust1.replace("-", "") == seq1
    assert rust2.replace("-", "") == seq2



def test_pairwise_affine_gap_implementation():
    protein1 = btmm.Protein("", "AAQPSFCGKLPRMIHLQWTEVFYFLMES")
    protein2 = btmm.Protein("", "AAQAFPSFCGRMIHLQVFYLMLIQKMES")

    align1, align2 = btpa.pairwise_align(protein1, protein2)

    assert "".join(align1).replace("-", "") == protein1.sequence
    assert "".join(align2).replace("-", "") == protein2.sequence
    assert (align1, align2) == (
        ["A", "A", "Q", "-", "-", "P", "S", "F", "C", "G", "K", "L", "P", "R", "M", "I", "H", "L", "Q", "W", "T", "E", "V", "F", "Y", "F", "L", "-", "-", "-", "-", "M", "E", "S"], 
        ["A", "A", "Q", "A", "F", "P", "S", "F", "C", "G", "-", "-", "-", "R", "M", "I", "H", "L", "Q", "-", "-", "-", "V", "F", "Y", "L", "M", "L", "I", "Q", "K", "M", "E", "S"]
        )


def test_pairwise_rust_implementation():
    protein1 = "AAQPSFCGKLPRMIHLQWTEVFYFLMES"
    protein2 = "AAQAFPSFCGRMIHLQVFYLMLIQKMES"

    align1, align2 = btrust._align(protein1, protein2, "protein")
    assert align1.replace("-", "") == protein1
    assert align2.replace("-", "") == protein2
    assert (align1, align2) == (
        "AAQ--PSFCGKLPRMIHLQWTEVFYFL----MES",
        "AAQAFPSFCG---RMIHLQ---VFYLMLIQKMES",
    )


def test_rust_matches_python_affine_gap():
    protein1 = btmm.Protein("Protein 1", "AAQPSFCGKLPRMIHLQWTEVFYFLMES")
    protein2 = btmm.Protein("Protein 2", "AAQAFPSFCGRMIHLQVFYLMLIQKMES")

    python1, python2 = btpa.pairwise_align(protein1, protein2)
    rust1, rust2 = btrust._align(protein1.sequence, protein2.sequence, "protein")

    assert (rust1, rust2) == ("".join(python1), "".join(python2))



proteins = btparse.parse_fasta(str(Path(__file__).parent / "data" / "synthetic_proteins.fasta"), allow_x = True)
@pytest.mark.parametrize("protein1, protein2", list(combinations(proteins, 2)))
def test_rust_matches_python_protein_parametrize(protein1, protein2):
    python1, python2 = btpa.pairwise_align(protein1, protein2)

    rust1, rust2 = btrust._align(protein1.sequence, protein2.sequence, "protein")

    assert rust1 == "".join(python1)
    assert rust2 == "".join(python2)


dna = btparse.parse_fasta(str(Path(__file__).parent / "data" / "synthetic_dna.fasta"), "dna")
@pytest.mark.parametrize("dna1, dna2", list(combinations(dna, 2)))
def test_rust_matches_python_dna_parametrize(dna1, dna2):
    python1, python2 = btpa.pairwise_align(dna1, dna2)

    rust1, rust2 = btrust._align(dna1.sequence, dna2.sequence, "dna")

    assert rust1 == "".join(python1)
    assert rust2 == "".join(python2)



def test_pairwise_align_does_not_mix_protein_dna():
    protein1 = btmm.Protein("","QAVMYERF")
    dna1 = btmm.DNA("","AGTGATTC")
    with pytest.raises(TypeError):
        btpa.pairwise_align(protein1, dna1)


def test_pairwise_align_python_rejects_memory_hog(monkeypatch):
    monkeypatch.setattr(btc, "MAX_ALIGNMENT_CELLS_PYTHON", 256)

    seq1 = btmm.Protein("1", "A" * 20)
    seq2 = btmm.Protein("2", "G" * 20)

    with pytest.raises(ValueError, match="exceeding the limit"):
        btpa.pairwise_align(seq1, seq2)