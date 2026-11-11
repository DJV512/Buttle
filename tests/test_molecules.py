import buttle.models.molecules as btmm
import pytest
from hypothesis import given, strategies as st
from collections import Counter

def test_charge():
    protein = btmm.Protein("", "ArFTyGPl")
    charge = protein.charge()
    assert charge == 1


def test_charge_no_charged_residues():
    protein = btmm.Protein("","AFTYGPL")
    charge = protein.charge()
    assert charge == 0


def test_charge_many_charged_residues():
    protein = btmm.Protein("","AFRTKKYDGEKLREPLDWD")
    charge = protein.charge()
    assert charge == 0


def test_protein_rejects_invalid_aas():
    with pytest.raises(ValueError):
        protein = btmm.Protein("","ASWPMCLRUFCQW")


@given(st.text(alphabet="ACDEFGHIKLMNPQRSTVWY",
               min_size = 1,
               max_size = 100))
def test_charge_hypothesis(protein):
    test_protein = btmm.Protein("", protein)
    expected = test_protein.sequence.count("K") + test_protein.sequence.count("R") - test_protein.sequence.count("D") - test_protein.sequence.count("E")
    assert test_protein.charge() == expected


def test_mw():
    protein = btmm.Protein("","AFRTKKYDGEKLREPLDWD")
    mw = protein.mw()
    assert mw == pytest.approx(2367.6321)


def test_ec():
    protein = btmm.Protein("","AFRTcKKYDgEKLREpLDwD")
    ec = protein.ec()
    assert ec == 7090


def test_glycosylation():
    protein = btmm.Protein("","YQSFRPMNLSQ")
    sites = protein.glycosylation()
    assert sites == [7]


def test_glycosylation_no_sites():
    protein = btmm.Protein("","YQNSFRPMNL")
    sites = protein.glycosylation()
    assert sites == []


def test_glycosylation_cant_have_pro_as_x():
    protein = btmm.Protein("","YQSFRPMNPS")
    sites = protein.glycosylation()
    assert sites == []


def test_glycosylation_overlapping():
    protein = btmm.Protein("","YQSFRPMNNST")
    sites = protein.glycosylation()
    assert sites == [7, 8]



def test_dna_class_rejects_invalid_bases():
    with pytest.raises(ValueError):
        dna = btmm.DNA("","acctagctagucgtagtcagtc")



def test_dna_class_rejeects_n_normally():
    with pytest.raises(ValueError):
        dna = btmm.DNA("", "acctagctagncgtagntcagtnc")


def test_dna_class_accepts_n_when_told_to_do_so():
    dna = btmm.DNA("", "acctagctagncgtagntcagtnc", allow_n=True)
    assert dna.sequence == "acctagctagncgtagntcagtnc".upper()


def test_transcribe():
    dna = btmm.DNA("","actgctcctagctatactca")
    transcription = dna.transcribe()
    assert transcription == "ACUGCUCCUAGCUAUACUCA"


def test_transcribe_with_n():
    dna = btmm.DNA("","actgctcctagctatacntca", allow_n=True)
    transcription = dna.transcribe()
    assert transcription == "ACUGCUCCUAGCUAUACNUCA"
   

def test_translate():
    dna = btmm.DNA("","gtatcatgtctacgtgcgcggactgtcatag")
    t1, t2, t3 = dna.translate()
    assert t1 == "VSCLRARTVI"
    assert t2 == "YHVYVRGLS*"
    assert t3 == "IMSTCADCH"


def test_translate_with_n():
    dna = btmm.DNA("", "acctagctagncgtagntcagtnc", allow_n=True)
    t1, t2, t3 = dna.translate()
    assert t1 == "T*LXVXQX"
    assert t2 == "PS*X*XS"
    assert t3 == "LAXRXSX"


def test_reverse_comp():
    dna = btmm.DNA("", "gtatcatgtctacgtgcgcggactgtcatag")
    rcomp = dna.reverse_comp()
    assert rcomp == "CTATGACAGTCCGCGCACGTAGACATGATAC"


def test_reverse_comp_with_n():
    dna = btmm.DNA("", "gtatcatgtctacgtgcgcggactgtcanag", allow_n=True)
    rcomp = dna.reverse_comp()
    assert rcomp == "CTNTGACAGTCCGCGCACGTAGACATGATAC"


def test_gc():
    dna = btmm.DNA("", "gtatcatgtctacgtgcgcggactgtcatag")
    length, contents, at_content, gc_content = dna.gc()
    assert length == 31
    assert contents == Counter({'G': 9, 'T': 9, 'C': 7, 'A': 6})
    assert at_content == pytest.approx(48.3871)
    assert gc_content == pytest.approx(51.6129)


def test_gc_with_n():
    dna = btmm.DNA("", "gtatcatgntctacgtgcgcggactgtcatag", allow_n=True)
    length, contents, at_content, gc_content = dna.gc()
    assert length == 32
    assert contents == Counter({'G': 9, 'T': 9, 'C': 7, 'A': 6, 'N': 1})
    assert at_content == pytest.approx(46.875)
    assert gc_content == pytest.approx(50)
