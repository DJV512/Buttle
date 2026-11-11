from hypothesis import given, Verbosity, settings, strategies as st
import buttle.multiple_alignment as btmsa
import buttle.models.molecules as btmm
import buttle.models.alignments as btma
import pytest

protein_strategy = st.text(alphabet="ACDEFGHIKLMNPQRSTVWY",
               min_size = 80,
               max_size = 120)
proteins_strategy = st.lists(protein_strategy, min_size = 10, max_size = 10)
@given(proteins_strategy)
@settings(verbosity=Verbosity.verbose, max_examples=10, deadline=None)
def test_msa_protein_hypothesis(proteins):

    originals = {}
    sequences = ""
    for n, protein in enumerate(proteins):
        sequences += f">protein{n}\n{protein}\n"
        originals[f"protein{n}"] = protein

    profile = btmsa.msa(sequences, -10, -3, backend="rust")

    lengths = set()
    for name, align in zip(profile.name, profile.alignment, strict=True):
        lengths.add(len(align))
        assert originals[name] == "".join(align).replace("-","")
    assert len(lengths) == 1


dna_strategy = st.text(alphabet="ACGT",
               min_size = 80,
               max_size = 120)
dnas_strategy = st.lists(dna_strategy, min_size = 10, max_size = 10)
@given(dnas_strategy)
@settings(verbosity=Verbosity.verbose, max_examples=10, deadline=None)
def test_msa_dna_hypothesis(dnas):

    originals = {}
    sequences = ""
    for n, dna in enumerate(dnas):
        sequences += f">dna{n}\n{dna}\n"
        originals[f"dna{n}"] = dna

    profile = btmsa.msa(sequences, -10, -3, backend="rust")

    lengths = set()
    for name, align in zip(profile.name, profile.alignment, strict=True):
        lengths.add(len(align))
        assert originals[name] == "".join(align).replace("-","")
    assert len(lengths) == 1



def test_msa_correctly_runs_with_string_fasta_input():
    fasta = ">seq1\nAPLPERDFGNMTREEWQ\n>seq2\nAPIDKDFGQMPTREDFQ\n>seq3\nGPVDREYFQMSKDEYN\n>seq4\nWRTDPYMNCQWSLKCF\n>seq5\nWRSDPFMNCQWTIRCY\n"
    profile = btmsa.msa(fasta, -10, -3, backend="python")

    originals = {
        "seq1": "APLPERDFGNMTREEWQ",
        "seq2": "APIDKDFGQMPTREDFQ",
        "seq3": "GPVDREYFQMSKDEYN",
        "seq4": "WRTDPYMNCQWSLKCF",
        "seq5": "WRSDPFMNCQWTIRCY"
    }

    lengths = set()
    for name, align in zip(profile.name, profile.alignment, strict=True):
        lengths.add(len(align))
        assert originals[name] == "".join(align).replace("-","")
    assert len(lengths) == 1


def test_msa_correctly_runs_with_file_fasta_input(tmp_path):
    fasta = tmp_path / "test.fasta"

    fasta.write_text(
        ">seq1\n"
        "APLPERDFGN\n"
        "MTREEWQ\n"
        "\n"
        ">seq2\n"
        "APIDK\n"
        "DFGQMPTREDFQ\n"
        "\n"
        "\n"
        ">seq3\n"
        "GPVDREYFQMSKD\n"
        "EYN\n"
        "\n"
        ">seq4\n"
        "WRTDPY\n"
        "MNCQWSLKCF\n"
        "\n"
        "\n"
        "\n"
        ">seq5\n"
        "WRS\n"
        "DPFMNCQWTIRCY\n"
        "\n"
        "\n"
    )

    profile = btmsa.msa(str(fasta), -10, -3, backend="python")

    originals = {
        "seq1": "APLPERDFGNMTREEWQ",
        "seq2": "APIDKDFGQMPTREDFQ",
        "seq3": "GPVDREYFQMSKDEYN",
        "seq4": "WRTDPYMNCQWSLKCF",
        "seq5": "WRSDPFMNCQWTIRCY"
    }

    lengths = set()
    for name, align in zip(profile.name, profile.alignment, strict=True):
        lengths.add(len(align))
        assert originals[name] == "".join(align).replace("-","")
    assert len(lengths) == 1


def test_msa_rust_vs_python_static_example():
    fasta = ">seq1\nAPLPERDFGNMTREEWQ\n>seq2\nAPIDKDFGQMPTREDFQ\n>seq3\nGPVDREYFQMSKDEYN\n>seq4\nWRTDPYMNCQWSLKCF\n>seq5\nWRSDPFMNCQWTIRCY\n"
    python_profile = btmsa.msa(fasta, -10, -3, backend="python")
    rust_profile = btmsa.msa(fasta, -10, -3, backend="rust")

    originals = {
        "seq1": "APLPERDFGNMTREEWQ",
        "seq2": "APIDKDFGQMPTREDFQ",
        "seq3": "GPVDREYFQMSKDEYN",
        "seq4": "WRTDPYMNCQWSLKCF",
        "seq5": "WRSDPFMNCQWTIRCY"
    }

    python_lengths = set()
    for name, align in zip(python_profile.name, python_profile.alignment, strict=True):
        python_lengths.add(len(align))
        assert originals[name] == "".join(align).replace("-","")
    assert len(python_lengths) == 1

    rust_lengths = set()
    for name, align in zip(rust_profile.name, rust_profile.alignment, strict=True):
        rust_lengths.add(len(align))
        assert originals[name] == "".join(align).replace("-","")
    assert len(rust_lengths) == 1

    assert python_profile == rust_profile


AMINO_ACIDS = "ACDEFGHIKLMNPQRSTVWY"
protein_sequences = st.text(
    alphabet=AMINO_ACIDS,
    min_size=1,
    max_size=100,
)
@given(seq1=protein_sequences, seq2=protein_sequences, seq3=protein_sequences, seq4=protein_sequences, seq5=protein_sequences, seq6=protein_sequences)
@settings(max_examples=500, deadline=None)
def test_msa_rust_vs_python_protein_hypothesis(seq1, seq2, seq3, seq4, seq5, seq6):
    fasta = f">seq1\n{seq1}\n>seq2\n{seq2}\n>seq3\n{seq3}\n>seq4\n{seq4}\n>seq5\n{seq5}\n>seq6\n{seq6}"

    python_profile = btmsa.msa(fasta, -10, -3, backend="python", molecule_type = "protein")
    rust_profile = btmsa.msa(fasta, -10, -3, backend="rust", molecule_type="protein")

    assert python_profile == rust_profile


DNA_BASES = "ACGT"
dna_sequences = st.text(
    alphabet=DNA_BASES,
    min_size=1,
    max_size=100,
)
@given(seq1=dna_sequences, seq2=dna_sequences, seq3=dna_sequences, seq4=dna_sequences, seq5=dna_sequences, seq6=dna_sequences)
@settings(max_examples=500, deadline=None)
def test_msa_rust_vs_python_dna_hypothesis(seq1, seq2, seq3, seq4, seq5, seq6):
    fasta = f">seq1\n{seq1}\n>seq2\n{seq2}\n>seq3\n{seq3}\n>seq4\n{seq4}\n>seq5\n{seq5}\n>seq6\n{seq6}"

    python_profile = btmsa.msa(fasta, -10, -3, backend="python", molecule_type="dna")
    rust_profile = btmsa.msa(fasta, -10, -3, backend="rust", molecule_type="dna")

    assert python_profile == rust_profile

