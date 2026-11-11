import buttle.parsers as btp
import buttle.models.molecules as btmm
import pytest

def test_fasta_reader_minimal_fasta(tmp_path):
    fasta = tmp_path / "test.fasta"

    fasta.write_text(
        ">seq1\n"
        "APLPERDFGNMTREEWQ\n"
        ">seq2\n"
        "APIDKDFGQMPTREDFQ\n"
        ">seq3\n"
        "GPVDREYFQMSKDEYN\n"
        ">seq4\n"
        "WRTDPYMNCQWSLKCF\n"
        ">seq5\n"
        "WRSDPFMNCQWTIRCY\n"
    )

    expected = [
        btmm.Protein(name='seq1', sequence='APLPERDFGNMTREEWQ'), 
        btmm.Protein(name='seq2', sequence='APIDKDFGQMPTREDFQ'), 
        btmm.Protein(name='seq3', sequence='GPVDREYFQMSKDEYN'), 
        btmm.Protein(name='seq4', sequence='WRTDPYMNCQWSLKCF'), 
        btmm.Protein(name='seq5', sequence='WRSDPFMNCQWTIRCY')
    ]

    output = btp.parse_fasta(str(fasta), "protein")
    assert output == expected


def test_fasta_reader_blank_lines(tmp_path):
    fasta = tmp_path / "test.fasta"

    fasta.write_text(
        ">seq1\n"
        "APLPERDFGNMTREEWQ\n"
        "\n"
        ">seq2\n"
        "APIDKDFGQMPTREDFQ\n"
        "\n"
        "\n"
        ">seq3\n"
        "GPVDREYFQMSKDEYN\n"
        "\n"
        ">seq4\n"
        "WRTDPYMNCQWSLKCF\n"
        "\n"
        ">seq5\n"
        "WRSDPFMNCQWTIRCY\n"
    )

    expected = [
        btmm.Protein(name='seq1', sequence='APLPERDFGNMTREEWQ'), 
        btmm.Protein(name='seq2', sequence='APIDKDFGQMPTREDFQ'), 
        btmm.Protein(name='seq3', sequence='GPVDREYFQMSKDEYN'), 
        btmm.Protein(name='seq4', sequence='WRTDPYMNCQWSLKCF'), 
        btmm.Protein(name='seq5', sequence='WRSDPFMNCQWTIRCY')
    ]
    
    output = btp.parse_fasta(str(fasta), "protein")
    assert output == expected


def test_fasta_reader_multiline(tmp_path):
    fasta = tmp_path / "test.fasta"

    fasta.write_text(
        ">seq1\n"
        "APLPERDFGN\n"
        "MTREEWQ\n"
        "\n"
        ">seq2\n"
        "APIDK\n"
        "DFGQMP\n"
        "TREDFQ\n"
        "\n"
        ">seq3\n"
        "G\n"
        "PVDREYFQMSKD\n"
        "EYN\n"
        "\n"
        ">seq4\n"
        "WRTDPYMNCQWSLKCF\n"
        "\n"
        ">seq5\n"
        "WRS\n"
        "DPFMNCQWTIRCY\n"
    )

    expected = [
        btmm.Protein(name='seq1', sequence='APLPERDFGNMTREEWQ'), 
        btmm.Protein(name='seq2', sequence='APIDKDFGQMPTREDFQ'), 
        btmm.Protein(name='seq3', sequence='GPVDREYFQMSKDEYN'), 
        btmm.Protein(name='seq4', sequence='WRTDPYMNCQWSLKCF'), 
        btmm.Protein(name='seq5', sequence='WRSDPFMNCQWTIRCY')
    ]

    output = btp.parse_fasta(str(fasta), "protein")
    assert output == expected


def test_fasta_reader_stress_test(tmp_path):
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

    expected = [
        btmm.Protein(name='seq1', sequence='APLPERDFGNMTREEWQ'), 
        btmm.Protein(name='seq2', sequence='APIDKDFGQMPTREDFQ'), 
        btmm.Protein(name='seq3', sequence='GPVDREYFQMSKDEYN'), 
        btmm.Protein(name='seq4', sequence='WRTDPYMNCQWSLKCF'), 
        btmm.Protein(name='seq5', sequence='WRSDPFMNCQWTIRCY')
    ]
    
    output = btp.parse_fasta(str(fasta), "protein")
    assert output == expected


def test_fasta_reader_one_sequence_allowed(tmp_path):
    fasta = tmp_path / "test.fasta"
    
    fasta.write_text(
        ">seq1\n"
        "APLPERDFGN\n"
        "MTREEWQ\n"
        "\n"
    )

    output = btp.parse_fasta(str(fasta), "protein")
    assert len(output) == 1


def test_fasta_reader_duplicate_names_rejected(tmp_path):
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
        ">seq3\n"
        "WRS\n"
        "DPFMNCQWTIRCY\n"
        "\n"
        "\n"
    )

    with pytest.raises(ValueError):
        output = btp.parse_fasta(str(fasta))


def test_fasta_reader_works_with_dna(tmp_path):
    fasta = tmp_path / "test.fasta"

    fasta.write_text(
    ">seq1\n"
    "ACGTCTAGCTGAGTCGCGCTAGTC\n"
    ">seq2\n"
    "CGCGCGTCTAGTGCGCGTAC\n"
    ">seq3\n"
    "CGTTATCGATTACTGAGCGCGTCGTGTCA\n"
    ">seq4\n"
    "ACTAGCTAGCGGCGCGCGCGTGTCA\n"
    ">seq5\n"
    "CGTACGTACGTCAGTACGTCAGTCGT\n"
    )

    expected = [
        btmm.DNA(name='seq1', sequence='ACGTCTAGCTGAGTCGCGCTAGTC'), 
        btmm.DNA(name='seq2', sequence='CGCGCGTCTAGTGCGCGTAC'), 
        btmm.DNA(name='seq3', sequence='CGTTATCGATTACTGAGCGCGTCGTGTCA'), 
        btmm.DNA(name='seq4', sequence='ACTAGCTAGCGGCGCGCGCGTGTCA'), 
        btmm.DNA(name='seq5', sequence='CGTACGTACGTCAGTACGTCAGTCGT')
    ]

    output = btp.parse_fasta(str(fasta), "dna")
    assert output == expected




def test_parse_fasta_string_works_with_proper_input():
    fasta = ">seq1\nAPLPERDFGNMTREEWQ\n>seq2\nAPIDKDFGQMPTREDFQ\n>seq3\nGPVDREYFQMSKDEYN\n>seq4\nWRTDPYMNCQWSLKCF\n>seq5\nWRSDPFMNCQWTIRCY\n"

    expected = [
        btmm.Protein(name='seq1', sequence='APLPERDFGNMTREEWQ'), 
        btmm.Protein(name='seq2', sequence='APIDKDFGQMPTREDFQ'), 
        btmm.Protein(name='seq3', sequence='GPVDREYFQMSKDEYN'), 
        btmm.Protein(name='seq4', sequence='WRTDPYMNCQWSLKCF'), 
        btmm.Protein(name='seq5', sequence='WRSDPFMNCQWTIRCY')
    ]

    output = btp.parse_fasta(fasta)
    assert output == expected


def test_parse_sequences_rejects_regular_string():
    fasta = "WRSDPFMNCQWTIRCY"

    with pytest.raises(ValueError):
        output = btp.parse_fasta(fasta)


def test_parse_sequences_rejects_missing_carrot():
    fasta = ">seq1\nAPLPERDFGNMTREEWQ\n>seq2\nAPIDKDFGQMPTREDFQ\nseq3\nGPVDREYFQMSKDEYN\n>seq4\nWRTDPYMNCQWSLKCF\n>seq5\nWRSDPFMNCQWTIRCY\n"

    with pytest.raises(ValueError):
        output = btp.parse_fasta(fasta)


def test_parse_sequences_rejects_missing_sequence():
    fasta = ">seq1\nAPLPERDFGNMTREEWQ\n>seq2\n>seq3\nGPVDREYFQMSKDEYN\n>seq4\nWRTDPYMNCQWSLKCF\n>seq5\nWRSDPFMNCQWTIRCY\n"

    with pytest.raises(ValueError):
        output = btp.parse_fasta(fasta)


def test_parse_sequences_rejects_duplicate_names():
    fasta = ">seq1\nAPLPERDFGNMTREEWQ\n>seq2\nAPIDKDFGQMPTREDFQ\n>seq3\nGPVDREYFQMSKDEYN\n>seq4\nWRTDPYMNCQWSLKCF\n>seq2\nWRSDPFMNCQWTIRCY\n"

    with pytest.raises(ValueError):
        output = btp.parse_fasta(fasta)


def test_parse_fasta_string_empty_input():
    with pytest.raises(ValueError, match="No sequences found"):
        btp.parse_fasta("")


def test_parse_fasta_stop_codon_error():
    fasta = ">seq1\nAPLPERDF*NMTREEWQ\n>seq2\nAPIDKDFGQMPTREDFQ\n>seq3\nGPVDREYFQMSKDEYN\n>seq4\nWRTDPYMNCQWSLKCF\n>seq5\nWRSDPFMNCQWTIRCY\n"

    with pytest.raises(ValueError, match="Protein sequences contain stop codons."):
        output = btp.parse_fasta(fasta, stop_action="error")



def test_parse_fasta_stop_codon_truncate():
    fasta = ">seq1\nAPLPERDF*NMTREEWQ\n>seq2\nAPIDKDFGQMPTREDFQ\n>seq3\nGPVDREYFQMSKDEYN\n>seq4\nWRTDPYMNCQWSLKCF\n>seq5\nWRSDPFMNCQWTIRCY\n"

    expected = [
        btmm.Protein(name='seq1', sequence='APLPERDF'), 
        btmm.Protein(name='seq2', sequence='APIDKDFGQMPTREDFQ'), 
        btmm.Protein(name='seq3', sequence='GPVDREYFQMSKDEYN'), 
        btmm.Protein(name='seq4', sequence='WRTDPYMNCQWSLKCF'), 
        btmm.Protein(name='seq5', sequence='WRSDPFMNCQWTIRCY')
    ]

    output = btp.parse_fasta(fasta, stop_action="truncate")
    assert output == expected


def test_parse_fasta_stop_codon_exclude():
    fasta = ">seq1\nAPLPERDF*NMTREEWQ\n>seq2\nAPIDKDFGQMPTREDFQ\n>seq3\nGPVDREYFQMSKDEYN\n>seq4\nWRTDPYMNCQWSLKCF\n>seq5\nWRSDPFMNCQWTIRCY\n"

    expected = [
        btmm.Protein(name='seq2', sequence='APIDKDFGQMPTREDFQ'), 
        btmm.Protein(name='seq3', sequence='GPVDREYFQMSKDEYN'), 
        btmm.Protein(name='seq4', sequence='WRTDPYMNCQWSLKCF'), 
        btmm.Protein(name='seq5', sequence='WRSDPFMNCQWTIRCY')
    ]

    output = btp.parse_fasta(fasta, stop_action="exclude")
    assert output == expected


def test_parse_fasta_auto_detect_protein_without_x():
    fasta = ">seq1\nAPLPERDFGNMTREEWQ\n>seq2\nAPIDKDFGQMPTREDFQ\n>seq3\nGPVDREYFQMSKDEYN\n>seq4\nWRTDPYMNCQWSLKCF\n>seq5\nWRSDPFMNCQWTIRCY\n"
    
    expected = [
        btmm.Protein(name='seq1', sequence='APLPERDFGNMTREEWQ'), 
        btmm.Protein(name='seq2', sequence='APIDKDFGQMPTREDFQ'), 
        btmm.Protein(name='seq3', sequence='GPVDREYFQMSKDEYN'), 
        btmm.Protein(name='seq4', sequence='WRTDPYMNCQWSLKCF'), 
        btmm.Protein(name='seq5', sequence='WRSDPFMNCQWTIRCY')
    ]

    output = btp.parse_fasta(fasta)
    assert output == expected
    

def test_parse_fasta_auto_detect_protein_with_x():
    fasta = ">seq1\nAPLPERDFGNMTREEWQ\n>seq2\nAPIDKDFGQMPTREDFQ\n>seq3\nGPVDREYFQXSKDEYN\n>seq4\nWRTDPYMNCQWSLKCF\n>seq5\nWRSDPFMNCQWTIRCY\n"
    
    expected = [
        btmm.Protein(name='seq1', sequence='APLPERDFGNMTREEWQ'), 
        btmm.Protein(name='seq2', sequence='APIDKDFGQMPTREDFQ'), 
        btmm.Protein(name='seq3', sequence='GPVDREYFQXSKDEYN', allow_x=True), 
        btmm.Protein(name='seq4', sequence='WRTDPYMNCQWSLKCF'), 
        btmm.Protein(name='seq5', sequence='WRSDPFMNCQWTIRCY')
    ]

    output = btp.parse_fasta(fasta)
    assert output == expected


def test_parse_fasta_protein_detect_protein_with_x_allow_x_false():
    fasta = ">seq1\nAPLPERDFGNMTREEWQ\n>seq2\nAPIDKDFGQMPTREDFQ\n>seq3\nGPVDREYFQXSKDEYN\n>seq4\nWRTDPYMNCQWSLKCF\n>seq5\nWRSDPFMNCQWTIRCY\n"

    with pytest.raises(ValueError, match="GPVDREYFQXSKDEYN is not a valid Protein sequence."):
        output = btp.parse_fasta(fasta, molecule_type="protein", allow_x=False)


def test_parse_fasta_protein_detect_protein_with_x_allow_x_true():
    fasta = ">seq1\nAPLPERDFGNMTREEWQ\n>seq2\nAPIDKDFGQMPTREDFQ\n>seq3\nGPVDREYFQXSKDEYN\n>seq4\nWRTDPYMNCQWSLKCF\n>seq5\nWRSDPFMNCQWTIRCY\n"
    
    expected = [
        btmm.Protein(name='seq1', sequence='APLPERDFGNMTREEWQ'), 
        btmm.Protein(name='seq2', sequence='APIDKDFGQMPTREDFQ'), 
        btmm.Protein(name='seq3', sequence='GPVDREYFQXSKDEYN', allow_x=True), 
        btmm.Protein(name='seq4', sequence='WRTDPYMNCQWSLKCF'), 
        btmm.Protein(name='seq5', sequence='WRSDPFMNCQWTIRCY')
    ]

    output = btp.parse_fasta(fasta, molecule_type="protein", allow_x=True)
    assert output == expected


def test_parse_fatsa_dna_with_n_allow_n_false(tmp_path):
    fasta = tmp_path / "test.fasta"
    
    fasta.write_text(
    ">seq1\n"
    "ACGTCTAGCTGAGTCGCGCTAGTC\n"
    ">seq2\n"
    "CGCGCGTCTAGTGCGCGTAC\n"
    ">seq3\n"
    "CGTTATCGATTACTGAGCGCGTCGTGTCA\n"
    ">seq4\n"
    "ACTAGCTAGCGGCGCGCGCGTGTCA\n"
    ">seq5\n"
    "CGTACGTACGTCAGTACGTCAGTCNT\n"
    )
    with pytest.raises(ValueError, match = "CGTACGTACGTCAGTACGTCAGTCNT is not a valid DNA sequence."):
        output = btp.parse_fasta(str(fasta), "dna", allow_n=False)


def test_parse_fatsa_dna_with_n_allow_n_true(tmp_path):
    fasta = tmp_path / "test.fasta"
    
    fasta.write_text(
    ">seq1\n"
    "ACGTCTAGCTGAGTCGCGCTAGTC\n"
    ">seq2\n"
    "CGCGCGTCTAGTGCGCGTAC\n"
    ">seq3\n"
    "CGTTATCGATTACTGAGCGCGTCGTGTCA\n"
    ">seq4\n"
    "ACTAGCTAGCGGCGCGCGCGTGTCA\n"
    ">seq5\n"
    "CGTACGTACGTCAGTACGTCAGTCNT\n"
    )

    expected = [
        btmm.DNA(name='seq1', sequence='ACGTCTAGCTGAGTCGCGCTAGTC'), 
        btmm.DNA(name='seq2', sequence='CGCGCGTCTAGTGCGCGTAC'), 
        btmm.DNA(name='seq3', sequence='CGTTATCGATTACTGAGCGCGTCGTGTCA'), 
        btmm.DNA(name='seq4', sequence='ACTAGCTAGCGGCGCGCGCGTGTCA'), 
        btmm.DNA(name='seq5', sequence='CGTACGTACGTCAGTACGTCAGTCNT', allow_n=True)
    ]

    output = btp.parse_fasta(str(fasta), "dna", allow_n=True)
    assert output == expected


@pytest.mark.parametrize(
    "contents",
    [
        pytest.param("", id="empty-file"),
        pytest.param(" \n\t\n", id="whitespace-only"),
        pytest.param(">a\n", id="header-only"),
        pytest.param(">a\n>b\nAC\n", id="empty-first-record"),
        pytest.param(
            ">a\nAC\n>b\n>c\nAT\n",
            id="empty-middle-record",
        ),
        pytest.param(">a\nAC\n>b\n", id="empty-last-record"),
        pytest.param(">a\n \t\n>b\nAC\n", id="blank-sequence"),
    ],
)
def test_fasta_file_rejects_empty_records(tmp_path, contents):
    fasta = tmp_path / "input.fasta"
    fasta.write_text(contents)

    with pytest.raises(ValueError, match="No FASTA records|Missing sequence"):
        btp.parse_fasta(str(fasta))


def test_fasta_file_rejects_sequence_before_header(tmp_path):
    fasta = tmp_path / "input.fasta"
    fasta.write_text("ACGT\n>a\nATGC\n")

    with pytest.raises(ValueError, match="before first header"):
        btp.parse_fasta(str(fasta))


def test_fasta_file_preserves_wrapped_sequences_and_names(tmp_path):
    fasta = tmp_path / "input.fasta"
    fasta.write_text(
        "\n"
        ">alpha\n"
        "ac\n"
        "GT\n"
        "\n"
        ">beta\n"
        "At\n"
        "gc\n"
        "\n"
    )

    sequences = btp.parse_fasta(str(fasta), molecule_type="dna")

    assert [(seq.name, seq.sequence) for seq in sequences] == [
        ("alpha", "ACGT"),
        ("beta", "ATGC"),
    ]


@pytest.fixture(params=["file", "string"])
def fasta_source(request, tmp_path):
    def make_source(contents):
        if request.param == "string":
            return contents

        fasta = tmp_path / "input.fasta"
        fasta.write_text(contents)
        return str(fasta)

    return make_source


def test_fasta_rejects_empty_sequence_after_truncation(fasta_source):
    source = fasta_source(">empty\n*WE\n>valid\nWF\n")

    with pytest.raises(ValueError, match="empty after stop-codon truncation"):
        btp.parse_fasta(source, stop_action="truncate")


def test_fasta_rejects_exclusion_of_every_record(fasta_source):
    source = fasta_source(">a\nWE*\n>b\nW*F\n")

    with pytest.raises(ValueError, match="No sequences remain"):
        btp.parse_fasta(source, stop_action="exclude")


def test_fasta_truncates_at_first_stop_codon(fasta_source):
    source = fasta_source(">a\nWE*AA*GG\n>b\nWF\n")

    sequences = btp.parse_fasta(source, stop_action="truncate")

    assert [(seq.name, seq.sequence) for seq in sequences] == [
        ("a", "WE"),
        ("b", "WF"),
    ]


def test_fasta_exclusion_preserves_remaining_records(fasta_source):
    source = fasta_source(">keep1\nWE\n>remove\nW*F\n>keep2\nWY\n")

    sequences = btp.parse_fasta(source, stop_action="exclude")

    assert [(seq.name, seq.sequence) for seq in sequences] == [
        ("keep1", "WE"),
        ("keep2", "WY"),
    ]


def test_fasta_rejects_unknown_stop_action(fasta_source):
    source = fasta_source(">a\nWE\n>b\nWF\n")

    with pytest.raises(ValueError, match="Unsupported stop action"):
        btp.parse_fasta(source, stop_action="ignore")