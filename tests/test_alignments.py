from buttle.models import alignments as btma

def test_add_consensus():
    profile = btma.Profile(
        name=['seq1', 'seq2', 'seq3'], 
        alignment=[
            ['A', 'P', 'L', 'P', 'E', 'R', 'D', 'F', 'G', 'N', 'M', '-', 'T', 'R', 'E', 'E', 'W', 'Q'], 
            ['A', 'P', 'I', '-', 'D', 'K', 'D', 'F', 'G', 'Q', 'M', 'P', 'T', 'R', 'E', 'D', 'F', 'Q'], 
            ['G', 'P', 'V', '-', 'D', 'R', 'E', 'Y', 'F', 'Q', 'M', '-', 'S', 'K', 'D', 'E', 'Y', 'N']
        ],
        molecule_type="protein"
    )

    expected = btma.Profile(
        name=['seq1', 'seq2', 'seq3', 'Consensus'], 
        alignment=[
            ['A', 'P', 'L', 'P', 'E', 'R', 'D', 'F', 'G', 'N', 'M', '-', 'T', 'R', 'E', 'E', 'W', 'Q'], 
            ['A', 'P', 'I', '-', 'D', 'K', 'D', 'F', 'G', 'Q', 'M', 'P', 'T', 'R', 'E', 'D', 'F', 'Q'], 
            ['G', 'P', 'V', '-', 'D', 'R', 'E', 'Y', 'F', 'Q', 'M', '-', 'S', 'K', 'D', 'E', 'Y', 'N'],
            ['A', 'P', 'L', '-', 'D', 'R', 'D', 'F', 'G', 'Q', 'M', '-', 'T', 'R', 'E', 'E', 'W', 'Q'], 
        ],
        molecule_type="protein"
    )

    profile.add_consensus()
    assert profile == expected