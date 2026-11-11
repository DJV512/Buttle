import buttle.sequences as bts


def test_length():
    sequence = "adkfhaksdhfkasf"
    length = bts.length(sequence)
    assert length == 15


def test_length_strip():
    sequence = "  adkfhaksdhfkasf "
    length = bts.length(sequence)
    assert length == 15


def test_find_pos():
    sequence = "asdfasdnaspjnfasdpknfasksjnfs"
    motif = "FAS"
    hits = bts.find_position(sequence, motif)
    assert hits == [4, 14, 21]


def test_find_pos_overlapping():
    sequence = "asdfasdfanaspjnfasdpknfasksjnfs"
    motif = "AsDfA"
    hits = bts.find_position(sequence, motif)
    assert hits == [1, 5]


def test_find_pos_motif_missing():
    sequence = "slkfjasngoigwrvnosanfalsd"
    motif = "dsafsdf"
    hits = bts.find_position(sequence, motif)
    assert hits == []