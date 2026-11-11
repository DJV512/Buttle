from .models import molecules as btmm
from pathlib import Path
from . import constants as btc



def parse_fasta(
        source: str,
        molecule_type: str = "auto",
        allow_x: bool = False,
        allow_n: bool = False,
        stop_action: str = "error",
    ) -> list[btmm.Protein] | list[btmm.DNA]:
    """
    Parse a FASTA file or formatted string into Protein or DNA objects.

    When molecule_type="auto", infer the molecule type from the
    combined sequence alphabet. Explicit molecule types override detection.
    """

    # If given an empty string, or no sequence input, raise an error.
    if not source:
        raise ValueError("No sequences found.")

    # Validate stop-codon handling.
    if stop_action not in {"error", "truncate", "exclude"}:
        raise ValueError(f"Unsupported stop action: {stop_action}")

    # Determine whether the source is FASTA-formatted text or a file.
    if source.lstrip().startswith(">"):
        check_fasta_string(source)
        data = source.splitlines()

    elif Path(source).is_file():
        path = Path(source)
        check_fasta_file(path)

        with path.open("r") as f:
            data = f.readlines()

    else:
        raise ValueError(f"Unreadable source {source}.")

    # Parse FASTA records.
    # Sequences may span multiple lines, with blank lines permitted.
    sequences = []
    name = None
    sequence = ""

    for line in data:
        line = line.strip()

        if line.startswith(">"):

            # Finalize the previous record before starting another.
            if name is not None:
                if not sequence:
                    raise ValueError(f"Missing sequence for header {name!r}")
                sequences.append(btmm.Sequence(name=name, sequence=sequence))
            name = line[1:].strip()
            sequence = ""

        elif line:
            if name is None:
                raise ValueError("FASTA sequence found before first header.")
            sequence += line.upper()

    # Validate and append the final record.
    if name is None:
        raise ValueError("No FASTA records found.")
    if not sequence:
        raise ValueError(f"Missing sequence for header {name!r}")
    sequences.append(btmm.Sequence(name=name, sequence=sequence))

    # Check to make sure there are no duplicate names.
    names = [seq.name for seq in sequences]
    if len(names) != len(set(names)):
        raise ValueError("FASTA sequence names must be unique.")

    # Construct the initial alphabet.
    base_or_aa = set()
    amino_acids = set(btc.AMINO_ACIDS)

    for seq in sequences:
        base_or_aa.update(seq.sequence)

    # Process stop codons according to the requested behavior.
    if "*" in base_or_aa:
        if stop_action == "error":
            raise ValueError("Protein sequences contain stop codons.")

        elif stop_action == "truncate":
            for seq in sequences:
                seq.sequence = seq.sequence.split("*")[0]
                if not seq.sequence:
                    raise ValueError(f"Sequence {seq.name!r} is empty after stop-codon truncation.")

        elif stop_action == "exclude":
            sequences = [seq for seq in sequences if "*" not in seq.sequence]

    if not sequences:
        raise ValueError("No sequences remain after stop-codon processing.")

    # Rebuild the alphabet after truncation or exclusion.
    base_or_aa = set()
    for seq in sequences:
        base_or_aa.update(seq.sequence)

    # Construct Protein or DNA objects.
    if molecule_type == "protein":
        final_sequences = [btmm.Protein(seq.name, seq.sequence, allow_x=allow_x) for seq in sequences]

    elif molecule_type == "dna":
        final_sequences = [btmm.DNA(seq.name, seq.sequence, allow_n=allow_n) for seq in sequences]

    elif molecule_type == "auto":
        if base_or_aa.issubset(set("ACGT")):
            final_sequences = [btmm.DNA(seq.name, seq.sequence, allow_n=False) for seq in sequences]

        elif base_or_aa.issubset(set("ACGNT")):
            final_sequences = [btmm.DNA(seq.name, seq.sequence, allow_n=True) for seq in sequences]

        elif base_or_aa.issubset(amino_acids):
            final_sequences = [btmm.Protein(seq.name, seq.sequence, allow_x=False) for seq in sequences]

        elif base_or_aa.issubset(amino_acids | {"X"}):
            final_sequences = [btmm.Protein(seq.name, seq.sequence, allow_x=True) for seq in sequences]

        else:
            raise ValueError("Auto detection could not resolve molecule type.")

    else:
        raise ValueError(f"Unsupported molecule type: {molecule_type}")

    if not final_sequences:
        raise ValueError("No valid DNA or Protein sequences were found.")

    return final_sequences


def check_fasta_file(path: Path) -> None:
    if path.stat().st_size > btc.MAX_FASTA_BYTES:
        raise ValueError("FASTA file exceeds size limit")


def check_fasta_string(source: str) -> None:
    if len(source.encode("utf-8")) > btc.MAX_FASTA_BYTES:
        raise ValueError("FASTA string input exceeds size limit")