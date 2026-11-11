from . import constants as btc

def check_alignment_size(len1: int, len2: int, backend: str) -> None:
    cells = (len1 + 1) * (len2 + 1)

    if backend == "rust":
        limit = btc.MAX_ALIGNMENT_CELLS_RUST
    elif backend == "python":
        limit = btc.MAX_ALIGNMENT_CELLS_PYTHON
    else:
        raise ValueError(f"Unknown backend {backend}.")

    if cells > limit:
        raise ValueError(f"Alignment would require {cells:,} matrix cells, exceeding the limit of {limit:,}.")