# Changelog

All notable changes to Buttle will be listed here.

## [0.4.0] - 2026-11-xx

### Added
- New DNA class, in models/molecules.py
- FASTA parsing support for DNA and for handling stop codons in Proteins
- Imported public facing Functions and Classes into __init__.py for easy importing (from buttle import XX")
- Support for DNA Multiple Sequence Alignment
- Added `--backend` argument to dna_align to choose "rust" or "python" (Rust defaults)
- Added `--stop-action` argument to FASTA parsing to select option when stop codons ("*") are present ("error", "truncate", "exclude"). Defaults to "error"
- Added `--molecule-type` argument to FASTA parsing to declare "dna" or "protein" (defaults to "auto", inferring from bases/residues present)
- Added `--tree` argument parameter to MSA which optionally prints the final alignment Cluster tree to the terminal
- Added molecule_type parameter to "Profile", for use in image rendering
- New validation.py module
    - check_alignment_size() ensures that an alignment will not use excessive memory before running
    - Uses different limits for Python and Rust backends

### Changed 
- Renamed msa.py to multiple_alignment.py to avoid name collision with the function msa()
- Merged DNA and Protein pairwise align functions into one
- Modified the pairwise rust _align function to work with Protein and DNA
- Modified the multiple alignment rust _merge function to work with Protein and DNA
- Moved models.py into a models subpackage, and split into molecules.py and alignments.py
- Moved DNA related functions into the new DNA class as methods
- Deleted decorators.py (DNA validation moved into the class during instantiation)
- Merged DNA and Protein image saving and terminal printing
- Moved function add_consensus() from multiple_alignment.py to a method on the profile class (models/alignments.py)
- Improved Cluster tree printing

### Testing
- Refactored relevant tests to work with the new interfaces
- Added Hypothesis-based testing to verify Python DNA multiple alignment algorithm retains integrity of all sequences
- Added Hypothesis-based testing to verify consistency between Python and Rust DNA alignment algorithms
- Added several tests for the new parser logic


## [0.3.0] - 2026-11-02
- Rebranded to Buttle: A command line utility for biological sequence analysis. At your service.


## [0.2.0] - 2026-10-30

### Added
- Accelerated protein MSA using a Rust backend, substantially improving performance over the Python implementation (Python functions still available)
- Protein MSA can be tuned with custom affine `--gap-open` and `--gap-extend` penalties
- Introduced _non_positive_int to validate gap penalty arguments
- Argument `--backend` added to MSA and protein_align, with choices "rust" and "python" (Rust defaults) 

### Changed
- Bundled JetBrains Mono Regular Font for cross-platform image generation, replacing the macOS-specific Menlo font dependency
- Updated the README

### Fixed
- Type _positive_int was updated so that floating point arguments return a helpful error message
- Changed Rust alignment scoring from f32 to f64 to improve numerical consistency with Python
- Added the proper repo URL to the README install instructions

### Testing
- Added Hypothesis-based randomized testing to validate consistency between the Rust and Python Protein MSA implementation
- Started a testing module for the image renderer


## [0.1.0] - 2026-10-20

### Added

#### DNA sequence tools
- Transcribe to RNA
- Translate in three reading frames
- Generate reverse complement
- Calculate GC content
- Pairwise sequence alignment

#### Protein sequence tools
- Calculate net charge
- Calculate molecular weight
- Calculate extinction coefficient
- Search for and highlight N-linked glycosylation sites
- Compare two sequences and look for conserved motifs of varying lengths
- Pairwise sequence alignment (default implementation in Rust, Python functions available)
- Multiple sequence alignment (MSA) from FASTA files or FASTA-Formatted strings

#### General sequence tools
- Search for sequences and report their positions
- Returns the length of a sequence
- Reports the percent identity between all sequences from a FASTA file

#### Rendering
- Alignments output to terminal (when no `--save` argument provided)
- Generate publication-quality PNG and TIFF alignment images
    - Save Pairwise DNA and Protein alignments
    - Save Multiple Sequence DNA and Protein Alignments
    - Alignment row length can be tuned with `--wrap`
    - Argument `--real-ruler` displays amino acid residue number based on the first sequence, skipping gaps
- Argument `--diff` displays alignment differences relative to a reference sequence
- MSA optionally adds a consensus sequence with `--consensus`

#### Testing
- Comprehensive unit testing of all Python algorithm modules
- End to end testing of the MSA pipeline, verifying preservation of original sequences, uniform alignment lengths, and consistency between the Rust and Python pairwise_protein implementations