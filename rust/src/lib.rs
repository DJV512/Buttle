mod constants;
mod pairwise;
mod merge;

pub use pairwise::pairwise_align;
pub use merge::merge;
pub use constants::Profile;

#[cfg(feature = "python")]
use pyo3::prelude::*;

/// Globally align two sequences using affine gap penalties.
///
/// Inputs must be normalized sequences from validated Protein or DNA objects.
/// This internal function does not perform input validation.
///
/// Returns the two aligned sequences, with '-' characters marking gaps.
#[cfg(feature = "python")]
#[pyfunction]
fn _align(seq1: &str, seq2: &str, molecule_type: &str) -> (String, String) {
    pairwise_align(seq1, seq2, molecule_type)
}

/// Globally align two profiles, made up of already aligned sequences, using affine gap penalties.
///
/// Inputs are two Profiles to merge (separated out into their constituent lists of names and sequences to pass to rust),
/// as well as the gap open and gap extend penalties to use. This internal function does not perform input validation.
///
/// Returns a tuple of a list of sequence names, and a list of sequence alignments.
#[cfg(feature = "python")]
#[pyfunction]
fn _merge(
    left_names: Vec<String>,
    left_sequences: Vec<String>,
    right_names: Vec<String>,
    right_sequences: Vec<String>,
    molecule_type: &str,
    gap_open: f64,
    gap_extend: f64,
) -> (Vec<String>, Vec<String>) {
    let left = Profile {
        names: left_names,
        sequences: left_sequences
            .into_iter()
            .map(|s| s.chars().collect())
            .collect(),
    };

    let right = Profile {
        names: right_names,
        sequences: right_sequences
            .into_iter()
            .map(|s| s.chars().collect())
            .collect(),
    };

    let result = merge(&left, &right, molecule_type, gap_open, gap_extend);

    (
        result.names,
        result.sequences
            .into_iter()
            .map(|seq| seq.into_iter().collect())
            .collect(),
    )
}

#[cfg(feature = "python")]
#[pymodule]
fn buttle_rust(m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add_function(wrap_pyfunction!(_align, m)?)?;
    m.add_function(wrap_pyfunction!(_merge, m)?)?;
    Ok(())
}


