# Changelog

## phylotypr 0.2.0

### New features

- Add new function
  [`classify_sequences()`](https://mothur.org/phylotypr/reference/classify_sequences.md),
  which has confidence score filtering built-in as an option and is
  multithreaded by default. Can be run by giving the function the name
  of a data frame and the columns corresponding to the sequence id and
  sequences or as a [strollur](https://github.com/mothur/strollur)
  object and function will be output results as appropriate

### Improvements

- Update vignettes and documentation for new functions and features
- Core functions can now process more than one sequence at a time
- Added type checking and check for column names in
  [`classify_sequences()`](https://mothur.org/phylotypr/reference/classify_sequences.md),
  [`filter_taxonomy()`](https://mothur.org/phylotypr/reference/filter_taxonomy.md),
  [`print_taxonomy()`](https://mothur.org/phylotypr/reference/print_taxonomy.md)
- [`build_kmer_database()`](https://mothur.org/phylotypr/reference/build_kmer_database.md)
  now includes `reference_data` argument
- modified
  [`filter_taxonomy()`](https://mothur.org/phylotypr/reference/filter_taxonomy.md)
  to take a data frame rather than a list
- modified
  [`print_taxonomy()`](https://mothur.org/phylotypr/reference/print_taxonomy.md)
  to take a data frame rather than a list and then output mothur-like
  taxonomy string to user specified depth. If specified depth is shorter
  than observed, string gets truncated. If depth is longer than
  specified it adds “\_unclassified”. User can turn off confidence
  scores
- `rdp_trainset9_pds` and `rdp_trainset9` now include reference data as
  a data frame in the reference_data slot; will be incorporated in
  databases provided by `{phylotyprrefdata}`. These were renamed from
  trainset9_rdp and trainset9_pds
- Changed `id` to `sequence_name` across functions
- Modify read_fasta to convert non-ATGCU-. to N; force capitalization.
  This was also applied to built in reference data
- Remove use of [`set.seed()`](https://rdrr.io/r/base/Random.html) to be
  within call to
  [`classify_sequences()`](https://mothur.org/phylotypr/reference/classify_sequences.md)
- Optimized
  [`build_kmer_database()`](https://mothur.org/phylotypr/reference/build_kmer_database.md)
  and
  [`classify_sequences()`](https://mothur.org/phylotypr/reference/classify_sequences.md)
  for performance making both significantly faster

### Deprecated functions

- [`classify_sequence()`](https://mothur.org/phylotypr/reference/classify_sequence.md)

### Bug fixes

- Nothing to repair - if you find something let us know!

## phylotypr 0.1.1

CRAN release: 2025-02-21

This is a patch release to… \* update the citation \* change the license
to GPL v3 \* make a small change to the behavior of the code for when a
k-mer shows up more than once in a sequence

## phylotypr 0.1.0

CRAN release: 2024-10-14

This is the first release version of phylotypr! 🎉

- Functions:
  - [`read_fasta()`](https://mothur.org/phylotypr/reference/read_fasta.md)
  - [`read_taxonomy()`](https://mothur.org/phylotypr/reference/read_taxonomy.md)
  - [`build_kmer_database()`](https://mothur.org/phylotypr/reference/build_kmer_database.md)
  - [`classify_sequence()`](https://mothur.org/phylotypr/reference/classify_sequence.md)
  - [`filter_taxonomy()`](https://mothur.org/phylotypr/reference/filter_taxonomy.md)
  - [`print_taxonomy()`](https://mothur.org/phylotypr/reference/print_taxonomy.md)
- Example datasets:
  - `trainset9_rdp`
  - `trainset9_pds`
- Vignette:
  - [Getting
    started](https://mothur.org/phylotypr/articles/phylotypr.html)
