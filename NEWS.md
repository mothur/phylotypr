# phylotypr 0.X.X

* core functions can now process more than one sequence at a time
* change argument name to `classify_sequence()` from `unknown_sequence` to 
  `sequence` <- Revisit whether to remove `classify_sequence` from public facing
  functions? Should we leave a deprecation warning or just pull?
* can run `classify_sequences()` by giving it the name of a data frame and the
  columns corresponding to the sequence id and sequences
* add new function `classify_sequences()`, which has confidence score filtering 
  built-in as an option
* modified `filter_taxonomy()` to take a data frame rather than a list
* modified `print_taxonomy()` to take a data frame rather than a list and then 
  output mothur-like taxonomy string to user specified depth. If specified
  depth is shorter than observed, string gets truncated. If depth is longer than
  specified it adds "_unclassified". User can turn off confidence scores
* Added type checking and check for column names in `classify_sequences()`, 
 `filter_taxonomy()`, `print_taxonomy()`.

# phylotypr 0.1.1

This is a patch release to...
* update the citation
* change the license to GPL v3
* make a small change to the behavior of the code for when a k-mer shows up more than once in a sequence


# phylotypr 0.1.0

This is the first release version of phylotypr! 🎉

- Functions:
    - `read_fasta()`
    - `read_taxonomy()`
    - `build_kmer_database()`
    - `classify_sequence()`
    - `filter_taxonomy()`
    - `print_taxonomy()`
- Example datasets:
    - `trainset9_rdp`
    - `trainset9_pds`
- Vignette:
    - [Getting started](https://mothur.org/phylotypr/articles/phylotypr.html)
