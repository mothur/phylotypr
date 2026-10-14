# Read in taxonomy files

Read a [mothur-formatted taxonomy
file](https://mothur.org/wiki/taxonomy_file/) into R as a data frame

## Usage

``` r
read_taxonomy(file)
```

## Arguments

- file:

  Either a path to a file, a connection, or literal data (either a
  single string or a raw vector) containing the sequence name and the
  taxonomy information for each sequence.

  Files ending in .gz, .bz2, .xz, or .zip will be automatically
  uncompressed. Files starting with `http://`, `https://`, `ftp://`, or
  `ftps://` will be automatically downloaded. Remote gz files can also
  be autom downloaded and decompressed.

## Value

A data frame with two columns. The `sequence_name` column contains a
name for each sequence and the `taxonomy` column, which contains the
taxonomy for each sequence. The string in the `taxonomy` column is a
series of taxonomic names separated by semi-colons. The string does not
have a semi-colon at the end of the sequence

## Note

There are no checks to insure that each sequence has a unique name. It
is also assumed that each sequence has the same number of taxonomic
levels represented in the second column of the input file.

## Examples

``` r
temp <- tempfile()
write("seqA\tA;B;C;", file = temp)
write("seqB\tA;B; C;", file = temp, append = TRUE)
write("seqC\tA; B;C;", file = temp, append = TRUE)
write("seqD\tA;B;C", file = temp, append = TRUE)
write("seqE\tA;B; C", file = temp, append = TRUE)
write("seqF\tA; B;C", file = temp, append = TRUE)
write("seq G\tA;B;C;", file = temp, append = TRUE)

read_taxonomy(temp)
#>   sequence_name taxonomy
#> 1          seqA    A;B;C
#> 2          seqB    A;B;C
#> 3          seqC    A;B;C
#> 4          seqD    A;B;C
#> 5          seqE    A;B;C
#> 6          seqF    A;B;C
#> 7         seq G    A;B;C
```
