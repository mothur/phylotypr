
<!-- README.md is generated from README.Rmd. Please edit that file -->

# phylotypr

<!-- badges: start -->
<!-- badges: end -->

## Overview

phylotypr is a package for classification based analysis of DNA
sequences. This package primarily implements Naive Bayesian Classifier
from the Ribosomal Database Project. Although you can classify any type
of sequence (assuming you have the proper database), this algorithm is
mainly used to classify 16S rRNA gene sequences.

## Installation

You can install the development version of phylotypr from
[GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
devtools::install_github("riffomonas/phylotypr")
```

You can also get the official release version from CRAN

``` r
install.packages("phylotypr")
```

## Usage

Be sure to see the [Getting Started](articles/phylotypr.html) article to
see an example of how you would build the database and classify
individual and multiple sequences.

You will also need a reference database to classify your sequences. The
`{phylotypr}` package ships with the [RDP training set
(v9)](reference/trainset9_rdp.html). This is relatively small and old
(2010) relative to their latest versions. Also, there are databases
available from greengenes and SILVA. You are encouraged to install newer
databases from the packages on GitHub:

- [RDP’s `{trainset19}`](https://github.com/mothur/trainset19)
- More on the way!

## More information about `{phylotypr}`

You can learn more about the underlying algorithm in the paper that
originally described the algorithm that was published in [*Applied and
Environmental
Microbiology*](https://journals.asm.org/doi/10.1128/aem.00062-07). If
you want to learn more about how this package was created, be sure to
check out the Riffomonas YouTube channel where a [playlist is
available](https://www.youtube.com/watch?v=XjolVT16YNw&list=PLmNrK_nkqBpIZlWa3yGEc2-wX7An2kpCL)
showing every step.
