library(microbenchmark)
library(tidyverse)
library(data.table)
library(Matrix)

set.seed(19760620)

n_kmers <- 64
n_genera <- 15
n <- 10

df_grow_rbind <- function(n){
  df <- data.frame(kmer = NULL, genus = NULL, count = NULL)
  for(i in 1:n) {
    df <- rbind(df, list(kmer = sample(n_kmers, 1),
                      genus = sample(n_genera, 1),
                      count = sample(10, 1)))
  }
  df
}

tbl_grow_rbind <- function(n){
  df <- tibble(kmer = NULL, genus = NULL, count = NULL)
  for(i in 1:n) {
    df <- bind_rows(df, tibble_row(kmer = sample(n_kmers, 1),
                         genus = sample(n_genera, 1),
                         count = sample(10, 1)))
  }
  df
}

dt_grow_rbind <- function(n){
  df <- data.table(kmer = NULL, genus = NULL, count = NULL)
  for(i in 1:n) {
    df <- rbindlist(list(df, list(kmer = sample(n_kmers, 1),
                                  genus = sample(n_genera, 1),
                                  count = sample(10, 1))))
  }
  df
}

df_grow_index <- function(n){
  df <- data.frame(kmer = 0, genus = 0, count = 0)
  for(i in 1:n) {
    df[i, ] <- list(kmer = sample(n_kmers, 1),
                    genus = sample(n_genera, 1),
                    count = sample(10, 1))
  }
  df
}

tbl_grow_index <- function(n){
  df <- tibble(kmer = 0, genus = 0, count = 0)
  for(i in 1:n) {
    df[i, ] <- tibble_row(kmer = sample(n_kmers, 1),
                          genus = sample(n_genera, 1),
                          count = sample(10, 1))
  }
  df
}

df_prealloc_index <- function(n){
  df <- data.frame(kmer = rep(0, n), genus = rep(0, n), count = rep(0, n))
  for(i in 1:n) {
    df[i, ] <- list(kmer = sample(n_kmers, 1),
                    genus = sample(n_genera, 1),
                    count = sample(10, 1))
  }
  df
}

tbl_prealloc_index <- function(n){
  df <- tibble(kmer = rep(0, n), genus = rep(0, n), count = rep(0, n))
  for(i in 1:n) {
    df[i, ] <- tibble_row(kmer = sample(n_kmers, 1),
                          genus = sample(n_genera, 1),
                          count = sample(10, 1))
  }
  df
}

dt_prealloc_index <- function(n){
  df <- data.table(kmer = rep(0, n), genus = rep(0, n), count = rep(0, n))
  for(i in 1:n) {
    df[i, (colnames(df)) := list(sample(n_kmers, 1),
                                 sample(n_genera, 1),
                                 sample(10, 1))]
  }
  df
}

df_predefine <- function(n) {
  kmer <- sample(n_kmers, n, replace = TRUE)
  genus <- sample(n_genera, n, replace = TRUE)
  count <- sample(10, n, replace = TRUE)
  data.frame(kmer = kmer, genus = genus, count = count)
}

tbl_predefine <- function(n) {
  kmer <- sample(n_kmers, n, replace = TRUE)
  genus <- sample(n_genera, n, replace = TRUE)
  count <- sample(10, n, replace = TRUE)
  tibble(kmer = kmer, genus = genus, count = count)
}

dt_predefine <- function(n) {
  kmer <- sample(n_kmers, n, replace = TRUE)
  genus <- sample(n_genera, n, replace = TRUE)
  count <- sample(10, n, replace = TRUE)
  data.table(kmer = kmer, genus = genus, count = count)
}

list_predefine <- function(n) {
  kmer <- sample(n_kmers, n, replace = TRUE)
  genus <- sample(n_genera, n, replace = TRUE)
  count <- sample(10, n, replace = TRUE)
  list(kmer = kmer, genus = genus, count = count)
}

vector_predefine <- function(n) {
  kmer <- sample(n_kmers, n, replace = TRUE)
  genus <- sample(n_genera, n, replace = TRUE)
  count <- sample(10, n, replace = TRUE)
}

matrix_predefine <- function(n) {
  kmer <- sample(n_kmers, n, replace = TRUE)
  genus <- sample(n_genera, n, replace = TRUE)
  count <- sample(10, n, replace = TRUE)

  m <- matrix(0, nrow = n_kmers, ncol = n_genera)

  for(i in 1:n)  {
    m[kmer[i], genus[i]] <- count[i]
  }
  m
}

matrix_sparseC <- function(n){
  kmer <- sample(n_kmers, n, replace = TRUE)
  genus <- sample(n_genera, n, replace = TRUE)
  count <- sample(10, n, replace = TRUE)

  m <- Matrix(0, nrow = n_kmers, ncol = n_genera)
  m <- as(m, "CsparseMatrix")

  for(i in 1:n)  {
    m[kmer[i], genus[i]] <- count[i]
  }
  m

}

matrix_sparseR <- function(n){
  kmer <- sample(n_kmers, n, replace = TRUE)
  genus <- sample(n_genera, n, replace = TRUE)
  count <- sample(10, n, replace = TRUE)

  m <- Matrix(0, nrow = n_kmers, ncol = n_genera)
  m <- as(m, "RsparseMatrix")

  for(i in 1:n)  {
    m[kmer[i], genus[i]] <- count[i]
  }
  m

}

matrix_sparseT <- function(n){
  kmer <- sample(n_kmers, n, replace = TRUE)
  genus <- sample(n_genera, n, replace = TRUE)
  count <- sample(10, n, replace = TRUE)

  m <- Matrix(0, nrow = n_kmers, ncol = n_genera)
  m <- as(m, "TsparseMatrix")

  for(i in 1:n)  {
    m[kmer[i], genus[i]] <- count[i]
  }
  m

}

Rcpp::sourceCpp("benchmarking_df.cpp")

n <- 100

microbenchmark(df_grow_rbind(n),
               df_grow_index(n),
               df_prealloc_index(n),
               df_predefine(n),
               tbl_grow_rbind(n),
               tbl_grow_index(n),
               tbl_prealloc_index(n),
               tbl_predefine(n),
               dt_grow_rbind(n),
               dt_prealloc_index(n),
               dt_predefine(n),
               matrix_sparseC(n),
               matrix_sparseR(n),
               matrix_sparseT(n),
               list_predefine(n),
               vector_predefine(n),
               matrix_predefine(n),
               df_rcpp(n)

               ) %>%
  group_by(expr) %>%
  summarize(median_time = median(time)) %>%
  arrange(-median_time)
