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

n <- 1000
set.seed(19760620)
df <- df_predefine(n)
set.seed(19760620)
dt <- dt_predefine(n)
set.seed(19760620)
tbl <- tbl_predefine(n)
set.seed(19760620)
mfull <- matrix_predefine(n)
set.seed(19760620)
msparse <- matrix_sparseT(n)
set.seed(19760620)
lst <- list_predefine(n)

get_df_single <- function(){
  df[df$kmer == 20,]
}

get_df_three <- function(){
  df[df$kmer == 20 | df$kmer == 30 | df$kmer == 50,]
}

get_dt_single <- function(){
  dt[kmer == 20,]
}

get_dt_three <- function(){
  dt[kmer == 20 | kmer == 30 | kmer == 50,]
}

get_tbl_single <- function(){
  filter(tbl, kmer == 20)
}

get_tbl_three <- function(){
  filter(tbl, kmer == 20 | kmer == 30 | kmer == 50)
}

get_tbl_threeJ <- function(){
  filter_table <- tibble(kmer = c(20, 30, 50))
  inner_join(tbl, filter_table, by = "kmer")
}

get_mfull_single <- function(){
  mfull[20,]
}

get_mfull_three <- function(){
  mfull[c(20, 30, 50),]
}


get_msparse_single <- function(){
  msparse[20,]
}

get_msparse_three <- function(){
  msparse[c(20, 30, 50),]
}

kmer_vec <- lst$kmer
genus_vec <- lst$genus
count_vec <- lst$count

get_vector_single <- function() {
  kmer_present <- kmer_vec == 20
  list(kmer = kmer_vec[kmer_present],
       genus = genus_vec[kmer_present],
       count = count_vec[kmer_present])
}

get_vector_threeOR <- function() {
  kmer_present <- kmer_vec == 20 | kmer_vec == 30 | kmer_vec == 50
  list(kmer = kmer_vec[kmer_present],
       genus = genus_vec[kmer_present],
       count = count_vec[kmer_present])
}

get_vector_threeIN <- function() {
  kmer_present <- kmer_vec %in% c(20, 30, 50)

  list(kmer = kmer_vec[kmer_present],
       genus = genus_vec[kmer_present],
       count = count_vec[kmer_present])
}

get_vector_str2lang <- function() {
  kmer_present <- eval(str2lang(paste(paste("kmer_vec == ", c(20, 30, 50)),
                                      collapse = " | ")))
  kmer_vec[kmer_present]
}


get_which_single <- function(){
  kmer_present <- which(kmer_vec == 20)
  list(kmer = kmer_vec[kmer_present],
       genus = genus_vec[kmer_present],
       count = count_vec[kmer_present])
}

get_which_three <- function(){
  kmer_present <- which(kmer_vec == 20 | kmer_vec == 30 | kmer_vec == 50)
  list(kmer = kmer_vec[kmer_present],
       genus = genus_vec[kmer_present],
       count = count_vec[kmer_present])
}

microbenchmark(get_df_single(),
               get_dt_single(),
               get_tbl_single(),
               get_mfull_single(),
               get_msparse_single(),
               get_vector_single(),
               get_which_single(),

               get_df_three(),
               get_dt_three(),
               get_tbl_three(),
               get_tbl_threeJ(),
               get_mfull_three(),
               get_msparse_three(),
               get_vector_threeOR(),
               get_vector_threeIN(),
               get_which_three(),
               get_vector_str2lang()) %>%
  group_by(expr) %>%
  summarize(median_time = median(time)) %>%
  arrange(-median_time)
