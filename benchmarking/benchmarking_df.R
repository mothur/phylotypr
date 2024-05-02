library(microbenchmark)
library(tidyverse)
library(data.table)
library(Matrix)
library(duckdb)
library(duckplyr)
library(arrow)

set.seed(19760620)

n_kmers <- 4^8
n_genera <- 4000

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

dt_prealloc_index_orig <- function(n){
  df <- data.table(kmer = rep(0, n), genus = rep(0, n), count = rep(0, n))
  for(i in 1:n) {
    df[i, (colnames(df)) := list(sample(n_kmers, 1),
                                 sample(n_genera, 1),
                                 sample(10, 1))]
  }
  df
}

dt_prealloc_index_set <- function(n){
  df <- data.table(kmer = rep(0, n), genus = rep(0, n), count = rep(0, n))
  for(i in 1:n) {
    set(df, i, 1:3, list(sample(n_kmers, 1),
                         sample(n_genera, 1),
                         sample(10, 1)))
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

matrix_sparseC_old <- function(n){
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

matrix_sparseC <- function(n){
  kmer <- sample(n_kmers, n, replace = TRUE)
  genus <- sample(n_genera, n, replace = TRUE)
  count <- sample(10, n, replace = TRUE)
  sparseMatrix(i = kmer, j = genus, x = count, repr = "C")
}

matrix_sparseR_old <- function(n){
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

matrix_sparseR <- function(n){
  kmer <- sample(n_kmers, n, replace = TRUE)
  genus <- sample(n_genera, n, replace = TRUE)
  count <- sample(10, n, replace = TRUE)
  sparseMatrix(i = kmer, j = genus, x = count, repr = "R")
}


matrix_sparseT_old <- function(n){
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

matrix_sparseT <- function(n){
  kmer <- sample(n_kmers, n, replace = TRUE)
  genus <- sample(n_genera, n, replace = TRUE)
  count <- sample(10, n, replace = TRUE)
  sparseMatrix(i = kmer, j = genus, x = count, repr = "T")
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
               dt_prealloc_index_orig(n),
               dt_prealloc_index_set(n),
               dt_predefine(n),
               matrix_sparseC_old(n),
               matrix_sparseR_old(n),
               matrix_sparseT_old(n),
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

n <- 1e7
set.seed(19760620)
df <- df_predefine(n)
set.seed(19760620)
dt <- dt_predefine(n)
dt2 <- copy(dt)
setkey(dt2, kmer)

con <- dbConnect(duckdb())
dbWriteTable(con, "duck", df)

duck <- as_duckplyr_df(df)

set.seed(19760620)
tbl <- tbl_predefine(n)
set.seed(19760620)
mfull <- matrix_predefine(n)
set.seed(19760620)
msparseT <- matrix_sparseT(n)
set.seed(19760620)
msparseC <- matrix_sparseC(n)
set.seed(19760620)
msparseR <- matrix_sparseR(n)
set.seed(19760620)
lst <- list_predefine(n)
arr <- arrow_table(df)


get_df_single <- function(){
  df[which(df$kmer == 20),]
}

get_df_three <- function(){
  df[which(df$kmer == 20 | df$kmer == 30 | df$kmer == 50),]
}

get_dt_single <- function(){
  dt[which(kmer == 20),]
}

get_dt_three <- function(){
  dt[which(kmer == 20 | kmer == 30 | kmer == 50),]
}

get_dt_singlek <- function(){
  dt2[.(20)]
}

get_dt_threek <- function(){
  dt2[.(c(20, 50, 30)),]
}

get_tbl_single <- function(){
  filter(tbl, kmer == 20)
}

get_tbl_three <- function(){
  filter(tbl, kmer == 20 | kmer == 30 | kmer == 50)
}

get_duck_single <- function(){
  filter(duck, kmer == 20)
}

get_duck_three <- function(){
  filter(duck, kmer == 20 | kmer == 30 | kmer == 50)
}

get_arrow_single <- function(){
  filter(arr, kmer == 20)
}

get_arrow_three <- function(){
  filter(arr, kmer == 20 | kmer == 30 | kmer == 50)
}

get_dbi_single <- function(){
  dbGetQuery(con, "SELECT * FROM duck WHERE kmer == 20")
}

get_dbi_three <- function(){
  dbGetQuery(con, "SELECT * FROM duck WHERE kmer == 20 OR kmer == 30 OR kmer == 50")
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


get_msparseT_single <- function(){
  msparseT[20,]
}

get_msparseT_three <- function(){
  msparseT[c(20, 30, 50),]
}

get_msparseR_single <- function(){
  msparseR[20,]
}

get_msparseR_three <- function(){
  msparseR[c(20, 30, 50),]
}

get_msparseC_single <- function(){
  msparseC[20,]
}

get_msparseC_three <- function(){
  msparseC[c(20, 30, 50),]
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
               get_duck_single(),
               get_arrow_single(),
               get_dt_single(),
               get_dt_singlek(),
               get_tbl_single(),
               get_dbi_single(),
               get_mfull_single(),
               get_msparseT_single(),
               get_msparseR_single(),
               get_msparseC_single(),
               # get_vector_single(),
               get_which_single(),

               get_df_three(),
               get_dt_three(),
               get_dt_threek(),
               get_dbi_three(),
               get_tbl_three(),
               get_duck_three(),
               get_arrow_three(),
               # get_tbl_threeJ(),
               get_mfull_three(),
               get_msparseT_three(),
               get_msparseR_three(),
               get_msparseC_three(),
               #get_vector_threeOR(),
               #get_vector_threeIN(),
               get_which_three(),
               #get_vector_str2lang(),
               times = 10
               ) %>%
  group_by(expr) %>%
  summarize(median_time = median(time)) %>%
  arrange(-median_time) %>%
  print(n = Inf)

