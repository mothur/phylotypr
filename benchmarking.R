library(microbenchmark)
library(tidyverse)
library(Rcpp)

vector_grow_c <- function(x) {
  output <- numeric()
  for(i in 1:x) {
    output <- c(output, i^2)
  }
  output
}

vector_grow_br <- function(x) {
  output <- numeric()
  for(i in 1:x) {
    output[i] <- i^2
  }
  output
}

vector_prealloc <- function(x) {
  output <- numeric(x)
  for(i in 1:x) {
    output[i] <- i^2
  }
  output
}

vector_colon <- function(x) {
  (1:x)^2
}

vector_seq <- function(x) {
  (seq(1, x, by = 1))^2
}

vector_sapply <- function(x) {
  sapply(1:x, \(i) i^2)
}

vector_lapply <- function(x) {
  sapply(1:x, \(i) i^2)
}

n <- 1e4
sourceCpp("benchmarking.cpp")

microbenchmark(vector_grow_c(n),
               vector_grow_br(n),
               vector_prealloc(n),
               vector_colon(n),
               vector_seq(n),
               vector_sapply(n),
               vector_lapply(n),
               vector_rcpp(n),
               vector_grow_c(n/100),
               vector_grow_br(n/100),
               vector_prealloc(n/100),
               vector_colon(n/100),
               vector_seq(n/100),
               vector_sapply(n/100),
               vector_lapply(n/100),
               vector_rcpp(n/100)) %>%
  group_by(expr) %>%
  summarize(median_time = median(time)) %>%
  arrange(-median_time)



x <- (1:100)^2

x[5]
x[c(5, 10)]

vector_get_one <- function(x) {

  index <- c(5, 10, 15, 20, 25, 30, 35, 40, 45, 50)

  x[5]

}

vector_get_ten <- function(x) {

  index <- c(5, 10, 15, 20, 25, 30, 35, 40, 45, 50)

  x[5]
  x[10]
  x[15]
  x[20]
  x[25]
  x[30]
  x[35]
  x[40]
  x[45]
  x[50]

}

vector_get_c <- function(x) {

  index <- c(5, 10, 15, 20, 25, 30, 35, 40, 45, 50)

  c(x[5], x[10], x[15], x[20], x[25],
    x[30], x[35], x[40], x[45], x[50])

}

vector_get_index <- function(x) {

  index <- c(5, 10, 15, 20, 25, 30, 35, 40, 45, 50)
  x[index]

}


long <- (1:1e4)^2
short <- (1:1e2)^2

microbenchmark(vector_get_one(long),
               vector_get_ten(long),
               vector_get_c(long),
               vector_get_index(long),
               vector_get_one(short),
               vector_get_ten(short),
               vector_get_c(short),
               vector_get_index(short)) %>%
  group_by(expr) %>%
  summarize(median_time = median(time)) %>%
  arrange(-median_time)
