library(microbenchmark)
library(tidyverse)
library(Rcpp)

vector_grow_append <- function(x) {
  output <- numeric()
  for (i in 1:x) {
    output <- append(output, i^2)
  }
  output
}

list_grow_append <- function(x) {
  output <- vector(mode = "list")
  for (i in 1:x) {
    output <- append(output, i^2)
  }
  output
}

vector_grow_c <- function(x) {
  output <- numeric()
  for (i in 1:x) {
    output <- c(output, i^2)
  }
  output
}

list_grow_c <- function(x) {
  output <- vector(mode = "list")
  for (i in 1:x) {
    output <- c(output, i^2)
  }
  output
}

vector_grow_br <- function(x) {
  output <- numeric()
  for (i in 1:x) {
    output[i] <- i^2
  }
  output
}

list_grow_br <- function(x) {
  output <- vector(mode = "list")
  for (i in 1:x) {
    output[i] <- i^2
  }
  output
}

vector_prealloc_sng <- function(x) {
  output <- numeric(x)
  for (i in 1:x) {
    output[i] <- i^2
  }
  output
}

vector_prealloc_dbl <- function(x) {
  output <- numeric(x)
  for (i in 1:x) {
    output[[i]] <- i^2
  }
  output
}

list_prealloc_dbl <- function(x) {
  output <- vector(mode = "list", x)
  for (i in 1:x) {
    output[[i]] <- i^2
  }
  output
}

vector_colon <- function(x) {
  (1:x)^2
}

list_colon <- function(x) {
  as.list((1:x)^2)
}

vector_seq <- function(x) {
  (seq(1, x, by = 1))^2
}

list_seq <- function(x) {
  as.list((seq(1, x, by = 1))^2)
}

vector_xapply <- function(x) {
  sapply(1:x, \(i) i^2)
}

list_xapply <- function(x) {
  lapply(1:x, \(i) i^2)
}

vector_map <- function(x) {
  map_dbl(1:x, \(i) i^2)
}

list_map <- function(x) {
  map(1:x, \(i) i^2)
}

vector_magrittr <- function(x) {
  1:x %>% (\(i) i^2)()
}

list_magrittr <- function(x) {
  1:x %>%
    (\(i) i^2)() %>%
    as.list()
}

vector_base <- function(x) {
  1:x |> (\(i) i^2)()
}

list_base <- function(x) {
  1:x |>
    (\(i) i^2)() |>
    as.list()
}


n <- 1e4
sourceCpp("benchmarking.cpp")

mb_by_n <- function(n) {
  microbenchmark(
    vector_grow_append(n),
    vector_grow_c(n),
    vector_grow_br(n),
    # vector_prealloc_sng(n),
    vector_prealloc_dbl(n),
    vector_colon(n),
    vector_seq(n),
    vector_xapply(n),
    vector_rcpp(n),
    vector_base(n),
    vector_magrittr(n),
    vector_map(n),
    list_grow_append(n),
    list_grow_c(n),
    list_grow_br(n),
    # list_prealloc_sng(n),
    list_prealloc_dbl(n),
    list_colon(n),
    list_seq(n),
    list_xapply(n),
    list_rcpp(n),
    list_base(n),
    list_magrittr(n),
    list_map(n)
  ) %>%
    group_by(expr) %>%
    summarize(median_time = median(time)) %>%
    arrange(-median_time) %>%
    mutate(n = n)
}

ns <- c(1, 10, 100, 1000, 2500, 5000, 7500, 10000, 12500, 15000)
mb_data <- map_dfr(ns, mb_by_n)

mb_data %>%
  mutate(expr = str_replace(expr, "vector_(.*)\\(n\\)", "\\1")) %>%
  ggplot(aes(x = n, y = median_time, color = expr, shape = expr)) +
  geom_line() +
  geom_point() +
  labs(x = "size of vector", y = "median time (ns)") +
  scale_color_manual(values = rep(c("tomato", "dodgerblue", "darkgray", "orange"), 3)) +
  scale_shape_manual(values = rep(c(15, 17, 19), each = 4)) +
  coord_cartesian(ylim = c(0, 5e6)) +
  theme_classic()


mb_by_n(1e4) %>%
  mutate(
    expr = str_replace(expr, "\\(n\\)", ""),
    expr = str_replace(expr, "_", "-")
  ) %>%
  separate_wider_delim(expr, names = c("structure", "f"), delim = "-") %>%
  select(-n) %>%
  pivot_wider(names_from = structure, values_from = median_time)

vector_get_one <- function(x) {
  index <- c(5, 10, 15, 20, 25, 30, 35, 40, 45, 50)

  x[5]
}

vector_get_one_named <- function(x) {
  index <- c(
    "ccx", "grs", "lmh", "bls", "gjy",
    "eee", "arn", "kwp", "ffz", "kcn"
  )

  x["lmh"]
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

vector_get_ten_named <- function(x) {
  index <- c(
    "ccx", "grs", "lmh", "bls", "gjy",
    "eee", "arn", "kwp", "ffz", "kcn"
  )

  x["ccx"]
  x["grs"]
  x["lmh"]
  x["bls"]
  x["gjy"]
  x["eee"]
  x["arn"]
  x["kwp"]
  x["ffz"]
  x["kcn"]
}


vector_get_c <- function(x) {
  index <- c(5, 10, 15, 20, 25, 30, 35, 40, 45, 50)

  c(
    x[5], x[10], x[15], x[20], x[25],
    x[30], x[35], x[40], x[45], x[50]
  )
}

vector_get_c_named <- function(x) {
  index <- c(
    "ccx", "grs", "lmh", "bls", "gjy",
    "eee", "arn", "kwp", "ffz", "kcn"
  )

  c(
    x["ccx"], x["grs"], x["lmh"], x["bls"],
    x["gjy"], x["eee"], x["arn"], x["kwp"],
    x["ffz"], x["kcn"]
  )
}

vector_get_index <- function(x) {
  index <- c(5, 10, 15, 20, 25, 30, 35, 40, 45, 50)
  x[index]
}


vector_get_index_named <- function(x) {
  index <- c(
    "ccx", "grs", "lmh", "bls", "gjy",
    "eee", "arn", "kwp", "ffz", "kcn"
  )
  x[index]
}

long_vector <- (1:1e4)^2
long_list <- (1:1e4)^2 |> as.list()

n <- expand_grid(letters, letters, letters) |>
  apply(X = _, 1, paste0, collapse = "")
long_namedlist <- long_list
names(long_namedlist) <- n[1:length(long_list)]

microbenchmark(
  vector_get_one(long_vector),
  vector_get_ten(long_vector),
  vector_get_c(long_vector),
  vector_get_index(long_vector),
  vector_get_one(long_list),
  vector_get_ten(long_list),
  vector_get_c(long_list),
  vector_get_index_named(long_list),
  vector_get_one_named(long_namedlist),
  vector_get_ten_named(long_namedlist),
  vector_get_c_named(long_namedlist),
  vector_get_index_named(long_namedlist)
) %>%
  group_by(expr) %>%
  summarize(median_time = median(time)) %>%
  arrange(-median_time) %>%
  mutate(
    expr = str_replace(expr, "\\((.*)\\)", " \\1"),
    expr = str_replace(expr, "vector_", "")
  ) %>%
  separate_wider_delim(expr, names = c("f", "data"), delim = " ") %>%
  mutate(f = str_replace(f, "_named", "")) %>%
  pivot_wider(names_from = data, values_from = median_time)
