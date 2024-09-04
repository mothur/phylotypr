animal_legs <- data.frame(
  animal = c("cow", "fish", "chicken", "dog", "sheep"),
  n_legs = c(4, 0, 2, 4, 4)
)

animal_sounds <- data.frame(
  animal = c("cow", "chicken", "cat", "sheep", "dog"),
  sounds = c("mooo", "cluck", "meow", "baaa", "bark")
)

animal_pats_farm <- data.frame(
  animals = c("cow", "fish", "chicken", "dog", "cat", "sheep"),
  pat_has = c(T, F, T, T, T, T)
)
# base R
merge(animal_legs, animal_sounds, all = FALSE) # inner join
merge(animal_legs, animal_sounds, all = FALSE, by = "animal") # inner join
merge(animal_legs, animal_sounds, all = TRUE) # full join
merge(animal_legs, animal_sounds, all.x = TRUE) # left join
merge(animal_legs, animal_sounds, all.y = TRUE) # right join

merge(animal_legs, animal_pats_farm,
  all = FALSE,
  by.x = "animal", by.y = "animals"
) # inner join


# dplyr
inner_join(animal_legs, animal_sounds)
inner_join(animal_legs, animal_sounds, by = "animal")
inner_join(animal_legs, animal_sounds, by = join_by(animal))
full_join(animal_legs, animal_sounds, by = "animal")
full_join(animal_sounds, animal_legs, by = "animal")
left_join(animal_legs, animal_sounds, by = "animal")
right_join(animal_legs, animal_sounds, by = "animal")

# inner_join(animal_legs, animal_pats_farm) #Inner join errors
inner_join(animal_legs, animal_pats_farm, by = c("animal" = "animals"))
inner_join(animal_legs, animal_pats_farm, by = join_by(animal == animals))

inner_join(animal_legs, animal_pats_farm, by = join_by(animal == animals)) |>
  inner_join(x = _, animal_sounds, by = "animal")

# %>%
inner_join(animal_legs, animal_pats_farm, by = join_by(animal == animals)) %>%
  inner_join(., animal_sounds, by = "animal")


# data.table
animal_legs_dt <- data.table::data.table(animal_legs, key = "animal")
animal_sounds_dt <- data.table::data.table(animal_sounds, key = "animal")


# animal_legs[row, column, by, other arguments]
animal_legs_dt[animal_sounds_dt, on = .(animal)] # right join
animal_sounds_dt[animal_legs_dt, on = .(animal)] # left join
animal_sounds_dt[animal_legs_dt, nomatch = NULL, on = .(animal)] # inner join
animal_sounds_dt[animal_legs_dt, nomatch = NULL] # inner join

uniq_animals <- unique(c(animal_sounds_dt$animal, animal_legs_dt$animal))
animal_legs_dt[animal_sounds_dt[uniq_animals], on = .(animal)] # full join

animal_pats_farm_dt <- data.table::data.table(animal_pats_farm, key = "animals")

animal_legs_dt[animal_pats_farm_dt, nomatch = NULL]
animal_legs_dt[animal_pats_farm_dt, nomatch = NULL, on = "animal == animals"]




fasta_df <- read_fasta(fasta)
genera <- read_taxonomy(taxonomy)

fasta_dtA <- data.table::data.table(fasta_df, key = "id")
genera_dtA <- data.table::data.table(genera, key = "id")

microbenchmark::microbenchmark(
  base = merge(fasta_df, genera, by = "id", all = FALSE),
  ij = dplyr::inner_join(fasta_df, genera, by = "id"),
  dt = {
    fasta_dt <- data.table::data.table(fasta_df, key = "id")
    genera_dt <- data.table::data.table(genera, key = "id")
    fasta_dt[genera_dt, nomatch = NULL, on = .(id)] # inner join
  },
  dtA = fasta_dtA[genera_dtA, nomatch = NULL, on = .(id)]
)
