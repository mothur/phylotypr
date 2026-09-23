test_that("Accurately detect kmers from a sequence", {
  sequence <- "ATGCGCTAGTAGCATGC"
  indices <- c(
    14749, 58995, 39372, 26413, 40115, 29386, 52005, 11412, 45647,
    51514
  )
  detected <- detect_kmers_rcpp(sequence)
  expect_equal(list(sort(detected[[1]])), list(sort(indices)))


  sequence <- "ATGCGCTAGTAGCATGCN"
  indices <- c(
    14749, 58995, 39372, 26413, 40115, 29386, 52005, 11412, 45647,
    51514
  )
  detected <- detect_kmers_rcpp(sequence)
  expect_equal(list(sort(detected[[1]])), list(sort(indices)))

  sequence <- "NATGCGCTAGTAGCATGCN"
  indices <- c(
    14749, 58995, 39372, 26413, 40115, 29386, 52005, 11412, 45647,
    51514
  )
  detected <- detect_kmers_rcpp(sequence)
  expect_equal(list(sort(detected[[1]])), list(sort(indices)))

  sequence <- "NATGNGCTANGTAGCATGCN"
  indices <- c(45647, 51514)
  detected <- detect_kmers_rcpp(sequence)
  expect_equal(list(sort(detected[[1]])), list(sort(indices)))

  sequence <- "ACCGCTATCGGGGTGGGGCTAAGCCATGGGAGTNGTACGCCCCGCCGCTGCGGNNNNTGGCGGACGGCTGAGTAACACGTGGCTAACCTACCCTCGGGAGGGGGATAACACCGGGAAACTGGTGCTAATCCCCCATAGAGGAGGGGGCCTGGAAGGGTCCCTCCTCCGAAAGGGTGTGGCGGGGGTTAACGCTGCCGCACCGCCCGAGGATGGGGCTACGGCCCATTAGGTAGTTGGCGGGGTAACGGCCCGCCAAGCCGATAANGGGTAGGGGCCGTGAGAGCGGGAGCCCCCAGATGGGCACTGAGACAAGGGCCCAGGCCCTAGGGGGCGCACCAGGCGCGAAACCTCCGCAATGCGGGAAACCGTGACGGGGCCACCCCGAGTGGTGCCNTACGGGGGCNTTTCCCCGCTGTAGGAAGGCGGGGGAATAAGCGGGGGGCAAGTCNGGTGTCAGCCGCCGCGGTAATACCNNCCCCGCGAGTGGTCGGAACGATTATTGGGCCTAAAGCGCCCGTAGCCGGCTCGACAAGTCTCCTCCTAGATCCTNGGGCTCAACCCGAGGATTGGAGGAGATACTGTCGAGCTAGGGGGTGGGAGAGGCCGAGGGTACTCCCGGGGTAGGGGCGAAATCCTATAGTCCCGGGAGGACCACCAGTGGCGNNGGCGCTCGGCTGGAACACGCCCGACGGTGANGGGCGAAAGCCGGGGAAGCGAACCGGATTAGATACCCGGGTAGTCCCGGCTGTAAACGATGCGAGCTAGGTGCTATGGTGGANNTAGAGCCCACCCAGTGCCGTAGGGAANCCGTTAAGCCCGCCGCCTGGGGAGTACGGCCGNAANNCTGAAANNNANNNNNNTTGGCGGAGGNGCACCACAAGGGGTGGAGCCTGCGGTTCAATTGGAGTCAACGCCGGGAACTTACACCGGGGGNGACAGCAGGATGACGGCCAGGNTAANGACCTTGCCNGACGCGCTGAGGGGAGGTGCATGGCCGTCGCCAGCTCGTGTNGTGAGGTGTCCGGTTAAGTCCGGAAACGAGTGAGACCCCCNCCCCTAGTTGNGACCCGGGGCTACGGGTNCGGGGCACACTAGGGGGACTGCCGCCGTTTAATGCGGAGGAAGGATGGGGCCACGGCAGGTCAGCATGCCCCGAATCTCCCGGGCTACACGCGGGCTACAATGGCGGGGACAATGGGATCCGACCCCGAAAGGGGGAGGCAATCCCTCAAACCCCGCCGTAGTTGGGATCGAGGGCTGCAACTCGCCCTCGTGAANGAGGAATCCCTAGTAACCGCGTGTCAACATCGCGCGGTGAATACGTCCCTGCTCCTTGCACACACCGCCCGTCGCTCCACCCGAGGGGAGGGGGAGTGAGGCCTGGCCCTCCGGGTCGGGTCGAACT"
  detected <- detect_kmers_rcpp(sequence)
  indices <- c(343, 367, 374, 396, 399, 492, 603, 614, 683, 687, 861, 1115, 1126, 1135, 1242, 1370, 1379, 1436, 1444, 1465, 1478, 1495, 1583, 1594, 1597, 1627, 1658, 1686, 1894, 1967, 1989, 2383, 2394, 2404, 2411, 2434, 2454, 2475, 2619, 2667, 2710, 2729, 2735, 2742, 2748, 2907, 2934, 3111, 3182, 3414, 3421, 3422, 3444, 3542, 3689, 3691, 3739, 3748, 4003, 4266, 4267, 4280, 4330, 4331, 4375, 4454, 4459, 4502, 4507, 4539, 4555, 4683, 4967, 5187, 5203, 5287, 5307, 5423, 5473, 5475, 5478, 5515, 5547, 5548, 5595, 5719, 5743, 5748, 5776, 5793, 5803, 5858, 5910, 5978, 6118, 6329, 6375, 6388, 6487, 6489, 6507, 6559, 6570, 6630, 6731, 6739, 6741, 6743, 6777, 6822, 6826, 7000, 7080, 7339, 7515, 7574, 7714, 7770, 7866, 7897, 7954, 8459, 8534, 8789, 8811, 8843, 8855, 8982, 8991, 9130, 9380, 9446, 9531, 9542, 9557, 9574, 9613, 9623, 9640, 9643, 9703, 9734, 9815, 9891, 9899, 10027, 10028, 10095, 10281, 10282, 10294, 10322, 10381, 10411, 10466, 10475, 10491, 10510, 10589, 10595, 10619, 10649, 10667, 10837, 10874, 10891, 10903, 10905, 10914, 10915, 10916, 10918, 10919, 10927, 10939, 10952, 10966, 10991, 11056, 11082, 11156, 11165, 11190, 11282, 11287, 11370, 11527, 11611, 11625, 11736, 11810, 11818, 11868, 11959, 11962, 12199, 12201, 12358, 12443, 12635, 12726, 12782, 12841, 13014, 13653, 13683, 13685, 13702, 13773, 13867, 13927, 13996, 14167, 14442, 14679, 14730, 14755, 14761, 14940, 14955, 14988, 14990, 14994, 15014, 15016, 15081, 15501, 15533, 15611, 16011, 16012, 16038, 16470, 16695, 16729, 16791, 16858, 16985, 17062, 17068, 17118, 17240, 17307, 17319, 17321, 17385, 17451, 17478, 17498, 17523, 17681, 17706, 17711, 17740, 17753, 17763, 17814, 17833, 17835, 18007, 18027, 18067, 18154, 18219, 18313, 18667, 18729, 18746, 18790, 18908, 19001, 19032, 19047, 19155, 19351, 19367, 19595, 19866, 20054, 20119, 20131, 20268, 20631, 20747, 20812, 20819, 20823, 20825, 20901, 21051, 21111, 21142, 21146, 21222, 21226, 21283, 21417, 21451, 21589, 21647, 21670, 21690, 21705, 21747, 21796, 21811, 21833, 21837, 21889, 21892, 21900, 21911, 21913, 21919, 21964, 22019, 22030, 22043, 22057, 22059, 22063, 22097, 22106, 22108, 22115, 22140, 22175, 22179, 22184, 22186, 22188, 22189, 22218, 22234, 22315, 22317, 22320, 22337, 22363, 22366, 22379, 22383, 22430, 22539, 22584, 22614, 22636, 22692, 22699, 22700, 22715, 22799, 22807, 22851, 22873, 22876, 22886, 22887, 22888, 22893, 22896, 22907, 22924, 22957, 22972, 22990, 23015, 23021, 23042, 23101, 23159, 23164, 23169, 23170, 23179, 23197, 23201, 23208, 23213, 23219, 23223, 23281, 23334, 23339, 23344, 23398, 23400, 23431, 23433, 23491, 23537, 23562, 23640, 23694, 23723, 23729, 23743, 23756, 23810, 23905, 23909, 23915, 23923, 23926, 23977, 23993, 24172, 24182, 24195, 24214, 24233, 24466, 24600, 24614, 24619, 24630, 24667, 24798, 24844, 24919, 25007, 25203, 25231, 25232, 25256, 25257, 25261, 25315, 25324, 25497, 25552, 25658, 25683, 25690, 25866, 25896, 25946, 25954, 25955, 25965, 25966, 25975, 26008, 26011, 26015, 26025, 26035, 26048, 26091, 26114, 26159, 26220, 26233, 26280, 26289, 26297, 26350, 26423, 26450, 26474, 26507, 26519, 26523, 26547, 26631, 26649, 26730, 26785, 26867, 26924, 26955, 26964, 26970, 27097, 27107, 27113, 27117, 27138, 27144, 27174, 27177, 27179, 27250, 27267, 27269, 27282, 27285, 27293, 27297, 27306, 27312, 27313, 27315, 27323, 27340, 27353, 27355, 27405, 27524, 27587, 27601, 27799, 27817, 27839, 27999, 28053, 28062, 28187, 28195, 28317, 28369, 28426, 28612, 28711, 28765, 28821, 28886, 28943, 28954, 29022, 29094, 29100, 29238, 29355, 29370, 29378, 29486, 29547, 29615, 29702, 29718, 29974, 30059, 30081, 30097, 30124, 30153, 30167, 30168, 30202, 30225, 30296, 30367, 30371, 30433, 30853, 30891, 30897, 30984, 31077, 31078, 31152, 31192, 31237, 31243, 31320, 31395, 31464, 31490, 31529, 31587, 31814, 32325, 32860, 32862, 32867, 32891, 32919, 32939, 32940, 32984, 33050, 33129, 33168, 33266, 33377, 33423, 33435, 33454, 33546, 33564, 33624, 33654, 33835, 33838, 33851, 33939, 34069, 34137, 34155, 34298, 34408, 34453, 34463, 34474, 34489, 34711, 34883, 34902, 34971, 34982, 35016, 35154, 35158, 35194, 35241, 35275, 35339, 35342, 35349, 35364, 35371, 35387, 35391, 35396, 35417, 35423, 35487, 35491, 35497, 35498, 35506, 35557, 35566, 35589, 35611, 35650, 35721, 35723, 35758, 35759, 35858, 35927, 35964, 36194, 36235, 36379, 36451, 36517, 36522, 36644, 36671, 36771, 36983, 37078, 37095, 37138, 37149, 37189, 37195, 37222, 37347, 37519, 37557, 37782, 37798, 37926, 37974, 37994, 38046, 38123, 38166, 38186, 38205, 38227, 38241, 38246, 38279, 38283, 38293, 38295, 38323, 38327, 38347, 38359, 38364, 38443, 38449, 38470, 38490, 38492, 38495, 38508, 38522, 38558, 38561, 38569, 38603, 38604, 38618, 38627, 38653, 38659, 38811, 38817, 38822, 38827, 38918, 38922, 38926, 38935, 39069, 39099, 39189, 39260, 39297, 39343, 39387, 39395, 39451, 39465, 39553, 39562, 39581, 39586, 39593, 39595, 39596, 39597, 39620, 39649, 39669, 39861, 39960, 39974, 39990, 40004, 40007, 40042, 40043, 40107, 40111, 40155, 40172, 40198, 40262, 40319, 40325, 40360, 40380, 40491, 40493, 40514, 40538, 40578, 40641, 40651, 40983, 40985, 40991, 41031, 41060, 41085, 41113, 41124, 41127, 41132, 41155, 41174, 41231, 41286, 41384, 41446, 41514, 41522, 41558, 41567, 41603, 41606, 41609, 41617, 41643, 41658, 41671, 41681, 41699, 41733, 41817, 41827, 41863, 41899, 41929, 41961, 42030, 42038, 42056, 42105, 42158, 42262, 42267, 42315, 42320, 42342, 42355, 42358, 42379, 42423, 42425, 42433, 42473, 42474, 42499, 42500, 42566, 42593, 42615, 42631, 42635, 42665, 42667, 42668, 42758, 42762, 42769, 42770, 42779, 42818, 42850, 42892, 42897, 42913, 42929, 43014, 43016, 43040, 43047, 43057, 43076, 43130, 43147, 43158, 43170, 43173, 43179, 43183, 43186, 43193, 43202, 43223, 43225, 43276, 43282, 43295, 43334, 43335, 43347, 43375, 43377, 43387, 43393, 43410, 43459, 43461, 43463, 43473, 43493, 43530, 43533, 43537, 43551, 43562, 43563, 43564, 43565, 43567, 43569, 43587, 43589, 43602, 43612, 43615, 43617, 43621, 43633, 43634, 43652, 43656, 43659, 43660, 43661, 43665, 43672, 43674, 43685, 43707, 43709, 43714, 43723, 43753, 43755, 43761, 43783, 43806, 43819, 43822, 43864, 43873, 43883, 43939, 43945, 43947, 43963, 43970, 44059, 44082, 44150, 44203, 44214, 44223, 44325, 44382, 44418, 44449, 44460, 44557, 44623, 44657, 44660, 44682, 44707, 44714, 44755, 44759, 44778, 44807, 44812, 44868, 45081, 45127, 45146, 45162, 45254, 45462, 45478, 45526, 45659, 45699, 45729, 45738, 45782, 45819, 46100, 46106, 46228, 46231, 46442, 46443, 46454, 46458, 46497, 46512, 46600, 46632, 46675, 46710, 46722, 46766, 46942, 47154, 47211, 47238, 47242, 47270, 47279, 47419, 47469, 47556, 47567, 47654, 47723, 47729, 47753, 47784, 47835, 47846, 47938, 47946, 47963, 48039, 48154, 48166, 48174, 48400, 48795, 48804, 48911, 49252, 49306, 49431, 49436, 49511, 49522, 49567, 49574, 49748, 49751, 49771, 49879, 50006, 50075, 50235, 50267, 50279, 50539, 50551, 50774, 50838, 50839, 50859, 50902, 51031, 51127, 51350, 51363, 51398, 51416, 51562, 51723, 51878, 51879, 51881, 51882, 51884, 51916, 51944, 51974, 52055, 52202, 52203, 52406, 52651, 52923, 53162, 53270, 53326, 53335, 53350, 53499, 53839, 53850, 54359, 54612, 54632, 54696, 54697, 54698, 54699, 54732, 54737, 54744, 54760, 54787, 54806, 54852, 54913, 54958, 54973, 55076, 55091, 55129, 55133, 55269, 55363, 55453, 55466, 55626, 55646, 55707, 55765, 55815, 55931, 55947, 55983, 55991, 56667, 56694, 57543, 57766, 57771, 57873, 57878, 57895, 58008, 58025, 58044, 58050, 58398, 58437, 58602, 58713, 58770, 58775, 58803, 58920, 59019, 59041, 59070, 59150, 59195, 59232, 59410, 59435, 59544, 59555, 59573, 59742, 59758, 59810, 59811, 59819, 59842, 59939, 59958, 59959, 59976, 59997, 60044, 60053, 60061, 60265, 60317, 60423, 60577, 60677, 60710, 60780, 60810, 61083, 61544, 61590, 61622, 61671, 61719, 62002, 62131, 62443, 62527, 62810, 63762, 64041, 64046, 64105, 64107, 64142, 64152, 64570, 64855)
  expect_equal(list(sort(detected[[1]])), list(sort(indices)))

  sequence <- "ATGCGCTAGTAGCATGCN"
  indices <- c(
    3688, 14749, 9843, 6604, 10029, 7347, 13002, 2853, 11412, 12879,
    2362
  )
  detected <- detect_kmers_rcpp(sequence, kmer_size = 7)
  expect_equal(list(sort(detected[[1]])), list(sort(indices)))

  # make sure that duplicate kmers are scrubbed out
  sequence <- "ATGCGCTATGCGCT"
  indices <- c(13210, 9844, 10042, 7399, 6607, 14749, 3688)
  detected <- detect_kmers_rcpp(sequence, kmer_size = 7)
  expect_equal(list(sort(detected[[1]])), list(sort(indices)))
})

test_that("Accurately detect kmers across multiple sequences", {
  kmer_size <- 3
  sequences <- c("ATGCGCTA", "ATGCGCTC")
  base4_sequences <- c("03212130", "03212131")

  expected <- vector(mode = "list", length = 2)
  expected[[1]] <- c(15, 26, 29, 39, 40, 58)
  expected[[2]] <- c(15, 26, 30, 39, 40, 58)

  detect_matrix <- detect_kmers_rcpp(sequences, kmer_size) |>
    lapply(sort)

  expect_equal(detect_matrix, expected)
})

test_that("Calcuate word specific priors", {
  kmer_size <- 3
  sequences <- c("ATGCGCTA", "ATGCGCTC", "ATGCGCTC")
  detect_list <- detect_kmers_rcpp(sequences, kmer_size)

  # nolint start: commented_code_linter
  # 26 - all 3 = (3+0.5) / (1 + 3) =0.875
  # 29 - only 1 = 0.375
  # 30 - only 2 and 3 = 0.625
  # 64 - none = 0.125
  # nolint end: commented_code_linter

  priors <- calc_word_specific_priors(detect_list, kmer_size)

  expect_equal(priors[26], 0.875)
  expect_equal(priors[29], 0.375)
  expect_equal(priors[30], 0.625)
  expect_equal(priors[64], 0.125)
})


test_that("Calculate genus-specific conditional probabilities", {
  kmer_size <- 3
  sequences <- c("ATGCGCTA", "ATGCGCTC", "ATGCGCTC")
  genera <- c(1, 2, 2)

  detect_list <- detect_kmers_rcpp(sequences, kmer_size)
  priors <- calc_word_specific_priors(detect_list, kmer_size)

  # nolint start: commented_code_linter
  # (m(wi) + Pi) / (M + 1)
  # 26 - all 3 = (c(1, 2)+0.5) / (c(1, 2) + 1) = 0.74 & 0.8333333
  # 29 - only 1 = (c(1, 0)+0.5) / (c(1, 2) + 1) = 0.7500000 0.1666667
  # 30 - only 2 and 3 = (c(0, 2)+0.5) / (c(1, 2) + 1) = 0.2500000 0.8333333
  # 64 - none = 0.125 = (c(0, 2)+0.5) / (c(1, 2) + 1) = 0.2500000 0.1666667
  # nolint end

  conditional_prob <- calc_genus_conditional_prob(
    detect_list,
    genera,
    priors
  )
  # (column - 1) * nrow + row #column-wise
  expect_equal(conditional_prob[c(51, 52)], log((c(1, 2) + 0.875) / (c(1, 2) + 1)))
  expect_equal(conditional_prob[c(57, 58)], log((c(1, 0) + 0.375) / (c(1, 2) + 1)))
  expect_equal(conditional_prob[c(59, 60)], log((c(0, 2) + 0.625) / (c(1, 2) + 1)))
  expect_equal(conditional_prob[c(127, 128)], log((c(0, 0) + 0.125) / (c(1, 2) + 1)))
})

test_that("Convert back and forth between genus names and indices", {
  genera_str <- matrix(c("A", "B", "B"), ncol = 1)
  genera_index <- c(1, 2, 2)

  expect_equal(genera_str_to_index(genera_str), genera_index)
  expect_equal(get_unique_genera(genera_str), I(matrix(c("A", "B"), ncol = 1)))


  genera_str <- matrix(c("A", "B", "B", "C", "a", "b", "b", "b"), ncol = 2)
  genera_index <- c(1, 2, 2, 3)

  expect_equal(genera_str_to_index(genera_str), genera_index)
  expect_equal(get_unique_genera(genera_str),
               I(matrix(c("A", "B", "C", "a", "b", "b"), ncol = 2)))

})

test_that("Create kmer database from sequences, taxonomy, and kmer size", {
  kmer_size <- 3
  sequences <- c("ATGCGCTA", "ATGCGCTC", "ATGCGCTC")
  genera <- c("A", "B", "B")

  class_db <- create_database(
    sequence_name = c(1, 2, 3),
    sequence = sequences,
    taxonomy = genera
  )

  db <- build_kmer_database(class_db, kmer_size)

  expect_equal(
    db[["conditional_prob"]][c(51, 52)],
    log((c(1, 2) + 0.875) / (c(1, 2) + 1))
  )
  expect_equal(
    db[["conditional_prob"]][c(57, 58)],
    log((c(1, 0) + 0.375) / (c(1, 2) + 1))
  )
  expect_equal(
    db[["conditional_prob"]][c(59, 60)],
    log((c(0, 2) + 0.625) / (c(1, 2) + 1))
  )
  expect_equal(
    db[["conditional_prob"]][c(127, 128)],
    log((c(0, 0) + 0.125) / (c(1, 2) + 1))
  )

  expect_equal(db[["genera"]][[1]], "A")
  expect_equal(db[["genera"]][[2]], "B")
})

test_that("Bootstrap sample 1/kmer_size of kmers from unknowns", {
  kmers <- 1:100
  kmer_size <- 8
  expected_n_kmers <- as.integer(1 / 8 * 100)

  detected <- bootstrap_kmers(kmers, kmer_size)

  expect_length(detected, expected_n_kmers)
  expect_in(detected, kmers)
})


test_that("Classify boostrap sample works", {
  kmer_size <- 3
  sequences <- c("ATGCGCTA", "ATGCGCTC", "ATGCGCTC")
  genera <- c("A", "B", "B")

  classification_database <- create_database(
    sequence_name = 1:3,
    sequence = sequences,
    taxonomy = genera
  )
  db <- build_kmer_database(classification_database, kmer_size)
  kmers <- detect_kmers_rcpp("ATGCGCTC", kmer_size)[[1]]
  expected_classification <- 2

  detected_classification <- classify_bs_rcpp(kmers, db$conditional_prob, classification_database$seq_taxonomy_data$taxonomy |> unique() |> length())
  expect_equal(detected_classification, expected_classification)
})


test_that("Consensus classification of bootstrap subsamples", {
  db <- list()
  db[["genera"]] <- matrix(c("A", "a", "A",
                             "A", "a", "B",
                             "A", "a", "C",
                             "A", "b", "A",
                             "A", "b", "C"), nrow = 5, byrow = TRUE)

  bs_class <- list(c(1, 1, 1, 1, 4))

  expected <- list()
  expected[["level"]] <- c(1, 2, 3)
  expected[["taxonomy"]] <- c("A", "a", "A")
  expected[["confidence"]] <- c(100, 80, 80)

  observed <- consensus_bs_class(bs_class, db$genera)[[1]]

  expect_equal(observed, expected)
})


test_that("Return correct consensus taxonomy and confidence", {
  taxonomy <- matrix(c("A", "A", "A", "A", "A"), ncol = 1)
  expected <- list()
  expected[["level"]] <- 1
  expected[["taxonomy"]] <- "A"
  expected[["confidence"]] <- 100
  observed <- get_consensus(taxonomy)
  expect_equal(observed, expected)

  taxonomy <- stri_split_fixed(c("A;a", "A;a", "A;b", "A;b", "A;b"),
    pattern = ";"
  ) |> do.call(rbind, args = _)
  expected <- list()
  expected[["level"]] <- c(1, 2)
  expected[["taxonomy"]] <- c("A", "b")
  expected[["confidence"]] <- c(100, 60)
  observed <- get_consensus(taxonomy)
  expect_equal(observed, expected)

  taxonomy <- stri_split_fixed(c("A;a;A;c", "A;a;A;A_unclassified", "A;a;A;A_unclassified", "A;a;A;A_unclassified", "A;b;A;A_unclassified"),
    pattern = ";"
  ) |> do.call(rbind, args = _)
  expected <- list()
  expected[["level"]] <- c(1, 2, 3, 4)
  expected[["taxonomy"]] <- c("A", "a", "A", "A_unclassified")
  expected[["confidence"]] <- c(100, 80, 80, 60)
  observed <- get_consensus(taxonomy)
  expect_equal(observed, expected)
})
