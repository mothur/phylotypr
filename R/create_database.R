create_database <- function(id, sequence, taxonomy, reference_data = NULL) {
  structure(
    list(
      seq_taxonomy_data = data.frame(
        id = id,
        sequence = sequence,
        taxonomy = taxonomy
      ),
      reference_data = reference_data
    ),
    class = "classification_database"
  )
}
