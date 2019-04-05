protein_interactions <- read_delim("stringdb_mouse/10090__protein_links.tsv.gz", delim=" ")

protein_interactions_score <- protein_interactions[, c(1, 2, 16)]

print("test")
