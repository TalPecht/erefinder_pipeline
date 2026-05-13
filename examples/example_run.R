source("R/erefinder_pipeline_full.R")
genes_of_interest <- c(
  "ESR1",
  "ESR2",
  "PGR",
  "GREB1",
  "TFF1",
  "TFF3",
  "PDZK1",
  "XBP1",
  "FOXA1",
  "GATA3",
  "CCND1",
  "MYC",
  "BCL2",
  "KRT8",
  "KRT18",
  "KRT19",
  "MUC1",
  "AGR2",
  "CA12",
  "SCUBE2",
  "RARA",
  "ERBB2",
  "AURKA",
  "MKI67",
  "TOP2A",
  "BIRC5",
  "UBE2C",
  "CDC20",
  "CDK1",
  "MCM2"
) # breast-cancer / estrogen-response example

universe <- readLines("./universe_example_erefinder_pipeline.txt)

res <- run_ere_pipeline(
  genes_of_interest = genes_of_interest,
  universe = universe,
  erefinder_path = "/path/to/EREfinder",
  n_perm = 1000,
  upstream = 10000,
  downstream = 10000
)

res$summaries$genes_of_interest
res$summaries$empirical_p_anymatch
