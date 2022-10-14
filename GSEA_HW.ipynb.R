setwd("/Users/hongb/Dropbox/github/생물정보학특론")
getwd()
# Create the data folder if it doesn't exist
if (!dir.exists("data")) {
  dir.create("data")
}

# Define the file path to the plots directory
plots_dir <- "plots" # Can replace with path to desired output plots directory

# Create the plots folder if it doesn't exist
if (!dir.exists(plots_dir)) {
  dir.create(plots_dir)
}

# Define the file path to the results directory
results_dir <- "results" # Can replace with path to desired output results directory

# Create the results folder if it doesn't exist
if (!dir.exists(results_dir)) {
  dir.create(results_dir)
}

# clusterProfiler installation
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("clusterProfiler")
"clusterProfiler" %in% installed.packages()

# msigdbr installation
if (!("msigdbr" %in% installed.packages())) {
  # Install this package if it isn't installed yet
  BiocManager::install("msigdbr", update = FALSE)
}
"msigdbr" %in% installed.packages()

# "org.Mm.eg.db" installation
if (!("org.Mm.eg.db" %in% installed.packages())) {
  # Install this package if it isn't installed yet
  BiocManager::install("org.Mm.eg.db", update = FALSE)
}

"org.Mm.eg.db" %in% installed.packages()


# Attach the library
library(clusterProfiler)

# Package that contains MSigDB gene sets in tidy format
library(msigdbr)

# Human annotation package we'll use for gene identifier conversion
library(org.Mm.eg.db)

# We will need this so we can use the pipe: %>%
library(magrittr)

# Define the url to your differential expression results file
dge_url <- "https://refinebio-examples.s3.us-east-2.amazonaws.com/03-rnaseq/results/SRP123625/SRP123625_differential_expression_results.tsv"

dge_results_file <- file.path(
  results_dir,
  "SRP123625_differential_expression_results.tsv"
)

download.file(
  dge_url,
  # The file will be saved to this location and with this name
  destfile = dge_results_file
)

# Check if the results file exists
file.exists(dge_results_file)

install.packages("tidyverse")

# Read in the contents of the differential expression results file
dge_df <- readr::read_tsv(dge_results_file)

head(dge_df)

msigdbr_species()

mm_hallmark_sets <- msigdbr(
  species = "Mus musculus", # Replace with species name relevant to your data
  category = "H"
)
head(mm_hallmark_sets)
keytypes(org.Mm.eg.db)

# First let's create a mapped data frame we can join to the differential
# expression stats
dge_mapped_df <- data.frame(
  gene_symbol = mapIds(
    # Replace with annotation package for the organism relevant to your data
    org.Mm.eg.db,
    keys = dge_df$Gene,
    # Replace with the type of gene identifiers in your data
    keytype = "ENSEMBL",
    # Replace with the type of gene identifiers you would like to map to
    column = "SYMBOL",
    # This will keep only the first mapped value for each Ensembl ID
    multiVals = "first"
  )
) %>%
  # If an Ensembl gene identifier doesn't map to a gene symbol, drop that
  # from the data frame
  dplyr::filter(!is.na(gene_symbol)) %>%
  # Make an `Ensembl` column to store the rownames
  tibble::rownames_to_column("Ensembl") %>%
  # Now let's join the rest of the expression data
  dplyr::inner_join(dge_df, by = c("Ensembl" = "Gene"))

head(dge_mapped_df)

any(duplicated(dge_mapped_df$gene_symbol))

dup_gene_symbols <- dge_mapped_df %>%
  dplyr::filter(duplicated(gene_symbol)) %>%
  dplyr::pull(gene_symbol)

dge_mapped_df %>%
  dplyr::filter(gene_symbol %in% dup_gene_symbols) %>%
  dplyr::arrange(gene_symbol)

filtered_dge_mapped_df <- dge_mapped_df %>%
  # Sort so that the highest absolute values of the log2 fold change are at the
  # top
  dplyr::arrange(dplyr::desc(abs(log2FoldChange))) %>%
  # Filter out the duplicated rows using `dplyr::distinct()`
  dplyr::distinct(gene_symbol, .keep_all = TRUE)
any(duplicated(filtered_dge_mapped_df$gene_symbol))

# Let's create a named vector ranked based on the log2 fold change values
lfc_vector <- filtered_dge_mapped_df$log2FoldChange
names(lfc_vector) <- filtered_dge_mapped_df$gene_symbol

# We need to sort the log2 fold change values in descending order here
lfc_vector <- sort(lfc_vector, decreasing = TRUE)
lfc_vector
# Run GSEA using the GSEA() function
set.seed(2020)

gsea_results <- GSEA(
  geneList = lfc_vector, # Ordered ranked gene list
  minGSSize = 25, # Minimum gene set size
  maxGSSize = 500, # Maximum gene set set
  pvalueCutoff = 0.05, # p-value cutoff
  eps = 0, # Boundary for calculating the p value
  seed = TRUE, # Set seed to make results reproducible
  pAdjustMethod = "BH", # Benjamini-Hochberg correction
  TERM2GENE = dplyr::select(
    mm_hallmark_sets,
    gs_name,
    gene_symbol
  )
)
head(gsea_results@result)

gsea_result_df <- data.frame(gsea_results@result)

gsea_result_df %>%
  # This returns the 3 rows with the largest NES values
  dplyr::slice_max(NES, n = 3)

most_positive_nes_plot <- enrichplot::gseaplot(
  gsea_results,
  geneSetID = "HALLMARK_MYC_TARGETS_V2",
  title = "HALLMARK_MYC_TARGETS_V2",
  color.line = "#0d76ff"
)
most_positive_nes_plot

ggplot2::ggsave(file.path(plots_dir, "SRP123625_gsea_enrich_positive_plot.png"),
                plot = most_positive_nes_plot
)

gsea_result_df %>%
  # Return the 3 rows with the smallest (most negative) NES values
  dplyr::slice_min(NES, n = 3)

most_negative_nes_plot <- enrichplot::gseaplot(
  gsea_results,
  geneSetID = "HALLMARK_HYPOXIA",
  title = "HALLMARK_HYPOXIA",
  color.line = "#0d76ff"
)
most_negative_nes_plot

ggplot2::ggsave(file.path(plots_dir, "SRP123625_gsea_enrich_negative_plot.png"),
                plot = most_negative_nes_plot
)

readr::write_tsv(
  gsea_result_df,
  file.path(
    results_dir,
    "SRP123625_gsea_results.tsv"
  )
)
