# Create the data folder if it doesn't exist
if (!dir.exists("data")) {
  dir.create("data")
}

# Define the file path to the plots directory
plots_dir <- "plots"

# Create the plots folder if it doesn't exist
if (!dir.exists(plots_dir)) {
  dir.create(plots_dir)
}

# Define the file path to the results directory
results_dir <- "results"

# Create the results folder if it doesn't exist
if (!dir.exists(results_dir)) {
  dir.create(results_dir)
}


# Define the file path to the data directory
# Replace with the path of the folder the files will be in
data_dir <- file.path("data", "SRP140558")

# Declare the file path to the gene expression matrix file
# inside directory saved as `data_dir`
# Replace with the path to your dataset file
data_file <- file.path(data_dir, "SRP140558.tsv")

# Declare the file path to the metadata file
# inside the directory saved as `data_dir`
# Replace with the path to your metadata file
metadata_file <- file.path(data_dir, "metadata_SRP140558.tsv")

# Check if the gene expression matrix file is at the path stored in `data_file`
file.exists(data_file)

# Check if the metadata file is at the file path stored in `metadata_file`
file.exists(metadata_file)


if (!("DESeq2" %in% installed.packages())) {
  # Install this package if it isn't installed yet
  BiocManager::install("DESeq2", update = FALSE)
}

if (!("GSVA" %in% installed.packages())) {
  # Install this package if it isn't installed yet
  BiocManager::install("GSVA", update = FALSE)
}

if (!("qusage" %in% installed.packages())) {
  # Install this package if it isn't installed yet
  BiocManager::install("qusage", update = FALSE)
}

if (!("org.Hs.eg.db" %in% installed.packages())) {
  # Install this package if it isn't installed yet
  BiocManager::install("org.Hs.eg.db", update = FALSE)
}

if (!("pheatmap" %in% installed.packages())) {
  # Install pheatmap
  install.packages("pheatmap", update = FALSE)
}

# Attach the DESeq2 library
library(DESeq2)

# Attach the `qusage` library
library(qusage)

# Attach the `GSVA` library
library(GSVA)

# Human annotation package we'll use for gene identifier conversion
library(org.Hs.eg.db)

# We will need this so we can use the pipe: %>%
library(magrittr)

# Read in metadata TSV file
metadata <- readr::read_tsv(metadata_file)

# Read in data TSV file
expression_df <- readr::read_tsv(data_file) %>%
  # Here we are going to store the gene IDs as row names so that we can have a numeric matrix to perform calculations on later
  tibble::column_to_rownames("Gene")

# Make the data in the order of the metadata
expression_df <- expression_df %>%
  dplyr::select(metadata$refinebio_accession_code)

# Check if this is in the same order
all.equal(colnames(expression_df), metadata$refinebio_accession_code)

expression_df <- expression_df %>%
  # Only keep rows that have total counts above the cutoff
  dplyr::filter(rowSums(.) >= 50) %>%
  # The next DESeq2 functions need the values to be converted to integers
  round()

# 4.2.1 Prepare data for DESeq2

# Create a `DESeqDataSet` object
dds <- DESeqDataSetFromMatrix(
  countData = expression_df, # Our prepped data frame with counts
  colData = metadata, # Data frame with annotation for our samples
  design = ~1 # Here we are not specifying a model
)

# Perform DESeq2 normalization and transformation
# Normalize and transform the data in the `DESeqDataSet` object using the `vst()`
# function from the `DESEq2` R package
dds_norm <- vst(dds)

# Retrieve the normalized data from the `DESeqDataSet`
vst_df <- assay(dds_norm) %>%
  as.data.frame() %>% # Make into a data frame
  tibble::rownames_to_column("ensembl_id") # Make Gene IDs into their own column

## Import Gene Sets
hallmark_gene_sets <- msigdbr::msigdbr(
  species = "Homo sapiens", # Can change this to what species you need
  category = "H" # Only hallmark gene sets
)
head(hallmark_gene_sets)

hallmarks_list <- split(
  hallmark_gene_sets$entrez_gene, # The genes we want split into pathways
  hallmark_gene_sets$gs_name # The pathways made as the higher levels of the list
)

head(hallmarks_list, n = 2)

# 4.4.2 Gene identifier conversion
keytypes(org.Hs.eg.db)

# First let's create a mapped data frame we can join to the gene expression values
mapped_df <- data.frame(
  "entrez_id" = mapIds(
    # Replace with annotation package for the organism relevant to your data
    org.Hs.eg.db,
    keys = vst_df$ensembl_id,
    # Replace with the type of gene identifiers in your data
    keytype = "ENSEMBL",
    # Replace with the type of gene identifiers you would like to map to
    column = "ENTREZID",
    # This will keep only the first mapped value for each Ensembl ID
    multiVals = "first"
  )
) %>%
  # If an Ensembl gene identifier doesn't map to a Entrez gene identifier,
  # drop that from the data frame
  dplyr::filter(!is.na(entrez_id)) %>%
  # Make an `Ensembl` column to store the row names
  tibble::rownames_to_column("Ensembl") %>%
  # Now let's join the rest of the expression data
  dplyr::inner_join(vst_df, by = c("Ensembl" = "ensembl_id"))

head(mapped_df)
sum(duplicated(mapped_df$entrez_id))

### 4.4.2.1 Handling duplicate gene identifiers
# First let's determine the gene means
gene_means <- rowMeans(mapped_df %>% dplyr::select(-Ensembl, -entrez_id))

# Let's add this as a column in our `mapped_df`.
mapped_df <- mapped_df %>%
  # Add gene_means as a column called gene_means
  dplyr::mutate(gene_means) %>%
  # Reorder the columns so `gene_means` column is upfront
  dplyr::select(Ensembl, entrez_id, gene_means, dplyr::everything())

filtered_mapped_df <- mapped_df %>%
  # Sort so that the highest mean expression values are at the top
  dplyr::arrange(dplyr::desc(gene_means)) %>%
  # Filter out the duplicated rows using `dplyr::distinct()`
  dplyr::distinct(entrez_id, .keep_all = TRUE)

# Let’s do our check again to see if we still have duplicates.
sum(duplicated(filtered_mapped_df$entrez_id))

# Now we should prep this data so GSVA can use it.
filtered_mapped_matrix <- filtered_mapped_df %>%
  # GSVA can't the Ensembl IDs so we should drop this column as well as the means
  dplyr::select(-Ensembl, -gene_means) %>%
  # We need to store our gene identifiers as row names
  tibble::column_to_rownames("entrez_id") %>%
  # Now we can convert our object into a matrix
  as.matrix()

## 4.5 Gene Set Variation Analysis
gsva_results <- gsva(
  filtered_mapped_matrix,
  hallmarks_list,
  method = "gsva",
  # Appropriate for our vst transformed data
  kcdf = "Gaussian",
  # Minimum gene set size
  min.sz = 15,
  # Maximum gene set size
  max.sz = 500,
  # Compute Gaussian-distributed scores
  mx.diff = TRUE,
  # Don't print out the progress bar
  verbose = FALSE
)
# Print 6 rows,
head(gsva_results[, 1:10])

gsva_results %>%
  as.data.frame() %>%
  tibble::rownames_to_column("pathway") %>%
  readr::write_tsv(file.path(
    results_dir,
    "SRP140558_gsva_results.tsv"
  ))

head(metadata$refinebio_title)

annot_df <- metadata %>%
  # We need the sample IDs and the main column that contains the metadata info
  dplyr::select(
    refinebio_accession_code,
    refinebio_title
  ) %>%
  # Create our `time_point` variable based on `refinebio_title`
  dplyr::mutate(
    time_point = dplyr::case_when(
      # Create our new variable based whether the refinebio_title column
      # contains _AV_ or _CV_
      stringr::str_detect(refinebio_title, "_AV_") ~ "acute illness",
      stringr::str_detect(refinebio_title, "_CV_") ~ "recovering"
    )
  ) %>%
  # We don't need the older version of the variable anymore
  dplyr::select(-refinebio_title)

annot_df <- annot_df %>%
  # pheatmap will want our sample names that match our data to
  tibble::column_to_rownames("refinebio_accession_code")

pathway_heatmap <- pheatmap::pheatmap(gsva_results,
                                      annotation_col = annot_df, # Add metadata labels!
                                      show_colnames = FALSE, # Don't show sample labels
                                      fontsize_row = 6 # Shrink the pathway labels a tad
)

# Print out heatmap here
pathway_heatmap


length(intersect(
  hallmarks_list$HALLMARK_INTERFERON_ALPHA_RESPONSE,
  hallmarks_list$HALLMARK_INTERFERON_GAMMA_RESPONSE
))
# Replace file name with a relevant output plot name
heatmap_png_file <- file.path(plots_dir, "SRP140558_heatmap.png")

# Open a PNG file - width and height arguments control the size of the output
png(heatmap_png_file, width = 1000, height = 800)

# Print your heatmap
pathway_heatmap

# Close the PNG file:
dev.off()
