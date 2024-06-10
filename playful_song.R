# Playful Song Script
# Author: Brandon Polzin
# Purpose: The script is intended to run relative to the home directory where the script is located.
# Libraries: tidyverse, magrittr, DESeq2

# Load Required Libraries
library(tidyverse)
library(magrittr)
library(DESeq2)
library(clusterProfiler)

# Function Definitions

# Prepare DESeq2 Data for RRHO
# Args: deseq_data: Dataframe output from DESeq2 analysis
# Returns: A dataframe formatted for RRHO
prepare_deseq_for_rrho <- function(deseq_data) {
  deseq_data$og <- row.names(deseq_data)
  deseq_data %>%
    remove_na_from_deseq() %>%
    deseq_results_for_rrho()
}

# Filter for Matching Genes
# Args: main_data: The main dataframe
#       filter_data: The filtering dataframe
# Returns: A dataframe containing only matching genes
filter_matching_genes <- function(main_data, filter_data) {
  main_data %>%
    dplyr::filter(Gene.symbol %in% filter_data$Gene.symbol)
}

# Apply Median Split to a Variable
# Args: df: Dataframe containing the variable to be split
#       var_name: The name of the variable to be split
# Returns: Dataframe with a new column indicating categories based on median split
apply_median_split <- function(df, var_name) {
  median_value <- median(df[[var_name]], na.rm = TRUE)
  new_var_name <- paste0(var_name, "_cat")
  df[[new_var_name]] <- ifelse(df[[var_name]] > median_value, 1, 0)
  return(df)
}

# Remove NA Values from DESeq2 Data
# Args: df: Dataframe containing DESeq2 results
# Returns: Dataframe with rows containing NA values in 'pvalue' or 'padj' removed
remove_na_from_deseq <- function(df) {
  cols <- c("base_mean", "log2fold_change", "lfc_se", "stat", "pvalue", "padj")
  if (all(cols %in% colnames(df))) {
    df_no_na <- df %>%
      dplyr::filter(!is.na(padj) | !is.na(pvalue))
  } else {
    stop("Some expected DESeq2 columns are missing.")
  }
  return(df_no_na)
}

# Prepare Data for RRHO
# Args: df: Dataframe containing DESeq2 results
# Returns: Dataframe with columns formatted for RRHO
deseq_results_for_rrho <- function(df) {
  df %<>% janitor::clean_names()
  df$log10p <- log10(df$pvalue)  
  df$sign_log_fc <- sign(df$log2fold_change)  
  df$sign1 <- -df$log10p * df$sign_log_fc
  deseq2_sign1 <- df %>% 
    dplyr::select(og, sign1) %>% 
    as.data.frame()
  colnames(deseq2_sign1)[1] <- "Gene.symbol"
  return(deseq2_sign1)
}

#### STARLING ANALYSIS ####

# ---------------------
# 1. Load Data Section
# ---------------------

# Load Gene Expression Data
# Load gene expression data (RSEM) from starlings into a data frame and clean column names.
starling_rsem <- read.csv(file.path(".", "1_data", "starling_rsem.csv"), row.names = 1) %>%
  janitor::clean_names()

# Load Behavior Data
# Load behavioral attributes like sex, singing time, etc., of starlings.
starling_behav <- read.csv(file.path(".", "1_data", "starling_behavior.csv"), row.names = 1) %>%
  janitor::clean_names()

# ---------------------------
# 2. Data Preprocessing Section
# ---------------------------

# Filter by Sex
# Create subsets of the behavior data for male and female starlings.
starling_males_behav <- starling_behav %>% dplyr::filter(sex == "male")
starling_females_behav <- starling_behav %>% dplyr::filter(sex == "female")

# Apply Median Splits
# Classify starlings within each sex based on the median value of 'seconds_singing' and 'fly_to_new_perch'.
starling_males_behav <- apply_median_split(starling_males_behav, "seconds_singing")
starling_females_behav <- apply_median_split(starling_females_behav, "seconds_singing")
starling_males_behav <- apply_median_split(starling_males_behav, "fly_to_new_perch")
starling_females_behav <- apply_median_split(starling_females_behav, "fly_to_new_perch")

# Combine Male and Female Data
# Merge the male and female behavior data into a single data frame.
starling_behav <- rbind(starling_females_behav, starling_males_behav)


# ------------------------
# 3. DESeq2 Analysis Section
# ------------------------

# Data Conversion for DESeq2
# Convert the gene expression data to integers for DESeq2 analysis.
starling_genes <- row.names(starling_rsem)
starling_rsem_integer <- lapply(starling_rsem, function(x) as.integer(round(as.numeric(x)))) %>%
  as.data.frame()
row.names(starling_rsem_integer) <- starling_genes

# Create Condition Matrices
# Prepare the condition matrices for song and locomotion to be used in DESeq2 analysis.
# The matrices categorize starlings based on their behavior: high singing (hsing) vs. low singing (lsing), and high locomotion (hloco) vs. low locomotion (lloco).
# ... for song
hsing <- starling_behav %>% dplyr::filter(seconds_singing_cat == 1) %>% row.names()
lsing <- starling_behav %>% dplyr::filter(seconds_singing_cat == 0) %>% row.names()
# ... for locomotion
hloco <- starling_behav %>% dplyr::filter(fly_to_new_perch_cat == 1) %>% row.names()
lloco <- starling_behav %>% dplyr::filter(fly_to_new_perch_cat == 0) %>% row.names()

# ... Continuing DESeq2 Analysis
# Create condition matrix for song
hsing_category <- rep("hsing", length(hsing))
lsing_category <- rep("lsing", length(lsing))

starling_cond_matrix <- data.frame(
  sample_id = c(hsing, lsing),
  song = c(hsing_category, lsing_category)
)

# Relabel Condition Factors
# Relabel the factor levels for song and locomotion to set references for DESeq2 analysis.
starling_cond_matrix$song <- relevel(factor(starling_cond_matrix$song), ref = "lsing")

# Add locomotion to condition matrix
starling_cond_matrix$loco <- ifelse(starling_cond_matrix$sample_id %in% hloco, "hloco", "lloco")
starling_cond_matrix$loco <- relevel(factor(starling_cond_matrix$loco), ref = "lloco")

# Data Alignment
# Filter the gene expression data to match the samples present in the condition matrix.
starling_rsem_integer <- starling_rsem_integer[, starling_cond_matrix$sample_id]

# ------------------------
# 4. Execute DESeq2 Analysis
# ------------------------

# Run DESeq2
# Run DESeq2 to analyze differences in gene expression related to song and locomotion.
dds_starling <- DESeqDataSetFromMatrix(countData = starling_rsem_integer,
                                       colData = starling_cond_matrix,
                                       design = ~ song + loco)
dds_starling <- DESeq(dds_starling, betaPrior = FALSE)

# ----------------------
# 5. Extract Results
# ----------------------

# Song-related Gene Expression
# Extract DESeq2 results for genes that differ in expression between high singing and low singing starlings.
starling_song_res <- results(dds_starling, name = "song_hsing_vs_lsing") %>% as.data.frame()

# ------------------------
# 6. Save Results to CSV
# ------------------------

# Write DESeq2 Results
# Uncomment the following line to save the DESeq2 results for song-related gene expression to a CSV file.
write.csv(starling_song_res, "starling_deseq2_results.csv")

### RAT ANALYSIS ###

# -----------------------
# 1. Load Data
# -----------------------

# Load Gene Expression Data
# Import RSEM gene expression data for rats and clean column names.
rat_rsem <- read.csv(file.path(".", "1_data", "rat_rsem.csv"), row.names = 1) %>%
  janitor::clean_names()

# Load Behavior Data
# Import rat behavior data and clean column names.
rat_behav <- read.csv(file.path(".", "1_data", "rat_behavior.csv"), row.names = 1) %>%
  janitor::clean_names()

# ------------------------
# 2. Data Preprocessing
# ------------------------

# Apply Median Split
# Perform a median split on the variable 'total_bouts_sum_all_days'.
rat_behav <- apply_median_split(rat_behav, "total_bouts_sum_all_days")

# Remove outliers from behavior matrix
rat_behav <- rat_behav[-which(row.names(rat_behav) %in% c("x4_mpoa_sd_pair", "x12_mpoa_sd_group")),]





# ------------------------
# 3. Prepare for DESeq2
# ------------------------

# Convert Data to Integer
# Convert gene expression data to integers for DESeq2 analysis.
rat_genes <- row.names(rat_rsem)
rat_rsem_integer <- lapply(rat_rsem, function(x) as.integer(round(as.numeric(x)))) %>%
  as.data.frame()
row.names(rat_rsem_integer) <- rat_genes

# Remove outliers from rat_rsem
rat_rsem_integer <- rat_rsem_integer %>%
  dplyr::select(-x4_mpoa_sd_pair, -x12_mpoa_sd_group)

# Create Condition Matrix
# Build a condition matrix for housing and play groups.
group <- rat_behav %>% 
  dplyr::filter(housing_cat == 1) %>% 
  row.names()
pair <- rat_behav %>% 
  dplyr::filter(housing_cat == 0) %>%
  row.names()

group_category <- rep("group", length(group))
pair_category <- rep("pair", length(pair))

rat_cond_matrix <- data.frame(
  sample_id = c(group, pair),
  housing_group = c(group_category, pair_category)
)

# Relabel Condition Factors
# Relabel the factor levels for housing and play groups.
rat_cond_matrix$housing_group <- relevel(factor(rat_cond_matrix$housing_group), ref = "pair")

hplay <- rat_behav %>% dplyr::filter(total_bouts_sum_all_days_cat == 1) %>% row.names()
lplay <- rat_behav %>% dplyr::filter(total_bouts_sum_all_days_cat == 0) %>% row.names()

rat_cond_matrix$play_group <- ifelse(rat_cond_matrix$sample_id %in% hplay, "hplay", "lplay")
rat_cond_matrix$play_group <- relevel(factor(rat_cond_matrix$play_group), ref = "lplay")

# Data Alignment
# Align the gene expression data with the samples in the condition matrix.
rat_rsem_integer <- rat_rsem_integer[, rat_cond_matrix$sample_id]

# ------------------------
# 4. Execute DESeq2 Analysis
# ------------------------

# Perform DESeq2 analysis to find differentially expressed genes.
dds_rat <- DESeqDataSetFromMatrix(countData = rat_rsem_integer,
                                  colData = rat_cond_matrix,
                                  design = ~ play_group) # only doing it off of play group
dds_rat <- DESeq(dds_rat, betaPrior = FALSE)

# ----------------------
# 5. Extract Results
# ----------------------

# Extract results related to play behavior.
rat_play_res <- results(dds_rat, name = "play_group_hplay_vs_lplay") %>% as.data.frame()

# -----------------------
# 6. Add Gene Names
# -----------------------

# Import and clean gene key.
gene_key <- read.csv(file.path(".", "1_data", "gene_key.csv")) %>%
  janitor::clean_names()

# Validate and add gene names to results.
identical(gene_key$gene_id, row.names(rat_play_res))
rat_play_res$gene_name <- gene_key$symbol

# Write DESeq2 Results
# Uncomment the following line to save the DESeq2 results for song-related gene expression to a CSV file.
write.csv(rat_play_res, "rat_deseq2_results.csv")


### IDENTIFY ORTHOLOGS BETWEEN RATS AND STARLINGS ###
# ------------------------
# 1. Load Ortholog Data
# ------------------------

# Load the CSV file that maps orthologous genes between rats and starlings
ortholog_data <- read.csv(file.path(".", "1_data", "rat_starling_orthologs.csv"))

# Filter Data
# Filter the data to only include orthologous groups at the Tetrapoda level
ortholog_data <- ortholog_data %>%
  dplyr::filter(grepl("at32523", og))

# Write shared ortholog data
write.csv(ortholog_data, "shared_orthologs.csv")

# Prepare DESeq2 Data
# Clean the names and add gene names as a new column in DESeq2 results for starlings and rats
s_deseq2 <- starling_song_res %>%
  janitor::clean_names()
s_deseq2$gene_name <- row.names(s_deseq2)

r_deseq2 <- rat_play_res %>%
  janitor::clean_names()

# Reshape Ortholog Data
# Separate gene names into individual rows for easier merging
reshaped_ortholog_data <- ortholog_data %>%
  tidyr::separate_rows(gene_s_names, sep = ",") %>%
  dplyr::mutate(gene_s_names = stringr::str_trim(gene_s_names)) %>%
  tidyr::separate_rows(gene_r_names, sep = ",") %>%
  dplyr::mutate(gene_r_names = stringr::str_trim(gene_r_names))

# Merge DESeq2 Results with Ortholog Data
# Merge reshaped ortholog data with DESeq2 results to get fold change for each gene
starling_with_fc <- dplyr::left_join(reshaped_ortholog_data, s_deseq2, by = c("gene_s_names" = "gene_name"))

# Filter Highest Fold Change
# For each ortholog group, keep only the gene with the highest absolute fold change in starlings
filtered_starling <- starling_with_fc %>%
  dplyr::group_by(og) %>%
  dplyr::arrange(desc(abs(log2fold_change))) %>%
  dplyr::slice(1)

# Do the same for rat DESeq2 results
rat_with_fc <- dplyr::left_join(reshaped_ortholog_data, r_deseq2, by = c("gene_r_names" = "gene_name"))

# For each ortholog group, keep only the gene with the highest absolute fold change in rats
filtered_rat <- rat_with_fc %>%
  dplyr::group_by(og) %>%
  dplyr::arrange(abs(log2fold_change)) %>%
  dplyr::slice(1)

# Prepare Data for RRHO Analysis
# Remove NA values and prepare DESeq2 results for RRHO analysis
filtered_starling_no_na <- remove_na_from_deseq(filtered_starling)
starling_rrho_ready <- deseq_results_for_rrho(filtered_starling_no_na)

filtered_rat_no_na <- remove_na_from_deseq(filtered_rat)
rat_rrho_ready <- deseq_results_for_rrho(filtered_rat_no_na)

# Make sure that the gene lists match between starling and rat
starling_rrho_genes <- starling_rrho_ready$Gene.symbol
rat_rrho_genes <- rat_rrho_ready$Gene.symbol

starling_rrho_ready_ident <- starling_rrho_ready %>%
  dplyr::filter(Gene.symbol %in% rat_rrho_genes)

rat_rrho_ready_ident <- rat_rrho_ready %>%
  dplyr::filter(Gene.symbol %in% starling_rrho_genes)

# Perform RRHO Analysis
# Initialize and execute the Rank-Rank Hypergeometric Overlap (RRHO) analysis
rrho_obj <- RRHO2::RRHO2_initialize(starling_rrho_ready_ident, rat_rrho_ready_ident, 
                                    labels = c("Starling OGs", "Rat OGs"),
                                    method = "hyper",
                                    stepsize = 100,
                                    multipleTesting = "BY")

starling_rat_heatmap <- RRHO2::RRHO2_heatmap(rrho_obj)
starling_rat_heatmap

# --------------------------------------
# 2. Perform Male Only Starling Analysis
# -------------------------------------
# Grab sample IDs of starling males
starling_males <- row.names(starling_males_behav)

# Filter condition matrix so it's only males
starling_male_cond_matrix <- starling_cond_matrix %>%
  dplyr::filter(starling_cond_matrix$sample_id %in% starling_males)

# Filter the gene expression data to match the samples in the condition matrix
starling_male_rsem_integer <- starling_rsem_integer[, starling_male_cond_matrix$sample_id]

# DESeq2 Analysis
dds_starling_male <- DESeqDataSetFromMatrix(countData = starling_male_rsem_integer,
                                            colData = starling_male_cond_matrix,
                                            design = ~ song + loco)
dds_starling_male <- DESeq(dds_starling_male, betaPrior = FALSE)
resultsNames(dds_starling_male)

# Extract Results
starling_male_song_res <- results(dds_starling_male, name = "song_hsing_vs_lsing") %>% as.data.frame()
write.csv(starling_male_song_res, "starling_male_deseq2_results.csv")

# Starling 
s_male_deseq2 <- starling_male_song_res %>%
  janitor::clean_names()
s_male_deseq2$gene_name <- row.names(s_male_deseq2)

starling_male_with_fc <- dplyr::left_join(reshaped_ortholog_data, s_male_deseq2, by = c("gene_s_names" = "gene_name"))

# Remove outliers with high FC and low counts across gene
star_high_fc_genes <- starling_male_with_fc %>%
  dplyr::filter(log2fold_change>4 | log2fold_change< (-4)) %>%
  .$gene_s_names

star_low_count <- starling_male_rsem_integer %>%
  mutate(count_over_zero = rowSums(. > 0)) %>%
  # Filter rows where the count of cells > 0 is at least 2 and row names are in 'star_high_fc_genes'
  filter(count_over_zero <= 2, rownames(.) %in% star_high_fc_genes) %>%
  rownames(.)

print(paste("Genes removed from starlings:", star_low_count))

# Finalized data
starling_male_with_fc <- starling_male_with_fc %>%
  dplyr::filter(!gene_s_names %in% star_low_count)

# ------------------------
# 3. Perform RRHO Visualization for Male Only
# ------------------------

# Starling 
s_male_deseq2 <- starling_male_song_res %>%
  janitor::clean_names()
s_male_deseq2$gene_name <- row.names(s_male_deseq2)

# For each ortholog group, keep only the gene with the highest absolute fold change in starlings
filtered_starling_male <- starling_male_with_fc %>%
  dplyr::group_by(og) %>%
  dplyr::arrange(desc(abs(log2fold_change))) %>%
  dplyr::slice(1)

# Add fold change from rats DESeq2 results to reshaped ortholog data
rat_with_fc <- dplyr::left_join(reshaped_ortholog_data, r_deseq2, by = c("gene_r_names" = "gene_name"))
 
# Remove outliers with high FC and low counts across gene
rat_high_fc_genes <- rat_with_fc %>%
  dplyr::filter(log2fold_change>4 | log2fold_change< (-4)) %>%
  .$gene_r_names

# Convert
converted_genes <- gene_key %>%
  dplyr::filter(symbol %in% rat_high_fc_genes) %>%
  dplyr::select(gene_id) %>%
  dplyr::pull(1)

rat_low_count <- rat_rsem_integer %>%
  mutate(count_over_zero = rowSums(. > 0)) %>%
  # Filter rows where the count of cells > 0 is at least 2 and row names are in 'star_high_fc_genes'
  filter(count_over_zero <= 2, rownames(.) %in% converted_genes) %>%
  rownames(.)

new_convert <- gene_key %>%
  dplyr::filter(gene_id %in% rat_low_count) %>%
  dplyr::select(symbol) %>%
  dplyr::pull(1)

print(paste("Genes removed from rats:", rat_low_count, new_convert))

# Finalized data
rat_with_fc <- rat_with_fc %>%
  dplyr::filter(!(gene_r_names %in% new_convert & abs(log2fold_change) > 4))

# For each ortholog group, keep only the gene with the highest absolute fold change in rats
filtered_rat <- rat_with_fc %>%
   dplyr::group_by(og) %>%
   dplyr::arrange(abs(log2fold_change)) %>%
   dplyr::slice(1)

# Remove any NAs & Prepare for RRHO
filtered_starling_male_no_na <- remove_na_from_deseq(filtered_starling_male)
starling_male_rrho_ready <- deseq_results_for_rrho(filtered_starling_male_no_na)

# Make sure they match each other
starling_male_rrho_genes <- starling_male_rrho_ready$Gene.symbol

starling_male_rrho_ready_ident <- starling_male_rrho_ready %>%
  dplyr::filter(Gene.symbol %in% rat_rrho_genes)

rat_rrho_ready_ident <- rat_rrho_ready %>%
  dplyr::filter(Gene.symbol %in% starling_male_rrho_genes)


rrho_obj <- RRHO2::RRHO2_initialize(starling_male_rrho_ready_ident, rat_rrho_ready_ident, 
                                    labels = c("Starling OGs", "Rat OGs"),
                                    method = "hyper",
                                    stepsize = 100)
starling_male_rat_heatmap <- RRHO2::RRHO2_heatmap(rrho_obj)
starling_male_rat_heatmap

# ------------------------
# 4. Perform Ortholog Boostrap RRHO for Male Only
# ------------------------

upregulated_gene_list <- list()
downregulated_gene_list <- list()
up_star_down_rat_gene_list <- list()
down_star_up_rat_gene_list <- list()

filtered_starling_male_list <- list()
filtered_rat_list <- list()

if (file.exists("./bootstrapping_gene_lists.Rdata") == FALSE) {
  for (i in 1:200) {
    
    # Set a different seed for each iteration
    set.seed(i)
    
    # For each ortholog group, grab a random gene
    # starling_males
    filtered_starling_male_rand <- starling_male_with_fc %>%
      dplyr::group_by(og) %>%
      dplyr::slice_sample(n = 1)
    
    filtered_starling_male_list[[i]] <- filtered_starling_male_rand
    
    # Rats
    filtered_rat_rand <- rat_with_fc %>%
      dplyr::group_by(og) %>%
      dplyr::slice_sample(n = 1)
    
    filtered_rat_list[[i]] <- filtered_rat_rand
    
    # Prepare the data as you would
    # Remove NAs
    filtered_rat_no_na <- remove_na_from_deseq(filtered_rat_rand)
    filtered_starling_male_no_na <- remove_na_from_deseq(filtered_starling_male_rand)
    
    # Prepare for RRHO
    rat_rrho_ready <- deseq_results_for_rrho(filtered_rat_no_na)
    starling_male_rrho_ready <- deseq_results_for_rrho(filtered_starling_male_no_na)
    
    # Remove genes that don't match
    starling_male_rrho_genes <- starling_male_rrho_ready$Gene.symbol
    rat_rrho_genes <- rat_rrho_ready$Gene.symbol
    
    starling_male_rrho_ready_ident <- starling_male_rrho_ready %>%
      dplyr::filter(Gene.symbol %in% rat_rrho_genes)
    
    rat_rrho_ready_ident <- rat_rrho_ready %>%
      dplyr::filter(Gene.symbol %in% starling_male_rrho_genes)
    
    # Creat RRHO object
    rrho_obj <- RRHO2::RRHO2_initialize(starling_male_rrho_ready_ident, rat_rrho_ready_ident,
                                        labels = c("starling_male OGs", "Rat OGs"),
                                        method = "hyper",
                                        stepsize = 100,
                                        multipleTesting = "BY")
    
    # Extract genes and to list
    upregulated_genes <- rrho_obj[["genelist_uu"]][["gene_list_overlap_uu"]]
    downregulated_genes <- rrho_obj[["genelist_dd"]][["gene_list_overlap_dd"]]
    up_star_down_rat_genes <- rrho_obj[["genelist_ud"]][["gene_list_overlap_ud"]]
    down_star_up_rat_genes <- rrho_obj[["genelist_du"]][["gene_list_overlap_du"]]
    
    # Append the genes each list
    upregulated_gene_list[[i]] <- upregulated_genes
    downregulated_gene_list[[i]] <- downregulated_genes
    up_star_down_rat_gene_list[[i]] <- up_star_down_rat_genes
    down_star_up_rat_gene_list[[i]] <- down_star_up_rat_genes
    
    print(paste("Iteration number", i, "is finished."))
  }
  
  save(upregulated_gene_list, downregulated_gene_list, up_star_down_rat_gene_list, 
       down_star_up_rat_gene_list, file = "./bootstrapping_gene_lists.Rdata")
  
} else {
  
  load("./bootstrapping_gene_lists.Rdata")
}

# ------------------------
# 5. Identify Genes in 95% of Iterations
# ------------------------

print("BOOSTRAPPED")

# Initialize a list to hold Orthologous Groups (OGs) from downregulated genes.
all_ogs <- list()

# Populate the list with OGs extracted from each entry in downregulated_gene_list.
for (i in 1:length(downregulated_gene_list)) {
  all_ogs <- c(all_ogs, downregulated_gene_list[[i]])
}

# Determine the frequency of occurrence for each OG.
og_freq <- table(unlist(all_ogs))

# Define a frequency threshold for selection.
threshold <- 190 # Adjust value as per requirements.

# Select OGs that meet or exceed the frequency threshold.
selected_ogs <- names(og_freq)[og_freq >= threshold]

print(paste("There are", length(selected_ogs), "consistently downregulated genes."))

# From the filtered_starling_male dataset, select and sort genes that belong to the identified OGs based on log fold change.
final_down_genes <- filtered_starling_male %>%
  dplyr::filter(og %in% selected_ogs) %>%
  dplyr::arrange(abs(log2fold_change))

# Retain only the top 500 entries.
final_down_genes <- final_down_genes[1:500, ]

# Export the results to a CSV file
write.csv(final_down_genes, "male_sig_downregulated_top_500.csv")

# The same process is repeated for upregulated genes and specific combinations of up and down regulation in different species.

# Identify top upregulated genes.
all_ogs <- list()

for (i in 1:length(upregulated_gene_list)) {
  all_ogs <- c(all_ogs, upregulated_gene_list[[i]])
}

og_freq <- table(unlist(all_ogs))
selected_ogs <- names(og_freq)[og_freq >= threshold]

print(paste("There are", length(selected_ogs), "consistently upregulated genes."))

final_genes <- filtered_starling_male %>%
  dplyr::filter(og %in% selected_ogs) %>%
  dplyr::arrange(desc(abs(log2fold_change)))

final_genes <- final_genes[1:500, ]

# Export the results to a CSV file
write.csv(final_genes, "male_sig_upregulated_top_500.csv")

# Identify genes upregulated in high-singing starlings and downregulated in high-playing rats.
all_ogs <- list()

for (i in 1:length(up_star_down_rat_gene_list)) {
  all_ogs <- c(all_ogs, up_star_down_rat_gene_list[[i]])
}

og_freq <- table(unlist(all_ogs))
selected_ogs <- names(og_freq)[og_freq >= threshold]

print(paste("There are", length(selected_ogs), "genes upregulated in high-singing starlings and downregulated in high_playing rats."))

final_genes <- filtered_starling_male %>%
  dplyr::filter(og %in% selected_ogs) %>%
  dplyr::arrange(desc(abs(log2fold_change)))

final_genes <- final_genes[1:500, ]

# Export the results to a CSV file
write.csv(final_genes, "male_up_high_singers_down_high_players_top_500.csv")

# Identify genes downregulated in high-singing starlings and upregulated in high-playing rats.
all_ogs <- list()

for (i in 1:length(down_star_up_rat_gene_list)) {
  all_ogs <- c(all_ogs, down_star_up_rat_gene_list[[i]])
}

og_freq <- table(unlist(all_ogs))
selected_ogs <- names(og_freq)[og_freq >= threshold]

print(paste("There are", length(selected_ogs), "genes downregulated in high-singing starlings and upregulated in high_playing rats."))


final_genes <- filtered_starling_male %>%
  dplyr::filter(og %in% selected_ogs) %>%
  dplyr::arrange(desc(abs(log2fold_change)))

final_genes <- final_genes[1:500, ]

# Export the results to a CSV file
write.csv(final_genes, "male_down_high_singers_up_high_players_top_500.csv")

### MOST DIFFERENTIALLY EXPRESSED GENES ###

# Starling 
s_male_deseq2 <- starling_male_song_res %>%
  janitor::clean_names()
s_male_deseq2$gene_name <- row.names(s_male_deseq2)

# For each ortholog group, keep only the gene with the highest absolute fold change in starlings
filtered_starling_male <- starling_male_with_fc %>%
  dplyr::group_by(og) %>%
  dplyr::arrange(desc(abs(log2fold_change))) %>%
  dplyr::slice(1)

# Add fold change from rats DESeq2 results to reshaped ortholog data
rat_with_fc <- dplyr::left_join(reshaped_ortholog_data, r_deseq2, by = c("gene_r_names" = "gene_name"))

# For each ortholog group, keep only the gene with the highest absolute fold change in rats
filtered_rat <- rat_with_fc %>%
  dplyr::group_by(og) %>%
  dplyr::arrange(abs(log2fold_change)) %>%
  dplyr::slice(1)

# Remove any NAs & Prepare for RRHO
filtered_starling_male_no_na <- remove_na_from_deseq(filtered_starling_male)
starling_male_rrho_ready <- deseq_results_for_rrho(filtered_starling_male_no_na)

# Make sure they match each other
starling_male_rrho_genes <- starling_male_rrho_ready$Gene.symbol


starling_male_rrho_ready_ident <- starling_male_rrho_ready %>%
  dplyr::filter(Gene.symbol %in% rat_rrho_genes)

rat_rrho_ready_ident <- rat_rrho_ready %>%
  dplyr::filter(Gene.symbol %in% starling_male_rrho_genes)


rrho_obj <- RRHO2::RRHO2_initialize(starling_male_rrho_ready_ident, rat_rrho_ready_ident, 
                                    labels = c("Starling OGs", "Rat OGs"),
                                    method = "hyper",
                                    stepsize = 100,
                                    multipleTesting = "BY")
starling_male_rat_heatmap <- RRHO2::RRHO2_heatmap(rrho_obj)
starling_male_rat_heatmap

upregulated_genes_nbs <- rrho_obj[["genelist_uu"]][["gene_list_overlap_uu"]]
downregulated_genes_nbs <- rrho_obj[["genelist_dd"]][["gene_list_overlap_dd"]]
up_star_down_rat_genes_nbs <- rrho_obj[["genelist_ud"]][["gene_list_overlap_ud"]]
down_star_up_rat_genes_nbs <- rrho_obj[["genelist_du"]][["gene_list_overlap_du"]]

print("MOST DIFFERENTIALLY EXPRESSED")
print(paste("There are", length(upregulated_genes_nbs), "consistently upregulated genes."))
print(paste("There are", length(downregulated_genes_nbs), "consistently downregulated genes."))
print(paste("There are", length(up_star_down_rat_genes_nbs), "genes upregulated in high-singing starlings and downregulated in high_playing rats."))
print(paste("There are", length(down_star_up_rat_genes_nbs), "genes downregulated in high-singing starlings and upregulated in high_playing rats."))

## Write out the genes and OGs
ogs_combined <- rat_rrho_ready_ident$Gene.symbol

ogs_combined <- filtered_rat %>%
  dplyr::filter(og %in% ogs_combined)

write.csv(ogs_combined, "ogs_in_rrho.csv")

### GSEA ###
## RAT ##

# reading in data from deseq2
# Remove the high FC genes
r_deseq2$ensembl <- rownames(r_deseq2)

r_deseq2 <- r_deseq2 %>%
  dplyr::filter(!ensembl %in% rat_low_count)

# we want the log2 fold change 
original_gene_list <- r_deseq2$log2fold_change
# name the vector
names(original_gene_list) <- r_deseq2$gene_name %>% toupper()

# omit any NA values 
gene_list<-na.omit(original_gene_list)

# sort the list in decreasing order (required for clusterProfiler)
gene_list = sort(gene_list, decreasing = TRUE)

# GSEA
library(org.Hs.eg.db)
gse <- gseGO(geneList=gene_list, 
             ont ="ALL", 
             keyType = "SYMBOL", 
             minGSSize = 3, 
             maxGSSize = 800, 
             pvalueCutoff = 0.05,
             verbose = TRUE, 
             OrgDb = "org.Hs.eg.db")

require(DOSE)
dotplot(gse, showCategory=10, split=".sign") + facet_grid(.~.sign)
write.csv(gse, "rat_gsea_results.csv")

positive_rat_enrichments <- gse@result %>%
  dplyr::filter(enrichmentScore > 0) %>%
  dplyr::pull(3)

negative_rat_enrichments <- gse@result %>%
  dplyr::filter(enrichmentScore < 0) %>%
  dplyr::pull(3)

## STARLING ##
s_deseq2 <- s_deseq2 %>%
  dplyr::filter(!gene_name %in% star_low_count)

# we want the log2 fold change 
original_gene_list <- s_deseq2$log2fold_change
# name the vector
names(original_gene_list) <- s_deseq2$gene_name %>% toupper()

# omit any NA values 
gene_list<-na.omit(original_gene_list)

# sort the list in decreasing order (required for clusterProfiler)
gene_list = sort(gene_list, decreasing = TRUE)

# GSEA'
library(org.Hs.eg.db)
gse <- gseGO(geneList=gene_list, 
             ont ="ALL", 
             keyType = "SYMBOL",
             minGSSize = 3, 
             maxGSSize = 800, 
             pvalueCutoff = 0.05,
             verbose = TRUE, 
             OrgDb = "org.Hs.eg.db")

require(DOSE)
dotplot(gse, showCategory=10, split=".sign") + facet_grid(.~.sign)
write.csv(gse, "starling_gsea_results.csv")

positive_starling_enrichments <- gse@result %>%
  dplyr::filter(enrichmentScore > 0) %>%
  dplyr::pull(3)

negative_starling_enrichments <- gse@result %>%
  dplyr::filter(enrichmentScore < 0) %>%
  dplyr::pull(3)

# Compare results
all_positive_enrichments <- positive_starling_enrichments %in% positive_rat_enrichments
all_positive_enrichments <- positive_starling_enrichments[all_positive_enrichments]

all_negative_enrichments <- negative_starling_enrichments %in% negative_rat_enrichments
all_negative_enrichments <- negative_starling_enrichments[all_negative_enrichments]

### COMPARE STARLING DATASETS ###
# Preprocessing
# Add Orthologous Group (OG) identifiers to each dataset
s_deseq2$og <- row.names(s_deseq2)  # Add OG identifiers to general Starling dataset
s_male_deseq2$og <- row.names(s_male_deseq2)  # Add OG identifiers to male-only Starling dataset

# Prepare data for RRHO analysis by removing NA values and formatting
s_deseq2_rrho_ready <- s_deseq2 %>%
  remove_na_from_deseq() %>%
  deseq_results_for_rrho()  # Prepare general Starling dataset

s_male_deseq2_rrho_ready <- s_male_deseq2 %>%
  remove_na_from_deseq() %>%
  deseq_results_for_rrho()  # Prepare male-only Starling dataset

# Gene Matching
# Extract gene symbols for alignment
s_genes <- s_deseq2_rrho_ready$Gene.symbol  # Gene symbols from general dataset
s_male_genes <- s_male_deseq2_rrho_ready$Gene.symbol  # Gene symbols from male-only dataset

# Align datasets by filtering for matching genes
s_deseq2_rrho_ready <- prepare_deseq_for_rrho(s_deseq2)  # Prepare general dataset for comparison
s_male_deseq2_rrho_ready <- prepare_deseq_for_rrho(s_male_deseq2)  # Prepare male-only dataset for comparison

s_deseq2_rrho_ready_ident <- filter_matching_genes(s_deseq2_rrho_ready, s_male_deseq2_rrho_ready)  # Filter general dataset for matching genes
s_male_deseq2_rrho_ready_ident <- filter_matching_genes(s_male_deseq2_rrho_ready, s_deseq2_rrho_ready)  # Filter male-only dataset for matching genes

# RRHO Analysis
# Initialize the RRHO2 object for comparative analysis
rrho_obj <- RRHO2::RRHO2_initialize(
  s_male_deseq2_rrho_ready_ident, 
  s_deseq2_rrho_ready_ident,
  method = "hyper",
  labels = c("Starlings (Males)", "Starlings (All Sexes)"),
  stepsize = 100)  # Configuration for RRHO analysis; stepsize changed to 10 for detailed visualization

# Generate heatmap visualization of the RRHO analysis
RRHO2::RRHO2_heatmap(rrho_obj)

### DESeq2 PCA Plot - sex and song starling ###
# ---------------------
# 1. Load Data Section
# ---------------------

# Load Gene Expression Data
# Load gene expression data (RSEM) from starlings into a data frame and clean column names.
starling_rsem <- read.csv(file.path(".", "1_data", "starling_rsem.csv"), row.names = 1) %>%
  janitor::clean_names()

# Load Behavior Data
# Load behavioral attributes like sex, singing time, etc., of starlings.
starling_behav <- read.csv(file.path(".", "1_data", "starling_behavior.csv"), row.names = 1) %>%
  janitor::clean_names()

# ---------------------------
# 2. Data Preprocessing Section
# ---------------------------

# Filter by Sex
# Create subsets of the behavior data for male and female starlings.
starling_males_behav <- starling_behav %>% dplyr::filter(sex == "male")
starling_females_behav <- starling_behav %>% dplyr::filter(sex == "female")

# Apply Median Splits
# Classify starlings within each sex based on the median value of 'seconds_singing' and 'fly_to_new_perch'.
starling_males_behav <- apply_median_split(starling_males_behav, "seconds_singing")
starling_females_behav <- apply_median_split(starling_females_behav, "seconds_singing")
starling_males_behav <- apply_median_split(starling_males_behav, "fly_to_new_perch")
starling_females_behav <- apply_median_split(starling_females_behav, "fly_to_new_perch")

# Combine Male and Female Data
# Merge the male and female behavior data into a single data frame.
starling_behav <- rbind(starling_females_behav, starling_males_behav)
starling_behav <- starling_behav %>%
  dplyr::mutate(sex = ifelse(sex == "male", 0, 1))


# ------------------------
# 3. DESeq2 Analysis Section
# ------------------------

# Data Conversion for DESeq2
# Convert the gene expression data to integers for DESeq2 analysis.
starling_genes <- row.names(starling_rsem)
starling_rsem_integer <- lapply(starling_rsem, function(x) as.integer(round(as.numeric(x)))) %>%
  as.data.frame()
row.names(starling_rsem_integer) <- starling_genes

# Create Condition Matrices
# Prepare the condition matrices for song and locomotion to be used in DESeq2 analysis.
# The matrices categorize starlings based on their behavior: high singing (hsing) vs. low singing (lsing), and high locomotion (hloco) vs. low locomotion (lloco).
# ... for song
hsing <- starling_behav %>% dplyr::filter(seconds_singing_cat == 1) %>% row.names()
lsing <- starling_behav %>% dplyr::filter(seconds_singing_cat == 0) %>% row.names()
# ... for locomotion
female <- starling_behav %>% dplyr::filter(sex == 1) %>% row.names()
male <- starling_behav %>% dplyr::filter(sex == 0) %>% row.names()

# ... Continuing DESeq2 Analysis
# Create condition matrix for song
hsing_category <- rep("hsing", length(hsing))
lsing_category <- rep("lsing", length(lsing))

starling_cond_matrix <- data.frame(
  sample_id = c(hsing, lsing),
  song = c(hsing_category, lsing_category)
)

# Relabel Condition Factors
# Relabel the factor levels for song and locomotion to set references for DESeq2 analysis.
starling_cond_matrix$song <- relevel(factor(starling_cond_matrix$song), ref = "lsing")

# Add locomotion to condition matrix
starling_cond_matrix$sex <- ifelse(starling_cond_matrix$sample_id %in% female, "female", "male")
starling_cond_matrix$sex <- relevel(factor(starling_cond_matrix$sex), ref = "male")

# Data Alignment
# Filter the gene expression data to match the samples present in the condition matrix.
starling_rsem_integer <- starling_rsem_integer[, starling_cond_matrix$sample_id]

# ------------------------
# 4. Execute DESeq2 Analysis
# ------------------------

# Run DESeq2
# Run DESeq2 to analyze differences in gene expression related to song and locomotion.
dds_starling <- DESeqDataSetFromMatrix(countData = starling_rsem_integer,
                                       colData = starling_cond_matrix,
                                       design = ~ song + sex)
dds_starling <- DESeq(dds_starling, betaPrior = FALSE)

vst <- vst(dds_starling, blind=FALSE)  # Setting blind to FALSE considers the experimental design in transformation

# Basic PCA plot
# plotPCA(vst, intgroup=c("song", "sex")) + ggtitle("PCA of Starling Data")


### Generate Plot ###
