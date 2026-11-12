# ============================================================
# BCalm MPRA Analysis - Script 2
#
# Purpose:
# 1. Load Script 1 BCalm activity results
# 2. Build analysis groups
# 3. Summarize activity relative to empirical inactive controls
# 4. Apply minimum barcode-support filter for formal WT/MUT analysis
# 5. Build WT/MUT comparison structure
# 6. Construct BCalm variant-style MPRASet objects
# 7. Run BCalm mpralm WT-vs-MUT differential analysis
# 8. Classify significant effects
# 9. Save per-cell-line and combined outputs
#
# ============================================================
# 1. CONFIGURATION
# ============================================================
INPUT_DIR <- "BCalm_activity_new_labels_troubleshooting"
# Keep this separate from any previous Script 2 analysis.
OUTPUT_DIR <- "BCalm_differential_analysis_min10_new_labels_troubleshooting"
CHECKPOINT_DIR <- file.path(
 OUTPUT_DIR,
 "checkpoints")
dir.create(
 OUTPUT_DIR,
 showWarnings = FALSE,
 recursive = TRUE)
dir.create(
 CHECKPOINT_DIR,
 showWarnings = FALSE,
 recursive = TRUE)
# Retained replicates from Script 1
CELL_LINES <- list(
 HEK = c(1, 3),
 U2OS = c(1, 2),
 U2987 = c(1, 2, 3),
 HepG2 = c(2, 3),
 `SK-MEL30` = c(1, 2, 3))
# Effect-size thresholds
WEAK_EFFECT_THRESHOLD <- 0.5
STRONG_EFFECT_THRESHOLD <- 1.0
# FDR threshold
FDR_THRESHOLD <- 0.05
# Minimum number of distinct barcodes required per oligo
#
# This threshold is applied ONLY to the formal WT/MUT
# differential analysis.
#
# Script 1 activity/QC results are not filtered.
MIN_BARCODES_PER_OLIGO <- 10
checkpoint_dir <- "BCalm_differential_analysis_min10_new_labels_troubleshooting/checkpoints"
dir.create(
 checkpoint_dir,
 recursive = TRUE,
 showWarnings = FALSE)
save_checkpoint(
 loaded_results,
 file.path(checkpoint_dir, "01_loaded_script1_workspaces.rds"))
# ============================================================
# 2. PACKAGES
# ============================================================
suppressPackageStartupMessages({
 library(data.table)
 library(dplyr)
 library(tidyr)
 library(limma)
 library(BCalm)
})
# ============================================================
# 3. CHECKPOINT FUNCTIONS
# ============================================================
save_checkpoint <- function(object, filename) {
 path <- file.path(
 CHECKPOINT_DIR,
 filename)
 saveRDS(
 object,
 path)
 cat("\nCheckpoint saved:\n")
 cat(path, "\n")
 invisible(path)}
load_checkpoint <- function(filename) {
 path <- file.path(
 CHECKPOINT_DIR,
 filename)
 if (!file.exists(path)) {
 stop(
 "\nCheckpoint does not exist:\n",
 path,
 "\n")}
 readRDS(path)}
# ============================================================
# 4. HELPER FUNCTIONS
# ============================================================
get_pair_id <- function(oligo_name) {
 oligo_name <- as.character(oligo_name)
 oligo_name <- sub(
 "^wt_",
 "",
 oligo_name)
 oligo_name <- sub(
 "^mut_",
 "",
 oligo_name)
 oligo_name}
get_allele <- function(oligo_name) {
 oligo_name <- as.character(oligo_name)
 ifelse(
 grepl("^wt_", oligo_name),
 "WT",
 ifelse(
 grepl("^mut_", oligo_name),
 "MUT",
 NA_character_))}
safe_median <- function(x) {
 x <- x[is.finite(x)]
 if (length(x) == 0) {
 return(NA_real_)}
 median(x)}
# ============================================================
# 5. LOAD SCRIPT 1 WORKSPACES
# ============================================================
cat("\n")
cat("============================================================\n")
cat("Loading Script 1 BCalm results\n")
cat("============================================================\n")
loaded_results <- list()
for (cell in names(CELL_LINES)) {
 workspace_file <- file.path(
 INPUT_DIR,
 cell,
 paste0(
 cell,
 "_BCalm_activity_workspace.RData"))
 if (!file.exists(workspace_file)) {
 stop(
 "\nMissing Script 1 workspace for ",
 cell,
 ":\n",
 workspace_file,
 "\n")}
 cat(
 "\nLoading: ",
 workspace_file,
 "\n",
 sep = "")
 env <- new.env(
 parent = emptyenv())
 load(
 workspace_file,
 envir = env)
 required_objects <- c(
 "bc",
 "dna_elem",
 "rna_elem",
 "label_vector_cell",
 "result",
 "Bc")
 missing_objects <- setdiff(
 required_objects,
 ls(envir = env))
 if (length(missing_objects) > 0) {
 stop(
 "\nWorkspace for ",
 cell,
 " is missing required objects:\n",
 paste(
 missing_objects,
 collapse = ", "),
 "\n")}
 loaded_results[[cell]] <- env
 cat(
 "Loaded ",
 nrow(env$dna_elem),
 " oligos x ",
 ncol(env$dna_elem),
 " barcode columns\n",
 sep = "")
 cat(
 "Barcode-level rows in bc: ",
 nrow(env$bc),
 "\n",
 sep = "")}
save_checkpoint(
 loaded_results,
 "01_loaded_script1_workspaces.rds")
# ============================================================
# 6. VALIDATE REPLICATE STRUCTURE
# ============================================================
cat("\n")
cat("============================================================\n")
cat("Validating replicate structure\n")
cat("============================================================\n")
replicate_info <- list()
for (cell in names(CELL_LINES)) {
 env <- loaded_results[[cell]]
 dna_cols <- colnames(env$dna_elem)
 rna_cols <- colnames(env$rna_elem)
 if (!identical(dna_cols, rna_cols)) {
 stop(
 "\nDNA/RNA column names differ for ",
 cell,
 "\n")}
 # Script 1 constructs barcode columns as:
 # sample_count_X_bc1
 # sample_count_X_bc2
 # etc.
 column_sample <- sub(
 "_bc[0-9]+$",
 "",
 dna_cols)
 observed_samples <- unique(
 column_sample)
 expected_samples <- paste0(
 "sample_count_",
 CELL_LINES[[cell]])
 if (!identical(
 observed_samples,
 expected_samples
)) {
 stop(
 "\nReplicate structure mismatch for ",
 cell,
 "\nExpected:\n",
 paste(
 expected_samples,
 collapse = ", "),
 "\nObserved:\n",
 paste(
 observed_samples,
 collapse = ", "),
 "\n")}
 columns_per_replicate <- table(
 column_sample)
 if (length(unique(columns_per_replicate)) != 1) {
 stop(
 "\nDifferent numbers of barcode slots across replicates for ",
 cell,
 "\n")}
 barcode_slots <- as.integer(
 unique(columns_per_replicate))
 n_reps <- length(
 expected_samples)
 block_vector <- rep(
 seq_len(n_reps),
 each = barcode_slots)
 if (length(block_vector) != ncol(env$dna_elem)) {
 stop(
 "\nBlock vector length mismatch for ",
 cell,
 "\n")}
 replicate_info[[cell]] <- list(
 original_replicates =
 CELL_LINES[[cell]],
 sample_names =
 expected_samples,
 n_replicates =
 n_reps,
 barcode_slots_per_replicate =
 barcode_slots,
 block_vector =
 block_vector)
 cat(
 "\n",
 cell,
 "\n",
 " retained replicates: ",
 paste(
 CELL_LINES[[cell]],
 collapse = ", "),
 "\n",
 " barcode slots/replicate: ",
 barcode_slots,
 "\n",
 " total barcode columns: ",
 ncol(env$dna_elem),
 "\n",
 sep = "")}
save_checkpoint(
 replicate_info,
 "02_replicate_structure.rds")
# ============================================================
# 7. BUILD ANALYSIS GROUPS
# ============================================================
cat("\n")
cat("============================================================\n")
cat("Building analysis groups using Script 1 pair mapping\n")
cat("============================================================\n")
# ------------------------------------------------------------
# Load the exact WT-derived pair mapping created in Script 1
# ------------------------------------------------------------
PAIR_LOOKUP_FILE <- file.path(
 INPUT_DIR,
 "pair_lookup.rds")
if (!file.exists(PAIR_LOOKUP_FILE)) {
 stop(
 "\nScript 1 pair mapping was not found:\n",
 PAIR_LOOKUP_FILE,
 "\nRun Script 1 with the pair_lookup save step first.")}
pair_lookup_script1 <- readRDS(
 PAIR_LOOKUP_FILE)
# Validate the mapping structure
required_pair_columns <- c(
 "oligo_name",
 "pair_id")
missing_pair_columns <- setdiff(
 required_pair_columns,
 names(pair_lookup_script1))
if (length(missing_pair_columns) > 0) {
 stop(
 "The Script 1 pair lookup is missing required columns: ",
 paste(missing_pair_columns, collapse = ", "))}
pair_lookup_script1 <- pair_lookup_script1 %>%
 select(
 oligo_name,
 pair_id
) %>%
 distinct()
# Each oligo must have no more than one pair_id.
# A WT can still be associated with multiple different MUTs.
duplicate_oligo_mappings <- pair_lookup_script1 %>%
 dplyr::count(
 oligo_name,
 name = "n_mappings"
) %>%
 dplyr::filter(
 n_mappings > 1)
if (nrow(duplicate_oligo_mappings) > 0) {
 print(duplicate_oligo_mappings)
 stop(
 "Some oligos have multiple rows in the Script 1 pair lookup. ",
 "Resolve these before differential analysis.")}
# ------------------------------------------------------------
# Build per-cell-line analysis groups
# ------------------------------------------------------------
analysis_groups <- list()
for (cell in names(CELL_LINES)) {
 env <- loaded_results[[cell]]
 result_df <- as.data.frame(
 env$result)
 # Keep the original Script 1 construct labels.
 # Determine WT/MUT allele from the oligo name.
 result_df <- result_df %>%
 mutate(
 cell_line = cell,
 allele = get_allele(oligo_name))
 # Attach the exact pair_id from Script 1.
 # Do not reconstruct pair IDs from oligo names here.
 result_df <- result_df %>%
 select(
 -any_of("pair_id")
) %>%
 left_join(
 pair_lookup_script1,
 by = "oligo_name")
 # ----------------------------------------------------------
 # Check the pair mapping for experimental WT/MUT constructs
 # ----------------------------------------------------------
 unmapped_test_oligos <- result_df %>%
 filter(
 label %in% c(
 "Test_WildType",
 "Test_Mutant"),
 is.na(pair_id)
) %>%
 select(
 oligo_name,
 label,
 allele,
 pair_id)
 if (nrow(unmapped_test_oligos) > 0) {
 cat(
 "\nWARNING: Experimental WT/MUT oligos without a ",
 "Script 1 pair mapping in ",
 cell,
 ":\n",
 sep = "")
 print(
 head(
 unmapped_test_oligos,
 30))}
 # ----------------------------------------------------------
 # Assign analysis groups using the CURRENT Script 1 labels
 # ----------------------------------------------------------
 result_df <- result_df %>%
 mutate(
 analysis_group = case_when(
 label %in% c(
 "Test_WildType",
 "Test_Mutant"
) ~ "Test_WT_MUT",
 label == "WT_ncncm" ~ "WT_ncncm",
 label == "MUT_ncncm" ~ "MUT_ncncm",
 label == "Expression_controls" ~
 "Expression_controls",
 label == "Positive_control" ~
 "Positive_control",
 TRUE ~ "other"))
 # ----------------------------------------------------------
 # Keep empirical inactive candidates as separate flags.
 # These flags do not replace the original construct label.
 # ----------------------------------------------------------
 result_df <- result_df %>%
 mutate(
 empirical_inactive_5pct =
 !is.na(candidate_negative_5pct) &
 candidate_negative_5pct,
 empirical_inactive_10pct =
 !is.na(candidate_negative_10pct) &
 candidate_negative_10pct)
 # ----------------------------------------------------------
 # Save per-cell-line results
 # ----------------------------------------------------------
 cell_output_dir <- file.path(
 OUTPUT_DIR,
 cell)
 dir.create(
 cell_output_dir,
 showWarnings = FALSE,
 recursive = TRUE)
 analysis_groups[[cell]] <- result_df
 write.csv(
 result_df,
 file.path(
 cell_output_dir,
 paste0(
 cell,
 "_analysis_groups.csv")),
 row.names = FALSE)
 cat(
 "\n",
 cell,
 " — analysis groups:\n",
 sep = "")
 print(
 table(
 result_df$analysis_group,
 useNA = "ifany"))
 cat(
 "Experimental WT/MUT oligos without pair IDs: ",
 nrow(unmapped_test_oligos),
 "\n",
 sep = "")}
save_checkpoint(
 analysis_groups,
 "03_analysis_groups.rds")
# ============================================================
# 8. ANALYZE ACTIVITY GROUPS
# ============================================================
cat("\n")
cat("============================================================\n")
cat("Analyzing activity groups\n")
cat("============================================================\n")
activity_results <- list()
activity_summaries <- list()
for (cell in names(CELL_LINES)) {
 df <- analysis_groups[[cell]]
 inactive_values <- df %>%
 filter(
 empirical_inactive_5pct
) %>%
 pull(logFC)
 inactive_baseline <- safe_median(
 inactive_values)
 df <- df %>%
 mutate(
 empirical_inactive_median =
 inactive_baseline,
 activity_relative_to_inactive =
 logFC -
 inactive_baseline,
 activity_weak =
 abs(activity_relative_to_inactive) >=
 WEAK_EFFECT_THRESHOLD,
 activity_strong =
 abs(activity_relative_to_inactive) >=
 STRONG_EFFECT_THRESHOLD,
 activity_direction =
 case_when(
 activity_relative_to_inactive >=
 WEAK_EFFECT_THRESHOLD ~
 "increased",
 activity_relative_to_inactive <=
 -WEAK_EFFECT_THRESHOLD ~
 "decreased",
 TRUE ~
 "near_inactive"))
 activity_results[[cell]] <- df
 summary_df <- df %>%
 group_by(
 analysis_group
) %>%
 summarise(
 n_oligos = n(),
 median_logFC =
 safe_median(logFC),
 median_activity_relative_to_inactive =
 safe_median(
 activity_relative_to_inactive),
 n_weak =
 sum(
 activity_weak,
 na.rm = TRUE),
 n_strong =
 sum(
 activity_strong,
 na.rm = TRUE),
 n_increased =
 sum(
 activity_direction ==
 "increased",
 na.rm = TRUE),
 n_decreased =
 sum(
 activity_direction ==
 "decreased",
 na.rm = TRUE),
 .groups = "drop"
) %>%
 mutate(
 cell_line = cell,
 empirical_inactive_5pct_median =
 inactive_baseline)
 activity_summaries[[cell]] <-
 summary_df
 cell_output_dir <- file.path(
 OUTPUT_DIR,
 cell)
 write.csv(
 df,
 file.path(
 cell_output_dir,
 paste0(
 cell,
 "_activity_analysis.csv")),
 row.names = FALSE)
 write.csv(
 summary_df,
 file.path(
 cell_output_dir,
 paste0(
 cell,
 "_activity_summary.csv")),
 row.names = FALSE)}
save_checkpoint(
 activity_results,
 "04_activity_analysis.rds")
# ============================================================
# 9. APPLY MINIMUM BARCODE-SUPPORT FILTER
#
# This is the critical new step.
#
# Script 1's bc object contains barcode-level observations.
# We calculate the number of DISTINCT barcodes per oligo
# directly from that object.
#
# This filter does NOT modify:
# - Script 1 workspaces
# - activity_results
# - activity_summaries
#
# It creates separate analysis matrices used only for
# WT/MUT differential analysis.
#
# Therefore Script 1 remains the complete QC/activity record.
# ============================================================
cat("\n")
cat("============================================================\n")
cat("Applying minimum barcode-support filter\n")
cat("============================================================\n")
analysis_matrices <- list()
barcode_filter_summary <- list()
for (cell in names(CELL_LINES)) {
 env <- loaded_results[[cell]]
 cat(
 "\n------------------------------------------------------------\n")
 cat(
 "Cell line: ",
 cell,
 "\n",
 sep = "")
 # ----------------------------------------------------------
 # Calculate distinct barcode support directly from Script 1
 # barcode-level object.
 # ----------------------------------------------------------
 barcode_support <- env$bc %>%
 as.data.frame() %>%
 group_by(
 oligo_name
) %>%
 summarise(
 n_barcodes =
 n_distinct(barcode),
 .groups = "drop")
 # ----------------------------------------------------------
 # Confirm barcode support is available for all matrix oligos.
 # ----------------------------------------------------------
 matrix_oligos <- rownames(
 env$dna_elem)
 missing_support <- setdiff(
 matrix_oligos,
 barcode_support$oligo_name)
 if (length(missing_support) > 0) {
 stop(
 "\nBarcode-support information is missing for ",
 length(missing_support),
 " oligos in ",
 cell,
 ".\n",
 "Example missing oligos:\n",
 paste(
 head(missing_support, 10),
 collapse = ", "),
 "\n")}
 # ----------------------------------------------------------
 # Match support counts to matrix rows.
 # ----------------------------------------------------------
 support_lookup <- barcode_support$n_barcodes[
 match(
 matrix_oligos,
 barcode_support$oligo_name)
 ]
 if (any(is.na(support_lookup))) {
 stop(
 "\nNA barcode-support values encountered for ",
 cell,
 ".\n")}
 names(support_lookup) <- matrix_oligos
 # ----------------------------------------------------------
 # Determine which oligos pass the threshold.
 # ----------------------------------------------------------
 keep_oligos <- matrix_oligos[
 support_lookup >=
 MIN_BARCODES_PER_OLIGO
 ]
 remove_oligos <- matrix_oligos[
 support_lookup <
 MIN_BARCODES_PER_OLIGO
 ]
 # ----------------------------------------------------------
 # Restrict DNA/RNA matrices to the same retained oligos.
 # ----------------------------------------------------------
 dna_filtered <- env$dna_elem[
 keep_oligos,
,
 drop = FALSE
 ]
 rna_filtered <- env$rna_elem[
 keep_oligos,
,
 drop = FALSE
 ]
 # ----------------------------------------------------------
 # Safety checks.
 # ----------------------------------------------------------
 if (!identical(
 rownames(dna_filtered),
 rownames(rna_filtered)
)) {
 stop(
 "\nFiltered DNA/RNA rownames differ for ",
 cell,
 "\n")}
 if (
 nrow(dna_filtered) !=
 length(keep_oligos)
) {
 stop(
 "\nFiltered DNA matrix row count mismatch for ",
 cell,
 "\n")}
 # ----------------------------------------------------------
 # Restrict labels to retained oligos.
 # ----------------------------------------------------------
 label_lookup <- as.character(
 env$label_vector_cell)
 names(label_lookup) <- names(
 env$label_vector_cell)
 label_filtered <- label_lookup[
 rownames(dna_filtered)
 ]
 if (any(is.na(label_filtered))) {
 stop(
 "\nMissing labels detected after barcode filtering for ",
 cell,
 "\n")}
 # ----------------------------------------------------------
 # Calculate WT/MUT-specific retention.
 #
 # This describes individual oligo retention, not complete
 # WT -> MUT comparison retention.
 # ----------------------------------------------------------
 labels_df <- data.frame(
 oligo_name =
 matrix_oligos,
 label =
 label_lookup[
 matrix_oligos
 ],
 n_barcodes =
 as.numeric(
 support_lookup),
 stringsAsFactors = FALSE)
 labels_df <- labels_df %>%
 mutate(
 allele =
 get_allele(
 oligo_name),
 passes_threshold =
 n_barcodes >=
 MIN_BARCODES_PER_OLIGO)
 allele_summary <- labels_df %>%
 filter(
 label %in%
 c(
 "Test_WildType",
 "Test_Mutant")
) %>%
 mutate(
 allele_group =
 case_when(
 label ==
 "Test_WildType" ~
 "WT",
 label ==
 "Test_Mutant" ~
 "MUT",
 TRUE ~
 NA_character_)
) %>%
 group_by(
 allele_group
) %>%
 summarise(
 total_oligos =
 n(),
 retained_oligos =
 sum(
 passes_threshold),
 removed_oligos =
 sum(
 !passes_threshold),
 percent_retained =
 100 *
 retained_oligos /
 total_oligos,
 median_n_barcodes =
 safe_median(
 n_barcodes),
 .groups = "drop")
 # ----------------------------------------------------------
 # Overall retention summary.
 # ----------------------------------------------------------
 total_oligos <- nrow(
 labels_df)
 retained_oligos <- length(
 keep_oligos)
 removed_oligos_count <- length(
 remove_oligos)
 overall_summary <- data.frame(
 cell_line =
 cell,
 minimum_barcode_threshold =
 MIN_BARCODES_PER_OLIGO,
 total_oligos =
 total_oligos,
 retained_oligos =
 retained_oligos,
 removed_oligos =
 removed_oligos_count,
 percent_retained =
 100 *
 retained_oligos /
 total_oligos,
 percent_removed =
 100 *
 removed_oligos_count /
 total_oligos,
 total_WT =
 sum(
 labels_df$label ==
 "Test_WildType"),
 retained_WT =
 sum(
 labels_df$label ==
 "Test_WildType" &
 labels_df$passes_threshold),
 total_MUT =
 sum(
 labels_df$label ==
 "Test_Mutant"),
 retained_MUT =
 sum(
 labels_df$label ==
 "Test_Mutant" &
 labels_df$passes_threshold))
 # ----------------------------------------------------------
 # Store filtered matrices and metadata.
 # ----------------------------------------------------------
 analysis_matrices[[cell]] <- list(
 dna_elem =
 dna_filtered,
 rna_elem =
 rna_filtered,
 label_vector_cell =
 label_filtered,
 barcode_support =
 support_lookup,
 retained_oligos =
 keep_oligos,
 removed_oligos =
 remove_oligos,
 n_total_oligos =
 total_oligos,
 n_retained_oligos =
 retained_oligos,
 n_removed_oligos =
 removed_oligos_count)
 barcode_filter_summary[[cell]] <- list(
 overall =
 overall_summary,
 allele =
 allele_summary,
 barcode_support =
 labels_df)
 # ----------------------------------------------------------
 # Save per-cell-line filtering tables.
 # ----------------------------------------------------------
 cell_output_dir <- file.path(
 OUTPUT_DIR,
 cell)
 write.csv(
 overall_summary,
 file.path(
 cell_output_dir,
 paste0(
 cell,
 "_barcode_filter_summary.csv")),
 row.names = FALSE)
 write.csv(
 allele_summary,
 file.path(
 cell_output_dir,
 paste0(
 cell,
 "_barcode_filter_allele_summary.csv")),
 row.names = FALSE)
 write.csv(
 labels_df,
 file.path(
 cell_output_dir,
 paste0(
 cell,
 "_barcode_support_per_oligo.csv")),
 row.names = FALSE)
 # ----------------------------------------------------------
 # Console summary.
 # ----------------------------------------------------------
 cat(
 "\n",
 cell,
 "\n",
 " minimum barcode threshold: ",
 MIN_BARCODES_PER_OLIGO,
 "\n",
 " total oligos: ",
 total_oligos,
 "\n",
 " retained oligos: ",
 retained_oligos,
 " (",
 round(
 100 *
 retained_oligos /
 total_oligos,
 1),
 "%)\n",
 " removed oligos: ",
 removed_oligos_count,
 " (",
 round(
 100 *
 removed_oligos_count /
 total_oligos,
 1),
 "%)\n",
 sep = "")
 if (nrow(allele_summary) > 0) {
 for (i in seq_len(
 nrow(allele_summary)
)) {
 cat(
 " ",
 allele_summary$allele_group[i],
 ": ",
 allele_summary$retained_oligos[i],
 "/",
 allele_summary$total_oligos[i],
 " retained (",
 round(
 allele_summary$percent_retained[i],
 1),
 "%)\n",
 sep = "")}}}
# ------------------------------------------------------------
# Combine overall filtering summaries.
# ------------------------------------------------------------
barcode_filter_summary_df <-
 bind_rows(
 lapply(
 barcode_filter_summary,
 function(x) x$overall))
barcode_filter_allele_summary_df <-
 bind_rows(
 lapply(
 names(barcode_filter_summary),
 function(cell) {
 barcode_filter_summary[[cell]]$allele %>%
 mutate(
 cell_line = cell
) %>%
 select(
 cell_line,
 everything())}))
write.csv(
 barcode_filter_summary_df,
 file.path(
 OUTPUT_DIR,
 "ALL_CELL_LINES_barcode_filter_summary.csv"),
 row.names = FALSE)
write.csv(
 barcode_filter_allele_summary_df,
 file.path(
 OUTPUT_DIR,
 "ALL_CELL_LINES_barcode_filter_allele_summary.csv"),
 row.names = FALSE)
# ------------------------------------------------------------
# Save checkpoints.
# ------------------------------------------------------------
save_checkpoint(
 analysis_matrices,
 "04b_min10_filtered_analysis_matrices.rds")
save_checkpoint(
 list(
 overall =
 barcode_filter_summary_df,
 allele =
 barcode_filter_allele_summary_df,
 per_cell =
 barcode_filter_summary),
 "04c_min10_barcode_filter_summary.rds")
# ============================================================
# 10. BUILD WT/MUT COMPARISON STRUCTURE
#
# IMPORTANT:
#
# WT/MUT comparisons are now constructed from the FILTERED
# analysis matrices.
#
# Therefore both the WT and MUT oligos in every comparison
# must have >= MIN_BARCODES_PER_OLIGO distinct barcodes.
#
# Multiple MUTs per WT remain allowed.
# Each exact WT -> MUT combination remains a separate
# comparison.
# ============================================================
cat("\n")
cat("============================================================\n")
cat("Constructing WT/MUT comparison structure\n")
cat("============================================================\n")
comparison_structures <- list()
for (cell in names(CELL_LINES)) {
 matrices <-
 analysis_matrices[[cell]]
 labels <- data.frame(
 oligo_name =
 rownames(
 matrices$dna_elem),
 label =
 as.character(
 matrices$label_vector_cell),
 stringsAsFactors = FALSE)
 labels$allele <-
 get_allele(
 labels$oligo_name)
 labels$pair_id <-
 get_pair_id(
 labels$oligo_name)
 wt <- labels %>%
 filter(
 label ==
 "Test_WildType",
 allele ==
 "WT")
 mut <- labels %>%
 filter(
 label ==
 "Test_Mutant",
 allele ==
 "MUT")
 if (nrow(wt) == 0) {
 stop(
 "\nNo Test_WildType oligos with >= ",
 MIN_BARCODES_PER_OLIGO,
 " barcodes found for ",
 cell,
 "\n")}
 if (nrow(mut) == 0) {
 stop(
 "\nNo Test_Mutant oligos with >= ",
 MIN_BARCODES_PER_OLIGO,
 " barcodes found for ",
 cell,
 "\n")}
 # IMPORTANT:
 # Multiple MUTs per WT are allowed.
 #
 # The join intentionally creates one row for every
 # WT -> MUT combination sharing a pair_id.
 #
 # Because both WT and MUT sets have already been filtered,
 # every comparison contains two oligos that individually
 # meet the minimum barcode threshold.
 comparisons <- mut %>%
 select(
 mut_oligo = oligo_name,
 pair_id
) %>%
 inner_join(
 wt %>%
 select(
 wt_oligo = oligo_name,
 pair_id),
 by = "pair_id")
 comparisons <- comparisons %>%
 mutate(
 cell_line = cell,
 comparison_id =
 paste(
 cell,
 wt_oligo,
 mut_oligo,
 sep = "__")
) %>%
 select(
 cell_line,
 comparison_id,
 pair_id,
 wt_oligo,
 mut_oligo)
 if (nrow(comparisons) == 0) {
 stop(
 "\nNo WT/MUT comparisons found after >= ",
 MIN_BARCODES_PER_OLIGO,
 "-barcode filtering for ",
 cell,
 "\n")}
 comparison_structures[[cell]] <-
 comparisons
 cell_output_dir <- file.path(
 OUTPUT_DIR,
 cell)
 write.csv(
 comparisons,
 file.path(
 cell_output_dir,
 paste0(
 cell,
 "_WT_MUT_comparison_structure.csv")),
 row.names = FALSE)
 cat(
 "\n",
 cell,
 ": ",
 nrow(comparisons),
 " WT/MUT comparisons after >= ",
 MIN_BARCODES_PER_OLIGO,
 "-barcode filtering\n",
 sep = "")
 cat(
 " unique WT oligos: ",
 n_distinct(
 comparisons$wt_oligo),
 "\n",
 sep = "")
 cat(
 " unique MUT oligos: ",
 n_distinct(
 comparisons$mut_oligo),
 "\n",
 sep = "")
 cat(
 " unique pair IDs: ",
 n_distinct(
 comparisons$pair_id),
 "\n",
 sep = "")}
save_checkpoint(
 comparison_structures,
 "05_WT_MUT_comparison_structure.rds")
# ============================================================
# 11. BUILD BCalm VARIANT-STYLE MPRASet
# ============================================================
build_variant_MPRASet <- function(
 dna_elem,
 rna_elem,
 comparisons,
 block_vector
) {
 stopifnot(
 identical(
 rownames(dna_elem),
 rownames(rna_elem)))
 n_barcode_slots_total <-
 ncol(dna_elem)
 n_reps <-
 length(
 unique(block_vector))
 if (
 n_barcode_slots_total %% n_reps != 0
) {
 stop(
 "Number of barcode columns is not divisible by number of replicates.")}
 slots_per_rep <-
 n_barcode_slots_total /
 n_reps
 n_comparisons <-
 nrow(comparisons)
 total_variant_columns <-
 n_reps *
 slots_per_rep *
 2
 dna_variant <- matrix(
 NA_real_,
 nrow = n_comparisons,
 ncol = total_variant_columns)
 rna_variant <- matrix(
 NA_real_,
 nrow = n_comparisons,
 ncol = total_variant_columns)
 variant_colnames <-
 character(
 total_variant_columns)
 # ----------------------------------------------------------
 # Create columns in this exact order:
 #
 # rep1 WT barcodes
 # rep1 MUT barcodes
 # rep2 WT barcodes
 # rep2 MUT barcodes
 # ...
 #
 # This is required so the block vector correctly identifies
 # replicate membership.
 # ----------------------------------------------------------
 col_index <- 1
 for (rep_i in seq_len(n_reps)) {
 for (allele in c("WT", "MUT")) {
 for (bc_i in seq_len(slots_per_rep)) {
 variant_colnames[col_index] <-
 paste0(
 "rep",
 rep_i,
 "_",
 allele,
 "_bc",
 bc_i)
 col_index <-
 col_index + 1}}}
 # ----------------------------------------------------------
 # Populate matrices
 # ----------------------------------------------------------
 for (rep_i in seq_len(n_reps)) {
 rep_start <-
 ((rep_i - 1) *
 slots_per_rep) + 1
 rep_end <-
 rep_i *
 slots_per_rep
 rep_cols <-
 rep_start:rep_end
 wt_start <-
 ((rep_i - 1) *
 slots_per_rep *
 2) + 1
 wt_end <-
 wt_start +
 slots_per_rep -
 1
 mut_start <-
 wt_end + 1
 mut_end <-
 mut_start +
 slots_per_rep -
 1
 for (
 comparison_i in
 seq_len(n_comparisons)
) {
 wt_name <-
 comparisons$wt_oligo[
 comparison_i
 ]
 mut_name <-
 comparisons$mut_oligo[
 comparison_i
 ]
 wt_dna <-
 as.numeric(
 dna_elem[
 wt_name,
 rep_cols
 ])
 mut_dna <-
 as.numeric(
 dna_elem[
 mut_name,
 rep_cols
 ])
 wt_rna <-
 as.numeric(
 rna_elem[
 wt_name,
 rep_cols
 ])
 mut_rna <-
 as.numeric(
 rna_elem[
 mut_name,
 rep_cols
 ])
 dna_variant[
 comparison_i,
 wt_start:wt_end
 ] <-
 wt_dna
 dna_variant[
 comparison_i,
 mut_start:mut_end
 ] <-
 mut_dna
 rna_variant[
 comparison_i,
 wt_start:wt_end
 ] <-
 wt_rna
 rna_variant[
 comparison_i,
 mut_start:mut_end
 ] <-
 mut_rna}}
 colnames(dna_variant) <-
 variant_colnames
 colnames(rna_variant) <-
 variant_colnames
 rownames(dna_variant) <-
 comparisons$comparison_id
 rownames(rna_variant) <-
 comparisons$comparison_id
 # ----------------------------------------------------------
 # Variant factor
 #
 # IMPORTANT:
 # There are WT and MUT columns within EVERY replicate.
 #
 # Therefore:
 #
 # rep(c("WT","MUT"), times=n_reps, each=slots_per_rep)
 #
 # rather than only one WT/MUT pair.
 # ----------------------------------------------------------
 variant_factor <- factor(
 rep(
 c("WT", "MUT"),
 times = n_reps,
 each = slots_per_rep),
 levels =
 c(
 "WT",
 "MUT"))
 design <-
 model.matrix(
 ~ variant_factor)
 # ----------------------------------------------------------
 # Replicate blocking
 # ----------------------------------------------------------
 variant_block <-
 rep(
 seq_len(n_reps),
 each =
 slots_per_rep * 2)
 if (
 length(variant_block) !=
 ncol(dna_variant)
) {
 stop(
 "Variant block vector does not match variant matrix columns.")}
 if (
 length(variant_factor) !=
 ncol(dna_variant)
) {
 stop(
 "Variant factor does not match variant matrix columns.")}
 if (
 nrow(design) !=
 ncol(dna_variant)
) {
 stop(
 "Design matrix rows do not match variant matrix columns.")}
 BcVariant <-
 BCalm::MPRASet(
 label =
 comparisons$comparison_id,
 DNA =
 dna_variant,
 RNA =
 rna_variant,
 eid =
 rownames(dna_variant),
 barcode =
 NULL)
 list(
 BcVariant =
 BcVariant,
 dna_variant =
 dna_variant,
 rna_variant =
 rna_variant,
 design =
 design,
 variant_factor =
 variant_factor,
 block =
 variant_block,
 slots_per_rep =
 slots_per_rep,
 n_reps =
 n_reps)}
# ============================================================
# 12. BUILD VARIANT OBJECTS
# ============================================================
cat("\n")
cat("============================================================\n")
cat("Building BCalm WT/MUT variant objects\n")
cat("============================================================\n")
variant_objects <- list()
for (cell in names(CELL_LINES)) {
 # IMPORTANT:
 # Use the filtered analysis matrices here.
 #
 # Do NOT use env$dna_elem/env$rna_elem directly because
 # those contain the pre-filter Script 1 matrices.
 matrices <-
 analysis_matrices[[cell]]
 comparisons <-
 comparison_structures[[cell]]
 block_vector <-
 replicate_info[[cell]]$block_vector
 cat(
 "\nBuilding ",
 cell,
 "...\n",
 sep = "")
 variant_objects[[cell]] <-
 build_variant_MPRASet(
 dna_elem =
 matrices$dna_elem,
 rna_elem =
 matrices$rna_elem,
 comparisons =
 comparisons,
 block_vector =
 block_vector)
 obj <-
 variant_objects[[cell]]
 cat(
 " comparisons: ",
 nrow(obj$dna_variant),
 "\n",
 sep = "")
 cat(
 " barcode columns: ",
 ncol(obj$dna_variant),
 "\n",
 sep = "")
 cat(
 " design rows: ",
 nrow(obj$design),
 "\n",
 sep = "")
 cat(
 " replicate blocks: ",
 paste(
 table(obj$block),
 collapse = ", "),
 "\n",
 sep = "")}
save_checkpoint(
 variant_objects,
 "06_WT_MUT_variant_objects.rds")
# ============================================================
# 13. VALIDATE VARIANT OBJECTS
# ============================================================
cat("\n")
cat("============================================================\n")
cat("Validating WT/MUT variant objects\n")
cat("============================================================\n")
validation <- list()
for (cell in names(variant_objects)) {
 obj <-
 variant_objects[[cell]]
 comparisons <-
 comparison_structures[[cell]]
 checks <- list(
 n_comparisons =
 nrow(obj$dna_variant),
 n_dna_columns =
 ncol(obj$dna_variant),
 n_rna_columns =
 ncol(obj$rna_variant),
 dna_rna_same_dimensions =
 identical(
 dim(obj$dna_variant),
 dim(obj$rna_variant)),
 dna_rna_same_colnames =
 identical(
 colnames(obj$dna_variant),
 colnames(obj$rna_variant)),
 rownames_match_comparisons =
 identical(
 rownames(obj$dna_variant),
 comparisons$comparison_id),
 design_rows =
 nrow(obj$design),
 design_columns =
 ncol(obj$design),
 block_length =
 length(obj$block),
 variant_factor_length =
 length(obj$variant_factor),
 matrix_columns_match_block =
 length(obj$block) ==
 ncol(obj$dna_variant),
 matrix_columns_match_factor =
 length(obj$variant_factor) ==
 ncol(obj$dna_variant),
 design_rows_match_matrix_columns =
 nrow(obj$design) ==
 ncol(obj$dna_variant),
 n_WT =
 sum(
 obj$variant_factor ==
 "WT"),
 n_MUT =
 sum(
 obj$variant_factor ==
 "MUT"),
 n_replicates =
 obj$n_reps,
 slots_per_replicate =
 obj$slots_per_rep,
 NA_DNA =
 sum(
 is.na(obj$dna_variant)),
 NA_RNA =
 sum(
 is.na(obj$rna_variant)))
 validation[[cell]] <-
 checks
 cat(
 "\n",
 cell,
 "\n",
 " comparisons: ",
 checks$n_comparisons,
 "\n",
 " columns: ",
 checks$n_dna_columns,
 "\n",
 " design rows: ",
 checks$design_rows,
 "\n",
 " replicates: ",
 checks$n_replicates,
 "\n",
 " WT columns: ",
 checks$n_WT,
 "\n",
 " MUT columns: ",
 checks$n_MUT,
 "\n",
 " NA DNA values: ",
 checks$NA_DNA,
 "\n",
 " NA RNA values: ",
 checks$NA_RNA,
 "\n",
 sep = "")}
save_checkpoint(
 validation,
 "07_variant_validation.rds")
# ------------------------------------------------------------
# Stop if structural validation fails
# ------------------------------------------------------------
for (cell in names(validation)) {
 v <-
 validation[[cell]]
 structural_checks <- c(
 v$dna_rna_same_dimensions,
 v$dna_rna_same_colnames,
 v$rownames_match_comparisons,
 v$matrix_columns_match_block,
 v$matrix_columns_match_factor,
 v$design_rows_match_matrix_columns,
 v$n_WT ==
 v$n_MUT)
 if (!all(structural_checks)) {
 stop(
 "\nStructural validation failed for ",
 cell,
 ".\n",
 "No differential BCalm model was run.\n")}}
cat(
 "\nAll structural validation checks passed.\n")
# ============================================================
# 14. RUN BCalm WT/MUT DIFFERENTIAL ANALYSIS
# ============================================================
cat("\n")
cat("============================================================\n")
cat("Running BCalm WT/MUT differential analysis\n")
cat("============================================================\n")
differential_fits <- list()
for (cell in names(CELL_LINES)) {
 obj <-
 variant_objects[[cell]]
 cat(
 "\n------------------------------------------------------------\n")
 cat(
 "Cell line: ",
 cell,
 "\n",
 sep = "")
 cat(
 "Comparisons: ",
 nrow(obj$dna_variant),
 "\n",
 sep = "")
 cat(
 "Replicates: ",
 obj$n_reps,
 "\n",
 sep = "")
 if (obj$n_reps < 2) {
 warning(
 "Skipping ",
 cell,
 " because fewer than 2 replicates are available.")
 next}
 cat(
 "Running mpralm()...\n")
 fit <-
 BCalm::mpralm(
 object =
 obj$BcVariant,
 design =
 obj$design,
 aggregate =
 "none",
 normalize =
 TRUE,
 model_type =
 "corr_groups",
 block =
 obj$block)
 differential_fits[[cell]] <-
 fit
 save_checkpoint(
 differential_fits,
 "08_WT_MUT_differential_fits_partial.rds")
 cat(
 "Fit completed for ",
 cell,
 "\n",
 sep = "")}
save_checkpoint(
 differential_fits,
 "08_WT_MUT_differential_fits.rds")
# ============================================================
# 15. EXTRACT DIFFERENTIAL RESULTS
# ============================================================
cat("\n")
cat("============================================================\n")
cat("Extracting WT/MUT differential results\n")
cat("============================================================\n")
differential_results <- list()
for (cell in names(differential_fits)) {
 fit <-
 differential_fits[[cell]]
 comparisons <-
 comparison_structures[[cell]]
 top <-
 limma::topTable(
 fit,
 coef = 2,
 number = Inf,
 sort.by = "none")
 top <-
 as.data.frame(top)
 if (
 nrow(top) !=
 nrow(comparisons)
) {
 stop(
 "\nNumber of BCalm results does not match number of WT/MUT comparisons for ",
 cell,
 "\n")}
 top$comparison_id <-
 rownames(top)
 top <-
 top %>%
 left_join(
 comparisons,
 by = "comparison_id")
 # ----------------------------------------------------------
 # Add Script 1 activity information
 #
 # IMPORTANT:
 # These activity values are descriptive information from
 # Script 1 and retain the original, unfiltered activity data.
 #
 # The differential model itself was already restricted to
 # >= MIN_BARCODES_PER_OLIGO oligos.
 # ----------------------------------------------------------
 activity_df <-
 activity_results[[cell]] %>%
 select(
 oligo_name,
 label,
 logFC,
 n_barcodes,
 replicate_mean,
 replicate_sd)
 wt_activity <-
 activity_df %>%
 dplyr::filter(
 label ==
 "Test_WildType"
) %>%
 dplyr::rename(
 wt_activity =
 logFC,
 wt_n_barcodes =
 n_barcodes,
 wt_replicate_mean =
 replicate_mean,
 wt_replicate_sd =
 replicate_sd
) %>%
 dplyr::select(
 oligo_name,
 wt_activity,
 wt_n_barcodes,
 wt_replicate_mean,
 wt_replicate_sd)
 mut_activity <-
 activity_df %>%
 dplyr::filter(
 label ==
 "Test_Mutant"
) %>%
 dplyr::rename(
 mut_activity =
 logFC,
 mut_n_barcodes =
 n_barcodes,
 mut_replicate_mean =
 replicate_mean,
 mut_replicate_sd =
 replicate_sd
) %>%
 dplyr::select(
 oligo_name,
 mut_activity,
 mut_n_barcodes,
 mut_replicate_mean,
 mut_replicate_sd)
 top <-
 top %>%
 left_join(
 wt_activity,
 by =
 c(
 wt_oligo =
 "oligo_name")
) %>%
 left_join(
 mut_activity,
 by =
 c(
 mut_oligo =
 "oligo_name"))
 # ----------------------------------------------------------
 # Differential effect classification
 #
 # Positive logFC:
 # MUT > WT
 #
 # Negative logFC:
 # MUT < WT
 # ----------------------------------------------------------
 top <-
 top %>%
 mutate(
 cell_line =
 cell,
 direction =
 case_when(
 logFC > 0 ~
 "gain",
 logFC < 0 ~
 "loss",
 TRUE ~
 "no_change"),
 effect_abs =
 abs(logFC),
 effect_tier =
 case_when(
 effect_abs >=
 STRONG_EFFECT_THRESHOLD ~
 "strong",
 effect_abs >=
 WEAK_EFFECT_THRESHOLD ~
 "weak",
 TRUE ~
 "minimal"),
 FDR_significant =
 adj.P.Val <=
 FDR_THRESHOLD,
 FDR_effect_weak =
 FDR_significant &
 effect_abs >=
 WEAK_EFFECT_THRESHOLD,
 FDR_effect_strong =
 FDR_significant &
 effect_abs >=
 STRONG_EFFECT_THRESHOLD,
 gain_FDR =
 FDR_significant &
 logFC >=
 WEAK_EFFECT_THRESHOLD,
 loss_FDR =
 FDR_significant &
 logFC <=
 -WEAK_EFFECT_THRESHOLD)
 differential_results[[cell]] <-
 top
 cell_output_dir <-
 file.path(
 OUTPUT_DIR,
 cell)
 write.csv(
 top,
 file.path(
 cell_output_dir,
 paste0(
 cell,
 "_WT_MUT_differential_results.csv")),
 row.names = FALSE)
 cat(
 "Results saved for ",
 cell,
 "\n",
 sep = "")}
save_checkpoint(
 differential_results,
 "09_WT_MUT_differential_results.rds")
# ============================================================
# 16. SAVE SIGNIFICANT SUBSETS
# ============================================================
cat("\n")
cat("============================================================\n")
cat("Saving significant WT/MUT results\n")
cat("============================================================\n")
significant_results <- list()
strong_results <- list()
gain_results <- list()
loss_results <- list()
for (cell in names(differential_results)) {
 df <-
 differential_results[[cell]]
 sig <-
 df %>%
 filter(
 FDR_significant
) %>%
 arrange(
 adj.P.Val)
 strong <-
 df %>%
 filter(
 FDR_effect_strong
) %>%
 arrange(
 adj.P.Val)
 gain <-
 df %>%
 filter(
 gain_FDR
) %>%
 arrange(
 adj.P.Val)
 loss <-
 df %>%
 filter(
 loss_FDR
) %>%
 arrange(
 adj.P.Val)
 significant_results[[cell]] <-
 sig
 strong_results[[cell]] <-
 strong
 gain_results[[cell]] <-
 gain
 loss_results[[cell]] <-
 loss
 cell_output_dir <-
 file.path(
 OUTPUT_DIR,
 cell)
 write.csv(
 sig,
 file.path(
 cell_output_dir,
 paste0(
 cell,
 "_WT_MUT_significant_FDR.csv")),
 row.names = FALSE)
 write.csv(
 strong,
 file.path(
 cell_output_dir,
 paste0(
 cell,
 "_WT_MUT_strong_effect_FDR.csv")),
 row.names = FALSE)
 write.csv(
 gain,
 file.path(
 cell_output_dir,
 paste0(
 cell,
 "_WT_MUT_gain_FDR.csv")),
 row.names = FALSE)
 write.csv(
 loss,
 file.path(
 cell_output_dir,
 paste0(
 cell,
 "_WT_MUT_loss_FDR.csv")),
 row.names = FALSE)
 cat(
 "\n",
 cell,
 "\n",
 " FDR < ",
 FDR_THRESHOLD,
 ": ",
 nrow(sig),
 "\n",
 " FDR + |effect| >= ",
 STRONG_EFFECT_THRESHOLD,
 ": ",
 nrow(strong),
 "\n",
 " FDR + gain >= ",
 WEAK_EFFECT_THRESHOLD,
 ": ",
 nrow(gain),
 "\n",
 " FDR + loss <= -",
 WEAK_EFFECT_THRESHOLD,
 ": ",
 nrow(loss),
 "\n",
 sep = "")}
save_checkpoint(
 list(
 significant =
 significant_results,
 strong =
 strong_results,
 gain =
 gain_results,
 loss =
 loss_results),
 "10_significant_subsets.rds")
# ============================================================
# 17. DIFFERENTIAL SUMMARY
# ============================================================
cat("\n")
cat("============================================================\n")
cat("Creating differential summary\n")
cat("============================================================\n")
differential_summary <- list()
for (cell in names(differential_results)) {
 df <-
 differential_results[[cell]]
 summary_df <-
 data.frame(
 cell_line =
 cell,
 n_WT_MUT_comparisons =
 nrow(df),
 n_FDR_significant =
 sum(
 df$FDR_significant,
 na.rm = TRUE),
 n_FDR_effect_weak =
 sum(
 df$FDR_effect_weak,
 na.rm = TRUE),
 n_FDR_effect_strong =
 sum(
 df$FDR_effect_strong,
 na.rm = TRUE),
 n_gain_FDR =
 sum(
 df$gain_FDR,
 na.rm = TRUE),
 n_loss_FDR =
 sum(
 df$loss_FDR,
 na.rm = TRUE),
 median_logFC =
 safe_median(
 df$logFC),
 median_abs_logFC =
 safe_median(
 abs(df$logFC)),
 max_abs_logFC =
 ifelse(
 all(
 !is.finite(
 df$logFC)),
 NA_real_,
 max(
 abs(df$logFC),
 na.rm = TRUE)))
 differential_summary[[cell]] <-
 summary_df}
differential_summary_df <-
 bind_rows(
 differential_summary)
# Add the barcode threshold to the summary.
differential_summary_df <-
 differential_summary_df %>%
 mutate(
 minimum_barcode_threshold =
 MIN_BARCODES_PER_OLIGO
) %>%
 select(
 cell_line,
 minimum_barcode_threshold,
 everything())
write.csv(
 differential_summary_df,
 file.path(
 OUTPUT_DIR,
 "ALL_CELL_LINES_WT_MUT_differential_summary.csv"),
 row.names = FALSE)
save_checkpoint(
 differential_summary_df,
 "11_differential_summary.rds")
# ============================================================
# 18. COMBINED OUTPUTS
# ============================================================
cat("\n")
cat("============================================================\n")
cat("Creating combined outputs\n")
cat("============================================================\n")
combined_differential <-
 bind_rows(
 differential_results)
combined_significant <-
 bind_rows(
 significant_results)
combined_strong <-
 bind_rows(
 strong_results)
combined_gain <-
 bind_rows(
 gain_results)
combined_loss <-
 bind_rows(
 loss_results)
write.csv(
 combined_differential,
 file.path(
 OUTPUT_DIR,
 "ALL_CELL_LINES_WT_MUT_differential_results.csv"),
 row.names = FALSE)
write.csv(
 combined_significant,
 file.path(
 OUTPUT_DIR,
 "ALL_CELL_LINES_WT_MUT_significant_FDR.csv"),
 row.names = FALSE)
write.csv(
 combined_strong,
 file.path(
 OUTPUT_DIR,
 "ALL_CELL_LINES_WT_MUT_strong_effect_FDR.csv"),
 row.names = FALSE)
write.csv(
 combined_gain,
 file.path(
 OUTPUT_DIR,
 "ALL_CELL_LINES_WT_MUT_gain_FDR.csv"),
 row.names = FALSE)
write.csv(
 combined_loss,
 file.path(
 OUTPUT_DIR,
 "ALL_CELL_LINES_WT_MUT_loss_FDR.csv"),
 row.names = FALSE)
# ============================================================
# 19. FINAL CHECKPOINT
# ============================================================
save_checkpoint(
 list(
 minimum_barcode_threshold =
 MIN_BARCODES_PER_OLIGO,
 barcode_filter_summary =
 barcode_filter_summary,
 barcode_filter_summary_df =
 barcode_filter_summary_df,
 barcode_filter_allele_summary_df =
 barcode_filter_allele_summary_df,
 activity_results =
 activity_results,
 activity_summaries =
 activity_summaries,
 analysis_matrices =
 analysis_matrices,
 comparison_structures =
 comparison_structures,
 variant_objects =
 variant_objects,
 differential_results =
 differential_results,
 significant_results =
 significant_results,
 strong_results =
 strong_results,
 gain_results =
 gain_results,
 loss_results =
 loss_results,
 differential_summary =
 differential_summary_df),
 "12_SCRIPT2_complete.rds")
# ============================================================
# 20. FINAL SUMMARY
# ============================================================
cat("\n")
cat("============================================================\n")
cat("SCRIPT 2 COMPLETE\n")
cat("============================================================\n")
cat(
 "\nMinimum barcode threshold for formal WT/MUT analysis: ",
 MIN_BARCODES_PER_OLIGO,
 "\n",
 sep = "")
cat(
 "\nCell lines analyzed:\n")
cat(
 paste(
 names(differential_results),
 collapse = ", "),
 "\n")
cat(
 "\nK562 was excluded because only one retained replicate was available.\n")
cat(
 "\nBarcode-support filtering summary:\n")
for (cell in names(
 barcode_filter_summary
)) {
 s <-
 barcode_filter_summary[[cell]]$overall
 cat(
 " ",
 cell,
 ": ",
 s$retained_oligos,
 "/",
 s$total_oligos,
 " oligos retained (",
 round(
 s$percent_retained,
 1),
 "%)\n",
 sep = "")}
cat(
 "\nCombined WT/MUT comparisons after >= ",
 MIN_BARCODES_PER_OLIGO,
 "-barcode filtering: ",
 nrow(combined_differential),
 "\n",
 sep = "")
cat(
 "FDR-significant comparisons: ",
 nrow(combined_significant),
 "\n",
 sep = "")
cat(
 "FDR-significant with |effect| >= ",
 STRONG_EFFECT_THRESHOLD,
 ": ",
 nrow(combined_strong),
 "\n",
 sep = "")
cat(
 "FDR-significant gains with effect >= ",
 WEAK_EFFECT_THRESHOLD,
 ": ",
 nrow(combined_gain),
 "\n",
 sep = "")
cat(
 "FDR-significant losses with effect <= -",
 WEAK_EFFECT_THRESHOLD,
 ": ",
 nrow(combined_loss),
 "\n",
 sep = "")
cat(
 "\nOutput directory:\n",
 OUTPUT_DIR,
 "\n")
cat(
 "\nCheckpoint directory:\n",
 CHECKPOINT_DIR,
 "\n")
cat(
 "\nScript 2 finished successfully.\n")
####
####
####
####
#Reloading partial differential fits
for (cell in names(CELL_LINES)) {
 # Skip cell lines that were already successfully fitted
 if (cell %in% names(differential_fits)) {
 cat(
 "\nSkipping ",
 cell,
 " -- already present in partial checkpoint.\n",
 sep = "")
 next}
 obj <-
 variant_objects[[cell]]
 cat(
 "\n------------------------------------------------------------\n")
 cat(
 "Cell line: ",
 cell,
 "\n",
 sep = "")
 cat(
 "Comparisons: ",
 nrow(obj$dna_variant),
 "\n",
 sep = "")
 cat(
 "Replicates: ",
 obj$n_reps,
 "\n",
 sep = "")
 if (obj$n_reps < 2) {
 warning(
 "Skipping ",
 cell,
 " because fewer than 2 replicates are available.")
 next}
 cat(
 "Running mpralm()...\n")
 fit <-
 BCalm::mpralm(
 object =
 obj$BcVariant,
 design =
 obj$design,
 aggregate =
 "none",
 normalize =
 TRUE,
 model_type =
 "corr_groups",
 block =
 obj$block)
 differential_fits[[cell]] <-
 fit
 save_checkpoint(
 differential_fits,
 "08_WT_MUT_differential_fits_partial.rds")
 cat(
 "Fit completed for ",
 cell,
 "\n",
 sep = "")}
