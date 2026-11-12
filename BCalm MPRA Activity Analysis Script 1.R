############################################################
# BCalm MPRA ACTIVITY ANALYSIS — SCRIPT 1
#
# PURPOSE
#
# Run BCalm separately for each MPRA cell line using the
# MPRAsnakeflow reporter-experiment barcode count files.
#
############################################################
# Test_WildType = experimental WT oligos
# Test_Mutant = experimental mutant oligos
# WT_ncncm = WT NCNCM constructs
# MUT_ncncm = mutant NCNCM constructs
# Expression_controls = expression controls
# Positive_control = positive control
#
# Candidate empirical negatives are tracked separately from
# the six construct labels. They are identified using the
# specified activity-based criteria and barcode threshold.
#
############################################################
#
# Candidate inactive oligos are instead identified among
# TEST WT/MUT oligos based on activity CLOSEST TO ZERO:
#
# distance_from_zero = abs(logFC)
#
############################################################
# 0. PACKAGES
############################################################
library(BCalm)
library(dplyr)
library(ggplot2)
############################################################
# 1. GLOBAL SETTINGS
############################################################
experiment_dirs <- list(
 run_exp_260818 =
 "../MPRAsnakeflow/run_exp_260818",
 run_second_cell_lines =
 "../MPRAsnakeflow/run_second_cell_lines")
labels_file <- "labels_updated.tsv"
output_dir <- "BCalm_activity_new_labels_troubleshooting"
dir.create(
 output_dir,
 showWarnings = FALSE,
 recursive = TRUE)
############################################################
# BCalm downsampling
############################################################
# BCalm default = 0.95.
#
# Set to NULL to disable downsampling.
downsample_percentile <- 0.95
# Minimum number of eligible barcodes required per oligo
# before and after RNA-outlier removal.
minimum_bcalm_barcodes <- 10
############################################################
# Candidate inactive settings
############################################################
candidate_negative_percentiles <- c(
 0.05,
 0.10)
minimum_candidate_barcodes <- 10
############################################################
# 2. CELL-LINE CONFIGURATION
############################################################
cell_lines <- list(
 HEK = list(
 run = "run_exp_260818",
 replicates = c(1, 3)),
 U2OS = list(
 run = "run_exp_260818",
 replicates = c(1, 2)),
 K562 = list(
 run = "run_exp_260818",
 replicates = c(3)),
 U2987 = list(
 run = "run_second_cell_lines",
 replicates = c(1, 2, 3)),
 HepG2 = list(
 run = "run_second_cell_lines",
 replicates = c(2, 3)),
 `SK-MEL30` = list(
 run = "run_second_cell_lines",
 replicates = c(1, 2, 3)))
############################################################
# 3. INFORMATION
############################################################
cat("\n")
cat("========================================================\n")
cat("BCalm MPRA ACTIVITY ANALYSIS\n")
cat("========================================================\n")
cat("\nR version:\n")
print(R.version.string)
cat("\nBCalm version:\n")
print(packageVersion("BCalm"))
cat("\nLabels file:\n")
cat(labels_file, "\n")
cat("\nOutput directory:\n")
cat(output_dir, "\n")
cat("\nDownsampling:\n")
cat(
 ifelse(
 is.null(downsample_percentile),
 "disabled",
 paste0(
 downsample_percentile * 100,
 "%")),
 "\n")
############################################################
# 4. INPUT CHECKS
############################################################
if (!file.exists(labels_file)) {
 stop(
 "Labels file does not exist:\n",
 labels_file)}
for (run_name in names(experiment_dirs)) {
 run_dir <- experiment_dirs[[run_name]]
 if (!dir.exists(run_dir)) {
 stop(
 "Experiment directory does not exist:\n",
 run_dir)}}
############################################################
# 5. READ LABELS
############################################################
cat("\n")
cat("========================================================\n")
cat("Reading labels\n")
cat("========================================================\n")
labels <- read.delim(
 labels_file,
 header = TRUE,
 sep = "\t",
 stringsAsFactors = FALSE,
 check.names = FALSE)
if (ncol(labels) != 2) {
 stop("Expected exactly two columns in labels_updated.tsv.")}
# The current file has generic column names: Column1, Column2
names(labels) <- c("oligo_name", "label")
# Remove accidental leading/trailing whitespace
labels$oligo_name <- trimws(labels$oligo_name)
labels$label <- trimws(labels$label)
cat("\nLabel dimensions:\n")
print(dim(labels))
cat("\nLabel distribution:\n")
print(table(labels$label))
############################################################
# 6. REQUIRED LABEL CHECK
############################################################
expected_labels <- c(
 "Test_WildType",
 "Test_Mutant",
 "WT_ncncm",
 "MUT_ncncm",
 "Expression_controls",
 "Positive_control")
missing_labels <- setdiff(
 expected_labels,
 unique(labels$label))
unexpected_labels <- setdiff(
 unique(labels$label),
 expected_labels)
if (length(missing_labels) > 0) {
 stop(
 "Expected labels missing from labels file: ",
 paste(missing_labels, collapse = ", "))}
if (length(unexpected_labels) > 0) {
 stop(
 "Unexpected labels found in labels file: ",
 paste(unexpected_labels, collapse = ", "))}
cat("\nAll six expected construct labels are present.\n")
############################################################
# 7. DUPLICATE LABEL CHECK
############################################################
duplicate_labels <- labels %>%
 dplyr::count(
 oligo_name
) %>%
 filter(
 n > 1)
cat(
 "\nOligos with duplicate labels:",
 nrow(duplicate_labels),
 "\n")
if (nrow(duplicate_labels) > 0) {
 print(
 head(
 duplicate_labels,
 20))
 stop(
 "Duplicate oligo names found in labels file.")}
############################################################
# 8. NAMED LABEL VECTOR
############################################################
label_vector <- labels$label
names(label_vector) <- labels$oligo_name
############################################################
# 9. WT/MUT PAIR IDENTIFIER
############################################################
candidate_file <- "Final_mpra_candidates_pancancer_20240111_incl_headers.csv"
if (!file.exists(candidate_file)) {
 stop("Candidate file not found: ", candidate_file)}
mpra_candidates <- read.csv(
 candidate_file,
 stringsAsFactors = FALSE,
 check.names = FALSE)
# Check required columns
required_columns <- c(
 "chr",
 "genomic_position",
 "oligo_name")
missing_columns <- setdiff(
 required_columns,
 names(mpra_candidates))
if (length(missing_columns) > 0) {
 stop(
 "Missing columns in candidate file: ",
 paste(missing_columns, collapse = ", "))}
############################################################
# 9A. CHECK FOR MULTIPLE WT OLIGOS AT THE SAME
# CHROMOSOME-POSITION COMBINATION
############################################################
wt_candidates <- mpra_candidates %>%
 filter(
 grepl("^wt_", oligo_name, ignore.case = TRUE))
multiple_wt_positions <- wt_candidates %>%
 group_by(chr, genomic_position) %>%
 summarise(
 n_wt = n_distinct(oligo_name),
 wt_names = paste(
 sort(unique(oligo_name)),
 collapse = " | "),
 .groups = "drop"
) %>%
 filter(n_wt > 1) %>%
 arrange(desc(n_wt))
cat("\n========================================================\n")
cat("9A. MULTIPLE WT OLIGOS AT THE SAME CHROMOSOME-POSITION\n")
cat("========================================================\n")
cat(
 "Chromosome-position combinations with multiple WT oligos:",
 nrow(multiple_wt_positions),
 "\n")
if (nrow(multiple_wt_positions) > 0) {
 print(
 as.data.frame(multiple_wt_positions),
 row.names = FALSE)
 # Show full candidate records for ambiguous WT combinations
 multiple_wt_details <- wt_candidates %>%
 semi_join(
 multiple_wt_positions,
 by = c("chr", "genomic_position")
) %>%
 arrange(chr, genomic_position, oligo_name) %>%
 select(
 chr,
 genomic_position,
 oligo_name,
 allele_name,
 ref,
 alt,
 oligo_number)
 cat("\nFull details for these WT oligos:\n")
 print(
 as.data.frame(multiple_wt_details),
 row.names = FALSE)}
############################################################
# BUILD WT LOOKUP
############################################################
# One WT-derived pair_id per chromosome-position combination.
# Do not silently choose a WT if multiple distinct WTs exist.
wt_lookup <- wt_candidates %>%
 distinct(
 chr,
 genomic_position,
 oligo_name
) %>%
 group_by(
 chr,
 genomic_position
) %>%
 summarise(
 n_wt = n_distinct(oligo_name),
 pair_id = if (n_distinct(oligo_name) == 1) {
 sub(
 "^wt_",
 "",
 dplyr::first(oligo_name),
 ignore.case = TRUE)
 } else {
 NA_character_},
 .groups = "drop"
) %>%
 filter(n_wt == 1) %>%
 select(
 chr,
 genomic_position,
 pair_id)
############################################################
# ASSIGN WT-DERIVED PAIR IDs TO ALL OLIGOS
############################################################
pair_lookup <- mpra_candidates %>%
 select(
 oligo_name,
 chr,
 genomic_position
) %>%
 distinct() %>%
 left_join(
 wt_lookup,
 by = c("chr", "genomic_position"))
############################################################
# 9B. INSPECT OLIGOS WITHOUT A WT MAPPING
############################################################
unmapped_oligos <- pair_lookup %>%
 filter(is.na(pair_id))
cat("\n========================================================\n")
cat("9B. OLIGOS WITHOUT A WT MAPPING\n")
cat("========================================================\n")
cat(
 "Oligos without a WT pair_id:",
 nrow(unmapped_oligos),
 "\n")
if (nrow(unmapped_oligos) > 0) {
 unmapped_details <- unmapped_oligos %>%
 left_join(
 mpra_candidates %>%
 select(
 oligo_name,
 allele_name,
 ref,
 alt,
 oligo_number
) %>%
 distinct(),
 by = "oligo_name"
) %>%
 distinct()
 print(
 as.data.frame(unmapped_details),
 row.names = FALSE)}
############################################################
# FINAL LOOKUP FUNCTION
############################################################
# Keep original oligo names unchanged.
# Pair IDs are stored separately.
pair_id_vector <- pair_lookup$pair_id
names(pair_id_vector) <- pair_lookup$oligo_name
get_pair_id <- function(x) {
 unname(pair_id_vector[x])}
cat("\n========================================================\n")
cat("PAIR MAPPING SUMMARY\n")
cat("========================================================\n")
cat("Candidate oligos:", nrow(pair_lookup), "\n")
cat(
 "Oligos with a WT pair_id:",
 sum(!is.na(pair_lookup$pair_id)),
 "\n")
cat(
 "Oligos without a WT pair_id:",
 sum(is.na(pair_lookup$pair_id)),
 "\n")
cat(
 "Ambiguous chromosome-position combinations:",
 nrow(multiple_wt_positions),
 "\n")
############################################################
# 10. REPLICATE ACTIVITY FUNCTION
############################################################
calculate_replicate_activity <- function(
 bc,
 keep_reps
) {
 result_list <- list()
 for (rep in keep_reps) {
 dna_col <- paste0(
 "dna_count_",
 rep)
 rna_col <- paste0(
 "rna_count_",
 rep)
 dna_total <- sum(
 bc[[dna_col]],
 na.rm = TRUE)
 rna_total <- sum(
 bc[[rna_col]],
 na.rm = TRUE)
 if (
 dna_total <= 0 ||
 rna_total <= 0
) {
 stop(
 "DNA or RNA total is zero for replicate ",
 rep)}
 ########################################################
 # Library-size normalization
 ########################################################
 dna_norm <- (
 bc[[dna_col]] /
 dna_total
) * 1e9
 rna_norm <- (
 bc[[rna_col]] /
 rna_total
) * 1e9
 ########################################################
 # Barcode-level log2 RNA/DNA ratio
 ########################################################
 barcode_logratio <- log2(
 (rna_norm + 1) /
 (dna_norm + 1))
 temp <- data.frame(
 oligo_name = bc$oligo_name,
 barcode_logratio = barcode_logratio,
 stringsAsFactors = FALSE)
 ########################################################
 # Summarize barcode activity by oligo
 ########################################################
 replicate_summary <- temp %>%
 group_by(
 oligo_name
) %>%
 summarise(
 activity =
 mean(
 barcode_logratio,
 na.rm = TRUE),
 n_barcodes =
 sum(
 is.finite(
 barcode_logratio)),
 .groups = "drop")
 ########################################################
 # Rename replicate-specific columns
 ########################################################
 colnames(replicate_summary)[
 colnames(replicate_summary) == "activity"
 ] <- paste0(
 "rep",
 rep,
 "_activity")
 colnames(replicate_summary)[
 colnames(replicate_summary) == "n_barcodes"
 ] <- paste0(
 "rep",
 rep,
 "_barcodes")
 ########################################################
 # Store replicate result
 ########################################################
 result_list[[as.character(rep)]] <-
 replicate_summary}
 ##########################################################
 # Merge replicate summaries
 ##########################################################
 result <- Reduce(
 function(x, y) {
 full_join(
 x,
 y,
 by = "oligo_name")},
 result_list)
 return(result)}
############################################################
# 11. RUN ONE CELL LINE
############################################################
run_bcalm_cellline <- function(
 cell_line,
 run_name,
 keep_reps,
 experiment_dirs,
 label_vector,
 output_dir
) {
 cat("\n\n")
 cat("########################################################\n")
 cat(
 "# Cell line: ",
 cell_line,
 "\n",
 sep = "")
 cat(
 "# Sequencing run: ",
 run_name,
 "\n",
 sep = "")
 cat("########################################################\n")
 ##########################################################
 # Output directory
 ##########################################################
 cell_output_dir <- file.path(
 output_dir,
 cell_line)
 dir.create(
 cell_output_dir,
 showWarnings = FALSE,
 recursive = TRUE)
 ##########################################################
 # Barcode file
 ##########################################################
 barcode_dir <- file.path(
 experiment_dirs[[run_name]],
 "results",
 "experiments",
 "countBasic")
 barcode_file <- file.path(
 barcode_dir,
 paste0(
 "reporter_experiment.barcode.",
 cell_line,
 ".fromFile.default.all.tsv.gz"))
 cat("\nBarcode file:\n")
 cat(
 barcode_file,
 "\n")
 if (!file.exists(barcode_file)) {
 stop(
 "Barcode file does not exist:\n",
 barcode_file)}
 ##########################################################
 # Replicate configuration
 ##########################################################
 all_replicates <- 1:3
 if (
 length(keep_reps) < 1
) {
 stop(
 "No replicates selected for ",
 cell_line)}
 if (
 !all(
 keep_reps %in% all_replicates)
) {
 stop(
 "Invalid replicate selected for ",
 cell_line)}
 excluded_reps <- setdiff(
 all_replicates,
 keep_reps)
 cat(
 "\nReplicates retained: ",
 paste(
 keep_reps,
 collapse = ", "),
 "\n",
 sep = "")
 cat(
 "Replicates excluded: ",
 ifelse(
 length(excluded_reps) == 0,
 "none",
 paste(
 excluded_reps,
 collapse = ", ")),
 "\n",
 sep = "")
 ##########################################################
 # READ BARCODE DATA
 ##########################################################
 cat(
 "\nReading barcode counts...\n")
 bc <- read.delim(
 gzfile(
 barcode_file),
 header = TRUE,
 sep = "\t",
 stringsAsFactors = FALSE,
 check.names = FALSE)
 cat(
 "Barcode table dimensions: ",
 paste(
 dim(bc),
 collapse = " x "),
 "\n",
 sep = "")
 ##########################################################
 # COLUMN CHECK
 ##########################################################
 expected_columns <- c(
 "barcode",
 "oligo_name",
 "dna_count_1",
 "rna_count_1",
 "dna_count_2",
 "rna_count_2",
 "dna_count_3",
 "rna_count_3")
 missing_columns <- setdiff(
 expected_columns,
 colnames(bc))
 if (
 length(missing_columns) > 0
) {
 stop(
 cell_line,
 " barcode table is missing:\n",
 paste(
 missing_columns,
 collapse = ", "))}
 ##########################################################
 # KEEP RETAINED REPLICATES
 ##########################################################
 dna_columns <- paste0(
 "dna_count_",
 keep_reps)
 rna_columns <- paste0(
 "rna_count_",
 keep_reps)
 bc <- bc[
,
 c(
 "barcode",
 "oligo_name",
 dna_columns,
 rna_columns),
 drop = FALSE
 ]
 ##########################################################
 # BARCODE SUPPORT
 ##########################################################
 ##########################################################
 # BARCODE ELIGIBILITY AND RNA OUTLIER FILTERING
 ##########################################################
 # Preserve the selected-replicate barcode data before filtering.
 bc_raw <- bc
 # All selected DNA and RNA count columns must be numeric.
 count_columns <- c(dna_columns, rna_columns)
 if (!all(vapply(bc[count_columns], is.numeric, logical(1)))) {
 stop(
 "DNA/RNA count columns must be numeric for ",
 cell_line)}
 cat("\nInitial barcode rows: ", nrow(bc), "\n", sep = "")
 ##########################################################
 # 1. DNA/RNA ELIGIBILITY
 ##########################################################
 # Retain a barcode only if:
 # - all selected DNA/RNA counts are non-missing
 # - every selected DNA/RNA count is >= 1
 complete_rows <- complete.cases(
 bc[, count_columns, drop = FALSE])
 positive_rows <- rep(FALSE, nrow(bc))
 positive_rows[complete_rows] <-
 rowSums(
 as.matrix(
 bc[complete_rows, count_columns, drop = FALSE]
) >= 1
) == length(count_columns)
 eligible_rows <- complete_rows & positive_rows
 cat(
 "Barcodes failing DNA/RNA eligibility: ",
 sum(!eligible_rows),
 "\n",
 sep = "")
 bc <- bc[eligible_rows,, drop = FALSE]
 cat(
 "Barcodes remaining after DNA/RNA eligibility: ",
 nrow(bc),
 "\n",
 sep = "")
 ##########################################################
 # 2. MINIMUM BARCODE SUPPORT BEFORE OUTLIER REMOVAL
 ##########################################################
 support_before_outliers <- bc %>%
 group_by(oligo_name) %>%
 summarise(
 n_barcodes = n_distinct(barcode),
 .groups = "drop")
 oligos_before_outliers <- support_before_outliers %>%
 filter(
 n_barcodes >= minimum_bcalm_barcodes
) %>%
 pull(oligo_name)
 cat(
 "Oligos before outlier filtering: ",
 n_distinct(bc$oligo_name),
 "\n",
 sep = "")
 cat(
 "Oligos with at least ",
 minimum_bcalm_barcodes,
 " eligible barcodes: ",
 length(oligos_before_outliers),
 "\n",
 sep = "")
 bc <- bc %>%
 filter(
 oligo_name %in% oligos_before_outliers)
 if (nrow(bc) == 0) {
 stop(
 "No barcodes remain after minimum barcode-support filtering for ",
 cell_line)}
 ##########################################################
 # 3. RNA OUTLIER DETECTION WITHIN EACH OLIGO
 ##########################################################
 # Evaluate each retained replicate separately.
 #
 # An RNA barcode is an outlier if its RNA count is:
 # greater than the oligo-specific mean + 3 SD.
 #
 # This is a two-tailed RNA-count filter.
 outlier_columns <- character(0)
 for (rep in keep_reps) {
 rna_col <- paste0(
 "rna_count_",
 rep)
 flag_col <- paste0(
 "outlier_rep_",
 rep)
 outlier_columns <- c(
 outlier_columns,
 flag_col)
 rep_outliers <- bc %>%
 group_by(oligo_name) %>%
 mutate(
 mean_RNA = mean(.data[[rna_col]]),
 sd_RNA = sd(.data[[rna_col]]),
 RNA_outlier =
 if_else(
 is.finite(sd_RNA) & sd_RNA > 0,
 .data[[rna_col]] > mean_RNA + 3 * sd_RNA |
 .data[[rna_col]] < mean_RNA - 3 * sd_RNA,
 FALSE)
) %>%
 ungroup() %>%
 pull(RNA_outlier)
 bc[[flag_col]] <- rep_outliers
 cat(
 "RNA outliers in replicate ",
 rep,
 ": ",
 sum(rep_outliers),
 "\n",
 sep = "")}
 ##########################################################
 # 4. COMBINE OUTLIER FLAGS ACROSS REPLICATES
 ##########################################################
 # Conservative paired-barcode approach:
 # If a barcode is an RNA outlier in ANY retained replicate,
 # remove that barcode from the entire barcode table.
 #
 # This preserves the same barcode set across retained
 # replicates for downstream BCalm analysis.
 bc$outlier_any_replicate <-
 rowSums(
 as.data.frame(
 bc[, outlier_columns, drop = FALSE])
) > 0
 outlier_flags_combined <- bc %>%
 select(
 oligo_name,
 barcode,
 all_of(outlier_columns),
 outlier_any_replicate)
 outlier_barcodes <- sum(
 bc$outlier_any_replicate)
 cat(
 "Barcodes flagged in at least one replicate: ",
 outlier_barcodes,
 "\n",
 sep = "")
 bc <- bc %>%
 filter(
 !outlier_any_replicate
) %>%
 select(
 -all_of(outlier_columns),
 -outlier_any_replicate)
 cat(
 "Barcodes remaining after RNA-outlier removal: ",
 nrow(bc),
 "\n",
 sep = "")
 ##########################################################
 # 5. REAPPLY MINIMUM BARCODE SUPPORT
 ##########################################################
 support_after_outliers <- bc %>%
 group_by(oligo_name) %>%
 summarise(
 n_barcodes = n_distinct(barcode),
 .groups = "drop")
 oligos_after_outliers <- support_after_outliers %>%
 filter(
 n_barcodes >= minimum_bcalm_barcodes
) %>%
 pull(oligo_name)
 cat(
 "Oligos remaining after outlier filtering: ",
 n_distinct(bc$oligo_name),
 "\n",
 sep = "")
 cat(
 "Oligos retaining at least ",
 minimum_bcalm_barcodes,
 " barcodes: ",
 length(oligos_after_outliers),
 "\n",
 sep = "")
 bc <- bc %>%
 filter(
 oligo_name %in% oligos_after_outliers)
 if (nrow(bc) == 0) {
 stop(
 "No oligos remain after RNA-outlier and minimum-support filtering for ",
 cell_line)}
 ##########################################################
 # 6. FINAL BARCODE SUPPORT BEFORE DOWNSAMPLING
 ##########################################################
 barcode_support <- bc %>%
 group_by(oligo_name) %>%
 summarise(
 n_barcodes = n_distinct(barcode),
 .groups = "drop")
 cat(
 "Final eligible barcode rows before downsampling: ",
 nrow(bc),
 "\n",
 sep = "")
 cat(
 "Final oligos before downsampling: ",
 nrow(barcode_support),
 "\n",
 sep = "")
 cat(
 "Median eligible barcode support: ",
 median(barcode_support$n_barcodes),
 "\n",
 sep = "")
 ##########################################################
 # LABEL MATCHING
 ##########################################################
 barcode_oligos <- unique(
 bc$oligo_name)
 label_oligos <- names(
 label_vector)
 oligos_without_labels <- setdiff(
 barcode_oligos,
 label_oligos)
 labels_without_counts <- setdiff(
 label_oligos,
 barcode_oligos)
 cat(
 "\nUnique oligos in barcode table: ",
 length(barcode_oligos),
 "\n",
 sep = "")
 cat(
 "Unique oligos in labels: ",
 length(label_oligos),
 "\n",
 sep = "")
 cat(
 "Barcode oligos without labels: ",
 length(oligos_without_labels),
 "\n",
 sep = "")
 cat(
 "Labels without barcode counts: ",
 length(labels_without_counts),
 "\n",
 sep = "")
 ##########################################################
 # COMMON OLIGOS
 ##########################################################
 common_oligos <- intersect(
 barcode_oligos,
 label_oligos)
 if (
 length(common_oligos) == 0
) {
 stop(
 "No common oligos between barcode counts and labels for ",
 cell_line)}
 ##########################################################
 # REPLICATE ACTIVITY
 ##########################################################
 ##########################################################
 # DOWNSAMPLE BARCODES
 ##########################################################
 cat(
 "\nDownsampling barcodes...\n")
 if (
 is.null(
 downsample_percentile)
) {
 bc_down <- bc
 } else {
 set.seed(123)
 bc_down <- downsample_barcodes(
 bc,
 id_column_name = "oligo_name",
 percentile =
 downsample_percentile)}
 cat(
 "Rows before downsampling: ",
 nrow(bc),
 "\n",
 sep = "")
 cat(
 "Rows after downsampling: ",
 nrow(bc_down),
 "\n",
 sep = "")
 ##########################################################
 # BARCODE SUPPORT AFTER DOWNSAMPLING
 ##########################################################
 barcode_support_downsampled <- bc_down %>%
 group_by(oligo_name) %>%
 summarise(
 n_barcodes_downsampled = n_distinct(barcode),
 .groups = "drop")
 ##########################################################
 # REPLICATE ACTIVITY
 ##########################################################
 cat(
 "\nCalculating replicate-level activity...\n")
 replicate_activity <-
 calculate_replicate_activity(
 bc = bc_down,
 keep_reps = keep_reps)
 ##########################################################
 # CREATE DNA/RNA MATRICES
 ##########################################################
 cat(
 "\nCreating DNA dataframe...\n")
 dna_elem <- create_dna_df(
 bc_down,
 id_column_name = "oligo_name")
 cat(
 "Creating RNA dataframe...\n")
 rna_elem <- create_rna_df(
 bc_down,
 id_column_name = "oligo_name")
 ##########################################################
 # COMMON OLIGOS AFTER BCalm MATRIX CREATION
 ##########################################################
 ##########################################################
 # COMMON OLIGOS AFTER BCalm MATRIX CREATION
 ##########################################################
 common_oligos <- Reduce(
 intersect,
 list(
 rownames(dna_elem),
 rownames(rna_elem),
 names(label_vector)))
 if (length(common_oligos) == 0) {
 stop(
 "No common oligos across DNA, RNA, and labels for ",
 cell_line)}
 dna_elem <- dna_elem[
 common_oligos,
,
 drop = FALSE
 ]
 rna_elem <- rna_elem[
 common_oligos,
,
 drop = FALSE
 ]
 label_vector_cell <- label_vector[common_oligos]
 # Verify identical row order
 stopifnot(
 identical(rownames(dna_elem), rownames(rna_elem)),
 identical(rownames(dna_elem), names(label_vector_cell)))
 ##########################################################
 # LABEL DISTRIBUTION
 ##########################################################
 cat(
 "\nLabels represented in this cell line:\n")
 print(
 table(
 label_vector_cell))
 ##########################################################
 # PAIR INFORMATION
 ##########################################################
 result_pair_info <- data.frame(
 oligo_name =
 names(label_vector_cell),
 label =
 as.character(
 label_vector_cell),
 pair_id =
 get_pair_id(
 names(label_vector_cell)),
 stringsAsFactors = FALSE)
 ##########################################################
 # WT/MUT PAIR STRUCTURE
 ##########################################################
 pair_summary <- result_pair_info %>%
 filter(
 label %in%
 c(
 "Test_WildType",
 "Test_Mutant"),
 !is.na(pair_id)
) %>%
 group_by(
 pair_id
) %>%
 summarise(
 n_WT =
 sum(
 label == "Test_WildType"),
 n_MUT =
 sum(
 label == "Test_Mutant"),
 .groups = "drop")
 cat(
 "\nWT/MUT construct structure:\n")
 cat(
 "Total WT/MUT construct groups: ",
 nrow(pair_summary),
 "\n",
 sep = "")
 cat(
 "Groups with 1 WT: ",
 sum(
 pair_summary$n_WT == 1),
 "\n",
 sep = "")
 cat(
 "Groups with >1 WT: ",
 sum(
 pair_summary$n_WT > 1),
 "\n",
 sep = "")
 cat(
 "Groups with 1 MUT: ",
 sum(
 pair_summary$n_MUT == 1),
 "\n",
 sep = "")
 cat(
 "Groups with >1 MUT: ",
 sum(
 pair_summary$n_MUT > 1),
 "\n",
 sep = "")
 ##########################################################
 # CREATE MPRASET
 ##########################################################
 cat(
 "\nCreating BCalm MPRASet...\n")
 Bc <- BCalm::MPRASet(
 label =
 label_vector_cell,
 DNA =
 dna_elem,
 RNA =
 rna_elem,
 eid =
 row.names(
 dna_elem),
 barcode =
 NULL)
 cat(
 "\nMPRASet created successfully.\n")
 print(Bc)
 ##########################################################
 # DETERMINE BLOCK STRUCTURE
 ##########################################################
 n_reps <- length(
 keep_reps)
 n_dna_columns <- ncol(
 dna_elem)
 if (
 n_dna_columns %% n_reps != 0
) {
 stop(
 "DNA column number is not divisible by number of ",
 "retained replicates for ",
 cell_line)}
 barcodes_per_replicate <-
 n_dna_columns / n_reps
 block_vector <- rep(
 seq_len(n_reps),
 each =
 barcodes_per_replicate)
 cat(
 "\nNumber of retained replicates: ",
 n_reps,
 "\n",
 sep = "")
 cat(
 "DNA columns: ",
 n_dna_columns,
 "\n",
 sep = "")
 cat(
 "Estimated barcodes per replicate: ",
 barcodes_per_replicate,
 "\n",
 sep = "")
 ##########################################################
 # BCalm FIT
 ##########################################################
 cat("\n")
 cat(
 "========================================================\n")
 cat(
 "Running BCalm activity analysis for ",
 cell_line,
 "\n",
 sep = "")
 cat(
 "========================================================\n")
 ##########################################################
 # MULTI-REPLICATE
 ##########################################################
 if (
 n_reps >= 2
) {
 cat(
 "\nUsing replicated BCalm model with ",
 n_reps,
 " biological replicates.\n",
 sep = "")
 fit <- fit_elements(
 object = Bc,
 normalize = TRUE,
 block = block_vector)
 analysis_mode <- "replicated"
 cat(
 "\nfit_elements() completed for ",
 cell_line,
 ".\n",
 sep = "")}
 ##########################################################
 # SINGLE REPLICATE
 ##########################################################
 if (
 n_reps == 1
) {
 cat(
 "\nWARNING: ",
 cell_line,
 " has only one retained biological replicate.\n",
 sep = "")
 cat(
 "No replicate-based BCalm statistical model will be ",
 "used.\n")
 analysis_mode <-
 "single_replicate"
 fit <- NULL}
 ##########################################################
 # EXTRACT RESULTS
 ##########################################################
 cat(
 "\nExtracting activity results...\n")
 if (
 analysis_mode == "replicated"
) {
 result <- topTable(
 fit,
 coef = 1,
 number = Inf)}
 if (
 analysis_mode == "single_replicate"
) {
 result <- replicate_activity %>%
 select(
 oligo_name,
 matches(
 "^rep[0-9]+_activity$"))
 activity_column <- grep(
 "^rep[0-9]+_activity$",
 colnames(result),
 value = TRUE)
 result$logFC <-
 result[[activity_column]]
 result$AveExpr <-
 NA_real_
 result$t <-
 NA_real_
 result$P.Value <-
 NA_real_
 result$adj.P.Val <-
 NA_real_
 result$B <-
 NA_real_
 rownames(result) <-
 result$oligo_name}
 ##########################################################
 # ADD IDENTIFIERS
 ##########################################################
 result$oligo_name <-
 rownames(result)
 result <- result %>%
 mutate(
 cell_line =
 cell_line,
 analysis =
 "BCalm_element_activity",
 analysis_mode =
 analysis_mode,
 label =
 label_vector_cell[
 oligo_name
 ],
 allele =
 case_when(
 grepl(
 "^wt_",
 oligo_name,
 ignore.case = TRUE)
 ~ "WT",
 grepl(
 "^mut_",
 oligo_name,
 ignore.case = TRUE)
 ~ "MUT",
 TRUE
 ~ NA_character_),
 pair_id =
 get_pair_id(
 oligo_name))
 ##########################################################
 # BARCODE SUPPORT
 ##########################################################
 result <- result %>%
 left_join(
 barcode_support,
 by = "oligo_name"
) %>%
 left_join(
 barcode_support_downsampled,
 by = "oligo_name")
 ##########################################################
 # REPLICATE ACTIVITY
 ##########################################################
 result <- result %>%
 left_join(
 replicate_activity,
 by = "oligo_name")
 ##########################################################
 # REPLICATE CONSISTENCY
 ##########################################################
 replicate_activity_columns <- grep(
 "^rep[0-9]+_activity$",
 colnames(result),
 value = TRUE)
 if (
 length(
 replicate_activity_columns
) >= 2
) {
 result$replicate_mean <-
 rowMeans(
 result[
,
 replicate_activity_columns,
 drop = FALSE
 ],
 na.rm = TRUE)
 result$replicate_sd <-
 apply(
 result[
,
 replicate_activity_columns,
 drop = FALSE
 ],
 1,
 sd,
 na.rm = TRUE)
 } else {
 result$replicate_mean <-
 result[[
 replicate_activity_columns[1]
 ]]
 result$replicate_sd <-
 NA_real_}
 ##########################################################
 # ACTIVITY RANK / PERCENTILE
 ##########################################################
 result <- result %>%
 mutate(
 activity_rank =
 rank(
 -logFC,
 ties.method = "average",
 na.last = "keep"),
 activity_percentile =
 percent_rank(
 logFC),
 activity_direction =
 case_when(
 logFC > 0
 ~ "positive",
 logFC < 0
 ~ "negative",
 TRUE
 ~ "zero"))
 ##########################################################
 # CANDIDATE INACTIVE ELEMENTS
 ##########################################################
 # IMPORTANT:
 #
 # Candidate inactive elements are NOT the lowest activity.
 #
 # They are the TEST WT/MUT elements whose activity is
 # closest to zero.
 #
 # abs(logFC) is therefore used.
 #
 # Strongly negative values may represent genuine repression
 # and should not automatically be treated as inactive.
 test_results <- result %>%
 filter(
 label %in%
 c(
 "Test_WildType",
 "Test_Mutant")
) %>%
 mutate(
 distance_from_zero =
 abs(logFC))
 ##########################################################
 # THRESHOLDS
 ##########################################################
 threshold_5 <- quantile(
 test_results$distance_from_zero,
 probs = 0.05,
 na.rm = TRUE)
 threshold_10 <- quantile(
 test_results$distance_from_zero,
 probs = 0.10,
 na.rm = TRUE)
 cat("\n")
 cat(
 "Empirical activity-distance-from-zero thresholds:\n")
 cat(
 "Closest 5% threshold: ",
 threshold_5,
 "\n",
 sep = "")
 cat(
 "Closest 10% threshold: ",
 threshold_10,
 "\n",
 sep = "")
 ##########################################################
 # CANDIDATE NEGATIVE FLAGS
 ##########################################################
 result$candidate_negative_5pct <-
 result$label %in%
 c(
 "Test_WildType",
 "Test_Mutant"
) &
 abs(
 result$logFC
) <= threshold_5 &
 result$n_barcodes >=
 minimum_candidate_barcodes
 result$candidate_negative_10pct <-
 result$label %in%
 c(
 "Test_WildType",
 "Test_Mutant"
) &
 abs(
 result$logFC
) <= threshold_10 &
 result$n_barcodes >=
 minimum_candidate_barcodes
 ##########################################################
 # REORDER RESULT COLUMNS
 ##########################################################
 preferred_columns <- c(
 "cell_line",
 "analysis",
 "oligo_name",
 "label",
 "allele",
 "pair_id",
 "n_barcodes",
 "logFC",
 "AveExpr",
 "t",
 "P.Value",
 "adj.P.Val",
 "B",
 "replicate_mean",
 "replicate_sd")
 remaining_columns <- setdiff(
 colnames(result),
 preferred_columns)
 result <- result[
,
 c(
 preferred_columns,
 remaining_columns),
 drop = FALSE
 ]
 ##########################################################
 # SAVE COMPLETE RESULTS
 ##########################################################
 result_file <- file.path(
 cell_output_dir,
 paste0(
 cell_line,
 "_BCalm_activity_results.csv"))
 write.csv(
 result,
 result_file,
 row.names = FALSE)
 ##########################################################
 # SAVE TEST WT/MUT ACTIVITY
 ##########################################################
 test_file <- file.path(
 cell_output_dir,
 paste0(
 cell_line,
 "_Test_WT_MUT_activity.csv"))
 write.csv(
 test_results,
 test_file,
 row.names = FALSE)
 ##########################################################
 # SAVE CANDIDATE INACTIVES
 ##########################################################
 candidate_negatives <- result %>%
 filter(
 candidate_negative_10pct
) %>%
 arrange(
 abs(logFC))
 write.csv(
 candidate_negatives,
 file.path(
 cell_output_dir,
 paste0(
 cell_line,
 "_candidate_negative_oligos.csv")),
 row.names = FALSE)
 ##########################################################
 # SAVE PAIR STRUCTURE
 ##########################################################
 write.csv(
 pair_summary,
 file.path(
 cell_output_dir,
 paste0(
 cell_line,
 "_WT_MUT_pair_structure.csv")),
 row.names = FALSE)
 ##########################################################
 # SAVE REPLICATE ACTIVITY
 ##########################################################
 write.csv(
 replicate_activity,
 file.path(
 cell_output_dir,
 paste0(
 cell_line,
 "_replicate_activity.csv")),
 row.names = FALSE)
 ##########################################################
 # SAVE BCalm FIT
 ##########################################################
 if (
 analysis_mode == "replicated"
) {
 saveRDS(
 fit,
 file = file.path(
 cell_output_dir,
 paste0(
 cell_line,
 "_BCalm_fit_elements.rds")))}
 ##########################################################
 # SAVE WORKSPACE
 ##########################################################
 save(
 bc,
 bc_raw,
 outlier_flags_combined,
 bc_down,
 dna_elem,
 rna_elem,
 label_vector_cell,
 result_pair_info,
 pair_summary,
 barcode_support,
 barcode_support_downsampled,
 replicate_activity,
 Bc,
 fit,
 analysis_mode,
 result,
 file = file.path(
 cell_output_dir,
 paste0(
 cell_line,
 "_BCalm_activity_workspace.RData")))
 ##########################################################
 # SUMMARY
 ##########################################################
 cat("\n")
 cat("--------------------------------------------------------\n")
 cat(cell_line, " summary\n")
 cat("--------------------------------------------------------\n")
 cat(
 "Total oligos analyzed: ",
 nrow(result),
 "\n",
 sep = "")
 cat(
 "Test WT: ",
 sum(result$label == "Test_WildType", na.rm = TRUE),
 "\n",
 sep = "")
 cat(
 "Test MUT: ",
 sum(result$label == "Test_Mutant", na.rm = TRUE),
 "\n",
 sep = "")
 cat(
 "Positive controls: ",
 sum(result$label == "Positive_control", na.rm = TRUE),
 "\n",
 sep = "")
 cat(
 "Expression controls: ",
 sum(result$label == "Expression_controls", na.rm = TRUE),
 "\n",
 sep = "")
 cat(
 "WT NCNCM controls: ",
 sum(result$label == "WT_ncncm", na.rm = TRUE),
 "\n",
 sep = "")
 cat(
 "MUT NCNCM controls: ",
 sum(result$label == "MUT_ncncm", na.rm = TRUE),
 "\n",
 sep = "")
 cat(
 "Total NCNCM/background controls: ",
 sum(
 result$label %in% c("WT_ncncm", "MUT_ncncm"),
 na.rm = TRUE),
 "\n",
 sep = "")
 cat(
 "Construct groups: ",
 nrow(pair_summary),
 "\n",
 sep = "")
 cat(
 "Candidate inactive, closest 5%: ",
 sum(result$candidate_negative_5pct, na.rm = TRUE),
 "\n",
 sep = "")
 cat(
 "Candidate inactive, closest 10%: ",
 sum(result$candidate_negative_10pct, na.rm = TRUE),
 "\n",
 sep = "")
 cat(
 "Median BCalm activity: ",
 median(result$logFC, na.rm = TRUE),
 "\n",
 sep = "")
 cat(
 "Median barcode support: ",
 median(result$n_barcodes, na.rm = TRUE),
 "\n",
 sep = "")
 ##########################################################
 # RETURN
 ##########################################################
 return(result)}
############################################################
# 12. RUN ALL CELL LINES
############################################################
all_results <- list()
failed_cell_lines <- list()
for (
 cell_line in names(cell_lines)
) {
 config <- cell_lines[[cell_line]]
 cat("\n")
 cat(
 "========================================================\n")
 cat(
 "STARTING CELL LINE: ",
 cell_line,
 "\n",
 sep = "")
 cat(
 "========================================================\n")
 tryCatch(
 {
 all_results[[cell_line]] <-
 run_bcalm_cellline(
 cell_line =
 cell_line,
 run_name =
 config$run,
 keep_reps =
 config$replicates,
 experiment_dirs =
 experiment_dirs,
 label_vector =
 label_vector,
 output_dir =
 output_dir)
 cat(
 "\nSUCCESS: ",
 cell_line,
 " completed.\n",
 sep = "")},
 error = function(e) {
 cat(
 "\nERROR in ",
 cell_line,
 ":\n",
 sep = "")
 cat(
 conditionMessage(e),
 "\n")
 failed_cell_lines[[cell_line]] <<-
 conditionMessage(e)})}
############################################################
# 13. COMBINE RESULTS
############################################################
cat("\n")
cat(
 "========================================================\n")
cat(
 "Combining cell-line activity results\n")
cat(
 "========================================================\n")
if (
 length(all_results) == 0
) {
 stop(
 "No cell lines completed successfully.")}
combined_results <- bind_rows(
 all_results)
cat("\nSuccessful cell lines:\n")
cat(
 paste(
 names(all_results),
 collapse = ", "),
 "\n")
if (
 length(failed_cell_lines) > 0
) {
 cat(
 "\nFailed cell lines:\n")
 for (
 x in names(failed_cell_lines)
) {
 cat(
 x,
 ": ",
 failed_cell_lines[[x]],
 "\n",
 sep = "")}
} else {
 cat(
 "\nNo cell lines failed.\n")}
############################################################
# 13A. WT-RELATIVE MUTANT ACTIVITY
############################################################
cat("\n")
cat("========================================================\n")
cat("WT-RELATIVE MUTANT ACTIVITY COMPARISON\n")
cat("========================================================\n")
# Identify replicate activity columns
rep_cols <- grep(
 "^rep[0-9]+_activity$",
 names(combined_results),
 value = TRUE)
# WT rows with valid pair IDs
wt_activity <- combined_results %>%
 filter(
 label == "Test_WildType",
 !is.na(pair_id)
) %>%
 select(
 cell_line,
 pair_id,
 wt_oligo_name = oligo_name,
 wt_logFC = logFC,
 all_of(rep_cols)
) %>%
 rename_with(
 ~ paste0("wt_", .x),
 all_of(rep_cols))
# Validate that each cell-line/pair_id has at most one WT
wt_counts <- wt_activity %>%
 count(
 cell_line,
 pair_id,
 name = "n_WT")
ambiguous_wt <- wt_counts %>%
 filter(n_WT > 1)
if (nrow(ambiguous_wt) > 0) {
 cat("\nERROR: Multiple WT results for a pair ID:\n")
 print(as.data.frame(ambiguous_wt), row.names = FALSE)
 stop(
 "Resolve multiple WT results before calculating ",
 "WT-relative activity.")}
# All labelled mutant constructs, including MUT, INS and DEL,
# provided they have the Test_Mutant label and a valid pair_id.
mutant_activity <- combined_results %>%
 filter(
 label == "Test_Mutant",
 !is.na(pair_id)
) %>%
 select(
 cell_line,
 pair_id,
 mutant_oligo_name = oligo_name,
 mutant_allele = allele,
 mutant_logFC = logFC,
 all_of(rep_cols)
) %>%
 rename_with(
 ~ paste0("mutant_", .x),
 all_of(rep_cols))
# Identify mutant groups without a WT result in that cell line
unpaired_mutants <- mutant_activity %>%
 anti_join(
 wt_activity %>%
 select(cell_line, pair_id),
 by = c("cell_line", "pair_id"))
cat(
 "\nMutant oligos with no WT result in the same cell line: ",
 nrow(unpaired_mutants),
 "\n",
 sep = "")
if (nrow(unpaired_mutants) > 0) {
 print(
 as.data.frame(
 unpaired_mutants %>%
 select(
 cell_line,
 pair_id,
 mutant_oligo_name,
 mutant_allele)),
 row.names = FALSE)}
# Join every mutant to its own WT.
# A single WT can match multiple mutant rows.
mutant_vs_wt <- mutant_activity %>%
 inner_join(
 wt_activity,
 by = c("cell_line", "pair_id")
) %>%
 mutate(
 # Difference in BCalm activity estimates
 delta_logFC = mutant_logFC - wt_logFC,
 # Ratio of estimated fold changes, assuming logFC is log2
 relative_fold_change = 2^delta_logFC)
# Calculate replicate-specific differences when both
# WT and mutant replicate estimates are available
for (rep_col in rep_cols) {
 wt_col <- paste0("wt_", rep_col)
 mutant_col <- paste0("mutant_", rep_col)
 delta_col <- paste0("delta_", rep_col)
 mutant_vs_wt[[delta_col]] <-
 mutant_vs_wt[[mutant_col]] -
 mutant_vs_wt[[wt_col]]}
# Summarize the replicate-level differences
delta_rep_cols <- paste0("delta_", rep_cols)
if (length(delta_rep_cols) > 0) {
 delta_matrix <- as.matrix(
 mutant_vs_wt[, delta_rep_cols, drop = FALSE])
 mutant_vs_wt$n_replicates_compared <-
 rowSums(is.finite(delta_matrix))
 mutant_vs_wt$replicate_delta_mean <-
 rowMeans(
 delta_matrix,
 na.rm = TRUE)
 mutant_vs_wt$replicate_delta_mean[
 mutant_vs_wt$n_replicates_compared == 0
 ] <- NA_real_
 mutant_vs_wt$replicate_delta_sd <- apply(
 delta_matrix,
 1,
 function(x) {
 x <- x[is.finite(x)]
 if (length(x) >= 2) sd(x) else NA_real_})}
# Order the results for easier inspection
mutant_vs_wt <- mutant_vs_wt %>%
 arrange(
 cell_line,
 pair_id,
 mutant_oligo_name)
# Save a separate WT-relative activity table
wt_relative_file <- file.path(
 output_dir,
 "ALL_CELL_LINES_WT_relative_mutant_activity.csv")
write.csv(
 mutant_vs_wt,
 wt_relative_file,
 row.names = FALSE)
cat("\nWT-relative comparisons created: ", nrow(mutant_vs_wt), "\n", sep = "")
cat("Saved to:\n", wt_relative_file, "\n", sep = "")
# Preview the main effect-size columns
print(
 head(
 mutant_vs_wt %>%
 select(
 cell_line,
 pair_id,
 wt_oligo_name,
 mutant_oligo_name,
 mutant_allele,
 wt_logFC,
 mutant_logFC,
 delta_logFC,
 relative_fold_change,
 replicate_delta_mean,
 replicate_delta_sd,
 n_replicates_compared),
 20))
############################################################
# 13A. OVERALL ACTIVITY SUMMARY BY CONSTRUCT GROUP
# AND CANDIDATE EMPIRICAL NEGATIVES
############################################################
# Six original construct groups
activity_summary_input <- combined_results %>%
 mutate(
 summary_group = label)
# Candidate empirical negatives are a separate analysis
# designation. Their original construct labels are preserved.
#
# The 10% definition matches the existing combined candidate
# negative output generated in Section 15.
candidate_empirical_negatives <- combined_results %>%
 filter(
 candidate_negative_10pct %in% TRUE
) %>%
 mutate(
 summary_group = "Candidate_empirical_negatives")
# Add the candidate-negative designation as a seventh group.
# These rows intentionally overlap the original six groups.
activity_summary_input <- bind_rows(
 activity_summary_input,
 candidate_empirical_negatives)
# Summarize activity for each cell line and analysis group.
activity_by_group <- activity_summary_input %>%
 group_by(
 cell_line,
 summary_group
) %>%
 summarise(
 n_oligos = n_distinct(oligo_name),
 median_logFC = median(
 logFC,
 na.rm = TRUE),
 mean_logFC = mean(
 logFC,
 na.rm = TRUE),
 sd_logFC = sd(
 logFC,
 na.rm = TRUE),
 median_abs_logFC = median(
 abs(logFC),
 na.rm = TRUE),
 median_barcodes = median(
 n_barcodes,
 na.rm = TRUE),
 n_positive = sum(
 logFC > 0,
 na.rm = TRUE),
 n_negative = sum(
 logFC < 0,
 na.rm = TRUE),
 n_zero = sum(
 logFC == 0,
 na.rm = TRUE),
 .groups = "drop"
) %>%
 arrange(
 cell_line,
 summary_group)
write.csv(
 activity_by_group,
 file.path(
 output_dir,
 "ALL_CELL_LINES_activity_summary_by_group.csv"),
 row.names = FALSE)
cat("\nActivity summary by cell line and analysis group:\n")
print(
 as.data.frame(activity_by_group),
 row.names = FALSE)
############################################################
# 14. SAVE COMBINED RESULTS
############################################################
combined_file <- file.path(
 output_dir,
 "ALL_CELL_LINES_BCalm_activity_results.csv")
write.csv(
 combined_results,
 combined_file,
 row.names = FALSE)
cat(
 "\nCombined results saved:\n",
 combined_file,
 "\n")
############################################################
# 15. COMBINED CANDIDATE INACTIVES
############################################################
combined_candidate_negatives <-
 combined_results %>%
 filter(
 candidate_negative_10pct
) %>%
 arrange(
 cell_line,
 abs(logFC))
write.csv(
 combined_candidate_negatives,
 file.path(
 output_dir,
 "ALL_CELL_LINES_candidate_negative_oligos.csv"),
 row.names = FALSE)
############################################################
# 16. CELL-LINE SUMMARY
############################################################
summary_table <- combined_results %>%
 group_by(cell_line) %>%
 summarise(
 oligos = n(),
 Test_WildType = sum(
 label == "Test_WildType",
 na.rm = TRUE),
 Test_Mutant = sum(
 label == "Test_Mutant",
 na.rm = TRUE),
 WT_ncncm = sum(
 label == "WT_ncncm",
 na.rm = TRUE),
 MUT_ncncm = sum(
 label == "MUT_ncncm",
 na.rm = TRUE),
 Expression_controls = sum(
 label == "Expression_controls",
 na.rm = TRUE),
 Positive_control = sum(
 label == "Positive_control",
 na.rm = TRUE),
 Candidate_empirical_negatives = sum(
 candidate_negative_10pct %in% TRUE,
 na.rm = TRUE),
 construct_groups = n_distinct(
 pair_id[
 !is.na(pair_id) &
 label %in% c(
 "Test_WildType",
 "Test_Mutant")
 ]),
 median_activity = median(
 logFC,
 na.rm = TRUE),
 median_barcodes = median(
 n_barcodes,
 na.rm = TRUE),
 candidate_negative_5pct = sum(
 candidate_negative_5pct %in% TRUE,
 na.rm = TRUE),
 candidate_negative_10pct = sum(
 candidate_negative_10pct %in% TRUE,
 na.rm = TRUE),
 .groups = "drop")
write.csv(
 summary_table,
 file.path(
 output_dir,
 "BCalm_activity_cell_line_summary.csv"),
 row.names = FALSE)
############################################################
# 17. FINAL OUTPUT
############################################################
cat("\n")
cat(
 "========================================================\n")
cat(
 "BCalm ACTIVITY ANALYSIS COMPLETE\n")
cat(
 "========================================================\n")
cat(
 "\nCell lines successfully analyzed: ",
 length(all_results),
 "\n",
 sep = "")
cat(
 "Total result rows: ",
 nrow(combined_results),
 "\n",
 sep = "")
cat(
 "Total candidate inactive oligos: ",
 nrow(combined_candidate_negatives),
 "\n",
 sep = "")
cat(
 "\nResults by cell line:\n")
print(
 summary_table)
cat(
 "\nOutput directory:\n")
cat(
 output_dir,
 "\n")
cat(
 "\nDone.\n")
#######
############################################################
# SAVE WT/MUT PAIR LOOKUP FOR DOWNSTREAM SCRIPTS
############################################################
saveRDS(
 pair_lookup,
 file = file.path(
 output_dir,
 "pair_lookup.rds"))
write.csv(
 pair_lookup,
 file = file.path(
 output_dir,
 "pair_lookup.csv"),
 row.names = FALSE)
cat(
 "\nSaved shared WT/MUT pair mapping:\n",
 file.path(output_dir, "pair_lookup.rds"),
 "\n",
 sep = "")
save(
 bc_raw,
 bc,
 bc_down,
 outlier_flags_combined,
 dna_elem,
 rna_elem,
 label_vector_cell,
 result_pair_info,
 pair_summary,
 barcode_support,
 barcode_support_downsampled,
 replicate_activity,
 Bc,
 fit,
 analysis_mode,
 result,
 file = file.path(
 cell_output_dir,
 paste0(
 cell_line,
 "_BCalm_activity_workspace.RData")))
