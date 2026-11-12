############################################################
# SCRIPT 3
# BCalm final visualization and cross-cell-line analysis
#
# PURPOSE
# Final visualization and biological candidate analysis
# using completed Script 1 and Script 2 results.
#
# IMPORTANT
# This script DOES NOT rerun BCalm or mpralm.
#
# PRIMARY BIOLOGICAL ANALYSIS
# WT -> MUT differential results from Script 2.
#
# STRONG SIGNIFICANT EFFECT
# FDR <= 0.05
# abs(logFC) >= 1
#
# WEAK SIGNIFICANT EFFECT
# FDR <= 0.05
# abs(logFC) >= 0.5
#
############################################################
# 1. LIBRARIES
############################################################
suppressPackageStartupMessages({
 library(dplyr)
 library(tidyr)
 library(ggplot2)
 library(readr)
 library(stringr)
})
############################################################
# 2. GLOBAL SETTINGS
############################################################
CELL_LINES <- c(
 "HEK",
 "U2OS",
 "U2987",
 "HepG2",
 "SK-MEL30")
# K562 deliberately excluded from replicated analysis.
WEAK_EFFECT_THRESHOLD <- 0.5
STRONG_EFFECT_THRESHOLD <- 1.0
FDR_THRESHOLD <- 0.05
# Number of candidates for detailed visualization
N_TOP_SHARED <- 10
N_TOP_SPECIFIC <- 10
# Reproducible random examples
RANDOM_SEED <- 12345
N_RANDOM_EXAMPLES <- 5
############################################################
# 3. DIRECTORIES
############################################################
SCRIPT2_DIR <- "BCalm_differential_analysis_min10_new_labels_troubleshooting"
CHECKPOINT_DIR <- file.path(
 SCRIPT2_DIR,
 "checkpoints")
OUTPUT_DIR <- file.path(
 SCRIPT2_DIR,
 "Script3_visualizations")
PLOT_DIR <- file.path(
 OUTPUT_DIR,
 "plots")
TABLE_DIR <- file.path(
 OUTPUT_DIR,
 "tables")
DATA_DIR <- file.path(
 OUTPUT_DIR,
 "plotting_data")
QC_DIR <- file.path(
 PLOT_DIR,
 "QC")
ACTIVITY_DIR <- file.path(
 PLOT_DIR,
 "element_activity")
DIFF_DIR <- file.path(
 PLOT_DIR,
 "WT_MUT_differential")
CROSS_DIR <- file.path(
 PLOT_DIR,
 "cross_cell_line")
CANDIDATE_DIR <- file.path(
 PLOT_DIR,
 "candidates")
dir.create(OUTPUT_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(PLOT_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(TABLE_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(DATA_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(QC_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(ACTIVITY_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(DIFF_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(CROSS_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(CANDIDATE_DIR, recursive = TRUE, showWarnings = FALSE)
############################################################
# 4. BASIC INFORMATION
############################################################
cat("\n============================================================\n")
cat("BCalm SCRIPT 3\n")
cat("Final visualization + cross-cell-line NCCM analysis\n")
cat("============================================================\n")
cat("\nCell lines:\n")
print(CELL_LINES)
cat("\nStrong effect threshold:", STRONG_EFFECT_THRESHOLD, "\n")
cat("Weak effect threshold:", WEAK_EFFECT_THRESHOLD, "\n")
cat("FDR threshold:", FDR_THRESHOLD, "\n")
cat("\nOutput directory:\n")
cat(OUTPUT_DIR, "\n")
############################################################
# 5. HELPER FUNCTIONS
############################################################
save_csv <- function(x, filename) {
 write.csv(
 x,
 file.path(TABLE_DIR, filename),
 row.names = FALSE)}
save_plot <- function(
 plot_object,
 filename,
 width = 8,
 height = 6
) {
 ggsave(
 filename = filename,
 plot = plot_object,
 width = width,
 height = height,
 dpi = 300)}
safe_cor <- function(x, y) {
 ok <- is.finite(x) & is.finite(y)
 if (sum(ok) < 3) {
 return(NA_real_)}
 suppressWarnings(
 cor(
 x[ok],
 y[ok],
 method = "spearman"))}
get_pair_id <- function(x) {
 x <- sub("^wt_", "", x)
 x <- sub("^mut_", "", x)
 x}
get_allele <- function(x) {
 ifelse(
 grepl("^wt_", x),
 "WT",
 ifelse(
 grepl("^mut_", x),
 "MUT",
 NA_character_))}
############################################################
# 6. LOAD SCRIPT 2 RESULTS
############################################################
cat("\n============================================================\n")
cat("Loading Script 2 results\n")
cat("============================================================\n")
differential_results_file <- file.path(
 CHECKPOINT_DIR,
 "09_WT_MUT_differential_results.rds")
if (!file.exists(differential_results_file)) {
 stop(
 "Could not find Script 2 differential results:\n",
 differential_results_file)}
differential_results <- readRDS(
 differential_results_file)
cat("\nDifferential results loaded.\n")
cat("\nCell lines present:\n")
print(names(differential_results))
missing_cells <- setdiff(
 CELL_LINES,
 names(differential_results))
if (length(missing_cells) > 0) {
 stop(
 "Missing differential results for:\n",
 paste(missing_cells, collapse = ", "))}
############################################################
# 7. COMBINE DIFFERENTIAL RESULTS
############################################################
cat("\n============================================================\n")
cat("Combining WT/MUT differential results\n")
cat("============================================================\n")
differential_all <- bind_rows(
 differential_results[CELL_LINES],
 .id = "cell_line")
cat(
 "\nTotal WT/MUT comparisons:",
 nrow(differential_all),
 "\n")
############################################################
# 8. VALIDATE REQUIRED COLUMNS
############################################################
required_diff_columns <- c(
 "cell_line",
 "pair_id",
 "wt_oligo",
 "mut_oligo",
 "logFC",
 "P.Value",
 "adj.P.Val")
missing_diff_columns <- setdiff(
 required_diff_columns,
 colnames(differential_all))
if (length(missing_diff_columns) > 0) {
 stop(
 "Missing required Script 2 columns:\n",
 paste(missing_diff_columns, collapse = ", "))}
############################################################
# 9. DEFINE EFFECT CLASSES
############################################################
cat("\n============================================================\n")
cat("Defining biological effect classes\n")
cat("============================================================\n")
differential_all <- differential_all %>%
 mutate(
 effect_abs = abs(logFC),
 strong_significant =
 !is.na(adj.P.Val) &
 adj.P.Val <= FDR_THRESHOLD &
 effect_abs >= STRONG_EFFECT_THRESHOLD,
 weak_significant =
 !is.na(adj.P.Val) &
 adj.P.Val <= FDR_THRESHOLD &
 effect_abs >= WEAK_EFFECT_THRESHOLD,
 strong_gain =
 strong_significant &
 logFC >= STRONG_EFFECT_THRESHOLD,
 strong_loss =
 strong_significant &
 logFC <= -STRONG_EFFECT_THRESHOLD,
 weak_gain =
 weak_significant &
 logFC >= WEAK_EFFECT_THRESHOLD,
 weak_loss =
 weak_significant &
 logFC <= -WEAK_EFFECT_THRESHOLD,
 effect_direction = case_when(
 logFC >= STRONG_EFFECT_THRESHOLD ~ "strong_gain",
 logFC <= -STRONG_EFFECT_THRESHOLD ~ "strong_loss",
 logFC >= WEAK_EFFECT_THRESHOLD ~ "weak_gain",
 logFC <= -WEAK_EFFECT_THRESHOLD ~ "weak_loss",
 TRUE ~ "minimal"),
 significance_class = case_when(
 strong_significant ~ "strong_significant",
 weak_significant ~ "weak_significant",
 TRUE ~ "not_significant"),
 # Preserve each WT/MUT comparison explicitly.
 comparison_key = paste(
 wt_oligo,
 mut_oligo,
 sep = "__"))
############################################################
# 10. SAVE COMBINED DIFFERENTIAL DATA
############################################################
save_csv(
 differential_all,
 "ALL_WT_MUT_differential_results.csv")
write.csv(
 differential_all,
 file.path(
 DATA_DIR,
 "ALL_WT_MUT_differential_results.csv"),
 row.names = FALSE)
############################################################
# 11. SUMMARY BY CELL LINE
############################################################
cat("\n============================================================\n")
cat("Summary of strong significant effects\n")
cat("============================================================\n")
cell_summary <- differential_all %>%
 group_by(cell_line) %>%
 summarise(
 total_comparisons = n(),
 strong_significant =
 sum(strong_significant, na.rm = TRUE),
 strong_gain =
 sum(strong_gain, na.rm = TRUE),
 strong_loss =
 sum(strong_loss, na.rm = TRUE),
 weak_significant =
 sum(weak_significant, na.rm = TRUE),
 .groups = "drop")
print(cell_summary)
save_csv(
 cell_summary,
 "cell_line_differential_summary.csv")
############################################################
# 12. STRONG SIGNIFICANT RESULTS BY CELL LINE
############################################################
strong_all <- differential_all %>%
 filter(
 strong_significant
) %>%
 arrange(
 cell_line,
 adj.P.Val)
save_csv(
 strong_all,
 "strong_significant_WT_MUT_comparisons.csv")
for (cell in CELL_LINES) {
 cell_results <- strong_all %>%
 filter(
 cell_line == cell)
 save_csv(
 cell_results,
 paste0(
 cell,
 "_strong_significant_WT_MUT.csv"))}
############################################################
# 13. ELEMENT-LEVEL ACTIVITY RESULTS
#
# These are separate from WT/MUT differential effects.
############################################################
cat("\n============================================================\n")
cat("Loading element-level activity results\n")
cat("============================================================\n")
activity_results_file <- file.path(
 CHECKPOINT_DIR,
 "04_activity_analysis.rds")
if (file.exists(activity_results_file)) {
 activity_results <- readRDS(
 activity_results_file)
 activity_all <- bind_rows(
 activity_results[CELL_LINES],
 .id = "cell_line")
 cat(
 "\nElement-level activity results loaded:",
 nrow(activity_all),
 "rows\n")
} else {
 cat(
 "\nWARNING: Element-level activity checkpoint not found.\n")
 activity_results <- NULL
 activity_all <- NULL}
############################################################
# 14. ELEMENT-LEVEL ACTIVITY DISTRIBUTION
############################################################
if (!is.null(activity_all)) {
 cat("\n============================================================\n")
 cat("Element-level activity visualizations\n")
 cat("============================================================\n")
 p_activity_distribution <- ggplot(
 activity_all %>%
 filter(
 is.finite(logFC)),
 aes(
 x = logFC)
) +
 geom_histogram(
 bins = 80
) +
 geom_vline(
 xintercept = 0,
 linetype = "dashed"
) +
 theme_bw() +
 facet_wrap(
 ~ cell_line,
 scales = "free_y"
) +
 labs(
 title = "Element-level BCalm activity",
 x = "BCalm log2 activity",
 y = "Number of oligos")
 save_plot(
 p_activity_distribution,
 file.path(
 ACTIVITY_DIR,
 "01_element_level_activity_distribution.png"),
 width = 12,
 height = 8)
 ##########################################################
 # Activity by analysis group
 ##########################################################
 if ("analysis_group" %in% colnames(activity_all)) {
 p_activity_groups <- ggplot(
 activity_all %>%
 filter(
 !is.na(logFC),
 !is.na(analysis_group)),
 aes(
 x = analysis_group,
 y = logFC)
) +
 geom_boxplot(
 outlier.alpha = 0.15
) +
 theme_bw() +
 facet_wrap(
 ~ cell_line,
 scales = "free_y"
) +
 theme(
 axis.text.x = element_text(
 angle = 45,
 hjust = 1)
) +
 labs(
 title = "Element-level activity by analysis group",
 x = "Analysis group",
 y = "BCalm log2 activity")
 save_plot(
 p_activity_groups,
 file.path(
 ACTIVITY_DIR,
 "02_element_activity_by_group.png"),
 width = 14,
 height = 8)}
 ##########################################################
 # Save activity data
 ##########################################################
 write.csv(
 activity_all,
 file.path(
 DATA_DIR,
 "ALL_element_level_activity.csv"),
 row.names = FALSE)}
############################################################
# 15. WT/MUT VOLCANO PLOTS
############################################################
cat("\n============================================================\n")
cat("WT/MUT volcano plots\n")
cat("============================================================\n")
for (cell in CELL_LINES) {
 cell_data <- differential_all %>%
 filter(
 cell_line == cell,
 !is.na(logFC),
 !is.na(P.Value),
 P.Value > 0
) %>%
 mutate(
 plot_class = case_when(
 strong_gain ~ "Strong gain",
 strong_loss ~ "Strong loss",
 weak_gain ~ "Weak gain",
 weak_loss ~ "Weak loss",
 TRUE ~ "Other"))
 p <- ggplot(
 cell_data,
 aes(
 x = logFC,
 y = -log10(P.Value))
) +
 geom_point(
 alpha = 0.25,
 size = 1
) +
 geom_vline(
 xintercept = c(
 -STRONG_EFFECT_THRESHOLD,
 STRONG_EFFECT_THRESHOLD),
 linetype = "dashed"
) +
 geom_hline(
 yintercept = -log10(FDR_THRESHOLD),
 linetype = "dashed"
) +
 theme_bw() +
 labs(
 title = paste(
 cell,
 "WT vs MUT differential analysis"),
 subtitle = "Dashed vertical lines = strong effect threshold; horizontal line = nominal FDR threshold",
 x = "WT → MUT log2 fold change",
 y = "-log10(p-value)")
 save_plot(
 p,
 file.path(
 DIFF_DIR,
 paste0(
 cell,
 "_WT_MUT_volcano.png")))}
############################################################
# 16. WT/MUT MA PLOTS
############################################################
cat("\n============================================================\n")
cat("WT/MUT MA plots\n")
cat("============================================================\n")
for (cell in CELL_LINES) {
 cell_data <- differential_all %>%
 filter(
 cell_line == cell,
 is.finite(logFC),
 is.finite(AveExpr)
) %>%
 mutate(
 strong_significant =
 !is.na(adj.P.Val) &
 adj.P.Val <= FDR_THRESHOLD &
 abs(logFC) >= STRONG_EFFECT_THRESHOLD)
 p <- ggplot(
 cell_data,
 aes(
 x = AveExpr,
 y = logFC)
) +
 geom_point(
 aes(
 alpha = strong_significant),
 size = 1
) +
 geom_hline(
 yintercept = 0,
 linetype = "dashed"
) +
 theme_bw() +
 labs(
 title = paste(
 cell,
 "WT vs MUT MA plot"),
 x = "Average expression",
 y = "WT → MUT log2 fold change",
 alpha = "Strong significant")
 save_plot(
 p,
 file.path(
 DIFF_DIR,
 paste0(
 cell,
 "_WT_MUT_MA_plot.png")))}
############################################################
# 17. WT/MUT LOGFC DISTRIBUTIONS
############################################################
for (cell in CELL_LINES) {
 cell_data <- differential_all %>%
 filter(
 cell_line == cell,
 is.finite(logFC))
 p <- ggplot(
 cell_data,
 aes(
 x = logFC)
) +
 geom_histogram(
 bins = 80
) +
 geom_vline(
 xintercept = 0,
 linetype = "dashed"
) +
 geom_vline(
 xintercept = c(
 -STRONG_EFFECT_THRESHOLD,
 STRONG_EFFECT_THRESHOLD),
 linetype = "dotted"
) +
 theme_bw() +
 labs(
 title = paste(
 cell,
 "WT → MUT effect distribution"),
 x = "WT → MUT log2 fold change",
 y = "Number of comparisons")
 save_plot(
 p,
 file.path(
 DIFF_DIR,
 paste0(
 cell,
 "_WT_MUT_logFC_distribution.png")))}
############################################################
# 18. TOP DIFFERENTIAL COMPARISONS
############################################################
for (cell in CELL_LINES) {
 top_cell <- differential_all %>%
 filter(
 cell_line == cell,
 !is.na(adj.P.Val)
) %>%
 arrange(
 adj.P.Val
) %>%
 slice_head(
 n = 30)
 if (nrow(top_cell) == 0) {
 next}
 top_cell$comparison_label <- paste(
 top_cell$wt_oligo,
 "→",
 top_cell$mut_oligo)
 top_cell$comparison_label <- factor(
 top_cell$comparison_label,
 levels = rev(
 top_cell$comparison_label))
 p <- ggplot(
 top_cell,
 aes(
 x = logFC,
 y = comparison_label)
) +
 geom_point() +
 geom_vline(
 xintercept = 0,
 linetype = "dashed"
) +
 theme_bw() +
 labs(
 title = paste(
 cell,
 "top WT/MUT differential comparisons"),
 x = "WT → MUT log2 fold change",
 y = "WT → MUT comparison")
 save_plot(
 p,
 file.path(
 DIFF_DIR,
 paste0(
 cell,
 "_top_WT_MUT_comparisons.png")),
 width = 10,
 height = 10)}
############################################################
# SECTION 18A: TOP 30 STRONG WT/MUT ACTIVITY HEATMAPS
#
# PURPOSE
#
# Create WT vs MUT activity heatmaps for the top 30
# strong significant gains and losses in each replicated
# cell line.
#
# Ranking is based on the formal Script 2 differential
# analysis:
#
# FDR <= 0.05
# abs(logFC) >= 1
#
# Activity values are taken directly from the Script 1
# replicate_activity object.
#
# For replicated cell lines, the heatmap displays the
# mean replicate activity across the retained biological
# replicates.
#
# IMPORTANT
#
# The heatmap does NOT recalculate BCalm activity from
# raw DNA/RNA counts.
#
# The formal differential ranking remains:
#
# logFC = Script 2 BCalm WT -> MUT effect
#
# The heatmap itself shows:
#
# WT activity
# MUT activity
#
# using the existing Script 1 replicate_activity values.
#
############################################################
cat(
 "\n============================================================\n")
cat(
 "SECTION 18A: TOP 30 WT/MUT ACTIVITY HEATMAPS\n")
cat(
 "============================================================\n")
############################################################
# Heatmap output directory
############################################################
WT_MUT_HEATMAP_DIR <- file.path(
 DIFF_DIR,
 "WT_MUT_activity_heatmaps")
dir.create(
 WT_MUT_HEATMAP_DIR,
 recursive = TRUE,
 showWarnings = FALSE)
############################################################
# Helper function:
# Extract mean BCalm replicate activity for WT/MUT oligos
############################################################
get_mean_replicate_activity <- function(
 env,
 oligo_names
) {
 if (!exists(
 "replicate_activity",
 envir = env
)) {
 stop(
 "replicate_activity not found in Script 1 workspace.")}
 activity_data <- get(
 "replicate_activity",
 envir = env)
 # ----------------------------------------------------------
 # Identify replicate activity columns
 # ----------------------------------------------------------
 activity_columns <- grep(
 "^rep[0-9]+_activity$",
 colnames(activity_data),
 value = TRUE)
 if (length(activity_columns) == 0) {
 stop(
 "No rep*_activity columns found in replicate_activity.")}
 # ----------------------------------------------------------
 # Keep only requested WT/MUT oligos
 # ----------------------------------------------------------
 activity_data <- activity_data %>%
 filter(
 oligo_name %in% oligo_names)
 # ----------------------------------------------------------
 # Calculate mean across retained biological replicates
 #
 # This uses the existing Script 1 activity values.
 # It does NOT recalculate activity from DNA/RNA counts.
 # ----------------------------------------------------------
 activity_data %>%
 mutate(
 mean_activity = rowMeans(
 across(
 all_of(activity_columns)),
 na.rm = TRUE)
) %>%
 select(
 oligo_name,
 mean_activity)}
SCRIPT1_BASE <- "BCalm_activity"
############################################################
# Loop through replicated cell lines
############################################################
# ------------------------------------------------------------
# 18A.1 Load Script 1 workspaces
# ------------------------------------------------------------
load_script1_workspace_heatmap <- function(cell) {
 workspace_file <- switch(
 cell,
 HEK = file.path(
 SCRIPT1_BASE,
 "HEK",
 "HEK_BCalm_activity_workspace.RData"),
 U2OS = file.path(
 SCRIPT1_BASE,
 "U2OS",
 "U2OS_BCalm_activity_workspace.RData"),
 U2987 = file.path(
 SCRIPT1_BASE,
 "U2987",
 "U2987_BCalm_activity_workspace.RData"),
 HepG2 = file.path(
 SCRIPT1_BASE,
 "HepG2",
 "HepG2_BCalm_activity_workspace.RData"),
 `SK-MEL30` = file.path(
 SCRIPT1_BASE,
 "SK-MEL30",
 "SK-MEL30_BCalm_activity_workspace.RData"),
 stop("Unknown cell line: ", cell))
 if (!file.exists(workspace_file)) {
 stop(
 "Script 1 workspace not found for ",
 cell,
 ": ",
 workspace_file)}
 env <- new.env(parent = emptyenv())
 load(workspace_file, envir = env)
 env}
for (cell in CELL_LINES) {
 cat(
 "\n------------------------------------------------------------\n")
 cat(
 "Creating WT/MUT activity heatmaps:",
 cell,
 "\n")
 cat(
 "------------------------------------------------------------\n")
 ##########################################################
 # Load Script 1 workspace
 ##########################################################
 env <- load_script1_workspace_heatmap(
 cell)
 ##########################################################
 # Identify top strong gains
 ##########################################################
 top30_gain <- differential_all %>%
 filter(
 cell_line == cell,
 strong_gain
) %>%
 arrange(
 desc(logFC),
 adj.P.Val
) %>%
 slice_head(
 n = 30)
 ##########################################################
 # Identify top strong losses
 ##########################################################
 top30_loss <- differential_all %>%
 filter(
 cell_line == cell,
 strong_loss
) %>%
 arrange(
 logFC,
 adj.P.Val
) %>%
 slice_head(
 n = 30)
 ##########################################################
 # Function to create one heatmap
 ##########################################################
 create_activity_heatmap <- function(
 candidate_table,
 direction_label,
 filename
) {
 if (
 nrow(candidate_table) == 0
) {
 cat(
 "\nNo strong significant",
 direction_label,
 "candidates for",
 cell,
 "\n")
 return(
 invisible(NULL))}
 ########################################################
 # Get WT and MUT oligo names
 ########################################################
 oligo_names <- unique(
 c(
 candidate_table$wt_oligo,
 candidate_table$mut_oligo))
 ########################################################
 # Get existing Script 1 BCalm activity
 ########################################################
 activity_data <- get_mean_replicate_activity(
 env,
 oligo_names)
 ########################################################
 # Add WT and MUT activity
 ########################################################
 heatmap_data <- candidate_table %>%
 select(
 pair_id,
 wt_oligo,
 mut_oligo,
 logFC,
 adj.P.Val
) %>%
 left_join(
 activity_data %>%
 dplyr::rename(
 wt_oligo = oligo_name,
 wt_activity = mean_activity),
 by = "wt_oligo"
) %>%
 left_join(
 activity_data %>%
 dplyr::rename(
 mut_oligo = oligo_name,
 mut_activity = mean_activity),
 by = "mut_oligo")
 ########################################################
 # Check for missing activity
 ########################################################
 missing_activity <- heatmap_data %>%
 filter(
 is.na(wt_activity) |
 is.na(mut_activity))
 if (nrow(missing_activity) > 0) {
 warning(
 "\n",
 cell,
 " | ",
 direction_label,
 ": ",
 nrow(missing_activity),
 " WT/MUT comparisons have missing activity values.")}
 ########################################################
 # Calculate observed activity difference
 #
 # This is descriptive:
 #
 # MUT activity - WT activity
 #
 # The formal statistical effect remains Script 2 logFC.
 ########################################################
 heatmap_data <- heatmap_data %>%
 mutate(
 observed_activity_difference =
 mut_activity -
 wt_activity)
 ########################################################
 # Save heatmap data
 ########################################################
 data_filename <- sub(
 "\\.png$",
 "_data.csv",
 filename)
 write.csv(
 heatmap_data,
 file.path(
 WT_MUT_HEATMAP_DIR,
 data_filename),
 row.names = FALSE)
 ########################################################
 # Convert to long format
 ########################################################
 plot_data <- heatmap_data %>%
 select(
 pair_id,
 wt_activity,
 mut_activity
) %>%
 pivot_longer(
 cols = c(
 wt_activity,
 mut_activity),
 names_to = "allele",
 values_to = "activity"
) %>%
 mutate(
 allele = recode(
 allele,
 wt_activity = "WT",
 mut_activity = "MUT"),
 pair_id = factor(
 pair_id,
 levels = rev(
 heatmap_data$pair_id)))
 ########################################################
 # Create heatmap
 ########################################################
 p <- ggplot(
 plot_data,
 aes(
 x = allele,
 y = pair_id,
 fill = activity)
) +
 geom_tile(
 color = "white",
 linewidth = 0.4
) +
 geom_text(
 aes(
 label = ifelse(
 is.na(activity),
 "",
 sprintf(
 "%.2f",
 activity))),
 size = 3
) +
 scale_fill_gradient2(
 midpoint = 0,
 name = "Activity"
) +
 labs(
 title = paste(
 cell,
 "Top 30 Strong",
 direction_label),
 subtitle = paste(
 "BCalm differential analysis: FDR <= 0.05, |logFC| >= 1;",
 "minimum 10 assigned barcodes"),
 x = NULL,
 y = "WT/MUT pair"
) +
 theme_minimal(
 base_size = 11
) +
 theme(
 panel.grid = element_blank(),
 axis.text.y = element_text(
 size = 8),
 plot.title = element_text(
 face = "bold"))
 ########################################################
 # Save plot
 ########################################################
 save_plot(
 p,
 file.path(
 WT_MUT_HEATMAP_DIR,
 filename),
 width = 7,
 height = 10)
 ########################################################
 # Print summary
 ########################################################
 cat(
 "\nCreated:",
 cell,
 "|",
 direction_label,
 "\n")
 cat(
 " Candidates:",
 nrow(candidate_table),
 "\n")
 cat(
 " Heatmap:",
 file.path(
 WT_MUT_HEATMAP_DIR,
 filename),
 "\n")
 invisible(
 heatmap_data)}
 ##########################################################
 # Create gain heatmap
 ##########################################################
 create_activity_heatmap(
 candidate_table = top30_gain,
 direction_label = "Mutation Gains",
 filename = paste0(
 cell,
 "_top30_strong_gains_activity_heatmap_min10.png"))
 ##########################################################
 # Create loss heatmap
 ##########################################################
 create_activity_heatmap(
 candidate_table = top30_loss,
 direction_label = "Mutation Losses",
 filename = paste0(
 cell,
 "_top30_strong_losses_activity_heatmap_min10.png"))}
############################################################
# Summary
############################################################
cat(
 "\n============================================================\n")
cat(
 "WT/MUT activity heatmaps completed.\n")
cat(
 "Output directory:\n")
cat(
 WT_MUT_HEATMAP_DIR,
 "\n")
cat(
 "============================================================\n")
# ============================================================
# 18A. WT/MUT DIFFERENTIAL EFFECT HEATMAPS
#
# Uses the FORMAL BCalm logFC directly.
#
# Strong gain:
# FDR <= 0.05
# logFC >= 1
#
# Strong loss:
# FDR <= 0.05
# logFC <= -1
#
# IMPORTANT:
# These heatmaps show BCalm logFC, NOT Script 1
# descriptive replicate_activity.
# ============================================================
WT_MUT_HEATMAP_DIR <- file.path(
 DIFF_DIR,
 "WT_MUT_BCalm_logFC_heatmaps")
dir.create(
 WT_MUT_HEATMAP_DIR,
 recursive = TRUE,
 showWarnings = FALSE)
# ------------------------------------------------------------
# Function to create BCalm logFC heatmap
# ------------------------------------------------------------
create_bcalm_logfc_heatmap <- function(
 data,
 cell_line,
 effect_type,
 output_file
) {
 if (nrow(data) == 0) {
 message(
 "No results available for ",
 cell_line,
 " ",
 effect_type)
 return(NULL)}
 # Order strongest effect at top
 if (effect_type == "strong_gain") {
 data <- data %>%
 arrange(desc(logFC))
 } else {
 data <- data %>%
 arrange(logFC)}
 # Create readable row labels
 data <- data %>%
 mutate(
 label = paste0(
 pair_id,
 "\nFDR = ",
 formatC(
 adj.P.Val,
 format = "e",
 digits = 2)),
 label = factor(
 label,
 levels = rev(label)))
 # ----------------------------------------------------------
 # Heatmap
 #
 # Tile color = FORMAL BCalm logFC
 # Text = BCalm logFC
 # ----------------------------------------------------------
 p <- ggplot(
 data,
 aes(
 x = "BCalm logFC",
 y = label,
 fill = logFC)
) +
 geom_tile(
 color = "white",
 linewidth = 0.4
) +
 geom_text(
 aes(
 label = sprintf("%.2f", logFC)),
 size = 3.5
) +
 scale_fill_gradient2(
 midpoint = 0,
 name = "BCalm logFC"
) +
 labs(
 title = paste0(
 cell_line,
 ": ",
 ifelse(
 effect_type == "strong_gain",
 "Top 30 strong gains",
 "Top 30 strong losses")),
 subtitle = paste(
 "Formal BCalm WT→MUT differential effect",
 "| FDR ≤ 0.05 | |logFC| ≥ 1"),
 x = NULL,
 y = NULL
) +
 theme_minimal(base_size = 11) +
 theme(
 axis.text.x = element_text(
 size = 10),
 axis.text.y = element_text(
 size = 8),
 panel.grid = element_blank(),
 plot.title = element_text(
 face = "bold"))
 ggsave(
 output_file,
 p,
 width = 7,
 height = max(8, nrow(data) * 0.28),
 dpi = 300,
 bg = "white")
 return(p)}
# ------------------------------------------------------------
# Generate heatmaps for every replicated cell line
# ------------------------------------------------------------
for (cell in CELL_LINES) {
 message(
 "\nCreating BCalm logFC heatmaps: ",
 cell)
 # ----------------------------------------------------------
 # Strong gains
 # ----------------------------------------------------------
 top30_gain <- differential_all %>%
 filter(
 cell_line == cell,
 strong_gain
) %>%
 arrange(
 desc(logFC),
 adj.P.Val
) %>%
 slice_head(n = 30)
 # Save underlying data
 write.csv(
 top30_gain,
 file.path(
 WT_MUT_HEATMAP_DIR,
 paste0(
 cell,
 "_top30_strong_gains_BCalm_logFC_data.csv")),
 row.names = FALSE)
 create_bcalm_logfc_heatmap(
 data = top30_gain,
 cell_line = cell,
 effect_type = "strong_gain",
 output_file = file.path(
 WT_MUT_HEATMAP_DIR,
 paste0(
 cell,
 "_top30_strong_gains_BCalm_logFC_heatmap.png")))
 # ----------------------------------------------------------
 # Strong losses
 # ----------------------------------------------------------
 top30_loss <- differential_all %>%
 filter(
 cell_line == cell,
 strong_loss
) %>%
 arrange(
 logFC,
 adj.P.Val
) %>%
 slice_head(n = 30)
 # Save underlying data
 write.csv(
 top30_loss,
 file.path(
 WT_MUT_HEATMAP_DIR,
 paste0(
 cell,
 "_top30_strong_losses_BCalm_logFC_data.csv")),
 row.names = FALSE)
 create_bcalm_logfc_heatmap(
 data = top30_loss,
 cell_line = cell,
 effect_type = "strong_loss",
 output_file = file.path(
 WT_MUT_HEATMAP_DIR,
 paste0(
 cell,
 "_top30_strong_losses_BCalm_logFC_heatmap.png")))}
message(
 "\nBCalm logFC heatmaps complete.\n",
 "Output directory:\n",
 WT_MUT_HEATMAP_DIR)
############################################################
# 19. CROSS-CELL-LINE SUMMARY
############################################################
cat("\n============================================================\n")
cat("Cross-cell-line NCCM analysis\n")
cat("============================================================\n")
############################################################
# IMPORTANT:
#
# Each WT/MUT comparison is preserved separately.
#
# We use comparison_key rather than pair_id alone so that
# multiple MUT alleles associated with the same WT are not
# accidentally collapsed.
############################################################
cross_summary <- differential_all %>%
 group_by(
 comparison_key,
 pair_id,
 wt_oligo,
 mut_oligo
) %>%
 summarise(
 n_cell_lines_tested =
 n_distinct(cell_line),
 n_cell_lines_strong =
 sum(
 strong_significant,
 na.rm = TRUE),
 n_cell_lines_weak =
 sum(
 weak_significant,
 na.rm = TRUE),
 strong_cell_lines =
 paste(
 cell_line[
 strong_significant
 ],
 collapse = ", "),
 strong_directions =
 paste(
 paste(
 cell_line[
 strong_significant
 ],
 ifelse(
 logFC[
 strong_significant
 ] > 0,
 "gain",
 "loss"),
 sep = ":"),
 collapse = ", "),
 mean_strong_logFC =
 ifelse(
 any(strong_significant),
 mean(
 logFC[
 strong_significant
 ],
 na.rm = TRUE),
 NA_real_),
 max_abs_strong_logFC =
 ifelse(
 any(strong_significant),
 max(
 abs(
 logFC[
 strong_significant
 ]),
 na.rm = TRUE),
 NA_real_),
 min_strong_FDR =
 ifelse(
 any(strong_significant),
 min(
 adj.P.Val[
 strong_significant
 ],
 na.rm = TRUE),
 NA_real_),
 .groups = "drop")
############################################################
# 20. SHARED STRONG NCCMs
############################################################
shared_strong <- cross_summary %>%
 filter(
 n_cell_lines_strong >= 2
) %>%
 arrange(
 desc(n_cell_lines_strong),
 min_strong_FDR,
 desc(max_abs_strong_logFC))
save_csv(
 shared_strong,
 "shared_strong_significant_NCCMs_2plus_cell_lines.csv")
############################################################
# 21. STRONG NCCMs IN 3+ CELL LINES
############################################################
shared_3plus <- cross_summary %>%
 filter(
 n_cell_lines_strong >= 3
) %>%
 arrange(
 desc(n_cell_lines_strong),
 min_strong_FDR)
save_csv(
 shared_3plus,
 "shared_strong_significant_NCCMs_3plus_cell_lines.csv")
############################################################
# 22. STRONG NCCMs IN 4+ CELL LINES
############################################################
shared_4plus <- cross_summary %>%
 filter(
 n_cell_lines_strong >= 4
) %>%
 arrange(
 desc(n_cell_lines_strong),
 min_strong_FDR)
save_csv(
 shared_4plus,
 "shared_strong_significant_NCCMs_4plus_cell_lines.csv")
############################################################
# 23. STRONG NCCMs IN ALL FIVE CELL LINES
############################################################
shared_all <- cross_summary %>%
 filter(
 n_cell_lines_strong == length(CELL_LINES)
) %>%
 arrange(
 min_strong_FDR,
 desc(max_abs_strong_logFC))
save_csv(
 shared_all,
 "strong_significant_NCCMs_ALL_5_cell_lines.csv")
############################################################
# 24. CELL-LINE-SPECIFIC STRONG NCCMs
############################################################
specific_summary <- cross_summary %>%
 filter(
 n_cell_lines_strong == 1
) %>%
 arrange(
 desc(max_abs_strong_logFC),
 min_strong_FDR)
save_csv(
 specific_summary,
 "cell_line_specific_strong_NCCMs.csv")
############################################################
# 25. CELL-LINE-SPECIFIC GAIN/LOSS TABLES
############################################################
specific_comparisons <- differential_all %>%
 filter(
 strong_significant
) %>%
 group_by(
 comparison_key
) %>%
 mutate(
 n_cell_lines_strong =
 sum(strong_significant, na.rm = TRUE)
) %>%
 ungroup() %>%
 filter(
 n_cell_lines_strong == 1)
specific_gains <- specific_comparisons %>%
 filter(
 strong_gain
) %>%
 arrange(
 cell_line,
 adj.P.Val)
specific_losses <- specific_comparisons %>%
 filter(
 strong_loss
) %>%
 arrange(
 cell_line,
 adj.P.Val)
save_csv(
 specific_gains,
 "cell_line_specific_strong_gains.csv")
save_csv(
 specific_losses,
 "cell_line_specific_strong_losses.csv")
############################################################
# 26. CONSISTENT-DIRECTION SHARED NCCMs
############################################################
cat("\n============================================================\n")
cat("Direction consistency\n")
cat("============================================================\n")
direction_summary <- differential_all %>%
 filter(
 strong_significant
) %>%
 group_by(
 comparison_key,
 pair_id,
 wt_oligo,
 mut_oligo
) %>%
 summarise(
 n_strong =
 n(),
 n_gain =
 sum(
 logFC > 0,
 na.rm = TRUE),
 n_loss =
 sum(
 logFC < 0,
 na.rm = TRUE),
 cell_lines =
 paste(
 cell_line,
 collapse = ", "),
 .groups = "drop"
) %>%
 mutate(
 direction_consistency = case_when(
 n_gain == n_strong ~ "consistent_gain",
 n_loss == n_strong ~ "consistent_loss",
 TRUE ~ "inconsistent"))
consistent_gain <- direction_summary %>%
 filter(
 direction_consistency == "consistent_gain",
 n_strong >= 2)
consistent_loss <- direction_summary %>%
 filter(
 direction_consistency == "consistent_loss",
 n_strong >= 2)
inconsistent_direction <- direction_summary %>%
 filter(
 direction_consistency == "inconsistent",
 n_strong >= 2)
save_csv(
 consistent_gain,
 "shared_consistent_strong_gains.csv")
save_csv(
 consistent_loss,
 "shared_consistent_strong_losses.csv")
save_csv(
 inconsistent_direction,
 "shared_strong_directionally_inconsistent.csv")
############################################################
# 27. CROSS-CELL-LINE SUMMARY TABLE
############################################################
cross_count_summary <- cross_summary %>%
 dplyr::count(
 n_cell_lines_strong,
 name = "n_NCCMs"
) %>%
 arrange(
 n_cell_lines_strong)
save_csv(
 cross_count_summary,
 "number_of_cell_lines_per_strong_NCCM.csv")
print(cross_count_summary)
############################################################
# 28. CROSS-CELL-LINE COUNT PLOT
############################################################
p_shared_counts <- ggplot(
 cross_count_summary,
 aes(
 x = factor(n_cell_lines_strong),
 y = n_NCCMs)
) +
 geom_col() +
 theme_bw() +
 labs(
 title = "Number of strong significant NCCMs by cell-line recurrence",
 x = "Number of cell lines with strong significant effect",
 y = "Number of NCCMs")
save_plot(
 p_shared_counts,
 file.path(
 CROSS_DIR,
 "01_NCCMs_by_number_of_cell_lines.png"))
############################################################
# 29. CROSS-CELL-LINE EFFECT-SIZE MATRIX
############################################################
cat("\nCreating cross-cell-line effect matrix...\n")
strong_candidates <- cross_summary %>%
 filter(
 n_cell_lines_strong >= 2)
if (nrow(strong_candidates) > 0) {
 effect_long <- differential_all %>%
 semi_join(
 strong_candidates,
 by = "comparison_key"
) %>%
 select(
 comparison_key,
 wt_oligo,
 mut_oligo,
 cell_line,
 logFC)
 effect_wide <- effect_long %>%
 select(
 comparison_key,
 wt_oligo,
 mut_oligo,
 cell_line,
 logFC
) %>%
 pivot_wider(
 names_from = cell_line,
 values_from = logFC)
 write.csv(
 effect_wide,
 file.path(
 TABLE_DIR,
 "shared_NCCM_effect_size_matrix.csv"),
 row.names = FALSE)
 ##########################################################
 # Heatmap
 ##########################################################
 effect_plot_data <- effect_long %>%
 mutate(
 comparison_label = paste(
 wt_oligo,
 "→",
 mut_oligo))
 # Order rows by recurrence and maximum effect.
 candidate_order <- strong_candidates %>%
 arrange(
 desc(n_cell_lines_strong),
 desc(max_abs_strong_logFC)
) %>%
 mutate(
 comparison_label = paste(
 wt_oligo,
 "→",
 mut_oligo)
) %>%
 pull(
 comparison_label)
 effect_plot_data$comparison_label <- factor(
 effect_plot_data$comparison_label,
 levels = rev(candidate_order))
 p_effect_heatmap <- ggplot(
 effect_plot_data,
 aes(
 x = cell_line,
 y = comparison_label,
 fill = logFC)
) +
 geom_tile(
 na.rm = FALSE
) +
 geom_text(
 aes(
 label = ifelse(
 is.na(logFC),
 "",
 sprintf("%.2f", logFC))),
 size = 2.5
) +
 theme_bw() +
 labs(
 title = "Strong significant NCCMs across cell lines",
 subtitle = "WT → MUT log2 fold change",
 x = "Cell line",
 y = "WT → MUT comparison",
 fill = "log2 fold change"
) +
 theme(
 axis.text.y = element_text(
 size = 6),
 axis.text.x = element_text(
 angle = 45,
 hjust = 1))
 save_plot(
 p_effect_heatmap,
 file.path(
 CROSS_DIR,
 "02_shared_NCCM_effect_size_heatmap.png"),
 width = 12,
 height = max(
 8,
 min(
 25,
 length(candidate_order) * 0.18)))}
############################################################
# 30. SIGNIFICANCE / EFFECT CATEGORY HEATMAP
############################################################
category_long <- differential_all %>%
 mutate(
 category = case_when(
 strong_gain ~ "Strong gain",
 strong_loss ~ "Strong loss",
 weak_gain ~ "Weak gain",
 weak_loss ~ "Weak loss",
 TRUE ~ "Not significant / minimal"))
category_candidates <- cross_summary %>%
 filter(
 n_cell_lines_strong >= 2)
if (nrow(category_candidates) > 0) {
 category_plot_data <- category_long %>%
 semi_join(
 category_candidates,
 by = "comparison_key"
) %>%
 mutate(
 comparison_label = paste(
 wt_oligo,
 "→",
 mut_oligo))
 category_order <- category_candidates %>%
 arrange(
 desc(n_cell_lines_strong),
 desc(max_abs_strong_logFC)
) %>%
 mutate(
 comparison_label = paste(
 wt_oligo,
 "→",
 mut_oligo)
) %>%
 pull(
 comparison_label)
 category_plot_data$comparison_label <- factor(
 category_plot_data$comparison_label,
 levels = rev(category_order))
 # Numeric representation for plotting.
 category_plot_data$category_value <- case_when(
 category_plot_data$category == "Strong gain" ~ 2,
 category_plot_data$category == "Weak gain" ~ 1,
 category_plot_data$category == "Not significant / minimal" ~ 0,
 category_plot_data$category == "Weak loss" ~ -1,
 category_plot_data$category == "Strong loss" ~ -2,
 TRUE ~ NA_real_)
 p_category_heatmap <- ggplot(
 category_plot_data,
 aes(
 x = cell_line,
 y = comparison_label,
 fill = category_value)
) +
 geom_tile() +
 theme_bw() +
 scale_fill_gradient2(
 midpoint = 0,
 breaks = c(
 -2,
 -1,
 0,
 1,
 2),
 labels = c(
 "Strong loss",
 "Weak loss",
 "Not significant / minimal",
 "Weak gain",
 "Strong gain")
) +
 labs(
 title = "Cross-cell-line strong/weak effect classification",
 x = "Cell line",
 y = "WT → MUT comparison",
 fill = "Effect class"
) +
 theme(
 axis.text.y = element_text(
 size = 6),
 axis.text.x = element_text(
 angle = 45,
 hjust = 1))
 save_plot(
 p_category_heatmap,
 file.path(
 CROSS_DIR,
 "03_shared_NCCM_effect_category_heatmap.png"),
 width = 12,
 height = max(
 8,
 min(
 25,
 length(category_order) * 0.18)))}
############################################################
# 31. ALL-CELL-LINE NCCM HEATMAP
############################################################
if (nrow(shared_all) > 0) {
 all5_data <- differential_all %>%
 semi_join(
 shared_all,
 by = "comparison_key"
) %>%
 mutate(
 comparison_label = paste(
 wt_oligo,
 "→",
 mut_oligo))
 all5_order <- shared_all %>%
 mutate(
 comparison_label = paste(
 wt_oligo,
 "→",
 mut_oligo)
) %>%
 arrange(
 desc(abs(max_abs_strong_logFC))
) %>%
 pull(
 comparison_label)
 all5_data$comparison_label <- factor(
 all5_data$comparison_label,
 levels = rev(all5_order))
 p_all5 <- ggplot(
 all5_data,
 aes(
 x = cell_line,
 y = comparison_label,
 fill = logFC)
) +
 geom_tile() +
 geom_text(
 aes(
 label = sprintf(
 "%.2f",
 logFC)),
 size = 2.5
) +
 geom_hline(
 yintercept = seq(
 0.5,
 length(all5_order) - 0.5,
 by = 1),
 alpha = 0.1
) +
 theme_bw() +
 labs(
 title = "Strong significant NCCMs in all five replicated cell lines",
 x = "Cell line",
 y = "WT → MUT comparison",
 fill = "log2 fold change"
) +
 theme(
 axis.text.y = element_text(
 size = 7),
 axis.text.x = element_text(
 angle = 45,
 hjust = 1))
 save_plot(
 p_all5,
 file.path(
 CROSS_DIR,
 "04_ALL_5_cell_lines_strong_NCCMs_heatmap.png"),
 width = 10,
 height = max(
 6,
 min(
 20,
 length(all5_order) * 0.3)))}
############################################################
# 32. CELL-LINE-SPECIFIC BAR PLOT
############################################################
specific_counts <- specific_comparisons %>%
 dplyr::count(
 cell_line,
 name = "n_specific_strong")
save_csv(
 specific_counts,
 "cell_line_specific_strong_NCCM_counts.csv")
p_specific_counts <- ggplot(
 specific_counts,
 aes(
 x = cell_line,
 y = n_specific_strong)
) +
 geom_col() +
 theme_bw() +
 labs(
 title = "Cell-line-specific strong significant NCCMs",
 x = "Cell line",
 y = "Number of cell-line-specific strong NCCMs")
save_plot(
 p_specific_counts,
 file.path(
 CROSS_DIR,
 "05_cell_line_specific_strong_NCCM_counts.png"))
############################################################
# 33. TOP CELL-LINE-SPECIFIC CANDIDATES
############################################################
for (cell in CELL_LINES) {
 cell_specific <- specific_comparisons %>%
 filter(
 cell_line == cell
) %>%
 arrange(
 adj.P.Val)
 save_csv(
 cell_specific,
 paste0(
 cell,
 "_cell_line_specific_strong_NCCMs.csv"))
 if (nrow(cell_specific) == 0) {
 next}
 top_specific <- cell_specific %>%
 slice_head(
 n = N_TOP_SPECIFIC
) %>%
 mutate(
 comparison_label = paste(
 wt_oligo,
 "→",
 mut_oligo))
 top_specific$comparison_label <- factor(
 top_specific$comparison_label,
 levels = rev(
 top_specific$comparison_label))
 p <- ggplot(
 top_specific,
 aes(
 x = logFC,
 y = comparison_label)
) +
 geom_point() +
 geom_vline(
 xintercept = 0,
 linetype = "dashed"
) +
 theme_bw() +
 labs(
 title = paste(
 cell,
 "top cell-line-specific strong NCCMs"),
 x = "WT → MUT log2 fold change",
 y = "WT → MUT comparison")
 save_plot(
 p,
 file.path(
 CANDIDATE_DIR,
 paste0(
 cell,
 "_top_cell_line_specific_NCCMs.png")),
 width = 10,
 height = 7)}
############################################################
# 34. TOP SHARED NCCM PLOT
############################################################
if (nrow(shared_strong) > 0) {
 top_shared <- shared_strong %>%
 slice_head(
 n = N_TOP_SHARED)
 top_shared_data <- differential_all %>%
 semi_join(
 top_shared,
 by = "comparison_key"
) %>%
 mutate(
 comparison_label = paste(
 wt_oligo,
 "→",
 mut_oligo))
 shared_order <- top_shared %>%
 mutate(
 comparison_label = paste(
 wt_oligo,
 "→",
 mut_oligo)
) %>%
 pull(
 comparison_label)
 top_shared_data$comparison_label <- factor(
 top_shared_data$comparison_label,
 levels = rev(shared_order))
 p_top_shared <- ggplot(
 top_shared_data,
 aes(
 x = cell_line,
 y = comparison_label,
 fill = logFC)
) +
 geom_tile() +
 geom_text(
 aes(
 label = sprintf(
 "%.2f",
 logFC)),
 size = 3
) +
 theme_bw() +
 labs(
 title = "Top shared strong significant NCCMs",
 x = "Cell line",
 y = "WT → MUT comparison",
 fill = "log2 fold change"
) +
 theme(
 axis.text.x = element_text(
 angle = 45,
 hjust = 1))
 save_plot(
 p_top_shared,
 file.path(
 CANDIDATE_DIR,
 "01_top_shared_strong_NCCMs.png"),
 width = 10,
 height = 7)}
############################################################
# 35. RANDOM GAIN EXAMPLES
############################################################
cat("\n============================================================\n")
cat("Selecting random strong gain examples\n")
cat("============================================================\n")
set.seed(RANDOM_SEED)
n_gain <- min(
 N_RANDOM_EXAMPLES,
 nrow(differential_all %>% filter(strong_gain)))
random_gains <- differential_all %>%
 filter(strong_gain) %>%
 slice_sample(n = n_gain)
save_csv(
 random_gains,
 "random_5_strong_gain_examples.csv")
############################################################
# 36. RANDOM LOSS EXAMPLES
############################################################
set.seed(RANDOM_SEED)
n_loss <- min(
 N_RANDOM_EXAMPLES,
 nrow(differential_all %>% filter(strong_loss)))
random_losses <- differential_all %>%
 filter(strong_loss) %>%
 slice_sample(n = n_loss)
save_csv(
 random_losses,
 "random_5_strong_loss_examples.csv")
# ============================================================
# SECTION 37: THRESHOLD SETTINGS AND HELPER FUNCTIONS
# ============================================================
SCRIPT1_BASE <- "BCalm_activity"
BARCODE_THRESHOLDS <- c(
 ALL = 0,
 min10 = 10,
 min20 = 20)
THRESHOLD_LABELS <- c(
 ALL = "ALL",
 min10 = ">=10",
 min20 = ">=20")
# ------------------------------------------------------------
# Load Script 1 workspace
# ------------------------------------------------------------
load_script1_workspace <- function(cell) {
 workspace_file <- file.path(
 SCRIPT1_BASE,
 cell,
 paste0(
 cell,
 "_BCalm_activity_workspace.RData"))
 env <- new.env()
 load(
 workspace_file,
 envir = env)
 env}
# ------------------------------------------------------------
# Filter barcode table by barcode threshold
#
# The barcode threshold is applied to the number of distinct
# barcodes supporting each oligo.
#
# ALL = no threshold
# >=10 = at least 10 distinct barcodes
# >=20 = at least 20 distinct barcodes
# ------------------------------------------------------------
filter_bc_threshold <- function(
 bc,
 threshold
) {
 barcode_support <- bc %>%
 group_by(
 oligo_name
) %>%
 summarise(
 n_barcodes = n_distinct(barcode),
 .groups = "drop")
 if (threshold == 0) {
 keep_oligos <- unique(
 bc$oligo_name)
 } else {
 keep_oligos <- barcode_support %>%
 filter(
 n_barcodes >= threshold
) %>%
 pull(
 oligo_name)}
 bc_filtered <- bc %>%
 filter(
 oligo_name %in% keep_oligos)
 list(
 bc = bc_filtered,
 barcode_support = barcode_support,
 keep_oligos = keep_oligos)}
# ------------------------------------------------------------
# Calculate threshold summary
# ------------------------------------------------------------
calculate_threshold_summary <- function(
 threshold_counts
) {
 all_n <- as.numeric(
 threshold_counts["ALL"])
 min10_n <- as.numeric(
 threshold_counts["min10"])
 min20_n <- as.numeric(
 threshold_counts["min20"])
 data.frame(
 all_n = all_n,
 min10_n = min10_n,
 min20_n = min20_n,
 lost_from_all = all_n - min10_n,
 retained_percent_of_all =
 ifelse(
 all_n > 0,
 100 * min10_n / all_n,
 NA_real_),
 lost_percent_of_all =
 ifelse(
 all_n > 0,
 100 * (all_n - min10_n) / all_n,
 NA_real_),
 lost_from_previous =
 min10_n - min20_n,
 lost_percent_from_previous =
 ifelse(
 min10_n > 0,
 100 * (min10_n - min20_n) / min10_n,
 NA_real_),
 retained_percent_from_min10 =
 ifelse(
 min10_n > 0,
 100 * min20_n / min10_n,
 NA_real_))}
# ------------------------------------------------------------
# Safe filename helper
# ------------------------------------------------------------
safe_filename <- function(x) {
 x <- gsub(
 "[^A-Za-z0-9_.-]",
 "_",
 x)
 x}
# ------------------------------------------------------------
# Calculate mean of upper triangle of correlation matrix
# ------------------------------------------------------------
mean_upper_triangle <- function(
 correlation_matrix
) {
 if (
 is.null(correlation_matrix) ||
 !is.matrix(correlation_matrix)
) {
 return(NA_real_)}
 values <- correlation_matrix[
 upper.tri(
 correlation_matrix)
 ]
 if (length(values) == 0) {
 return(NA_real_)}
 mean(
 values,
 na.rm = TRUE)}
# ------------------------------------------------------------
# Get BCalm oligo-level replicate activity
#
# IMPORTANT:
#
# Activity is NOT recalculated from the raw DNA/RNA counts.
#
# Script 1 / BCalm has already calculated the activity values
# in replicate_activity.
#
# This function retrieves those existing BCalm values and
# applies the threshold only as an oligo-level filter.
# ------------------------------------------------------------
get_bcalm_activity_by_threshold <- function(
 env,
 keep_oligos
) {
 if (!exists(
 "replicate_activity",
 envir = env
)) {
 stop(
 "replicate_activity not found in Script 1 workspace.")}
 activity_data <- get(
 "replicate_activity",
 envir = env)
 activity_data %>%
 filter(
 oligo_name %in% keep_oligos)}
# ============================================================
# SECTION 38: THRESHOLD QC AND CORRELATION ANALYSIS
# ============================================================
cat(
 "\n============================================================\n")
cat(
 "SECTION 38: THRESHOLD QC AND CORRELATION ANALYSIS\n")
cat(
 "============================================================\n")
script3_threshold_results <- list()
threshold_retention_results <- list()
for (cell in CELL_LINES) {
 cat(
 "\n\n============================================================\n")
 cat(
 "CELL LINE:",
 cell,
 "\n")
 cat(
 "============================================================\n")
 # ----------------------------------------------------------
 # Load Script 1 workspace
 # ----------------------------------------------------------
 env <- load_script1_workspace(
 cell)
 bc_original <- get(
 "bc",
 envir = env)
 # ----------------------------------------------------------
 # Find DNA and RNA columns
 # ----------------------------------------------------------
 dna_columns <- grep(
 "^dna_count_[0-9]+$",
 colnames(bc_original),
 value = TRUE)
 rna_columns <- grep(
 "^rna_count_[0-9]+$",
 colnames(bc_original),
 value = TRUE)
 # ----------------------------------------------------------
 # Determine actual replicate numbers
 #
 # Example:
 # HEK -> dna_count_1, dna_count_3
 #
 # Do not assume the replicates are sequential.
 # ----------------------------------------------------------
 dna_replicates <- sub(
 "^dna_count_",
 "",
 dna_columns)
 rna_replicates <- sub(
 "^rna_count_",
 "",
 rna_columns)
 common_replicates <- intersect(
 dna_replicates,
 rna_replicates)
 if (length(common_replicates) == 0) {
 warning(
 "No matching DNA/RNA replicate columns found for ",
 cell)
 next}
 # Keep the actual replicate ordering from the DNA/RNA data.
 common_replicates <- sort(
 as.integer(
 common_replicates))
 dna_columns <- paste0(
 "dna_count_",
 common_replicates)
 rna_columns <- paste0(
 "rna_count_",
 common_replicates)
 n_activity_reps <- length(
 common_replicates)
 cat(
 "\nReplicates detected:",
 paste(
 common_replicates,
 collapse = ", "),
 "\n")
 # ----------------------------------------------------------
 # Check BCalm activity object
 # ----------------------------------------------------------
 if (!exists(
 "replicate_activity",
 envir = env
)) {
 warning(
 "replicate_activity not found for ",
 cell,
 ". Skipping cell line.")
 next}
 bcalm_activity <- get(
 "replicate_activity",
 envir = env)
 cat(
 "BCalm activity rows:",
 nrow(bcalm_activity),
 "\n")
 # ----------------------------------------------------------
 # Check expected activity columns
 # ----------------------------------------------------------
 activity_columns <- paste0(
 "rep",
 common_replicates,
 "_activity")
 barcode_activity_columns <- paste0(
 "rep",
 common_replicates,
 "_barcodes")
 missing_activity_columns <- setdiff(
 activity_columns,
 colnames(bcalm_activity))
 if (length(missing_activity_columns) > 0) {
 warning(
 "Missing BCalm activity columns for ",
 cell,
 ": ",
 paste(
 missing_activity_columns,
 collapse = ", "))
 next}
 # ----------------------------------------------------------
 # Output directory
 # ----------------------------------------------------------
 cell_output_dir <- file.path(
 SCRIPT2_DIR,
 "Script3_visualizations",
 cell)
 dir.create(
 cell_output_dir,
 recursive = TRUE,
 showWarnings = FALSE)
 # ----------------------------------------------------------
 # Store results for this cell line
 # ----------------------------------------------------------
 cell_threshold_results <- list()
 threshold_counts <- c(
 ALL = NA_integer_,
 min10 = NA_integer_,
 min20 = NA_integer_)
 # ==========================================================
 # LOOP THROUGH ALL / >=10 / >=20
 # ==========================================================
 for (
 threshold_name in names(BARCODE_THRESHOLDS)
) {
 threshold_value <- BARCODE_THRESHOLDS[
 threshold_name
 ]
 threshold_label <- THRESHOLD_LABELS[
 threshold_name
 ]
 cat(
 "\n------------------------------------------------------------\n")
 cat(
 cell,
 "|",
 threshold_label,
 "\n")
 cat(
 "------------------------------------------------------------\n")
 # --------------------------------------------------------
 # Apply barcode threshold
 # --------------------------------------------------------
 threshold_result <- filter_bc_threshold(
 bc_original,
 threshold_value)
 bc <- threshold_result$bc
 keep_oligos <- threshold_result$keep_oligos
 threshold_counts[
 threshold_name
 ] <- length(
 keep_oligos)
 # --------------------------------------------------------
 # QC input check
 # --------------------------------------------------------
 cat(
 "\nQC input check:",
 cell,
 "|",
 threshold_label,
 "\n")
 cat(
 "Barcode rows:",
 nrow(bc),
 "\n")
 cat(
 "Distinct oligos:",
 dplyr::n_distinct(
 bc$oligo_name),
 "\n")
 cat(
 "Distinct barcodes:",
 dplyr::n_distinct(
 bc$barcode),
 "\n")
 cat(
 "Oligos retained:",
 length(keep_oligos),
 "\n")
 # --------------------------------------------------------
 # Barcode count distribution
 # --------------------------------------------------------
 barcode_support_filtered <- bc %>%
 group_by(
 oligo_name
) %>%
 summarise(
 n_barcodes = n_distinct(
 barcode),
 .groups = "drop")
 barcode_count_summary <- data.frame(
 cell_line = cell,
 threshold = threshold_label,
 n_oligos = nrow(
 barcode_support_filtered),
 median_barcodes = median(
 barcode_support_filtered$n_barcodes),
 mean_barcodes = mean(
 barcode_support_filtered$n_barcodes),
 min_barcodes = min(
 barcode_support_filtered$n_barcodes),
 max_barcodes = max(
 barcode_support_filtered$n_barcodes))
 # --------------------------------------------------------
 # BCalm activity values
 #
 # IMPORTANT:
 #
 # These are the values already calculated by BCalm.
 # We do NOT reconstruct activity from DNA/RNA counts.
 # --------------------------------------------------------
 activity_oligo <- get_bcalm_activity_by_threshold(
 env,
 keep_oligos)
 # --------------------------------------------------------
 # Store basic threshold results
 # --------------------------------------------------------
 cell_threshold_results[[threshold_name]] <- list(
 threshold = threshold_label,
 threshold_value = threshold_value,
 n_oligos = length(keep_oligos),
 n_barcode_rows = nrow(bc),
 barcode_support = barcode_support_filtered,
 barcode_count_summary = barcode_count_summary,
 bcalm_activity = activity_oligo)
 # ========================================================
 # BARCODE-LEVEL CORRELATIONS
 # ========================================================
 # --------------------------------------------------------
 # DNA barcode-level correlation
 # --------------------------------------------------------
 dna_test <- as.data.frame(
 bc[
 dna_columns
 ])
 dna_cor_spearman_barcode <- cor(
 dna_test,
 method = "spearman",
 use = "pairwise.complete.obs")
 dna_cor_pearson_barcode <- cor(
 dna_test,
 method = "pearson",
 use = "pairwise.complete.obs")
 # --------------------------------------------------------
 # RNA barcode-level correlation
 # --------------------------------------------------------
 rna_test <- as.data.frame(
 bc[
 rna_columns
 ])
 rna_cor_spearman_barcode <- cor(
 rna_test,
 method = "spearman",
 use = "pairwise.complete.obs")
 rna_cor_pearson_barcode <- cor(
 rna_test,
 method = "pearson",
 use = "pairwise.complete.obs")
 # ========================================================
 # AGGREGATE DNA/RNA COUNTS TO OLIGO LEVEL
 # ========================================================
 # Aggregate DNA and RNA counts to oligo level
 oligo_counts <- bc %>%
 group_by(oligo_name) %>%
 summarise(
 across(
 all_of(c(dna_columns, rna_columns)),
 sum,
 na.rm = TRUE),
 .groups = "drop")
 # --------------------------------------------------------
 # Oligo-level DNA correlations
 # --------------------------------------------------------
 dna_oligo_test <- as.data.frame(
 oligo_counts[
 dna_columns
 ])
 dna_cor_spearman_oligo <- cor(
 dna_oligo_test,
 method = "spearman",
 use = "pairwise.complete.obs")
 dna_cor_pearson_oligo <- cor(
 dna_oligo_test,
 method = "pearson",
 use = "pairwise.complete.obs")
 # --------------------------------------------------------
 # Oligo-level RNA correlations
 # --------------------------------------------------------
 rna_oligo_test <- as.data.frame(
 oligo_counts[
 rna_columns
 ])
 rna_cor_spearman_oligo <- cor(
 rna_oligo_test,
 method = "spearman",
 use = "pairwise.complete.obs")
 rna_cor_pearson_oligo <- cor(
 rna_oligo_test,
 method = "pearson",
 use = "pairwise.complete.obs")
 # ========================================================
 # BCALM ACTIVITY CORRELATIONS
 #
 # These now use the activity values directly generated by
 # BCalm in Script 1.
 # ========================================================
 activity_oligo_test <- as.data.frame(
 activity_oligo[
 activity_columns
 ])
 colnames(
 activity_oligo_test
) <- activity_columns
 activity_cor_spearman_oligo <- cor(
 activity_oligo_test,
 method = "spearman",
 use = "pairwise.complete.obs")
 activity_cor_pearson_oligo <- cor(
 activity_oligo_test,
 method = "pearson",
 use = "pairwise.complete.obs")
 # --------------------------------------------------------
 # Barcode-level activity correlation
 #
 # This is a descriptive QC measure calculated directly
 # from the barcode-level DNA and RNA counts.
 #
 # It is separate from the BCalm oligo-level activity
 # values used below.
 # --------------------------------------------------------
 activity_barcode <- data.frame(
 barcode = bc$barcode,
 oligo_name = bc$oligo_name)
 for (i in seq_along(common_replicates)) {
 activity_name <- paste0(
 "activity_rep_",
 common_replicates[i])
 activity_barcode[[activity_name]] <- log2(
 (bc[[rna_columns[i]]] + 1) /
 (bc[[dna_columns[i]]] + 1))}
 activity_barcode_columns <- paste0(
 "activity_rep_",
 common_replicates)
 activity_barcode_test <- as.data.frame(
 activity_barcode[
 activity_barcode_columns
 ])
 activity_cor_spearman_barcode <- cor(
 activity_barcode_test,
 method = "spearman",
 use = "pairwise.complete.obs")
 activity_cor_pearson_barcode <- cor(
 activity_barcode_test,
 method = "pearson",
 use = "pairwise.complete.obs")
 # ========================================================
 # SAVE CORRELATION MATRICES
 # ========================================================
 threshold_file_label <- safe_filename(
 threshold_label)
 # --------------------------------------------------------
 # Barcode-level correlation files
 # --------------------------------------------------------
 write.csv(
 dna_cor_spearman_barcode,
 file.path(
 cell_output_dir,
 paste0(
 threshold_file_label,
 "_DNA_correlation_spearman_barcode_level.csv")))
 write.csv(
 dna_cor_pearson_barcode,
 file.path(
 cell_output_dir,
 paste0(
 threshold_file_label,
 "_DNA_correlation_pearson_barcode_level.csv")))
 write.csv(
 rna_cor_spearman_barcode,
 file.path(
 cell_output_dir,
 paste0(
 threshold_file_label,
 "_RNA_correlation_spearman_barcode_level.csv")))
 write.csv(
 rna_cor_pearson_barcode,
 file.path(
 cell_output_dir,
 paste0(
 threshold_file_label,
 "_RNA_correlation_pearson_barcode_level.csv")))
 write.csv(
 activity_cor_spearman_barcode,
 file.path(
 cell_output_dir,
 paste0(
 threshold_file_label,
 "_activity_correlation_spearman_barcode_level.csv")))
 write.csv(
 activity_cor_pearson_barcode,
 file.path(
 cell_output_dir,
 paste0(
 threshold_file_label,
 "_activity_correlation_pearson_barcode_level.csv")))
 # --------------------------------------------------------
 # Oligo-level correlation files
 # --------------------------------------------------------
 write.csv(
 dna_cor_spearman_oligo,
 file.path(
 cell_output_dir,
 paste0(
 threshold_file_label,
 "_DNA_correlation_spearman_oligo_level.csv")))
 write.csv(
 dna_cor_pearson_oligo,
 file.path(
 cell_output_dir,
 paste0(
 threshold_file_label,
 "_DNA_correlation_pearson_oligo_level.csv")))
 write.csv(
 rna_cor_spearman_oligo,
 file.path(
 cell_output_dir,
 paste0(
 threshold_file_label,
 "_RNA_correlation_spearman_oligo_level.csv")))
 write.csv(
 rna_cor_pearson_oligo,
 file.path(
 cell_output_dir,
 paste0(
 threshold_file_label,
 "_RNA_correlation_pearson_oligo_level.csv")))
 write.csv(
 activity_cor_spearman_oligo,
 file.path(
 cell_output_dir,
 paste0(
 threshold_file_label,
 "_activity_correlation_spearman_oligo_level.csv")))
 write.csv(
 activity_cor_pearson_oligo,
 file.path(
 cell_output_dir,
 paste0(
 threshold_file_label,
 "_activity_correlation_pearson_oligo_level.csv")))
 # ========================================================
 # SAVE BCALM OLIGO-LEVEL ACTIVITY DATA
 # ========================================================
 write.csv(
 activity_oligo,
 file.path(
 cell_output_dir,
 paste0(
 threshold_file_label,
 "_BCalm_activity_oligo_level.csv")),
 row.names = FALSE)
 # ========================================================
 # STORE CORRELATION RESULTS
 # ========================================================
 cell_threshold_results[[threshold_name]]$correlations <- list(
 dna = list(
 barcode_level = list(
 pearson = dna_cor_pearson_barcode,
 spearman = dna_cor_spearman_barcode),
 oligo_level = list(
 pearson = dna_cor_pearson_oligo,
 spearman = dna_cor_spearman_oligo)),
 rna = list(
 barcode_level = list(
 pearson = rna_cor_pearson_barcode,
 spearman = rna_cor_spearman_barcode),
 oligo_level = list(
 pearson = rna_cor_pearson_oligo,
 spearman = rna_cor_spearman_oligo)),
 activity = list(
 barcode_level = list(
 pearson = activity_cor_pearson_barcode,
 spearman = activity_cor_spearman_barcode),
 oligo_level = list(
 pearson = activity_cor_pearson_oligo,
 spearman = activity_cor_spearman_oligo)))
 # ========================================================
 # PRINT CORRELATION SUMMARY
 # ========================================================
 cat(
 "\nCorrelation summary:",
 cell,
 "|",
 threshold_label,
 "\n")
 cat(
 "DNA barcode-level Spearman:",
 round(
 mean_upper_triangle(
 dna_cor_spearman_barcode),
 3),
 "\n")
 cat(
 "DNA oligo-level Spearman:",
 round(
 mean_upper_triangle(
 dna_cor_spearman_oligo),
 3),
 "\n")
 cat(
 "DNA oligo-level Pearson:",
 round(
 mean_upper_triangle(
 dna_cor_pearson_oligo),
 3),
 "\n")
 cat(
 "RNA barcode-level Spearman:",
 round(
 mean_upper_triangle(
 rna_cor_spearman_barcode),
 3),
 "\n")
 cat(
 "RNA oligo-level Spearman:",
 round(
 mean_upper_triangle(
 rna_cor_spearman_oligo),
 3),
 "\n")
 cat(
 "RNA oligo-level Pearson:",
 round(
 mean_upper_triangle(
 rna_cor_pearson_oligo),
 3),
 "\n")
 cat(
 "BCalm activity barcode-level Spearman:",
 round(
 mean_upper_triangle(
 activity_cor_spearman_barcode),
 3),
 "\n")
 cat(
 "BCalm activity oligo-level Spearman:",
 round(
 mean_upper_triangle(
 activity_cor_spearman_oligo),
 3),
 "\n")
 cat(
 "BCalm activity oligo-level Pearson:",
 round(
 mean_upper_triangle(
 activity_cor_pearson_oligo),
 3),
 "\n")}
 # ==========================================================
 # THRESHOLD RETENTION SUMMARY
 # ==========================================================
 cell_threshold_summary <- calculate_threshold_summary(
 threshold_counts)
 cell_threshold_summary$cell_line <- cell
 threshold_retention_results[[cell]] <- cell_threshold_summary
 # ----------------------------------------------------------
 # Print retention summary
 # ----------------------------------------------------------
 cat(
 "\n\nThreshold retention:",
 cell,
 "\n")
 print(
 cell_threshold_summary)
 script3_threshold_results[[cell]] <- cell_threshold_results}
# ============================================================
# SECTION 39: THRESHOLD-AWARE CANDIDATE PLOTS
# ============================================================
cat(
 "\n============================================================\n")
cat(
 "SECTION 39: THRESHOLD-AWARE CANDIDATE PLOTS\n")
cat(
 "============================================================\n")
# Formal differential results remain unchanged.
#
# Script 3 only determines whether each WT/MUT candidate pair
# passes the selected barcode threshold in BOTH alleles.
candidate_plot_results <- list()
for (
 candidate_type in c(
 "shared",
 "random_gain",
 "random_loss")
) {
 if (!exists(
 candidate_type
)) {
 next}
 candidates <- get(
 candidate_type)
 if (
 !is.data.frame(candidates) ||
 nrow(candidates) == 0
) {
 next}
 candidate_plot_results[[candidate_type]] <- list()
 for (
 cell in CELL_LINES
) {
 cell_candidates <- candidates %>%
 filter(
 cell_line == cell)
 if (nrow(cell_candidates) == 0) {
 next}
 env <- load_script1_workspace(
 cell)
 bc_original <- get(
 "bc",
 envir = env)
 barcode_support <- bc_original %>%
 group_by(
 oligo_name
) %>%
 summarise(
 n_barcodes = n_distinct(
 barcode),
 .groups = "drop")
 candidate_plot_results[[candidate_type]][[cell]] <- list()
 for (
 threshold_name in names(BARCODE_THRESHOLDS)
) {
 threshold_value <- BARCODE_THRESHOLDS[
 threshold_name
 ]
 threshold_label <- THRESHOLD_LABELS[
 threshold_name
 ]
 # ------------------------------------------------------
 # Check WT and MUT barcode support separately
 # ------------------------------------------------------
 candidate_threshold <- cell_candidates %>%
 left_join(
 barcode_support %>%
 dplyr::rename(
 wt_n_barcodes_threshold =
 n_barcodes),
 by = c(
 "wt_oligo" = "oligo_name")
) %>%
 left_join(
 barcode_support %>%
 dplyr::rename(
 mut_n_barcodes_threshold =
 n_barcodes),
 by = c(
 "mut_oligo" = "oligo_name")
) %>%
 mutate(
 threshold = threshold_label,
 threshold_value = threshold_value,
 wt_passes_threshold =
 ifelse(
 threshold_value == 0,
 TRUE,
 !is.na(
 wt_n_barcodes_threshold
) &
 wt_n_barcodes_threshold >=
 threshold_value),
 mut_passes_threshold =
 ifelse(
 threshold_value == 0,
 TRUE,
 !is.na(
 mut_n_barcodes_threshold
) &
 mut_n_barcodes_threshold >=
 threshold_value),
 passes_threshold =
 wt_passes_threshold &
 mut_passes_threshold)
 # ------------------------------------------------------
 # Only retain candidates where BOTH alleles pass
 # ------------------------------------------------------
 candidate_threshold_pass <- candidate_threshold %>%
 filter(
 passes_threshold)
 # ------------------------------------------------------
 # Output directory
 # ------------------------------------------------------
 threshold_output_dir <- file.path(
 SCRIPT2_DIR,
 "Script3_visualizations",
 cell,
 paste0(
 "threshold_",
 safe_filename(
 threshold_label),
 "____",
 candidate_type))
 dir.create(
 threshold_output_dir,
 recursive = TRUE,
 showWarnings = FALSE)
 candidate_plot_results[[candidate_type]][[cell]][[threshold_name]] <- candidate_threshold_pass
 if (nrow(candidate_threshold_pass) == 0) {
 cat(
 "\nNo candidates passed:",
 cell,
 "|",
 candidate_type,
 "|",
 threshold_label,
 "\n")
 next}
 # ------------------------------------------------------
 # Export threshold-filtered candidate table
 # ------------------------------------------------------
 write.csv(
 candidate_threshold_pass,
 file.path(
 threshold_output_dir,
 paste0(
 safe_filename(cell),
 "_",
 safe_filename(candidate_type),
 "_",
 safe_filename(threshold_label),
 "_candidates.csv")),
 row.names = FALSE)
 # ------------------------------------------------------
 # Candidate activity plots
 #
 # Use the existing Script 1 BCalm activity values.
 #
 # Do NOT recalculate activity here.
 # ------------------------------------------------------
 if (!exists(
 "result",
 envir = env
)) {
 warning(
 "Script 1 result object not found for ",
 cell)
 next}
 bcalm_result <- get(
 "result",
 envir = env)
 for (
 candidate_index in seq_len(
 nrow(candidate_threshold_pass))
) {
 candidate <- candidate_threshold_pass[
 candidate_index,
,
 drop = FALSE
 ]
 wt_oligo <- candidate$wt_oligo
 mut_oligo <- candidate$mut_oligo
 # ----------------------------------------------------
 # Get BCalm activity values from Script 1
 # ----------------------------------------------------
 activity_data <- bcalm_result %>%
 filter(
 oligo_name %in% c(
 wt_oligo,
 mut_oligo))
 if (nrow(activity_data) == 0) {
 next}
 # ----------------------------------------------------
 # Candidate identifier
 # ----------------------------------------------------
 candidate_id <- paste0(
 safe_filename(
 wt_oligo),
 "__",
 safe_filename(
 mut_oligo))
 candidate_plot_results[[candidate_type]] [[cell]] [[threshold_name]] [[candidate_index]] <- list(
 candidate_id = candidate_id,
 activity_data = activity_data)}}}}
# ============================================================
# SECTION 40: RANDOM CANDIDATE EXPORT
# ============================================================
cat(
 "\n============================================================\n")
cat(
 "SECTION 40: RANDOM CANDIDATE EXPORT\n")
cat(
 "============================================================\n")
for (
 candidate_type in c(
 "random_gain",
 "random_loss")
) {
 if (
 !candidate_type %in%
 names(candidate_plot_results)
) {
 next}
 for (
 cell in names(candidate_plot_results[[candidate_type]])
) {
 for (threshold_name in names(candidate_plot_results[[candidate_type]][[cell]])
) {
 candidate_table <- candidate_plot_results[[candidate_type]][[cell]][[threshold_name]]
 if (
 is.data.frame(candidate_table) &&
 nrow(candidate_table) > 0
) {
 output_file <- file.path(
 SCRIPT2_DIR,
 "Script3_visualizations",
 cell,
 paste0(
 candidate_type,
 "_",
 safe_filename(
 THRESHOLD_LABELS[
 threshold_name
 ]),
 "_export.csv"))
 write.csv(
 candidate_table,
 output_file,
 row.names = FALSE)}}}}
# ============================================================
# SECTION 41: THRESHOLD RETENTION SUMMARY
# ============================================================
cat(
 "\n============================================================\n")
cat(
 "SECTION 41: THRESHOLD RETENTION SUMMARY\n")
cat(
 "============================================================\n")
threshold_retention_table <- bind_rows(
 threshold_retention_results)
write.csv(
 threshold_retention_table,
 file.path(
 SCRIPT2_DIR,
 "Script3_threshold_retention_summary.csv"),
 row.names = FALSE)
cat(
 "\nThreshold retention summary:\n")
print(
 threshold_retention_table)
# ------------------------------------------------------------
# Compact cell-by-threshold table
# ------------------------------------------------------------
threshold_retention_long <- bind_rows(
 lapply(
 names(
 script3_threshold_results),
 function(cell) {
 cell_results <- script3_threshold_results[[cell]]
 bind_rows(
 lapply(
 names(cell_results),
 function(threshold_name) {
 data.frame(cell_line = cell,
 threshold =
 THRESHOLD_LABELS[
 threshold_name
 ],
 threshold_value =
 BARCODE_THRESHOLDS[
 threshold_name
 ],
 n_oligos =
 cell_results[[threshold_name]]$n_oligos,
 n_barcode_rows =
 cell_results[[threshold_name]]$n_barcode_rows)}))}))
write.csv(
 threshold_retention_long,
 file.path(
 SCRIPT2_DIR,
 "Script3_threshold_retention_long.csv"),
 row.names = FALSE)
# ============================================================
# SECTION 42: FINAL FORMAL DIFFERENTIAL SUMMARY
# ============================================================
cat(
 "\n============================================================\n")
cat(
 "SECTION 42: FINAL FORMAL DIFFERENTIAL SUMMARY\n")
cat(
 "============================================================\n")
# IMPORTANT:
#
# Formal differential results remain those generated by
# Script 2 using the >=10 barcode threshold.
#
# Script 3 does NOT recalculate mpralm or FDR statistics.
formal_differential_summary <- data.frame()
if (exists(
 "differential_results"
)) {
 formal_differential_summary <- differential_results
} else if (exists(
 "all_differential_results"
)) {
 formal_differential_summary <- all_differential_results
} else if (exists(
 "results"
)) {
 formal_differential_summary <- results}
if (
 is.data.frame(formal_differential_summary) &&
 nrow(formal_differential_summary) > 0
) {
 write.csv(
 formal_differential_summary,
 file.path(
 SCRIPT2_DIR,
 "Script3_formal_differential_summary.csv"),
 row.names = FALSE)
 cat(
 "\nFormal differential results:",
 nrow(
 formal_differential_summary),
 "rows\n")
} else {
 cat(
 "\nNo formal differential results object found.\n")
 cat(
 "Script 2 results remain in:",
 SCRIPT2_DIR,
 "\n")}
# ============================================================
# SECTION 43: MASTER SCRIPT 3 RESULTS OBJECT
# ============================================================
cat(
 "\n============================================================\n")
cat(
 "SECTION 43: MASTER SCRIPT 3 RESULTS OBJECT\n")
cat(
 "============================================================\n")
script3_results <- list(
 threshold_settings = list(
 thresholds = BARCODE_THRESHOLDS,
 labels = THRESHOLD_LABELS),
 threshold_qc = script3_threshold_results,
 threshold_retention =
 threshold_retention_table,
 threshold_retention_long =
 threshold_retention_long,
 candidate_plots =
 candidate_plot_results,
 formal_differential_results =
 formal_differential_summary)
save(
 script3_results,
 file = file.path(
 SCRIPT2_DIR,
 "Script3_results.RData"))
cat(
 "\nSaved master Script 3 results object.\n")
# ============================================================
# SECTION 44: OUTPUT MANIFEST AND COMPLETION SUMMARY
# ============================================================
cat(
 "\n============================================================\n")
cat(
 "SECTION 44: OUTPUT MANIFEST AND COMPLETION SUMMARY\n")
cat(
 "============================================================\n")
# ------------------------------------------------------------
# Output manifest
# ------------------------------------------------------------
output_files <- list.files(
 SCRIPT2_DIR,
 recursive = TRUE,
 full.names = TRUE)
output_manifest <- data.frame(
 file = output_files,
 stringsAsFactors = FALSE)
write.csv(
 output_manifest,
 file.path(
 SCRIPT2_DIR,
 "Script3_output_manifest.csv"),
 row.names = FALSE)
# ------------------------------------------------------------
# Final summary
# ------------------------------------------------------------
cat(
 "\n\n============================================================\n")
cat(
 "FINAL SCRIPT 3 SUMMARY\n")
cat(
 "============================================================\n")
cat(
 "\nThresholds analyzed:\n")
cat(
 " ALL = no barcode threshold\n")
cat(
 " >=10 = at least 10 distinct barcodes per oligo\n")
cat(
 " >=20 = at least 20 distinct barcodes per oligo\n")
cat(
 "\nCorrelation levels analyzed:\n")
cat(
 " 1. Barcode-level Pearson\n")
cat(
 " 2. Barcode-level Spearman\n")
cat(
 " 3. Oligo-level Pearson\n")
cat(
 " 4. Oligo-level Spearman\n")
cat(
 "\nMeasures analyzed:\n")
cat(
 " DNA\n")
cat(
 " RNA\n")
cat(
 " BCalm activity\n")
cat(
 "\nBCalm activity calculation:\n")
cat(
 " Activity values are taken directly from the Script 1\n")
cat(
 " BCalm replicate_activity object.\n")
cat(
 " Script 3 does not reconstruct BCalm activity from\n")
cat(
 " summed DNA/RNA counts.\n")
cat(
 "\nFormal differential statistics:\n")
cat(
 " Script 2 >=10 analysis retained unchanged.\n")
cat(
 " Script 3 does not recalculate mpralm or FDR.\n")
cat(
 "\nCandidate plots:\n")
cat(
 " Candidates are only retained at a threshold when BOTH WT\n")
cat(
 " and MUT pass that barcode threshold.\n")
cat(
 "\nOutput directory:\n")
cat(
 " ",
 SCRIPT2_DIR,
 "\n",
 sep = "")
cat(
 "\nScript 3 completed successfully.\n")
cat(
 "============================================================\n")
# ============================================================
# SECTION 45: THRESHOLD RETENTION AND CORRELATION SUMMARY
# ============================================================
cat(
 "\n============================================================\n")
cat(
 "SECTION 45: THRESHOLD RETENTION AND CORRELATION SUMMARY\n")
cat(
 "============================================================\n")
# ------------------------------------------------------------
# Helper function
# ------------------------------------------------------------
get_mean_upper_triangle <- function(
 correlation_matrix
) {
 if (
 is.null(correlation_matrix) ||
 length(correlation_matrix) == 0
) {
 return(NA_real_)}
 if (
 !is.matrix(correlation_matrix) ||
 nrow(correlation_matrix) < 2
) {
 return(NA_real_)}
 values <- correlation_matrix[
 upper.tri(correlation_matrix)
 ]
 if (length(values) == 0) {
 return(NA_real_)}
 if (all(is.na(values))) {
 return(NA_real_)}
 mean(
 values,
 na.rm = TRUE)}
# ------------------------------------------------------------
# Build summary table
# ------------------------------------------------------------
threshold_correlation_summary <- list()
for (cell in CELL_LINES) {
 if (
 !cell %in%
 names(script3_threshold_results)
) {
 next}
 cell_results <- script3_threshold_results[[cell]]
 # ----------------------------------------------------------
 # Get ALL oligo count for this cell
 # ----------------------------------------------------------
 if (!"ALL" %in% names(cell_results)) {
 next}
 all_n <- cell_results[["ALL"]]$n_oligos
 # ----------------------------------------------------------
 # Loop through thresholds
 # ----------------------------------------------------------
 for (
 threshold_name in names(BARCODE_THRESHOLDS)
) {
 if (
 !threshold_name %in%
 names(cell_results)
) {
 next}
 threshold_result <-
 cell_results[[threshold_name]]
 # --------------------------------------------------------
 # Oligo count
 # --------------------------------------------------------
 n_oligos <-
 threshold_result$n_oligos
 # --------------------------------------------------------
 # Retention from ALL
 # --------------------------------------------------------
 retained_percent_of_all <-
 100 * n_oligos / all_n
 # --------------------------------------------------------
 # Loss from previous threshold
 # --------------------------------------------------------
 if (
 threshold_name == "ALL"
) {
 lost_percent_from_previous <- 0
 } else if (
 threshold_name == "min10"
) {
 previous_n <-
 cell_results[["ALL"]]$n_oligos
 lost_percent_from_previous <-
 100 * (
 previous_n - n_oligos
) / previous_n
 } else if (
 threshold_name == "min20"
) {
 previous_n <-
 cell_results[["min10"]]$n_oligos
 lost_percent_from_previous <-
 100 * (
 previous_n - n_oligos
) / previous_n
 } else {
 lost_percent_from_previous <- NA_real_}
 # --------------------------------------------------------
 # Correlations
 # --------------------------------------------------------
 dna_spearman <-
 get_mean_upper_triangle(
 threshold_result$correlations$dna$oligo_level$spearman)
 dna_pearson <-
 get_mean_upper_triangle(
 threshold_result$correlations$dna$oligo_level$pearson)
 rna_spearman <-
 get_mean_upper_triangle(
 threshold_result$correlations$rna$oligo_level$spearman)
 rna_pearson <-
 get_mean_upper_triangle(
 threshold_result$correlations$rna$oligo_level$pearson)
 activity_spearman <-
 get_mean_upper_triangle(
 threshold_result$correlations$activity$oligo_level$spearman)
 activity_pearson <-
 get_mean_upper_triangle(
 threshold_result$correlations$activity$oligo_level$pearson)
 # --------------------------------------------------------
 # Store row
 # --------------------------------------------------------
 threshold_correlation_summary[[
 paste0(
 cell,
 "_",
 threshold_name)
 ]] <- data.frame(
 cell_line = cell,
 threshold =
 THRESHOLD_LABELS[
 threshold_name
 ],
 threshold_value =
 BARCODE_THRESHOLDS[
 threshold_name
 ],
 n_oligos =
 n_oligos,
 retained_percent_of_all =
 retained_percent_of_all,
 lost_percent_from_previous =
 lost_percent_from_previous,
 dna_spearman =
 dna_spearman,
 dna_pearson =
 dna_pearson,
 rna_spearman =
 rna_spearman,
 rna_pearson =
 rna_pearson,
 bcalm_activity_spearman =
 activity_spearman,
 bcalm_activity_pearson =
 activity_pearson,
 stringsAsFactors = FALSE)}}
# ------------------------------------------------------------
# Combine rows
# ------------------------------------------------------------
threshold_retention_correlation_summary <-
 bind_rows(
 threshold_correlation_summary)
# ------------------------------------------------------------
# Save summary
# ------------------------------------------------------------
summary_file <- file.path(
 SCRIPT2_DIR,
 "Script3_threshold_retention_correlation_summary.csv")
write.csv(
 threshold_retention_correlation_summary,
 summary_file,
 row.names = FALSE)
# ------------------------------------------------------------
# Print summary
# ------------------------------------------------------------
cat(
 "\nSaved threshold retention/correlation summary:\n")
cat(
 " ",
 summary_file,
 "\n",
 sep = "")
print(
 threshold_retention_correlation_summary)
############################################################
# 9. DEFINE EFFECT CLASSES
############################################################
cat("\n============================================================\n")
cat("Defining biological effect classes\n")
cat("============================================================\n")
differential_all <- differential_all %>%
 mutate(
 effect_abs = abs(logFC),
 strong_significant = is.finite(adj.P.Val) & adj.P.Val <= FDR_THRESHOLD &
 is.finite(logFC) & effect_abs >= STRONG_EFFECT_THRESHOLD,
 weak_significant = is.finite(adj.P.Val) & adj.P.Val <= FDR_THRESHOLD &
 is.finite(logFC) & effect_abs >= WEAK_EFFECT_THRESHOLD,
 strong_gain = strong_significant & logFC >= STRONG_EFFECT_THRESHOLD,
 strong_loss = strong_significant & logFC <= -STRONG_EFFECT_THRESHOLD,
 weak_gain = weak_significant & logFC >= WEAK_EFFECT_THRESHOLD,
 weak_loss = weak_significant & logFC <= -WEAK_EFFECT_THRESHOLD,
 effect_direction = case_when(
 logFC >= STRONG_EFFECT_THRESHOLD ~ "strong_gain",
 logFC <= -STRONG_EFFECT_THRESHOLD ~ "strong_loss",
 logFC >= WEAK_EFFECT_THRESHOLD ~ "weak_gain",
 logFC <= -WEAK_EFFECT_THRESHOLD ~ "weak_loss",
 TRUE ~ "minimal"),
 significance_class = case_when(
 strong_significant ~ "strong_significant",
 weak_significant ~ "weak_significant",
 TRUE ~ "not_significant"),
 comparison_key = paste(wt_oligo, mut_oligo, sep = "__"))
cat("Effect classes defined.\n")
############################################################
# 9A. STRONG AND WEAK FDR HEATMAPS
############################################################
# 9A.1 Output directories
FDR_HEATMAP_DIR <- file.path(PLOT_DIR, "WT_MUT_FDR_heatmaps")
STRONG_HEATMAP_DIR <- file.path(FDR_HEATMAP_DIR, "strong_FDR")
WEAK_HEATMAP_DIR <- file.path(FDR_HEATMAP_DIR, "weak_FDR")
dir.create(STRONG_HEATMAP_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(WEAK_HEATMAP_DIR, recursive = TRUE, showWarnings = FALSE)
# 9A.2 Separate strong and weak effects
strong_fdr_data <- differential_all %>%
 filter(is.finite(adj.P.Val), adj.P.Val <= FDR_THRESHOLD,
 is.finite(logFC), abs(logFC) >= STRONG_EFFECT_THRESHOLD)
weak_fdr_data <- differential_all %>%
 filter(is.finite(adj.P.Val), adj.P.Val <= FDR_THRESHOLD,
 is.finite(logFC), abs(logFC) >= WEAK_EFFECT_THRESHOLD,
 abs(logFC) < STRONG_EFFECT_THRESHOLD)
cat("\nStrong significant comparisons:", nrow(strong_fdr_data), "\n")
cat("Weak significant comparisons:", nrow(weak_fdr_data), "\n")
# 9A.3 Shared heatmap plotting function
plot_fdr_heatmap <- function(data, title, x_levels) {
 max_logFC <- max(1, abs(data$logFC), na.rm = TRUE)
 ggplot(data, aes(x = cell_line, y = comparison_label, fill = logFC)) +
 geom_tile(color = "white", linewidth = 0.6, width = 0.96, height = 0.96) +
 geom_text(aes(label = ifelse(is.finite(logFC), sprintf("%.2f", logFC), "")),
 size = 4.5, fontface = "bold") +
 scale_fill_gradient2(low = "#B2182B", mid = "white", high = "#2166AC",
 midpoint = 0, limits = c(-max_logFC, max_logFC),
 name = "WT → MUT\nlogFC", na.value = "grey92") +
 scale_x_discrete(limits = x_levels, drop = FALSE) +
 labs(title = title, subtitle = paste("FDR ≤", FDR_THRESHOLD),
 x = NULL, y = "WT → MUT comparison") +
 theme_minimal(base_size = 12) +
 theme(panel.grid = element_blank(),
 axis.text.x = element_text(size = 12, face = "bold"),
 axis.text.y = element_text(size = 9.5, lineheight = 0.95),
 axis.title.y = element_text(size = 11, face = "bold"),
 plot.title = element_text(size = 15, face = "bold"),
 plot.subtitle = element_text(size = 11),
 legend.position = "bottom")}
# 9A.4 Cross-cell-line heatmaps
# A comparison must be significant in at least 2 cell lines.
create_cross_cell_heatmap <- function(class_data, effect_label, output_dir,
 max_comparisons = 30) {
 sig_data <- class_data %>%
 filter(is.finite(logFC), is.finite(adj.P.Val),
 adj.P.Val <= FDR_THRESHOLD) %>%
 distinct(comparison_key, cell_line, .keep_all = TRUE)
 selected <- sig_data %>%
 group_by(comparison_key, pair_id, wt_oligo, mut_oligo) %>%
 summarise(n_cells = n_distinct(cell_line),
 best_fdr = min(adj.P.Val, na.rm = TRUE),
 max_effect = max(abs(logFC), na.rm = TRUE),
 .groups = "drop") %>%
 filter(n_cells >= 2) %>%
 arrange(desc(n_cells), best_fdr, desc(max_effect)) %>%
 slice_head(n = max_comparisons) %>%
 mutate(comparison_label = paste0(wt_oligo, "\n→ ", mut_oligo))
 if (nrow(selected) == 0) {
 message("No cross-cell-line ", effect_label, " comparisons found.")
 return(invisible(NULL))}
 plot_data <- tidyr::crossing(
 comparison_key = selected$comparison_key,
 cell_line = CELL_LINES
) %>%
 left_join(selected %>% select(comparison_key, comparison_label),
 by = "comparison_key") %>%
 left_join(sig_data %>% select(comparison_key, cell_line, logFC),
 by = c("comparison_key", "cell_line"))
 plot_data$comparison_label <- factor(
 plot_data$comparison_label,
 levels = rev(selected$comparison_label))
 p <- plot_fdr_heatmap(
 plot_data,
 paste("Cross-cell-line", effect_label, "FDR comparisons"),
 CELL_LINES)
 # Height is based on unique comparisons, not tile count.
 plot_height <- min(45, max(9, 4 + nrow(selected) * 0.48))
 save_plot(
 p,
 file.path(output_dir, paste0("cross_cell_", effect_label, "_FDR_heatmap.png")),
 width = 16,
 height = plot_height)
 cat("Created cross-cell-line", effect_label, "heatmap:",
 nrow(selected), "comparisons\n")
 invisible(plot_data)}
# 9A.5 Cell-line-specific heatmaps
# A comparison must be significant in exactly 1 cell line.
create_specific_cell_heatmaps <- function(class_data, effect_label, output_dir,
 max_comparisons = 30) {
 specific_data <- class_data %>%
 filter(is.finite(logFC), is.finite(adj.P.Val),
 adj.P.Val <= FDR_THRESHOLD) %>%
 distinct(comparison_key, cell_line, .keep_all = TRUE) %>%
 group_by(comparison_key) %>%
 filter(n_distinct(cell_line) == 1) %>%
 ungroup()
 for (cell in CELL_LINES) {
 cell_data <- specific_data %>%
 filter(cell_line == cell) %>%
 arrange(adj.P.Val, desc(abs(logFC))) %>%
 slice_head(n = max_comparisons) %>%
 mutate(comparison_label = paste0(wt_oligo, "\n→ ", mut_oligo))
 if (nrow(cell_data) == 0) {
 message("No cell-line-specific ", effect_label, " comparisons for ", cell)
 next}
 cell_data$comparison_label <- factor(
 cell_data$comparison_label,
 levels = rev(cell_data$comparison_label))
 p <- plot_fdr_heatmap(
 cell_data,
 paste(cell, "-", effect_label, "cell-line-specific comparisons"),
 cell)
 plot_height <- min(45, max(9, 4 + nrow(cell_data) * 0.48))
 save_plot(
 p,
 file.path(output_dir, paste0(cell, "_", effect_label, "_FDR_heatmap.png")),
 width = 9,
 height = plot_height)
 cat("Created", cell, effect_label, "specific heatmap:",
 nrow(cell_data), "comparisons\n")}
 invisible(specific_data)}
# 9A.6 Generate strong heatmaps
strong_cross_data <- create_cross_cell_heatmap(
 class_data = strong_fdr_data,
 effect_label = "strong",
 output_dir = STRONG_HEATMAP_DIR,
 max_comparisons = 30)
strong_specific_data <- create_specific_cell_heatmaps(
 class_data = strong_fdr_data,
 effect_label = "strong",
 output_dir = STRONG_HEATMAP_DIR,
 max_comparisons = 30)
# 9A.7 Generate weak heatmaps
weak_cross_data <- create_cross_cell_heatmap(
 class_data = weak_fdr_data,
 effect_label = "weak",
 output_dir = WEAK_HEATMAP_DIR,
 max_comparisons = 30)
weak_specific_data <- create_specific_cell_heatmaps(
 class_data = weak_fdr_data,
 effect_label = "weak",
 output_dir = WEAK_HEATMAP_DIR,
 max_comparisons = 30)
# 9A.8 Confirm saved files
heatmap_files <- list.files(
 FDR_HEATMAP_DIR,
 pattern = "\\.(png|pdf)$",
 recursive = TRUE,
 full.names = TRUE)
cat("\n============================================================\n")
cat("SECTION 9A: FDR HEATMAP GENERATION FINISHED\n")
cat("Saved heatmap files:", length(heatmap_files), "\n")
cat("Output directory:", FDR_HEATMAP_DIR, "\n")
cat("============================================================\n")
if (length(heatmap_files) == 0) {
 warning("No heatmap files were found. Check the plotting functions and output paths.")}
save_plot <- function(plot_object, filename, width = 8, height = 6) {
 ggsave(filename, plot_object, width = width, height = height, dpi = 300)
 pdf_filename <- sub("\\.png$", ".pdf", filename)
 ggsave(pdf_filename, plot_object, width = width, height = height, device = "pdf")}
############################################################
# 9. DEFINE EFFECT CLASSES
############################################################
cat("\n============================================================\n")
cat("Defining biological effect classes\n")
cat("============================================================\n")
differential_all <- differential_all %>%
 mutate(
 effect_abs = abs(logFC),
 strong_significant =
 !is.na(adj.P.Val) &
 adj.P.Val <= FDR_THRESHOLD &
 effect_abs >= STRONG_EFFECT_THRESHOLD,
 # Retain original inclusive weak definition for
 # downstream summaries and classifications.
 weak_significant =
 !is.na(adj.P.Val) &
 adj.P.Val <= FDR_THRESHOLD &
 effect_abs >= WEAK_EFFECT_THRESHOLD,
 strong_gain =
 strong_significant &
 logFC >= STRONG_EFFECT_THRESHOLD,
 strong_loss =
 strong_significant &
 logFC <= -STRONG_EFFECT_THRESHOLD,
 weak_gain =
 weak_significant &
 logFC >= WEAK_EFFECT_THRESHOLD,
 weak_loss =
 weak_significant &
 logFC <= -WEAK_EFFECT_THRESHOLD,
 effect_direction = case_when(
 logFC >= STRONG_EFFECT_THRESHOLD ~ "strong_gain",
 logFC <= -STRONG_EFFECT_THRESHOLD ~ "strong_loss",
 logFC >= WEAK_EFFECT_THRESHOLD ~ "weak_gain",
 logFC <= -WEAK_EFFECT_THRESHOLD ~ "weak_loss",
 TRUE ~ "minimal"),
 significance_class = case_when(
 strong_significant ~ "strong_significant",
 weak_significant ~ "weak_significant",
 TRUE ~ "not_significant"),
 comparison_key = paste(
 wt_oligo,
 mut_oligo,
 sep = "__"))
############################################################
# 9A. WT/MUT FDR HEATMAPS
############################################################
cat("\n============================================================\n")
cat("Creating WT/MUT FDR heatmaps\n")
cat("============================================================\n")
# ----------------------------------------------------------
# 1. OUTPUT DIRECTORIES AND SETTINGS
# ----------------------------------------------------------
FDR_HEATMAP_DIR <- file.path(
 PLOT_DIR,
 "WT_MUT_FDR_heatmaps")
STRONG_HEATMAP_DIR <- file.path(
 FDR_HEATMAP_DIR,
 "strong")
WEAK_HEATMAP_DIR <- file.path(
 FDR_HEATMAP_DIR,
 "weak")
dir.create(
 STRONG_HEATMAP_DIR,
 recursive = TRUE,
 showWarnings = FALSE)
dir.create(
 WEAK_HEATMAP_DIR,
 recursive = TRUE,
 showWarnings = FALSE)
MAX_HEATMAP_COMPARISONS <- 30
# Consistent scale across all strong and weak heatmaps.
# Zero is always white; negative is red; positive is blue.
all_fdr_logFC <- differential_all %>%
 filter(
 !is.na(adj.P.Val),
 adj.P.Val <= FDR_THRESHOLD,
 is.finite(logFC)
) %>%
 pull(logFC)
if (length(all_fdr_logFC) == 0) {
 stop("No FDR-significant finite logFC values are available for heatmaps.")}
HEATMAP_LIMIT <- max(
 1,
 abs(all_fdr_logFC),
 na.rm = TRUE)
cat(
 "\nHeatmap colour limits: ",
 -HEATMAP_LIMIT,
 " to ",
 HEATMAP_LIMIT,
 "\n",
 sep = "")
# ----------------------------------------------------------
# 2. DEFINE STRONG AND WEAK HEATMAP DATA
# ----------------------------------------------------------
strong_fdr_data <- differential_all %>%
 filter(
 !is.na(adj.P.Val),
 adj.P.Val <= FDR_THRESHOLD,
 is.finite(logFC),
 abs(logFC) >= STRONG_EFFECT_THRESHOLD)
# Weak is deliberately exclusive of strong:
# 0.5 <= abs(logFC) < 1.0
weak_fdr_data <- differential_all %>%
 filter(
 !is.na(adj.P.Val),
 adj.P.Val <= FDR_THRESHOLD,
 is.finite(logFC),
 abs(logFC) >= WEAK_EFFECT_THRESHOLD,
 abs(logFC) < STRONG_EFFECT_THRESHOLD)
cat("\nStrong FDR-significant comparisons:\n")
print(
 strong_fdr_data %>%
 count(cell_line, name = "n"))
cat("\nWeak-only FDR-significant comparisons:\n")
print(
 weak_fdr_data %>%
 count(cell_line, name = "n"))
# ----------------------------------------------------------
# 3. COMMON HEATMAP THEME AND PLOT FUNCTION
# ----------------------------------------------------------
plot_fdr_heatmap <- function(
 plot_data,
 title,
 subtitle,
 x_levels,
 specific = FALSE
) {
 plot_data <- plot_data %>%
 mutate(
 cell_line = factor(
 cell_line,
 levels = x_levels),
 comparison_label = factor(
 comparison_label,
 levels = rev(unique(comparison_label))))
 p <- ggplot(
 plot_data,
 aes(
 x = cell_line,
 y = comparison_label,
 fill = logFC)
) +
 geom_tile(
 color = "white",
 linewidth = 0.7,
 width = 0.95,
 height = 0.95
) +
 geom_text(
 aes(
 label = ifelse(
 is.finite(logFC),
 sprintf("%.2f", logFC),
 "")),
 size = 4.2,
 fontface = "bold",
 color = "black"
) +
 scale_fill_gradient2(
 low = "#B2182B",
 mid = "white",
 high = "#2166AC",
 midpoint = 0,
 limits = c(-HEATMAP_LIMIT, HEATMAP_LIMIT),
 oob = scales::squish,
 na.value = "grey88",
 name = "WT \u2192 MUT\nlog2 fold change"
) +
 scale_x_discrete(
 limits = x_levels,
 drop = FALSE
) +
 labs(
 title = title,
 subtitle = subtitle,
 x = "Cell line",
 y = "WT/MUT comparison"
) +
 theme_bw(base_size = 14) +
 theme(
 plot.title = element_text(
 size = 16,
 face = "bold",
 hjust = 0.5,
 margin = margin(b = 8)),
 plot.subtitle = element_text(
 size = 11,
 hjust = 0.5,
 margin = margin(b = 12)),
 axis.title.x = element_text(
 size = 13,
 face = "bold",
 margin = margin(t = 10)),
 axis.title.y = element_text(
 size = 13,
 face = "bold",
 margin = margin(r = 10)),
 axis.text.x = element_text(
 size = 12,
 face = "bold",
 color = "black",
 angle = ifelse(specific, 0, 0),
 hjust = 0.5),
 axis.text.y = element_text(
 size = 10.5,
 color = "black"),
 legend.title = element_text(
 size = 11,
 face = "bold"),
 legend.text = element_text(
 size = 10),
 legend.key.height = grid::unit(
 0.8,
 "cm"),
 panel.grid = element_blank(),
 panel.border = element_rect(
 color = "grey35",
 linewidth = 0.6),
 plot.margin = margin(
 t = 15,
 r = 20,
 b = 15,
 l = 15))
 p}
# ----------------------------------------------------------
# 4. SAVE HEATMAP AS PNG AND PDF
# ----------------------------------------------------------
save_heatmap <- function(
 plot_object,
 filename,
 width,
 height
) {
 # Keep dimensions safely below ggplot's 50-inch limit.
 width <- min(width, 30)
 height <- min(height, 45)
 png_file <- paste0(
 filename,
 ".png")
 pdf_file <- paste0(
 filename,
 ".pdf")
 ggsave(
 filename = png_file,
 plot = plot_object,
 width = width,
 height = height,
 units = "in",
 dpi = 300,
 bg = "white",
 limitsize = TRUE)
 ggsave(
 filename = pdf_file,
 plot = plot_object,
 width = width,
 height = height,
 units = "in",
 device = "pdf",
 bg = "white",
 limitsize = TRUE)
 cat(
 "\nSaved:\n",
 png_file,
 "\n",
 pdf_file,
 "\n")}
# ----------------------------------------------------------
# 5. CROSS-CELL-LINE HEATMAP FUNCTION
# ----------------------------------------------------------
create_cross_cell_heatmap <- function(
 class_data,
 effect_label,
 output_dir,
 max_comparisons = MAX_HEATMAP_COMPARISONS
) {
 cat(
 "\nCreating shared ",
 effect_label,
 " heatmap...\n",
 sep = "")
 if (nrow(class_data) == 0) {
 cat("No comparisons available.\n")
 return(invisible(NULL))}
 # Count the number of cell lines with significant
 # effects in this particular effect class.
 shared_summary <- class_data %>%
 group_by(comparison_key) %>%
 summarise(
 n_cell_lines = n_distinct(cell_line),
 best_fdr = min(adj.P.Val, na.rm = TRUE),
 max_effect = max(abs(logFC), na.rm = TRUE),
 pair_id = first(as.character(pair_id)),
 .groups = "drop"
) %>%
 filter(
 n_cell_lines >= 2
) %>%
 arrange(
 desc(n_cell_lines),
 best_fdr,
 desc(max_effect)
) %>%
 slice_head(
 n = max_comparisons)
 if (nrow(shared_summary) == 0) {
 cat("No shared comparisons found.\n")
 return(invisible(NULL))}
 # Use pair_id rather than long WT and MUT oligo names.
 shared_summary <- shared_summary %>%
 mutate(
 comparison_label = as.character(pair_id))
 # Protect against duplicate pair_id labels.
 shared_summary$comparison_label <- make.unique(
 shared_summary$comparison_label)
 # Create every comparison x cell-line combination.
 # A cell without a significant effect in this class
 # remains NA and is displayed grey.
 plot_data <- tidyr::expand_grid(
 comparison_key = shared_summary$comparison_key,
 cell_line = CELL_LINES
) %>%
 left_join(
 shared_summary %>%
 select(
 comparison_key,
 comparison_label),
 by = "comparison_key"
) %>%
 left_join(
 class_data %>%
 select(
 comparison_key,
 cell_line,
 logFC),
 by = c(
 "comparison_key",
 "cell_line"))
 plot_data$comparison_label <- factor(
 plot_data$comparison_label,
 levels = rev(shared_summary$comparison_label))
 p <- plot_fdr_heatmap(
 plot_data = plot_data,
 title = paste(
 "Shared",
 effect_label,
 "WT/MUT effects across cell lines"),
 subtitle = paste0(
 "FDR \u2264 ",
 FDR_THRESHOLD,
 " | Significant in at least 2 cell lines | Top ",
 nrow(shared_summary),
 " comparisons"),
 x_levels = CELL_LINES)
 plot_height <- min(
 45,
 max(
 8,
 3.5 + nrow(shared_summary) * 0.38))
 save_heatmap(
 p,
 file.path(
 output_dir,
 paste0(
 "shared_",
 effect_label,
 "_FDR_heatmap")),
 width = 15,
 height = plot_height)
 invisible(shared_summary)}
# ----------------------------------------------------------
# 6. CELL-LINE-SPECIFIC HEATMAP FUNCTION
# ----------------------------------------------------------
create_specific_cell_heatmaps <- function(
 class_data,
 effect_label,
 output_dir,
 max_comparisons = MAX_HEATMAP_COMPARISONS
) {
 cat(
 "\nCreating cell-line-specific ",
 effect_label,
 " heatmaps...\n",
 sep = "")
 if (nrow(class_data) == 0) {
 cat("No comparisons available.\n")
 return(invisible(NULL))}
 # Specific means significant in exactly one cell line
 # within the current effect class.
 specific_keys <- class_data %>%
 group_by(comparison_key) %>%
 summarise(
 n_cell_lines = n_distinct(cell_line),
 .groups = "drop"
) %>%
 filter(
 n_cell_lines == 1
) %>%
 pull(
 comparison_key)
 specific_data <- class_data %>%
 filter(
 comparison_key %in% specific_keys)
 if (nrow(specific_data) == 0) {
 cat("No cell-line-specific comparisons found.\n")
 return(invisible(NULL))}
 saved_counts <- setNames(
 integer(length(CELL_LINES)),
 CELL_LINES)
 for (cell in CELL_LINES) {
 cell_data <- specific_data %>%
 filter(
 cell_line == cell
) %>%
 arrange(
 adj.P.Val,
 desc(abs(logFC))
) %>%
 slice_head(
 n = max_comparisons)
 if (nrow(cell_data) == 0) {
 cat(
 "\n",
 cell,
 ": no specific ",
 effect_label,
 " comparisons.\n",
 sep = "")
 next}
 # Compact row labels.
 cell_data <- cell_data %>%
 mutate(
 comparison_label = as.character(pair_id))
 cell_data$comparison_label <- make.unique(
 cell_data$comparison_label)
 p <- plot_fdr_heatmap(
 plot_data = cell_data,
 title = paste(
 cell,
 "cell-line-specific",
 effect_label,
 "significant effects"),
 subtitle = paste0(
 "FDR \u2264 ",
 FDR_THRESHOLD,
 " | Significant in exactly one cell line | ",
 nrow(cell_data),
 " comparisons"),
 x_levels = cell,
 specific = TRUE)
 plot_height <- min(
 45,
 max(
 8,
 3.5 + nrow(cell_data) * 0.38))
 save_heatmap(
 p,
 file.path(
 output_dir,
 paste0(
 cell,
 "_specific_",
 effect_label,
 "_FDR_heatmap")),
 width = 9,
 height = plot_height)
 saved_counts[cell] <- nrow(cell_data)}
 invisible(saved_counts)}
# ----------------------------------------------------------
# 7. GENERATE STRONG HEATMAPS
# ----------------------------------------------------------
strong_shared_summary <- create_cross_cell_heatmap(
 class_data = strong_fdr_data,
 effect_label = "strong",
 output_dir = STRONG_HEATMAP_DIR)
strong_specific_counts <- create_specific_cell_heatmaps(
 class_data = strong_fdr_data,
 effect_label = "strong",
 output_dir = STRONG_HEATMAP_DIR)
# ----------------------------------------------------------
# 8. GENERATE WEAK HEATMAPS
# ----------------------------------------------------------
weak_shared_summary <- create_cross_cell_heatmap(
 class_data = weak_fdr_data,
 effect_label = "weak",
 output_dir = WEAK_HEATMAP_DIR)
weak_specific_counts <- create_specific_cell_heatmaps(
 class_data = weak_fdr_data,
 effect_label = "weak",
 output_dir = WEAK_HEATMAP_DIR)
# ----------------------------------------------------------
# 9. VERIFY OUTPUT FILES
# ----------------------------------------------------------
heatmap_files <- list.files(
 FDR_HEATMAP_DIR,
 pattern = "\\.(png|pdf)$",
 recursive = TRUE,
 full.names = TRUE)
cat("\n============================================================\n")
cat("FDR HEATMAP OUTPUT CHECK\n")
cat("============================================================\n")
cat(
 "\nTotal heatmap files saved:",
 length(heatmap_files),
 "\n")
if (length(heatmap_files) > 0) {
 print(
 data.frame(
 file = heatmap_files,
 size_MB = round(
 file.info(heatmap_files)$size / 1024^2,
 2)))
} else {
 warning("No heatmap PNG or PDF files were created.")}
cat("\nFDR heatmap section finished.\n")
