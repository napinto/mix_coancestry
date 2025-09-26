# Clear existing variables in the environment for a clean session
rm(list = ls())

# Load required packages
library(dplyr)
library(readr)

# Source auxiliary functions related to likelihood ratio calculations per marker
source("3cont_auxiliary_funcs.R")

# Read allele frequency data from CSV file
# The CSV should be formatted as a matrix:
# - First column contains allele identifiers (e.g., allele names or codes)
# - Subsequent columns represent different genetic markers
# - Cell [i,j] contains the frequency of allele i-1 for marker j-1
FREQ_ALELICAS <- read_csv("FREQ_ALELICAS.csv",show_col_types = FALSE)

# Number of profile simulations to run per marker
nsims <- 1000000

# Extract marker names from frequency data (excluding allele column)
markers <- FREQ_ALELICAS[-1] %>% colnames()

# Define LR interval cutoff points for grouping results
# it can be defined explicitly, e.g. cuts = c(10^-4,10^-3,10^-2,10^-1,10^0,10^1,10^2,10^3,10^4)
# or as a sequence e.g. cuts = 10^(-10:10)
cuts = seq(0,1000,0.1)
cuts_2 = 10^(-10:10) # c(-Inf, 1/2, 1, 2, Inf)

results_all_markers <-purrr::map(
  .x = markers,
  .f = ~{ 
    LR_per_marker(
      marker = .x,
      nsims = nsims,
      FREQ_ALELICAS = FREQ_ALELICAS,
      seed = NULL,
      cuts = cuts
    )
  }
) %>% 
  setNames(nm = markers)

results_all_markers %>% print()

summary_stats_all_markers <- extract_summary_and_counts(results_all_markers)
summary_stats_all_markers %>% print()
openxlsx::write.xlsx(x = summary_stats_all_markers,file = "3cont_summary_stats.xlsx")

full_profiles = concatenate_columns(results_all_markers)

results_per_mixture <- tibble(
  "LR[0]" = multiply_columns(results_all_markers,"LR[0]"),
  "LR[0.01]" = multiply_columns(results_all_markers,"LR[0.01]"),
  "LR[0]/LR[0.01]" = multiply_columns(results_all_markers,"LR[0]") / multiply_columns(results_all_markers,"LR[0.01]")
) %>%
  mutate(group = cut(`LR[0]/LR[0.01]`, cuts_2))
results_per_mixture %>% print()
#results_per_mixture[1:10,] %>% openxlsx::write.xlsx(x = .,file = "3cont_results_per_mixture_10.xlsx")
summary_stats_per_mixture <- extract_summary_and_counts(list(results_per_mixture = results_per_mixture))

bind_rows(summary_stats_per_mixture,summary_stats_all_markers)
openxlsx::write.xlsx(x = summary_stats_per_mixture,file = "3cont_summary_stats_per_mixture.xlsx")

# Results computed overall (combining all markers) and by marker
by_case_type_by_marker <- purrr::imap(
  c(list(overall = results_all_markers %>% bind_rows()), results_all_markers), function(df, marker) {
  df <- df %>%
    mutate(
      case_type = factor(paste0("( ",n_1," , " , n_mixture, " )")),
      group = factor(group)
    ) %>%
    group_by(case_type) 
  
  full_join(
    df %>%
      count(group) %>%
      mutate(prc = round(100 * n / sum(n), 1)) %>%
      ungroup()  %>% 
      tidyr::pivot_wider(names_from = "group",values_from = c("n","prc"), names_glue = "{group} {.value}"),
    df %>% 
      summarise(
        Min = min(`LR[0]/LR[0.01]`), 
        Max = max(`LR[0]/LR[0.01]`),
        Mean = mean(`LR[0]/LR[0.01]`),
        Median = median(`LR[0]/LR[0.01]`)
      ),
    by = join_by(case_type)
  ) %>%
    mutate(marker = marker,.before = 1) 
}) %>% 
  bind_rows() %>%
  select(marker, case_type, sort(names(.)))

by_case_type_by_marker %>% openxlsx::write.xlsx(x = .,file = "3cont_by_case_type_by_marker.xlsx")

