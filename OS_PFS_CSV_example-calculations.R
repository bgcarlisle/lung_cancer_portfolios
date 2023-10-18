library(tidyverse)
library(bayesmeta)
library(ggplot2)

pivotal <- read_csv("nccn_pilot_clean_formatted.csv")

nonpivotal <- tribble(
    ~nctid, ~endpoint, ~hr, ~ci_lower, ~ci_upper, ~p_value, ~studylabel,
    "NCT00000007", "os", 0.9, 0.80, 1.00, 0.04, "Jones 2001",
    "NCT00000008", "os", 0.9, 0.85, 0.95, 0.04, "Jones 2002",
    "NCT00000009", "os", 0.5, 0.40, 0.60, 0.05, "Jones 2003",
    "NCT00000010", "os", 0.6, 0.45, 0.85, 0.04, "Jones 2004",
    "NCT00000011", "os", 0.7, 0.55, 0.95, 0.04, "Jones 2005",
    "NCT00000012", "os", 0.9, 0.85, 0.95, 0.04, "Jones 2006",
    "NCT00000007", "pfs", 0.3, 0.20, 0.40, 0.05, "Jones 2001",
    "NCT00000008", "pfs", 0.9, 0.85, 0.95, 0.04, "Jones 2002",
    "NCT00000009", "pfs", 0.5, 0.40, 0.60, 0.05, "Jones 2003",
    "NCT00000010", "pfs", 0.6, 0.45, 0.85, 0.04, "Jones 2004",
    "NCT00000011", "pfs", 0.7, 0.55, 0.95, 0.04, "Jones 2005",
    "NCT00000012", "pfs", 0.9, 0.85, 0.95, 0.04, "Jones 2006"
)

pivotal <- pivotal %>%
    mutate(pivotal = TRUE)

nonpivotal <- nonpivotal %>%
    mutate(pivotal = FALSE)

trials <- pivotal %>%
    bind_rows(nonpivotal) %>%
    arrange(hr)

## Calculate SE
trials <- trials %>%
    mutate(se = (log(ci_upper) - log(ci_lower))/3.92)

## Split by endpoint
os_trials <- trials %>%
    filter(endpoint == "os")

pfs_trials <- trials %>%
    filter(endpoint == "pfs")

## Calculate meta-analyses and write to disk, unless it's already
## on-disk, in which case read the file and don't re-calculate

## (If you want to recalculate, you'll have to delete the .rds files)

if (! file.exists("os_meta.rds")) {

    os_meta <- bayesmeta(
        y = log(os_trials$hr),
        sigma = os_trials$se,
        labels = os_trials$studylabel,
        tau.prior = function(t){
            dhalfnormal(t, scale=0.5)
        }
    )
    
    pfs_meta <- bayesmeta(
        y = log(pfs_trials$hr),
        sigma = pfs_trials$se,
        labels = pfs_trials$studylabel,
        tau.prior = function(t){
            dhalfnormal(t, scale=0.5)
        }
    )
    
    saveRDS(os_meta, "os_meta.rds")
    saveRDS(os_trials, "os_trials.rds")
    
    saveRDS(os_meta, "pfs_meta.rds")
    saveRDS(os_trials, "pfs_trials.rds")
    
    message("Saved os_meta.rds, os_trials.rds, pfs_meta.rds and pfs_trial.rds")
} else {
    os_meta <- readRDS("os_meta.rds")
    os_trials <- readRDS("os_trials.rds")
    
    pfs_meta <- readRDS("pfs_meta.rds")
    pfs_trials <- readRDS("pfs_trials.rds")
    message("Read os_meta.rds, os_trials.rds, pfs_meta.rds and pfs_trial.rds from disk")
}

## OS
os_results <- tibble(
    study_label = os_meta$labels,
    original_hr = exp(os_meta$theta[1,1:os_meta$k]) %>% as.numeric(),
    original_ci_lower = os_trials$ci_lower,
    original_ci_upper = os_trials$ci_upper,
    shrinkage_hr = exp(os_meta$theta[5,1:os_meta$k]) %>% as.numeric(),
    shrinkage_ci_lower = exp(os_meta$theta[7,1:os_meta$k]) %>% as.numeric(),
    shrinkage_ci_upper = exp(os_meta$theta[8,1:os_meta$k]) %>% as.numeric()
) %>%
    mutate(rank = row_number())

os_results_original <- os_results %>%
    select(rank, study_label, original_hr, original_ci_lower, original_ci_upper) %>%
    rename(hr = original_hr) %>%
    rename(ci_lower = original_ci_lower) %>%
    rename(ci_upper = original_ci_upper) %>%
    mutate(estimate = "Original")

os_results_shrinkage <- os_results %>%
    select(rank, study_label, shrinkage_hr, shrinkage_ci_lower, shrinkage_ci_upper) %>%
    rename(hr = shrinkage_hr) %>%
    rename(ci_lower = shrinkage_ci_lower) %>%
    rename(ci_upper = shrinkage_ci_upper) %>%
    mutate(estimate = "Shrinkage")

os_plot_data <- os_results_original %>%
    bind_rows(os_results_shrinkage)

os_plot <- ggplot(
    aes(
        x = hr,
        y = rank,
        xmin = ci_lower,
        xmax = ci_upper,
        colour = estimate
    ),
    data = os_plot_data
) +
    geom_point() +
    geom_errorbar() +
    scale_y_discrete(
        limits = os_plot_data$rank,
        labels = os_plot_data$study_label
    ) +
    scale_x_continuous(
        breaks = seq(0, 3.5, 0.25),
        limits = c(0.25, 1.25)
    ) +
    labs(
        x = "Hazard ratio",
        y = "",
        colour = "Estimate",
        title = "Original and shrinkage-corrected OS effect size estimates for all trials of lung cancer therapy in our sample"
    ) +
    geom_vline(
        xintercept = 1
    )

pdf(
    "os-all.pdf",
    width = 10,
    height = 5
)
os_plot
dev.off()

## PFS
pfs_results <- tibble(
    study_label = pfs_meta$labels,
    original_hr = exp(pfs_meta$theta[1,1:pfs_meta$k]) %>% as.numeric(),
    original_ci_lower = pfs_trials$ci_lower,
    original_ci_upper = pfs_trials$ci_upper,
    shrinkage_hr = exp(pfs_meta$theta[5,1:pfs_meta$k]) %>% as.numeric(),
    shrinkage_ci_lower = exp(pfs_meta$theta[7,1:pfs_meta$k]) %>% as.numeric(),
    shrinkage_ci_upper = exp(pfs_meta$theta[8,1:pfs_meta$k]) %>% as.numeric()
) %>%
    mutate(rank = row_number())

pfs_results_original <- pfs_results %>%
    select(rank, study_label, original_hr, original_ci_lower, original_ci_upper) %>%
    rename(hr = original_hr) %>%
    rename(ci_lower = original_ci_lower) %>%
    rename(ci_upper = original_ci_upper) %>%
    mutate(estimate = "Original")

pfs_results_shrinkage <- pfs_results %>%
    select(rank, study_label, shrinkage_hr, shrinkage_ci_lower, shrinkage_ci_upper) %>%
    rename(hr = shrinkage_hr) %>%
    rename(ci_lower = shrinkage_ci_lower) %>%
    rename(ci_upper = shrinkage_ci_upper) %>%
    mutate(estimate = "Shrinkage")

pfs_plot_data <- pfs_results_original %>%
    bind_rows(pfs_results_shrinkage)

pfs_plot <- ggplot(
    aes(
        x = hr,
        y = rank,
        xmin = ci_lower,
        xmax = ci_upper,
        colour = estimate
    ),
    data = pfs_plot_data
) +
    geom_point() +
    geom_errorbar() +
    scale_y_discrete(
        limits = pfs_plot_data$rank,
        labels = pfs_plot_data$study_label
    ) +
    scale_x_continuous(
        breaks = seq(0, 3.5, 0.25),
        limits = c(0.25, 1.25)
    ) +
    labs(
        x = "Hazard ratio",
        y = "",
        colour = "Estimate",
        title = "Original and shrinkage-corrected PFS effect size estimates for all trials of lung cancer therapy in our sample"
    ) +
    geom_vline(
        xintercept = 1
    )

pdf(
    "pfs-all.pdf",
    width = 10,
    height = 5
)
pfs_plot
dev.off()

pdf(
    "combined.pdf",
    width = 10,
    height = 5
)
os_plot
pfs_plot
dev.off()
