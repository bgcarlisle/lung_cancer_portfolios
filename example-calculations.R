library(tidyverse)
library(bayesmeta)
library(ggplot2)

pivotal <- tribble(
    ~nctid, ~endpoint, ~hr, ~ci_lower, ~ci_upper, ~p_value, ~studylabel,
    "NCT00000001", "os", 0.9, 0.85, 0.95, 0.04, "Smith 2001",
    "NCT00000002", "os", 0.9, 0.75, 0.95, 0.04, "Smith 2002",
    "NCT00000003", "os", 0.5, 0.65, 0.95, 0.04, "Smith 2003",
    "NCT00000004", "os", 0.4, 0.75, 0.95, 0.04, "Smith 2004",
    "NCT00000005", "os", 0.7, 0.85, 0.95, 0.04, "Smith 2005",
    "NCT00000006", "os", 0.7, 0.85, 0.95, 0.04, "Smith 2006"
)

nonpivotal <- tribble(
    ~nctid, ~endpoint, ~hr, ~ci_lower, ~ci_upper, ~p_value, ~studylabel,
    "NCT00000007", "os", 0.9, 0.85, 0.95, 0.04, "Jones 2001",
    "NCT00000008", "os", 0.9, 0.75, 0.95, 0.04, "Jones 2002",
    "NCT00000009", "os", 0.5, 0.65, 0.95, 0.04, "Jones 2003",
    "NCT00000010", "os", 0.4, 0.75, 0.95, 0.04, "Jones 2004",
    "NCT00000011", "os", 0.7, 0.85, 0.95, 0.04, "Jones 2005",
    "NCT00000012", "os", 0.7, 0.85, 0.95, 0.04, "Jones 2006"
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
    saveRDS(os_meta, "os_meta.rds")
    saveRDS(os_trials, "os_trials.rds")
    message("Saved os_meta.rds and os_trials.rds")
} else {
    os_meta <- readRDS("os_meta.rds")
    os_trials <- readRDS("os_trials.rds")
    message("Read os_meta.rds and os_trials.rds from disk")
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
        limits = c(0, 1.25)
    ) +
    labs(
        x = "Hazard ratio",
        y = "",
        colour = "Estimate",
        title = "Original and shrinkage-corrected OS effect size estimates for all trials of lung cancer therapy in our sample"
    ) +
    geom_vline(
        xintercept = 1
    ) +
    geom_vline(
        xintercept = 0.75,
        linetype = "dashed"
    ) +
    geom_vline(
        xintercept = 0.7,
        linetype = "dashed"
    ) +
    geom_vline(
        xintercept = 0.65,
        linetype = "dashed"
    )

pdf(
    "os-all.pdf",
    width = 10,
    height = 5
)
os_plot
dev.off()
