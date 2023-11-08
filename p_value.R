library(tidyverse)
library(knitr)

pivotal <- read_csv("nccn_pilot_clean_formatted.csv")

nonpivotal <- read_csv("nccn_dummy_data_final.csv")

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

## FWER p-value recalculation
B = 1000
boot.p = matrix(NA, nrow = nrow(trials), ncol = B)
res.vec = matrix(NA, nrow = nrow(trials), ncol = B)

for(i in 1:B){
  p.bt <- rlnorm(
    n=nrow(trials),
    meanlog=0,
    sdlog=mean(trials$se)
  )
  boot.p[,i] = pmin(1, 2*(1-p.bt))
  res.vec[,i] = trials$hr > p.bt
}

trials$fwer_p <- pmin(1, 2*apply(res.vec, 1, sum)/B)

trials %>%
    arrange(desc(hr)) %>%
    filter(pivotal) %>%
    mutate(
        original_z = abs(log(hr))/se
    ) %>%
    mutate(
        original_p = 2 * (1-pnorm(original_z))
    ) %>%
    select(studylabel, original_p, fwer_p) %>%
    mutate(
        dangerzone = ifelse (
            original_p <= 0.05 & fwer_p > 0.05,
            "Danger",
            ""
        )
    ) %>%
    kable(col.names=c(
              "Pivotal trial and FDA label",
              "Original PFS p-value",
              "FWER PFS p-value",
              "Danger zone"))

trials %>%
    write_csv("fwer_pval.csv")
