## ---- include = FALSE---------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----setup--------------------------------------------------------------------
library(powerHaDeX)

## -----------------------------------------------------------------------------
spectrum1 <- simulate_theoretical_spectra("LVRKDLQN", protection_factor = 10, charge = 2, times = 60)
head(spectrum1)

## -----------------------------------------------------------------------------

spectrum2 <- simulate_theoretical_spectra("LVRKDLQN", protection_factor = 1000, charge = 2, times = 60)
paired_spectra <- rbind(spectrum1, spectrum2)

get_noisy_deuteration_curves(paired_spectra,
                             n_replicates = 4,
                             n_experiments = 1,
                             reference = "all")[[1]][[1]]


