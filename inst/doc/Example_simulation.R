## ---- include = FALSE---------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  fig.width = 5,
  fig.align = 'center'
)

## ----setup, warning=FALSE-----------------------------------------------------
library(powerHaDeX)

## -----------------------------------------------------------------------------
set.seed(17)

theo_spectrum <- simulate_theoretical_spectra(sequence = "CHERICHERILADY",
                                              charge = 4,
                                              protection_factor = 100,
                                              times = 0.167,
                                              pH = 7.5,
                                              temperature = 15,
                                              n_molecules = 500,
                                              time_step_const = 1,
                                              use_markov = TRUE)
theo_spectrum

## -----------------------------------------------------------------------------
plot_spectra(theo_spectrum)

## -----------------------------------------------------------------------------
plot_spectra(theo_spectrum, control_time = TRUE)

## -----------------------------------------------------------------------------
set.seed(17)

theo_spectra <- simulate_theoretical_spectra(sequence = "CHERICHERILADY",
                                             charge = c(3, 5),
                                             protection_factor = 100,
                                             times = c(0.167, 5),
                                             pH = 7.5,
                                             temperature = 15,
                                             n_molecules = 500,
                                             time_step_const = 1,
                                             use_markov = TRUE)

head(theo_spectra)

## -----------------------------------------------------------------------------
plot_spectra(theo_spectra)

## -----------------------------------------------------------------------------
undeuterated_mass = get_undeuterated_mass(theo_spectra)
spectra = get_spectra_list(theo_spectra)
replicated_spectra = add_noise_to_spectra(spectra,
                                          undeuterated_mass = undeuterated_mass,
                                          n_experiments = 2)
replicated_spectra

## -----------------------------------------------------------------------------

theo_spectra_pf_100 <- theo_spectra
theo_spectra_pf_200 <- simulate_theoretical_spectra(sequence = "CHERICHERILADY",
                                                    charge = c(3, 5),
                                                    protection_factor = 200,
                                                    times = c(0.167, 5),
                                                    pH = 7.5,
                                                    temperature = 15,
                                                    n_molecules = 500,
                                                    time_step_const = 1,
                                                    use_markov = TRUE)

theo_spectra_two_states <- rbind(theo_spectra_pf_100, theo_spectra_pf_200)


## -----------------------------------------------------------------------------

undeuterated_mass = get_undeuterated_mass(theo_spectra_two_states)
spectra = get_spectra_list(theo_spectra_two_states,
                           compare_pairs = TRUE,
                           reference = "all")
replicated_spectra_paired = add_noise_to_spectra(spectra,
                                                 undeuterated_mass = undeuterated_mass,
                                                 n_experiments = 2)
replicated_spectra_paired

## -----------------------------------------------------------------------------
deuteration_curves <- get_deuteration_curves_from_spectra(replicated_spectra)
deuteration_curves

## -----------------------------------------------------------------------------
get_noisy_deuteration_curves(theo_spectra,
                             n_replicates = 4,
                             n_experiments = 2,
                             compare_pairs = FALSE)

## -----------------------------------------------------------------------------
deuteration_curves_paired_states <- get_noisy_deuteration_curves(theo_spectra_two_states,
                                                                 n_replicates = 4,
                                                                 n_experiments = 2,
                                                                 compare_pairs = TRUE,
                                                                 reference = "all")
deuteration_curves_paired_states

## -----------------------------------------------------------------------------
calculate_hdx_power(deuteration_curves_paired_states,
                    tests = list(test_houde),
                    summarized = FALSE)

## -----------------------------------------------------------------------------
calculate_hdx_power(deuteration_curves_paired_states,
                    tests = list(test_houde),
                    summarized = TRUE)

## ---- eval = FALSE------------------------------------------------------------
#  example_test <- function(data, significance_level) {
#    States = unique(data$State)
#  
#    # testing procedure here
#  
#    return(data.table::data.table(Test = "Example test",
#                                  State_1 = States[1],
#                                  State_2 = States[2],
#                                  Test_statistic = NA,
#                                  P_value = NA,
#                                  Significant_difference = #TRUE or FALSE,
#                                    Time = NA,
#                                  Transformation = NA,
#                                  AIC = NA,
#                                  logLik = NA))
#  }
#  

