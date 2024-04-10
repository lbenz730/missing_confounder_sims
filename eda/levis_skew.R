### Script to investigate Skewness in Levis IWOR Estimator
library(tidyverse)
library(glue)
library(mgcv)
library(here)

source(here('generate_data.R'))
source(here('helpers.R'))
source(here('estimate_treatment_effects.R'))
source(here('levis_estimators.R'))
source(here('levis_helpers.R'))
source(here('compute_truth.R'))

################################################################################
###################### Simulation Set Up #######################################
################################################################################
### Set up Parallelization
library(furrr)
n_cores <- 32
plan(future::multicore(workers = n_cores))
options(future.fork.enable = T)
options(future.globals.maxSize = 8 * 1024^3)

### Simulation Parameters
sim_id <- 4
params <- read_rds(here(glue('inputs/sim_params_{sim_id}.rds')))
n_patients <- params$n_patients
n_sims <- params$n_sims
dgp <- params$data_generating_process

### Compute true ATE given shifts
df_truth <- compute_truth(params)
params$true_ate <- df_truth$ate
params$true_te1 <- df_truth$te_1
params$true_te0 <- df_truth$te_0

### Generate Simulated datasets w/ missing data
load_models(params)

df_estimates <- NULL

patients <- c(5000, 20000, 50000, 100000)
for(i in 1:length(patients)) {
  cat('Checking', i, '\n')
  n_patients <- patients[i]
  
  sim_datasets <- future_map(1:n_sims, ~generate_data(n_patients, 
                                                      dgp, 
                                                      multiple_lp = params$multiple_lp), 
                             .options = furrr_options(seed = 8796))
  
  
  ############ Levis Estimators #############
  df_levis <- 
    future_map_dfr(sim_datasets, ~compute_levis_estimators(.x, 
                                                           model_types = params$models, 
                                                           multiple_lp = params$multiple_lp)) %>% 
    mutate('dataset_id' = 1:n_sims) %>% 
    pivot_longer(cols = starts_with('I'),
                 names_sep = '_',
                 names_to = c('method', 'treatment'),
                 values_to = 'effect') %>% 
    pivot_wider(names_from = treatment, 
                names_glue = 'te{treatment}_hat',
                values_from = effect) %>% 
    mutate('ate_hat' = te1_hat - te0_hat) %>% 
    mutate('method' = paste('Levis', method))
  
  
  ### Combine Estimates and Save Outputs
  # Collect in each section
  df_estimates <- 
    df_estimates %>% 
    bind_rows(df_levis %>% 
                bind_cols(as_tibble(params[setdiff(names(params), c('shifts', 'models'))]) %>% dplyr::slice(1)) %>%  
                mutate('n_patients' = patients[i]))
  
}

write_csv(df_estimates, here('eda/levis_skew.csv'))
