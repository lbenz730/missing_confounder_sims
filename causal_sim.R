library(tidyverse)
library(glue)
library(mgcv)
library(ranger)
library(glmnet)

source('generate_data.R')
source('helpers.R')
source('estimate_treatment_effects.R')
source('levis_estimators.R')
source('levis_helpers.R')
source('compute_truth.R')

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
args <- commandArgs(trailingOnly = T)
sim_id <- as.numeric(args[1])

### Special Processing for Simulation 21
if(sim_id == 21) {
  plan(future::multisession(workers = n_cores))
  options(future.fork.enable = F)
  options(future.globals.maxSize = 8 * 1024^3)
}

t0 <- Sys.time()
cat('Running Simulation', sim_id, '\n')

params <- read_rds(glue('inputs/sim_params_{sim_id}.rds'))
n_patients <- params$n_patients
n_sims <- params$n_sims
dgp <- params$data_generating_process
mlp <- params$multiple_lp

### Process the data in batches
batch_size <- ifelse(n_sims == 1000, 1000, 500)
batches <- split(1:n_sims, ceiling(1:n_sims / batch_size))

cat('Compute true ATE given shifts\n')
### Compute true ATE given shifts
df_truth <- compute_truth(params)
params$true_ate <- df_truth$ate
params$true_te1 <- df_truth$te_1
params$true_te0 <- df_truth$te_0

### Generate Simulated datasets w/ missing data
cat('Generate Simulated datasets w/ missing data\n')
load_models(params)
sim_datasets <- 
  future_map(1:n_sims, ~generate_data(n_patients, dgp, mlp), 
             .options = furrr_options(seed = 8796))

### Impute missing confounder
###
### In case of |Lp| > 1, ols = linear regression for continuous, logistic regression
### for binary. ols_interact = same but w/ all pairwise interactions specified
cat('Impute Missing Confounder Data\n')
sim_datasets_imputed <- 
  list('true_dgp' = future_map(sim_datasets,
                               ~impute(.x, params$models, method = 'true_dgp', multiple_lp = mlp), 
                               .options = furrr_options(seed = 989)),
       'ols' = future_map(sim_datasets, 
                          ~impute(.x, params$models, method = 'lm', multiple_lp = mlp), 
                          .options = furrr_options(seed = 989)),
       'ols_interact' = future_map(sim_datasets, 
                                   ~impute(.x, params$models, method = 'lm_interact', multiple_lp = mlp), 
                                   .options = furrr_options(seed = 989)))

################################################################################
###################### Estimate Treatment Effects ##############################
################################################################################

cat('OLS Outcome Regressions\n')
######### Direct Modeling of Outcome + G-Formula ###########
### OLS Outcome Model (Various Imputation Models)
df_ols_cc <-
  batch_treatment_estimate(datasets = sim_datasets,
                           batches = batches,
                           method = 'OLS',
                           multiple_lp = mlp)  %>%
  mutate('dataset_id' = 1:n_sims) %>%
  mutate('imputation_model' = NA,
         'ipw_model' = NA,
         'interactions' = F,
         'outcome_model' = 'OLS',
         'interactions' = F)

df_ols_imputed_true <-
  batch_treatment_estimate(datasets = sim_datasets_imputed[['true_dgp']],
                           batches = batches,
                           method = 'OLS',
                           multiple_lp = mlp) %>%
  mutate('dataset_id' = 1:n_sims) %>%
  mutate('imputation_model' = 'True Gamma',
         'ipw_model' = NA,
         'outcome_model' = 'OLS',
         'interactions' = F)

df_ols_imputed_ols <-
  batch_treatment_estimate(datasets = sim_datasets_imputed[['ols']],
                           batches = batches,
                           method = 'OLS',
                           multiple_lp = mlp) %>%
  mutate('dataset_id' = 1:n_sims) %>%
  mutate('imputation_model' = 'OLS',
         'ipw_model' = NA,
         'outcome_model' = 'OLS',
         'interactions' = F)

df_ols_imputed_ols_interact <-
  batch_treatment_estimate(datasets = sim_datasets_imputed[['ols_interact']],
                           batches = batches,
                           method = 'OLS',
                           multiple_lp = mlp) %>%
  mutate('dataset_id' = 1:n_sims) %>%
  mutate('imputation_model' = 'OLS w/ Interactions',
         'ipw_model' = NA,
         'outcome_model' = 'OLS',
         'interactions' = F)

df_ols <-
  bind_rows(df_ols_cc,
            df_ols_imputed_true,
            df_ols_imputed_ols,
            df_ols_imputed_ols_interact)

cat('OLS Outcome Regressions w/ Interactions\n')
### OLS Outcome Model (w/ Interactions)
df_ols_interact_cc <-
  batch_treatment_estimate(datasets = sim_datasets,
                           batches = batches,
                           method = 'OLS',
                           interactions = T,
                           multiple_lp = mlp,
                           seed = 28732) %>%
  mutate('dataset_id' = 1:n_sims) %>%
  mutate('imputation_model' = NA,
         'ipw_model' = NA,
         'outcome_model' = 'OLS',
         'interactions' = T)

df_ols_interact_imputed_true <-
  batch_treatment_estimate(datasets = sim_datasets_imputed[['true_dgp']],
                           batches = batches,
                           method = 'OLS',
                           interactions = T,
                           multiple_lp = mlp,
                           seed = 91876) %>%
  mutate('dataset_id' = 1:n_sims) %>%
  mutate('imputation_model' = 'True Gamma',
         'ipw_model' = NA,
         'outcome_model' = 'OLS',
         'interactions' = T)

df_ols_interact_imputed_ols <-
  batch_treatment_estimate(datasets = sim_datasets_imputed[['ols']],
                           batches = batches,
                           method = 'OLS',
                           interactions = T,
                           multiple_lp = mlp,
                           seed = 91092) %>%
  mutate('dataset_id' = 1:n_sims) %>%
  mutate('imputation_model' = 'OLS',
         'ipw_model' = NA,
         'outcome_model' = 'OLS',
         'interactions' = T)

df_ols_interact_imputed_ols_interact <-
  batch_treatment_estimate(datasets = sim_datasets_imputed[['ols_interact']],
                           batches = batches,
                           method = 'OLS',
                           interactions = T,
                           multiple_lp = mlp,
                           seed = 88782) %>%
  mutate('dataset_id' = 1:n_sims) %>%
  mutate('imputation_model' = 'OLS w/ Interactions',
         'ipw_model' = NA,
         'outcome_model' = 'OLS',
         'interactions' = T)

df_ols_interact <-
  bind_rows(df_ols_interact_cc,
            df_ols_interact_imputed_true,
            df_ols_interact_imputed_ols,
            df_ols_interact_imputed_ols_interact)

cat('GAM Outcome Regressions\n')
### GAM Outcome Model (Various Imputation Models)
df_gam_cc <-
  batch_treatment_estimate(datasets = sim_datasets,
                           batches = batches,
                           method = 'gam',
                           multiple_lp = mlp) %>%
  mutate('dataset_id' = 1:n_sims) %>%
  mutate('imputation_model' = NA,
         'ipw_model' = NA,
         'outcome_model' = 'GAM',
         'interactions' = F)

df_gam_imputed_true <-
  batch_treatment_estimate(datasets = sim_datasets_imputed[['true_dgp']],
                           batches = batches,
                           method = 'gam',
                           multiple_lp = mlp) %>%
  mutate('dataset_id' = 1:n_sims) %>%
  mutate('imputation_model' = 'True Gamma',
         'ipw_model' = NA,
         'outcome_model' = 'GAM',
         'interactions' = F)

df_gam_imputed_ols <-
  batch_treatment_estimate(datasets = sim_datasets_imputed[['ols']],
                           batches = batches,
                           method = 'gam',
                           multiple_lp = mlp) %>%
  mutate('dataset_id' = 1:n_sims) %>%
  mutate('imputation_model' = 'OLS',
         'ipw_model' = NA,
         'outcome_model' = 'GAM',
         'interactions' = F)

df_gam_imputed_ols_interact <-
  batch_treatment_estimate(datasets = sim_datasets_imputed[['ols_interact']],
                           batches = batches,
                           method = 'gam',
                           multiple_lp = mlp) %>%
  mutate('dataset_id' = 1:n_sims) %>%
  mutate('imputation_model' = 'OLS w/ Interactions',
         'ipw_model' = NA,
         'outcome_model' = 'GAM',
         'interactions' = F)

df_gam <-
  bind_rows(df_gam_cc,
            df_gam_imputed_true,
            df_gam_imputed_ols,
            df_gam_imputed_ols_interact)

cat('GAM Outcome Regressions (w/ Interactions)\n')
### GAM Outcome Model (w/interactions)
df_gam_interact_cc <-
  batch_treatment_estimate(datasets = sim_datasets,
                           batches = batches,
                           method = 'gam',
                           multiple_lp = mlp,
                           interactions = T) %>%
  mutate('dataset_id' = 1:n_sims) %>%
  mutate('imputation_model' = NA,
         'ipw_model' = NA,
         'outcome_model' = 'GAM',
         'interactions' = T)

df_gam_interact_imputed_true <-
  batch_treatment_estimate(datasets = sim_datasets_imputed[['true_dgp']],
                           batches = batches,
                           method = 'gam',
                           multiple_lp = mlp,
                           interactions = T) %>%
  mutate('dataset_id' = 1:n_sims) %>%
  mutate('imputation_model' = 'True Gamma',
         'ipw_model' = NA,
         'outcome_model' = 'GAM',
         'interactions' = T)

df_gam_interact_imputed_ols <-
  batch_treatment_estimate(datasets = sim_datasets_imputed[['ols']],
                           batches = batches,
                           method = 'gam',
                           multiple_lp = mlp,
                           interactions = T) %>%
  mutate('dataset_id' = 1:n_sims) %>%
  mutate('imputation_model' = 'OLS',
         'ipw_model' = NA,
         'outcome_model' = 'GAM',
         'interactions' = T)

df_gam_interact_imputed_ols_interact <-
  batch_treatment_estimate(datasets = sim_datasets_imputed[['ols_interact']],
                           batches = batches,
                           method = 'gam',
                           multiple_lp = mlp,
                           interactions = T)%>%
  mutate('dataset_id' = 1:n_sims) %>%
  mutate('imputation_model' = 'OLS w/ Interactions',
         'ipw_model' = NA,
         'outcome_model' = 'GAM',
         'interactions' = T)

df_gam_interact <-
  bind_rows(df_gam_interact_cc,
            df_gam_interact_imputed_true,
            df_gam_interact_imputed_ols,
            df_gam_interact_imputed_ols_interact)

cat('Random Forest Outcome Regressions\n')
### Random Forest Outcome Model (Various Imputation Models)
df_rf_cc <-
  batch_treatment_estimate(datasets = sim_datasets,
                           batches = batches,
                           method = 'rf',
                           multiple_lp = mlp,
                           seed = 99991) %>%
  mutate('dataset_id' = 1:n_sims) %>%
  mutate('imputation_model' = NA,
         'ipw_model' = NA,
         'outcome_model' = 'Random Forest',
         'interactions' = NA)

df_rf_imputed_true <-
  batch_treatment_estimate(datasets = sim_datasets_imputed[['true_dgp']],
                           batches = batches,
                           method = 'rf',
                           multiple_lp = mlp,
                           seed = 134) %>%
  mutate('dataset_id' = 1:n_sims) %>%
  mutate('imputation_model' = 'True Gamma',
         'ipw_model' = NA,
         'outcome_model' = 'Random Forest',
         'interactions' = NA)

df_rf_imputed_ols <-
  batch_treatment_estimate(datasets = sim_datasets_imputed[['ols']],
                           batches = batches,
                           method = 'rf',
                           multiple_lp = mlp,
                           seed = 90901) %>%
  mutate('dataset_id' = 1:n_sims) %>%
  mutate('imputation_model' = 'OLS',
         'ipw_model' = NA,
         'outcome_model' = 'Random Forest',
         'interactions' = NA)

df_rf_imputed_ols_interact <-
  batch_treatment_estimate(datasets = sim_datasets_imputed[['ols_interact']],
                           batches = batches,
                           method = 'rf',
                           multiple_lp = mlp,
                           seed = 90912) %>%
  mutate('dataset_id' = 1:n_sims) %>%
  mutate('imputation_model' = 'OLS w/ Interactions',
         'ipw_model' = NA,
         'outcome_model' = 'Random Forest',
         'interactions' = NA)

df_rf <-
  bind_rows(df_rf_cc,
            df_rf_imputed_true,
            df_rf_imputed_ols,
            df_rf_imputed_ols_interact)

################## Inverse Probability Weighting ##############################
cat('IPW (Logistic Regressions)\n')
### Logistic Regressions
df_ipw_true_lr <-
  batch_treatment_estimate(datasets = sim_datasets_imputed[['true_dgp']],
                           batches = batches,
                           method = 'IPW',
                           weights = 'log_reg',
                           multiple_lp = mlp) %>%
  mutate('dataset_id' = 1:n_sims) %>%
  pivot_longer(cols = starts_with('H'),
               names_sep = '_',
               names_to = c('method', 'treatment'),
               values_to = 'effect') %>%
  pivot_wider(names_from = treatment,
              names_glue = 'te{treatment}_hat',
              values_from = effect) %>%
  mutate('ate_hat' = te1_hat - te0_hat) %>%
  mutate('imputation_model' = 'True Gamma',
         'ipw_model' = 'Logistic Regression',
         'outcome_model' = NA,
         'interactions' = F)

df_ipw_ols_lr <-
  batch_treatment_estimate(datasets = sim_datasets_imputed[['ols']],
                           batches = batches,
                           method = 'IPW',
                           weights = 'log_reg',
                           multiple_lp = mlp) %>%
  mutate('dataset_id' = 1:n_sims) %>%
  pivot_longer(cols = starts_with('H'),
               names_sep = '_',
               names_to = c('method', 'treatment'),
               values_to = 'effect') %>%
  pivot_wider(names_from = treatment,
              names_glue = 'te{treatment}_hat',
              values_from = effect) %>%
  mutate('ate_hat' = te1_hat - te0_hat) %>%
  mutate('imputation_model' = 'OLS',
         'ipw_model' = 'Logistic Regression',
         'outcome_model' = NA,
         'interactions' = F)

df_ipw_ols_interact_lr <-
  batch_treatment_estimate(datasets = sim_datasets_imputed[['ols_interact']],
                           batches = batches,
                           method = 'IPW',
                           weights = 'log_reg',
                           multiple_lp = mlp) %>%
  mutate('dataset_id' = 1:n_sims) %>%
  pivot_longer(cols = starts_with('H'),
               names_sep = '_',
               names_to = c('method', 'treatment'),
               values_to = 'effect') %>%
  pivot_wider(names_from = treatment,
              names_glue = 'te{treatment}_hat',
              values_from = effect) %>%
  mutate('ate_hat' = te1_hat - te0_hat) %>%
  mutate('imputation_model' = 'OLS w/ Interactions',
         'ipw_model' = 'Logistic Regression',
         'outcome_model' = NA,
         'interactions' = F)

df_ipw_cc_lr <-
  batch_treatment_estimate(datasets = sim_datasets,
                           batches = batches,
                           method = 'IPW',
                           weights = 'log_reg',
                           multiple_lp = mlp) %>%
  mutate('dataset_id' = 1:n_sims) %>%
  pivot_longer(cols = starts_with('H'),
               names_sep = '_',
               names_to = c('method', 'treatment'),
               values_to = 'effect') %>%
  pivot_wider(names_from = treatment,
              names_glue = 'te{treatment}_hat',
              values_from = effect) %>%
  mutate('ate_hat' = te1_hat - te0_hat) %>%
  mutate('imputation_model' = NA,
         'ipw_model' = 'Logistic Regression',
         'outcome_model' = NA,
         'interactions' = F)


cat('IPW (Logistic Regressions w/ Interactions)\n')
### Logistic Regressions w/ Interactions
df_ipw_true_lr_interact <-
  batch_treatment_estimate(datasets = sim_datasets_imputed[['true_dgp']],
                           batches = batches,
                           method = 'IPW',
                           weights = 'log_reg',
                           multiple_lp = mlp,
                           interactions = T,
                           seed = 101012) %>%
  mutate('dataset_id' = 1:n_sims) %>%
  pivot_longer(cols = starts_with('H'),
               names_sep = '_',
               names_to = c('method', 'treatment'),
               values_to = 'effect') %>%
  pivot_wider(names_from = treatment,
              names_glue = 'te{treatment}_hat',
              values_from = effect) %>%
  mutate('ate_hat' = te1_hat - te0_hat) %>%
  mutate('imputation_model' = 'True Gamma',
         'ipw_model' = 'Logistic Regression',
         'outcome_model' = NA,
         'interactions' = T)

df_ipw_ols_lr_interact <-
  batch_treatment_estimate(datasets = sim_datasets_imputed[['ols']],
                           batches = batches,
                           method = 'IPW',
                           weights = 'log_reg',
                           multiple_lp = mlp,
                           interactions = T,
                           seed = 12121) %>%
  mutate('dataset_id' = 1:n_sims) %>%
  pivot_longer(cols = starts_with('H'),
               names_sep = '_',
               names_to = c('method', 'treatment'),
               values_to = 'effect') %>%
  pivot_wider(names_from = treatment,
              names_glue = 'te{treatment}_hat',
              values_from = effect) %>%
  mutate('ate_hat' = te1_hat - te0_hat) %>%
  mutate('imputation_model' = 'OLS',
         'ipw_model' = 'Logistic Regression',
         'outcome_model' = NA,
         'interactions' = T)

df_ipw_ols_interact_lr_interact <-
  batch_treatment_estimate(datasets = sim_datasets_imputed[['ols_interact']],
                           batches = batches,
                           method = 'IPW',
                           weights = 'log_reg',
                           multiple_lp = mlp,
                           interactions = T,
                           seed = 911) %>%
  mutate('dataset_id' = 1:n_sims) %>%
  pivot_longer(cols = starts_with('H'),
               names_sep = '_',
               names_to = c('method', 'treatment'),
               values_to = 'effect') %>%
  pivot_wider(names_from = treatment,
              names_glue = 'te{treatment}_hat',
              values_from = effect) %>%
  mutate('ate_hat' = te1_hat - te0_hat) %>%
  mutate('imputation_model' = 'OLS w/ Interactions',
         'ipw_model' = 'Logistic Regression',
         'outcome_model' = NA,
         'interactions' = T)

df_ipw_cc_lr_interact <-
  batch_treatment_estimate(datasets = sim_datasets,
                           batches = batches,
                           method = 'IPW',
                           weights = 'log_reg',
                           multiple_lp = mlp,
                           interactions = T,
                           seed = 911) %>%
  mutate('dataset_id' = 1:n_sims) %>%
  pivot_longer(cols = starts_with('H'),
               names_sep = '_',
               names_to = c('method', 'treatment'),
               values_to = 'effect') %>%
  pivot_wider(names_from = treatment,
              names_glue = 'te{treatment}_hat',
              values_from = effect) %>%
  mutate('ate_hat' = te1_hat - te0_hat) %>%
  mutate('imputation_model' = NA,
         'ipw_model' = 'Logistic Regression',
         'outcome_model' = NA,
         'interactions' = T)


cat('IPW (GAM Logistic Regressions)\n')
### GAM Logistic Regressions
df_ipw_true_gam <-
  batch_treatment_estimate(datasets = sim_datasets_imputed[['true_dgp']],
                           batches = batches,
                           method = 'IPW',
                           weights = 'gam',
                           multiple_lp = mlp) %>%
  mutate('dataset_id' = 1:n_sims) %>%
  pivot_longer(cols = starts_with('H'),
               names_sep = '_',
               names_to = c('method', 'treatment'),
               values_to = 'effect') %>%
  pivot_wider(names_from = treatment,
              names_glue = 'te{treatment}_hat',
              values_from = effect) %>%
  mutate('ate_hat' = te1_hat - te0_hat) %>%
  mutate('imputation_model' = 'True Gamma',
         'ipw_model' = 'GAM',
         'outcome_model' = NA,
         'interactions' = F)

df_ipw_ols_gam <-
  batch_treatment_estimate(datasets = sim_datasets_imputed[['ols']],
                           batches = batches,
                           method = 'IPW',
                           weights = 'gam',
                           multiple_lp = mlp) %>%
  mutate('dataset_id' = 1:n_sims) %>%
  pivot_longer(cols = starts_with('H'),
               names_sep = '_',
               names_to = c('method', 'treatment'),
               values_to = 'effect') %>%
  pivot_wider(names_from = treatment,
              names_glue = 'te{treatment}_hat',
              values_from = effect) %>%
  mutate('ate_hat' = te1_hat - te0_hat) %>%
  mutate('imputation_model' = 'OLS',
         'ipw_model' = 'GAM',
         'outcome_model' = NA,
         'interactions' = F)

df_ipw_ols_interact_gam <-
  batch_treatment_estimate(datasets = sim_datasets_imputed[['ols_interact']],
                           batches = batches,
                           method = 'IPW',
                           weights = 'gam',
                           multiple_lp = mlp) %>%
  mutate('dataset_id' = 1:n_sims) %>%
  pivot_longer(cols = starts_with('H'),
               names_sep = '_',
               names_to = c('method', 'treatment'),
               values_to = 'effect') %>%
  pivot_wider(names_from = treatment,
              names_glue = 'te{treatment}_hat',
              values_from = effect) %>%
  mutate('ate_hat' = te1_hat - te0_hat) %>%
  mutate('imputation_model' = 'OLS w/ Interactions',
         'ipw_model' = 'GAM',
         'outcome_model' = NA,
         'interactions' = F)

df_ipw_cc_gam <-
  batch_treatment_estimate(datasets = sim_datasets,
                           batches = batches,
                           method = 'IPW',
                           weights = 'gam',
                           multiple_lp = mlp) %>%
  mutate('dataset_id' = 1:n_sims) %>%
  pivot_longer(cols = starts_with('H'),
               names_sep = '_',
               names_to = c('method', 'treatment'),
               values_to = 'effect') %>%
  pivot_wider(names_from = treatment,
              names_glue = 'te{treatment}_hat',
              values_from = effect) %>%
  mutate('ate_hat' = te1_hat - te0_hat) %>%
  mutate('imputation_model' = NA,
         'ipw_model' = 'GAM',
         'outcome_model' = NA,
         'interactions' = F)


cat('IPW (GAM Logistic Regressions w/ Interactions)\n')
### GAM Logistic Regressions
df_ipw_true_gam_interact <- 
  batch_treatment_estimate(datasets = sim_datasets_imputed[['true_dgp']],
                           batches = batches,
                           method = 'IPW',
                           weights = 'gam', 
                           interaction = T,
                           multiple_lp = mlp) %>% 
  mutate('dataset_id' = 1:n_sims) %>% 
  pivot_longer(cols = starts_with('H'),
               names_sep = '_',
               names_to = c('method', 'treatment'),
               values_to = 'effect') %>% 
  pivot_wider(names_from = treatment, 
              names_glue = 'te{treatment}_hat',
              values_from = effect) %>%   
  mutate('ate_hat' = te1_hat - te0_hat) %>% 
  mutate('imputation_model' = 'True Gamma',
         'ipw_model' = 'GAM',
         'outcome_model' = NA,
         'interactions' = T)

df_ipw_ols_gam_interact <- 
  batch_treatment_estimate(datasets = sim_datasets_imputed[['ols']],
                           batches = batches,
                           method = 'IPW',
                           weights = 'gam', 
                           interaction = T,
                           multiple_lp = mlp) %>% 
  mutate('dataset_id' = 1:n_sims) %>% 
  pivot_longer(cols = starts_with('H'),
               names_sep = '_',
               names_to = c('method', 'treatment'),
               values_to = 'effect') %>% 
  pivot_wider(names_from = treatment, 
              names_glue = 'te{treatment}_hat',
              values_from = effect) %>%   
  mutate('ate_hat' = te1_hat - te0_hat) %>% 
  mutate('imputation_model' = 'OLS',
         'ipw_model' = 'GAM',
         'outcome_model' = NA,
         'interactions' = T)


df_ipw_ols_interact_gam_interact <- 
  batch_treatment_estimate(datasets = sim_datasets_imputed[['ols_interact']],
                           batches = batches,
                           method = 'IPW',
                           weights = 'gam', 
                           interaction = T,
                           multiple_lp = mlp) %>% 
  mutate('dataset_id' = 1:n_sims) %>% 
  pivot_longer(cols = starts_with('H'),
               names_sep = '_',
               names_to = c('method', 'treatment'),
               values_to = 'effect') %>% 
  pivot_wider(names_from = treatment, 
              names_glue = 'te{treatment}_hat',
              values_from = effect) %>%   
  mutate('ate_hat' = te1_hat - te0_hat) %>% 
  mutate('imputation_model' = 'OLS w/ Interactions',
         'ipw_model' = 'GAM',
         'outcome_model' = NA,
         'interactions' = T)

df_ipw_cc_gam_interact <- 
  batch_treatment_estimate(datasets = sim_datasets,
                           batches = batches,
                           method = 'IPW',
                           weights = 'gam', 
                           interaction = T,
                           multiple_lp = mlp) %>% 
  mutate('dataset_id' = 1:n_sims) %>% 
  pivot_longer(cols = starts_with('H'),
               names_sep = '_',
               names_to = c('method', 'treatment'),
               values_to = 'effect') %>% 
  pivot_wider(names_from = treatment, 
              names_glue = 'te{treatment}_hat',
              values_from = effect) %>%   
  mutate('ate_hat' = te1_hat - te0_hat) %>% 
  mutate('imputation_model' = NA,
         'ipw_model' = 'GAM',
         'outcome_model' = NA,
         'interactions' = T)


cat('IPW (Random Forest Classifier)\n')
### Random Forest Classifier
df_ipw_true_rf <- 
  batch_treatment_estimate(datasets = sim_datasets_imputed[['true_dgp']],
                           batches = batches,
                           method = 'IPW',
                           weights = 'rf', 
                           multiple_lp = mlp,
                           seed = 12312) %>% 
  mutate('dataset_id' = 1:n_sims) %>% 
  pivot_longer(cols = starts_with('H'),
               names_sep = '_',
               names_to = c('method', 'treatment'),
               values_to = 'effect') %>% 
  pivot_wider(names_from = treatment, 
              names_glue = 'te{treatment}_hat',
              values_from = effect) %>%   
  mutate('ate_hat' = te1_hat - te0_hat) %>% 
  mutate('imputation_model' = 'True Gamma',
         'ipw_model' = 'Random Forest',
         'outcome_model' = NA,
         'interactions' = NA)

df_ipw_ols_rf <- 
  batch_treatment_estimate(datasets = sim_datasets_imputed[['ols']],
                           batches = batches,
                           method = 'IPW',
                           weights = 'rf', 
                           multiple_lp = mlp,
                           seed = 11) %>% 
  mutate('dataset_id' = 1:n_sims) %>% 
  pivot_longer(cols = starts_with('H'),
               names_sep = '_',
               names_to = c('method', 'treatment'),
               values_to = 'effect') %>% 
  pivot_wider(names_from = treatment, 
              names_glue = 'te{treatment}_hat',
              values_from = effect) %>%   
  mutate('ate_hat' = te1_hat - te0_hat) %>% 
  mutate('imputation_model' = 'OLS',
         'ipw_model' = 'Random Forest',
         'outcome_model' = NA,
         'interactions' = NA)

df_ipw_ols_interact_rf <- 
  batch_treatment_estimate(datasets = sim_datasets_imputed[['ols_interact']],
                           batches = batches,
                           method = 'IPW',
                           weights = 'rf', 
                           multiple_lp = mlp,
                           seed = 1231) %>% 
  mutate('dataset_id' = 1:n_sims) %>% 
  pivot_longer(cols = starts_with('H'),
               names_sep = '_',
               names_to = c('method', 'treatment'),
               values_to = 'effect') %>% 
  pivot_wider(names_from = treatment, 
              names_glue = 'te{treatment}_hat',
              values_from = effect) %>%   
  mutate('ate_hat' = te1_hat - te0_hat) %>% 
  mutate('imputation_model' = 'OLS w/ Interactions',
         'ipw_model' = 'Random Forest',
         'outcome_model' = NA,
         'interactions' = NA)

df_ipw_cc_rf <- 
  batch_treatment_estimate(datasets = sim_datasets,
                           batches = batches,
                           method = 'IPW',
                           weights = 'rf', 
                           multiple_lp = mlp,
                           seed = 1231) %>% 
  mutate('dataset_id' = 1:n_sims) %>% 
  pivot_longer(cols = starts_with('H'),
               names_sep = '_',
               names_to = c('method', 'treatment'),
               values_to = 'effect') %>% 
  pivot_wider(names_from = treatment, 
              names_glue = 'te{treatment}_hat',
              values_from = effect) %>%   
  mutate('ate_hat' = te1_hat - te0_hat) %>% 
  mutate('imputation_model' = NA,
         'ipw_model' = 'Random Forest',
         'outcome_model' = NA,
         'interactions' = NA)

df_ipw <- 
  bind_rows(df_ipw_true_lr, df_ipw_ols_lr, df_ipw_ols_interact_lr, df_ipw_cc_lr,
            df_ipw_true_lr_interact, df_ipw_ols_lr_interact, df_ipw_ols_interact_lr_interact, df_ipw_cc_lr_interact,
            df_ipw_true_gam, df_ipw_ols_gam, df_ipw_ols_interact_gam, df_ipw_cc_gam,
            df_ipw_true_gam_interact, df_ipw_ols_gam_interact, df_ipw_ols_interact_gam_interact, df_ipw_cc_gam_interact,
            df_ipw_true_rf, df_ipw_ols_rf, df_ipw_ols_interact_rf, df_ipw_cc_rf)


cat('Levis Estimators\n')
############ Levis Estimators #############
df_levis <- 
  batch_levis_estimators(datasets = sim_datasets,
                         batches = batches,
                         model_types = params$models,
                         multiple_lp = mlp,
                         seed = 1) %>% 
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


cat('Save Results\n')
### Combine Estimates and Save Outputs
# Collect in each section
df_estimates <- 
  bind_rows(df_ols, df_ols_interact, 
            df_gam, df_gam_interact,
            df_rf,
            df_ipw, 
            df_levis) %>% 
  bind_cols(as_tibble(params[setdiff(names(params), c('shifts', 'models'))]) %>% dplyr::slice(1)) 

write_rds(df_estimates, glue('outputs/sim_results_{sim_id}.rds'))

t1 <- Sys.time()
cat('\n----------------------\n')
cat('Simulation completed:', difftime(t1, t0, units='mins'), 'miutes\n')
cat('----------------------\n')
