library(dplyr)
library(purrr)
library(readr)

set.seed(158)

### Simulation parameters
n_subjects <- 5000
n_sims <- 1000
n_folds <- 2

n <- n_sims * n_subjects 
batch_size <- 2500
n_batches <- n/batch_size

### Set Up data
A <- rbinom(n, 1, 0.5)
Y <- rbeta(n, shape1 = ifelse(A == 0, 2, 4), shape2 = ifelse(A == 0, 4, 2)) 
S <- rbinom(n, 1, plogis(-0.35 + 0.5 * A + 0.18 * Y + 0.05 * A * Y))
L1 <- rbinom(n, 1, plogis(-0.6 + 0.5 * A + 0.25 * Y + 0.1 * A * Y))
L2 <- rnorm(n, 1.0 * A + Y + 2.5 * L1 * Y, 1.25)

df_sim <- 
  tibble('dataset_id' = rep(1:n_sims, each = n_subjects),
         'subject_id' = rep(1:n_subjects, n_sims),
         'fold_id' = rep(rep(1:n_folds,  each = n_subjects/2), n_sims),
         'batch_id' = rep(1:n_batches, each = batch_size),
         'A' = A,
         'Y' = Y,
         'S' = S,
         'L1' = L1,
         'L2' = L2)

### Save data out one fold at a time
folds <- 
  df_sim %>% 
  group_by(dataset_id, fold_id) %>% 
  group_split() 

df_ids <- 
  df_sim %>% 
  distinct(dataset_id, fold_id)

for(i in 1:nrow(df_ids)) {
  dataset_id <- df_ids$dataset_id[i]
  fold_id <- df_ids$fold_id[i]
  write_csv(folds[[i]], paste0('input_data/dataset_', dataset_id, '_fold_', fold_id, '.csv'))
}

### Key 
df_key <- 
  df_sim %>% 
  mutate('subject_uuid' = 1:n) %>% 
  distinct(dataset_id, fold_id, batch_id, subject_id, subject_uuid)

write_csv(df_key, 'subject_key.csv')
  