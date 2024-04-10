library(tidyverse)

args <- commandArgs(TRUE)
b_id <- as.numeric(args[1])

### Load the input data corresponding to this batch 
df_key <- read_csv('subject_key.csv')
subject_info <- 
  df_key %>% 
  filter(batch_id == b_id) 

dataset_id <- subject_info$dataset_id[1]
fold_id <- subject_info$fold_id[1]


df <- read_csv(paste0('input_data/dataset_', dataset_id, '_fold_', fold_id, '.csv'))
load(paste0('model_objects/objects_dataset_', dataset_id, '_fold_', fold_id, '_parametric.Rdata'))
contrib <- NULL

### Pre-Compute Before Any Subjects
delta.lambda <- 0.0555
lambda.grid <- seq(-3.3, 7.8, by = delta.lambda)

gamma.hyp.1.0 <- rep(NA, nrow = length(lambda.grid))
beta.hyp.1.0 <- rep(NA, nrow = length(lambda.grid))
gamma.hyp.0.0 <- rep(NA, nrow = length(lambda.grid))
beta.hyp.0.0 <- rep(NA, nrow = length(lambda.grid))
gamma.hyp.1.1 <- rep(NA, nrow = length(lambda.grid))
beta.hyp.1.1 <- rep(NA, nrow = length(lambda.grid))
gamma.hyp.0.1 <- rep(NA, nrow = length(lambda.grid))
beta.hyp.0.1 <- rep(NA, nrow = length(lambda.grid))

for (j in 1:length(lambda.grid)) {
  # L1 = 1, A = 1
  lambda.hyp.hat <-
    lambda1.pmf(a = 1, y = y.grid, l1 = 1) * 
    lambda2.pdf(a = 1, y = y.grid, l1 = 1, l2 = lambda.grid[j])
  gamma.hyp.1.1[j] <- delta.y * sum(lambda.hyp.hat * mu.pred.1)
  beta.hyp.1.1[j] <- delta.y * sum(y.grid * lambda.hyp.hat * mu.pred.1)
  
  # L1 = 1, A = 0
  lambda.hyp.hat <-
    lambda1.pmf(a = 0, y = y.grid, l1 = 1) * 
    lambda2.pdf(a = 0, y = y.grid, l1 = 1, l2 = lambda.grid[j])
  gamma.hyp.1.0[j] <- delta.y * sum(lambda.hyp.hat * mu.pred.0)
  beta.hyp.1.0[j] <- delta.y * sum(y.grid * lambda.hyp.hat * mu.pred.0)
  
  # L1 = 0, A = 1
  lambda.hyp.hat <-
    lambda1.pmf(a = 1, y = y.grid, l1 = 0) * 
    lambda2.pdf(a = 1, y = y.grid, l1 = 0, l2 = lambda.grid[j])
  gamma.hyp.0.1[j] <- delta.y * sum(lambda.hyp.hat * mu.pred.1)
  beta.hyp.0.1[j] <- delta.y * sum(y.grid * lambda.hyp.hat * mu.pred.1)
  
  # L1 = 0, A = 0
  lambda.hyp.hat <-
    lambda1.pmf(a = 0, y = y.grid, l1 = 0) * 
    lambda2.pdf(a = 0, y = y.grid, l1 = 0, l2 = lambda.grid[j])
  gamma.hyp.0.0[j] <- delta.y * sum(lambda.hyp.hat * mu.pred.0)
  beta.hyp.0.0[j] <- delta.y * sum(y.grid * lambda.hyp.hat * mu.pred.0)
}


for(i in 1:nrow(subject_info)) {
  subject_id <- subject_info$subject_id[i]
  cat('Subject', subject_id, '\n')
  subject_uuid <- subject_info$subject_id[i]
  ix <- which(df$subject_id == subject_id)
  
  b.12 <- b.02 <- 0
  lambda.obs.hat.0 <- 
    lambda2.pdf(a = df$A[ix], y = df$Y[ix], l1 = 0, l2 = lambda.grid)
  
  lambda.obs.hat.1 <- 
    lambda2.pdf(a = df$A[ix], y = df$Y[ix], l1 = 1, l2 = lambda.grid)
  
  b.01 <- delta.lambda * 
    sum((1 - lambda1.hat[ix]) * lambda.obs.hat.0 * beta.hyp.0.0 / 
          gamma.hyp.0.0 +
          lambda1.hat[ix] * lambda.obs.hat.1 * beta.hyp.1.0 / 
          gamma.hyp.1.0)
  b.11 <- delta.lambda * 
    sum((1 - lambda1.hat[ix]) * lambda.obs.hat.0 * beta.hyp.0.1 / 
          gamma.hyp.0.1 +
          lambda1.hat[i] * lambda.obs.hat.1 * beta.hyp.1.1 / 
          gamma.hyp.1.1)
  
  if (df$A[ix] == 0) {
    b.02 <- delta.lambda * 
      sum((1 - lambda1.hat[ix]) * lambda.obs.hat.0 * 
            ((1 - eta.hat) * (gamma.hyp.0.0) + 
               eta.hat * (gamma.hyp.0.1)) *
            (df$Y[ix] - beta.hyp.0.0 / gamma.hyp.0.0) /
            gamma.hyp.0.0 +
            lambda1.hat[ix] * lambda.obs.hat.1 *
            ((1 - eta.hat) * (gamma.hyp.1.0) + 
               eta.hat * (gamma.hyp.1.1)) *
            (df$Y[ix] - beta.hyp.1.0 / gamma.hyp.1.0) /
            gamma.hyp.1.0)
  } else {
    b.12 <- delta.lambda * 
      sum((1 - lambda1.hat[ix]) * lambda.obs.hat.0 * 
            ((1 - eta.hat) * (gamma.hyp.0.0) + 
               eta.hat * (gamma.hyp.0.1)) *
            (df$Y[ix] - beta.hyp.0.1 / gamma.hyp.0.1) /
            gamma.hyp.0.1 +
            lambda1.hat[ix] * lambda.obs.hat.1 *
            ((1 - eta.hat) * (gamma.hyp.1.0) + 
               eta.hat * (gamma.hyp.1.1)) *
            (df$Y[i] - beta.hyp.1.1 / gamma.hyp.1.1) /
            gamma.hyp.1.1)
  }
  
  
  
  
  
  export <- 
    data.frame(dataset_id = dataset_id,
               fold_id = fold_id,
               batch_id = b_id,
               subject_id = subject_id,
               subject_uuid = subject_uuid,
               b.01 = b.01,
               b.11 = b.11,
               b.02 = b.02, 
               b.12 = b.12)
  
  contrib <- bind_rows(contrib, export)
}

### Save Contributions
gamma.hat.1 <- sapply(1:n, function(i) {
  ifelse(df$S[i] == 0, 0,
         delta.y * sum(lambda.obs.mat.1[,i] * mu.pred.1))
}, simplify = 0)
gamma.hat.0 <- sapply(1:n, function(i) {
  ifelse(df$S[i] == 0, 0,
         delta.y * sum(lambda.obs.mat.0[,i] * mu.pred.0))
}, simplify = 0)

beta.hat.1 <- sapply(1:n, function(i) {
  ifelse(df$S[i] == 0, 0,
         delta.y * sum(y.grid * lambda.obs.mat.1[,i] * mu.pred.1))
}, simplify = 0)
beta.hat.0 <- sapply(1:n, function(i) {
  ifelse(df$S[i] == 0, 0,
         delta.y * sum(y.grid * lambda.obs.mat.0[,i] * mu.pred.0))
}, simplify = 0)



contrib$gamma.hat.0 <- gamma.hat.0
contrib$gamma.hat.1 <- gamma.hat.1
contrib$beta.hat.0 <- beta.hat.0
contrib$beta.hat.1 <- beta.hat.1
contrib$pi.hat <- pi.hat.1
contrib$eta.hat <- eta.hat


write_csv(contrib, paste0("subject_contributions/parametric_batch_", b_id, ".csv"))