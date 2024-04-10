library(haldensify)
library(SuperLearner)
library(statmod)
library(tidyverse)
library(data.table)

### Which dataset to work on
args <- commandArgs(trailingOnly = T)
dataset_id <- as.numeric(args[1])
n_folds <- 2

### Fit Models Seperately for Each Fold
for(fold in 1:n_folds) {
  df_train <- read_csv(paste0('input_data/dataset_', dataset_id, '_fold_', fold, '.csv'))
  df_test <- read_csv(paste0('input_data/dataset_', dataset_id, '_fold_', setdiff(c(1:2), fold), '.csv'))
  
  ### Eta
  eta.hat <- 0.5
  
  ### Pi
  pi.mod <- 
    SuperLearner(Y = df_train$S, 
                 X = as.matrix(select(df_train, A, Y)),
                 family = "binomial",
                 SL.library = c("SL.glm.interaction",
                                "SL.polymars"))
  
  ### Mu
  # HAL-based density estimate of Y|A
  hal_mu <- 
    haldensify(
      A = df_train$Y, 
      W = df_train$A,
      n_bins = 25, 
      grid_type = "equal_range",
      lambda_seq = exp(seq(-1, -10, length = 100)),
      # arguments passed to hal9001::fit_hal()
      max_degree = 3,
      reduce_basis = 1 / sqrt(nrow(df_train))
    )
  
  ### Lambda 1
  lambda1.mod <- 
    SuperLearner(Y = df_train$L1[df_train$S == 1], 
                 X = as.matrix(select(filter(df_train, S == 1), A, Y)),
                 family = "binomial",
                 SL.library = c("SL.glm.interaction",
                                "SL.polymars"))
  
  ### Lambda 2
  # HAL-based density estimate of L_p2 | L_p1, A, Y
  hal_lambda <- 
    haldensify(
      A = df_train$L2[df_train$S == 1], 
      W = as.matrix(select(filter(df_train, S == 1), L1, A, Y)),
      n_bins = 35, 
      grid_type = "equal_range",
      lambda_seq = exp(seq(-1, -10, length = 50)),
      # arguments passed to hal9001::fit_hal()
      max_degree = 3,
      reduce_basis = 1 / sqrt(nrow(df_train))
    )
  
  
  ### Make Predictions to Save Out
  ### On current fold
  delta.y <- 0.02
  new_y <- seq(0, 1, by = delta.y)
  new_dat <- 
    as.data.table(list(
      y = new_y,
      a_zero = rep(0, length(new_y)),
      a_one = rep(1, length(new_y))
    ))
  new_dat[, pred_a_zero := predict(hal_mu, 
                                   trim = F,
                                   new_A = new_dat$y,
                                   new_W = new_dat$a_zero)]
  new_dat[, pred_a_one := predict(hal_mu, 
                                  trim = F,
                                  new_A = new_dat$y,
                                  new_W = new_dat$a_one)]
  
  ### On the opposite fold
  pi.hat.1 <- 
    predict.SuperLearner(pi.mod, 
                         type = "response",
                         newdata = cbind.data.frame(A = df_test$A, Y = df_test$Y))$pred
  lambda1.hat <- 
    predict.SuperLearner(lambda1.mod, 
                         type = "response",
                         newdata = cbind.data.frame(A = df_test$A, Y = df_test$Y))$pred
  
  lambda.obs.mat.1 <- sapply(1:nrow(df_test), function(i) {
    if (df_test$S[i] == 0) { rep(NA, length(new_y)) }
    else { 
      ifelse(df_test$L1[i] == 1, lambda1.hat[i], 1 - lambda1.hat[i]) *
        predict(hal_lambda, trim = F,
                new_A = rep(df_test$L2[i], length(new_y)),
                new_W = cbind(L1 = rep(df_test$L1[i], length(new_y)), 
                              A = rep(1, length(new_y)),
                              Y = new_y))
    }
  })
  
  lambda.obs.mat.0 <- sapply(1:nrow(df_test), function(i) {
    if (df_test$S[i] == 0) { rep(NA, length(new_y)) }
    else { 
      ifelse(df_test$L1[i] == 1, lambda1.hat[i], 1 - lambda1.hat[i]) *
        predict(hal_lambda, trim = F,
                new_A = rep(df_test$L2[i], length(new_y)),
                new_W = cbind(L1 = rep(df_test$L1[i], length(new_y)), 
                              A = rep(0, length(new_y)),
                              Y = new_y))
    }
  })
  
  ### Save Data and Objects
  save(new_dat, delta.y, new_y,
       eta.hat, pi.mod, pi.hat.1,
       lambda1.mod, lambda1.hat,
       lambda.obs.mat.0, lambda.obs.mat.1,
       hal_mu, hal_lambda,
       file = paste0('model_objects/objects_dataset_', dataset_id, '_fold_', setdiff(c(1:2), fold), '.Rdata'))

}



