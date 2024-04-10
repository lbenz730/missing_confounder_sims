library(haldensify)
library(SuperLearner)
library(statmod)
library(tidyverse)
library(data.table)
library(EnvStats)

### Which dataset to work on
args <- commandArgs(trailingOnly = T)
dataset_id <- as.numeric(args[1])
n_folds <- 2

### Fit Models Seperately for Each Fold
for(fold in 1:n_folds) {
  df_train <- read_csv(paste0('input_data/dataset_', dataset_id, '_fold_', fold, '.csv'))
  df_test <- read_csv(paste0('input_data/dataset_', dataset_id, '_fold_', setdiff(c(1:2), fold), '.csv'))
  n <- nrow(df_train)

  delta.y <- 0.02
  y.grid <- seq(0, 1, by = delta.y)
  
  ### fit nuisance models
  
  ## A
  eta.hat <- 0.5
  
  ## Y | A
  ## use MLE for two shape parameters, stratified by A
  mu.pars.0 <- ebeta(df_train$Y[df_train$A == 0], method = "mle")$parameters
  mu.pars.1 <- ebeta(df_train$Y[df_train$A == 1], method = "mle")$parameters
  
  ## S | A, Y
  pi.mod <- 
    glm(S ~ A * Y, family = "binomial",
        data = df_train)
  
  ## L1 | A, Y, S = 1
  lambda1.mod <- 
    glm(L1 ~ A * Y,
        family = "binomial",
        data = filter(df_train, S == 1))
  
  ## L2 | A, Y, S = 1, L1
  lambda2.mod <- 
    lm(L2 ~ A + Y + Y:L1,
       data = filter(df_train, S == 1))
  
  
  ### estimate parameters
  pi.hat.1 <-
    predict.glm(pi.mod, type = "response",
                newdata = df_test)
  
  lambda1.pmf <- function(a, y, l1) {
    p <- predict.glm(lambda1.mod, type = "response",
                     newdata = cbind.data.frame(A = a, Y = y))
    ifelse(l1 == 1, p, 1 - p)
  }
  lambda1.pmf <- Vectorize(lambda1.pmf)
  lambda1.hat <- predict.glm(lambda1.mod, type = "response",
                             newdata = cbind.data.frame(A = df_test$A, Y = df_test$Y))
  
  mu.pdf <- function(a, y) {
    ifelse(a == 0, dbeta(y, mu.pars.0[1], mu.pars.0[2]),
           dbeta(y, mu.pars.1[1], mu.pars.1[2]))
  }
  mu.pdf <- Vectorize(mu.pdf)
  
  mu.pred.0 <- mu.pdf(a = 0, y = y.grid)
  mu.pred.1 <- mu.pdf(a = 1, y = y.grid)
  
  lambda2.pdf <- function(a, y, l1, l2) {
    dnorm(l2, predict.lm(lambda2.mod, newdata = cbind.data.frame(A = a, Y = y, L1 = l1)), sigma(lambda2.mod))
  }
  lambda2.pdf <- Vectorize(lambda2.pdf)
  
  lambda.obs.mat.0 <- sapply(1:n, function(i) {
    if (df_test$S[i] == 0) { rep(NA, length(y.grid)) }
    else { 
      ifelse(df_test$L1[i] == 1, lambda1.hat[i], 1 - lambda1.hat[i]) *
        lambda2.pdf(a = 0, y = y.grid, l1 = df_test$L1[i], l2 = df_test$L2[i])
    }
  })
  lambda.obs.mat.1 <- sapply(1:n, function(i) {
    if (df_test$S[i] == 0) { rep(NA, length(y.grid)) }
    else { 
      ifelse(df_test$L1[i] == 1, lambda1.hat[i], 1 - lambda1.hat[i]) *
        lambda2.pdf(a = 1, y = y.grid, l1 = df_test$L1[i], l2 = df_test$L2[i])
    }
  })
  
  ### Save Data and Objects
  save(n,delta.y, y.grid,
       eta.hat, pi.mod, pi.hat.1,
       lambda1.mod, lambda1.hat, lambda1.pmf,
       lambda.obs.mat.0, lambda.obs.mat.1,
       mu.pars.0, mu.pars.1, mu.pdf, mu.pred.0, mu.pred.1,
       lambda2.mod, lambda2.pdf,
       file = paste0('model_objects/objects_dataset_', dataset_id, '_fold_', setdiff(c(1:2), fold), '_parametric.Rdata'))
  
}



