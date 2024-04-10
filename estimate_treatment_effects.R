treatment_estimate <- function(df, method, interactions = F, multiple_lp = F, ...) {
  .l <- list(...)
  
  ### Either it's complete case or we've already imputed the missing confounder
  df <- drop_na(df)
  
  if(method == 'OLS') {
    if(!interactions) {
      df0 <- df1 <- df
      df0$A <- 0
      df1$A <- 1
      
      ### fit model
      if(!multiple_lp) {
        ols <- lm(ptwc ~ A + GENDER + bbl + raceeth_hispanic + cbl, data = df)
      } else {
        ols <- lm(ptwc ~ A + GENDER + bbl + raceeth_hispanic + cbl + smoker, data = df)
      }
      
      ### Estimate Individual Treatment Effects
      df <-
        df %>% 
        mutate('mu_hat0' = predict(ols, newdata = df0),
               'mu_hat1' = predict(ols, newdata = df1)) %>% 
        filter(!is.na(mu_hat0), !is.na(mu_hat1))
      
    } else {
      ### Specify Pairwise Interactions
      df <- 
        df %>% 
        mutate('A_bbl' = A*bbl,
               'A_cbl' = A*cbl,
               'A_GENDER' = A*GENDER,
               'A_race' = A*raceeth_hispanic,
               'GENDER_bbl' = GENDER*bbl,
               'GENDER_cbl' = GENDER*cbl,
               'GENDER_race' = GENDER*raceeth_hispanic,
               'race_bbl' = raceeth_hispanic*bbl,
               'race_cbl' = raceeth_hispanic*cbl,
               'bbl_cbl' = bbl*cbl)
      
      df1 <- 
        df %>% 
        mutate('A' = 1) %>% 
        mutate('A_bbl' = A*bbl,
               'A_cbl' = A*cbl,
               'A_GENDER' = A*GENDER,
               'A_race' = A*raceeth_hispanic)
      
      df0 <- 
        df %>% 
        mutate('A' = 0) %>% 
        mutate('A_bbl' = A*bbl,
               'A_cbl' = A*cbl,
               'A_GENDER' = A*GENDER,
               'A_race' = A*raceeth_hispanic)
      
      ### Add in Interactions for 2nd Lp if applicable
      if(multiple_lp) {
        df <- 
          df %>% 
          mutate('smoker_A' = smoker * A,
                 'smoker_bbl' = smoker * bbl,
                 'smoker_cbl' = smoker * cbl,
                 'smoker_GENDER' = smoker * GENDER,
                 'smoker_race' = smoker * raceeth_hispanic)
        
        df0 <- 
          df0 %>% 
          mutate('smoker_A' = smoker * A,
                 'smoker_bbl' = smoker * bbl,
                 'smoker_cbl' = smoker * cbl,
                 'smoker_GENDER' = smoker * GENDER,
                 'smoker_race' = smoker * raceeth_hispanic)
        
        df1 <- 
          df1 %>% 
          mutate('smoker_A' = smoker * A,
                 'smoker_bbl' = smoker * bbl,
                 'smoker_cbl' = smoker * cbl,
                 'smoker_GENDER' = smoker * GENDER,
                 'smoker_race' = smoker * raceeth_hispanic)
      }
      
      
      ### Fit Lasso Model
      Y <- as.matrix(df$ptwc)
      X <- as.matrix(select(df, -ptwc, -R))
      glmnet_cv <- cv.glmnet(y = Y, x = X, alpha = 1)
      model <- glmnet(y = Y, x = X, lambda = glmnet_cv$lambda.min, alpha = 1)
      
      ### Estimate Individual Treatment Effects
      df <-
        df %>% 
        mutate('mu_hat0' = predict(model, newx = as.matrix(select(df0, -ptwc, -R))),
               'mu_hat1' = predict(model, newx = as.matrix(select(df1, -ptwc, -R)))) %>% 
        filter(!is.na(mu_hat0), !is.na(mu_hat1))
    }
    
    
    ### Use G-Formula to estimate average counterfactual treatment effects
    estimates <- 
      tibble('te0_hat' = mean(df$mu_hat0),
             'te1_hat' = mean(df$mu_hat1),
             'ate_hat' = mean(df$mu_hat1 - df$mu_hat0))
    
  } else if(method == 'gam') {
    
    if(!interactions) {
      ### fit model
      if(!multiple_lp) {
        outcome_model <- 
          gam(ptwc ~ 
                A + GENDER + raceeth_hispanic + 
                s(bbl, bs = 'cr') + s(cbl, bs = 'cr'),
              data = df)
      } else {
        outcome_model <- 
          gam(ptwc ~ 
                A + GENDER + raceeth_hispanic + smoker + 
                s(bbl, bs = 'cr') + s(cbl, bs = 'cr'),
              data = df)
      }
    } else {
      if(!multiple_lp) {
        outcome_model <- 
          gam(ptwc ~ 
                ### Discrete terms + Interactions
                A + GENDER + raceeth_hispanic + 
                A:GENDER + A:raceeth_hispanic + GENDER:raceeth_hispanic +
                
                ### BBL Splines (by discrete interaction terms)
                s(bbl, bs = 'cs') + 
                s(bbl, bs = 'cs', by = GENDER) + 
                s(bbl, bs = 'cs', by = raceeth_hispanic) + 
                s(bbl, bs = 'cs', by = A) + 
                
                ### CBL Splines (by discrete interaction terms)
                s(cbl, bs = 'cs') + 
                s(cbl, bs = 'cs', by = GENDER) + 
                s(cbl, bs = 'cs', by = raceeth_hispanic) + 
                s(cbl, bs = 'cs', by = A) + 
                
                ### Multivariate Spline
                te(cbl, bbl, bs = 'cs'),
              
              data = df)
      } else {
        outcome_model <- 
          gam(ptwc ~ 
                ### Discrete terms + Interactions
                A + GENDER + raceeth_hispanic + smoker +
                A:GENDER + A:raceeth_hispanic + GENDER:raceeth_hispanic +
                A:smoker + raceeth_hispanic:smoker + GENDER:smoker +
                
                ### BBL Splines (by discrete interaction terms)
                s(bbl, bs = 'cs') + 
                s(bbl, bs = 'cs', by = GENDER) + 
                s(bbl, bs = 'cs', by = raceeth_hispanic) + 
                s(bbl, bs = 'cs', by = A) + 
                s(bbl, bs = 'cs', by = smoker) +
                
                ### CBL Splines (by discrete interaction terms)
                s(cbl, bs = 'cs') + 
                s(cbl, bs = 'cs', by = GENDER) + 
                s(cbl, bs = 'cs', by = raceeth_hispanic) + 
                s(cbl, bs = 'cs', by = A) + 
                s(cbl, bs = 'cs', by = smoker) +
                
                ### Multivariate Spline
                te(cbl, bbl, bs = 'cs'),
              
              data = df)
        
      }
    }
    
    df0 <- df1 <- df
    df0$A <- 0
    df1$A <- 1
    
    df <-
      df %>% 
      mutate('mu_hat0' = predict(outcome_model, newdata = df0),
             'mu_hat1' = predict(outcome_model, newdata = df1)) %>% 
      filter(!is.na(mu_hat0), !is.na(mu_hat1))
    
    ### Use G-Formula to estimate average counterfactual treatment effects
    estimates <- 
      tibble('te0_hat' = mean(df$mu_hat0),
             'te1_hat' = mean(df$mu_hat1),
             'ate_hat' = mean(df$mu_hat1 - df$mu_hat0))
    
    
  } else if(method == 'rf') {
    df0 <- df1 <- df
    df0$A <- 0
    df1$A <- 1
    
    ### Fit Random Forrest Using Defaults
    if(!multiple_lp) {
      rf <- 
        ranger(ptwc ~ A + bbl + GENDER + cbl + raceeth_hispanic, 
               mtry = 3,
               num.trees = 500,
               min.node.size = 10,
               data = df)
    } else {
      rf <- 
        ranger(ptwc ~ A + bbl + GENDER + cbl + raceeth_hispanic + smoker, 
               mtry = 3,
               num.trees = 500,
               min.node.size = 10,
               data = df)
    }
    
    df <-
      df %>% 
      mutate('mu_hat0' = predict(rf, data = df0)$predictions,
             'mu_hat1' = predict(rf, data = df1)$predictions) %>% 
      filter(!is.na(mu_hat0), !is.na(mu_hat1))
    
    ### Use G-Formula to estimate average counterfactual treatment effects
    estimates <- 
      tibble('te0_hat' = mean(df$mu_hat0),
             'te1_hat' = mean(df$mu_hat1),
             'ate_hat' = mean(df$mu_hat1 - df$mu_hat0))
    
  } else if(method == 'IPW') {
    if(is.null(.l$weights) | .l$weights == 'log_reg') {
      ### Estimate probability of treatment weights
      if(!interactions) {
        if(!multiple_lp) {
          weight_model <- 
            glm(A ~ GENDER + bbl + raceeth_hispanic + cbl, 
                data = df, 
                family = 'binomial')
        } else {
          weight_model <- 
            glm(A ~ GENDER + bbl + raceeth_hispanic + cbl + smoker, 
                data = df, 
                family = 'binomial')
        }
        
        probs <- predict(weight_model, newdata = df, type = 'response')
        
      } else {
        
        ### Specify Pairwise Interactions
        df <- 
          df %>% 
          mutate('GENDER_bbl' = GENDER*bbl,
                 'GENDER_cbl' = GENDER*cbl,
                 'GENDER_race' = GENDER*raceeth_hispanic,
                 'race_bbl' = raceeth_hispanic*bbl,
                 'race_cbl' = raceeth_hispanic*cbl,
                 'bbl_cbl' = bbl*cbl)
        
        if(multiple_lp) {
          df <- 
            df %>% 
            mutate('smoker_GENDER' = smoker*GENDER,
                   'smoker_bbl' = smoker*bbl,
                   'smoker_cbl' = smoker * cbl,
                   'smoker_race' = smoker * raceeth_hispanic)
          
          
        }
        
        Y <- as.factor(df$A)
        X <- as.matrix(select(df, -ptwc, -A, -R))
        glmnet_cv <- cv.glmnet(y = Y, x = X, alpha = 1, family = 'binomial')
        weight_model <- 
          glmnet(y = Y, 
                 x = X, 
                 lambda = glmnet_cv$lambda.min, 
                 alpha = 1,
                 family = 'binomial')
        
        probs <- predict(weight_model, newx = X, type = 'response')
        
      }
      
    } else if(.l$weights == 'gam') {
      if(!interactions) {
        if(!multiple_lp) {
          weight_model <- gam(A ~ GENDER + raceeth_hispanic + 
                                s(bbl, bs = 'cs') + s(cbl, bs = 'cs'),
                              data = df,
                              family = binomial(link = 'logit'))
        } else {
          weight_model <- gam(A ~ GENDER + raceeth_hispanic + smoker + 
                                s(bbl, bs = 'cs') + s(cbl, bs = 'cs'),
                              data = df,
                              family = binomial(link = 'logit'))
        }
      } else {
        if(!multiple_lp) {
          weight_model <- 
            gam(A ~ 
                  ### Discrete Interactions
                  GENDER + raceeth_hispanic + GENDER:raceeth_hispanic +
                  
                  ### BBL Splines (by discrete interaction terms)
                  s(bbl, bs = 'cs') + 
                  s(bbl, bs = 'cs', by = GENDER) + 
                  s(bbl, bs = 'cs', by = raceeth_hispanic) + 
                  
                  ### CBL Splins (by discrete interaction terms)
                  s(cbl, bs = 'cs') + 
                  s(cbl, bs = 'cs', by = GENDER) + 
                  s(cbl, bs = 'cs', by = raceeth_hispanic) + 
                  
                  ### Multivariate Spline
                  te(cbl, bbl, bs = 'cs'),
                
                data = df,
                family = binomial(link = 'logit'))
        } else {
          weight_model <- 
            gam(A ~ 
                  ### Discrete Interactions
                  GENDER + raceeth_hispanic + smoker + 
                  smoker:GENDER + smoker:raceeth_hispanic + GENDER:raceeth_hispanic +
                  
                  ### BBL Splines (by discrete interaction terms)
                  s(bbl, bs = 'cs') + 
                  s(bbl, bs = 'cs', by = GENDER) + 
                  s(bbl, bs = 'cs', by = raceeth_hispanic) + 
                  s(bbl, bs = 'cs', by = smoker) + 
                  
                  ### CBL Splins (by discrete interaction terms)
                  s(cbl, bs = 'cs') + 
                  s(cbl, bs = 'cs', by = GENDER) + 
                  s(cbl, bs = 'cs', by = raceeth_hispanic) + 
                  s(cbl, bs = 'cs', by = smoker) + 
                  
                  ### Multivariate Spline
                  te(cbl, bbl, bs = 'cs'),
                
                data = df,
                family = binomial(link = 'logit'))
          
        }
      }
      
      probs <- predict(weight_model, newdata = df, type = 'response')
      
    } else if(.l$weights == 'rf') {
      if(!multiple_lp) {
        weight_model <- 
          ranger(as.factor(A) ~ bbl + GENDER + cbl + raceeth_hispanic, 
                 # min.node.size = 10,
                 # max.depth = 5,
                 data = df,
                 probability = T,
                 classification = T)
      } else {
        weight_model <- 
          ranger(as.factor(A) ~ bbl + GENDER + cbl + raceeth_hispanic + smoker,
                 # min.node.size = 10,
                 # max.depth = 5,
                 data = df,
                 probability = T,
                 classification = T)
      }
      
      probs <- predict(weight_model, data = df)$predictions[,'1']
    }
    
    ### Threshold Probs in case some are exactly 0 or 1 
    probs <- 
      case_when(probs < 0.001 ~ 0.001,
                probs > 0.999 ~ 0.999,
                T ~ probs)
    
    df <- 
      df %>% 
      mutate('prob_A1' = probs) %>% 
      mutate('weight_1' = 1/prob_A1,
             'weight_0' = 1/(1-prob_A1),
             'weight' = weight_1 * A + weight_0 * (1-A))
    
    wls <- lm(ptwc ~ A, data = df, weights = weight)
    
    estimates <- 
      tibble('Hajek_0' = wls$coefficients[1],
             'Hajek_1' = wls$coefficients[1] + wls$coefficients[2],
             'Horvitz-Thompson_0' = mean(df$ptwc * (1-df$A) * df$weight_0),
             'Horvitz-Thompson_1' = mean(df$ptwc * df$A * df$weight_1))
    
  } 
  
  return(estimates)
}


### Function to call estimate_treatment_effects over batches (helps lower 
### amount of memory being used at any one time
batch_treatment_estimate <- function(datasets, batches, method, multiple_lp = F, interactions = F, ...) {
  .l <- list(...)
  df <- NULL
  
  if(!is.null(.l$seed)) {
    ### Need seeds to parallelization
    if(length(batches) > 1) {
      set.seed(.l$seed)
      seeds <- sample(1:1e7, length(batches))
    } else {
      seeds <- .l$seed
    }
    
    if(is.null(.l$weights)) {
      for(i in 1:length(batches)) {
        tmp <- 
          future_map_dfr(datasets[ batches[[i]] ], ~treatment_estimate(df = .x, 
                                                                       method = method, 
                                                                       multiple_lp = multiple_lp,
                                                                       interactions = interactions),
                         .options = furrr_options(seed = seeds[i])) 
        df <- bind_rows(df, tmp)
        garbage <- gc()
      }
    } else {
      for(i in 1:length(batches)) {
        tmp <- 
          future_map_dfr(datasets[ batches[[i]] ], ~treatment_estimate(df = .x, 
                                                                       method = method, 
                                                                       multiple_lp = multiple_lp,
                                                                       interactions = interactions, 
                                                                       weights = .l$weights),
                         .options = furrr_options(seed = seeds[i])) 
        df <- bind_rows(df, tmp)
        garbage <- gc()
      }
    }
  } else {
    if(is.null(.l$weights)) {
      for(i in 1:length(batches)) {
        tmp <- 
          future_map_dfr(datasets[ batches[[i]] ], ~treatment_estimate(df = .x, 
                                                              method = method, 
                                                              multiple_lp = multiple_lp,
                                                              interactions = interactions)) 
        df <- bind_rows(df, tmp)
        garbage <- gc()
      }
    } else {
      for(i in 1:length(batches)) {
        tmp <- 
          future_map_dfr(datasets[ batches[[i]] ], ~treatment_estimate(df = .x, 
                                                              method = method, 
                                                              multiple_lp = multiple_lp,
                                                              interactions = interactions, 
                                                              weights = .l$weights)) 
        df <- bind_rows(df, tmp)
        garbage <- gc()
      }
    }
  }
  
  return(df)
}