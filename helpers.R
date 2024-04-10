### Inverse Logit Function
expit <- function(x) {
  return(exp(x)/(1 + exp(x)))
}

### Logit Function
logit <- function(x) {
  return(log(x/(1-x)))
}

### Function to shift predictions by pre-specified amplification factors
shift_preds <- function(preds, shift, df) {
  for(i in 1:length(shift)) {
    coeff <- names(shift)[i]
    if(coeff == '(Intercept)') { ### Intercept Shift
      preds <- preds + shift[[coeff]]
    } else if(grepl(':', coeff)) { ### Interaction Shift
      coeffs <- unlist(strsplit(coeff, ':'))
      preds <- preds + df[[ coeffs[1] ]] * df[[ coeffs[2] ]] * shift[[coeff]]
    }  else if(grepl('\\^2\\)$', coeff)) { ### quadratic term
      variable <- gsub('\\^2\\)$', '', gsub('^I\\(', '', coeff))
      preds <- preds + df[[ variable ]]^2 * shift[[coeff]]
    } else { ### Single coeff shift
      preds <- preds + df[[coeff]] * shift[[coeff]]
    }
  }
  
  return(preds)
}

### Function to cap predictions
pred_cap <- function(preds, upper = 0.999, lower = 0.001) {
  preds <- 
    case_when(preds > upper ~ upper,
              preds < lower ~ lower,
              T ~ preds)
  return(preds)
}

### Function to impute the missing confounder with the true model structure
### Expand for multiple imputation [TODO]
impute <- function(df, model_types, method = 'true_dgp', multiple_lp = F) {
  if(method == 'true_dgp') {
    if(model_types$lambda == 'levis_paper') {
      lambda_model <- 
        glm(cbl ~ A + GENDER + bbl + raceeth_hispanic + ptwc + A:ptwc, 
            data = df[df$R == 1,],
            family = Gamma(link = "log"))
    } else if(model_types$lambda == 'interact_nonlinear') {
      lambda_model <- 
        glm(cbl ~ A + GENDER + bbl + raceeth_hispanic + ptwc + A:ptwc 
            + I(bbl^2) + bbl:ptwc, 
            data = df[df$R == 1,],
            family = Gamma(link = "log"))
    } else if(model_types$lambda == 'alternative' | 
              model_types$lambda == 'levis_flexible' | 
              model_types$lambda == 'levis_flexible_linear') {
      lambda_model <- 
        glm(cbl ~ GENDER + bbl + raceeth_hispanic,
            data = df[df$R == 1,],
            family = Gamma(link = "log")) 
    }
    
    lambda_shape_hat <- 1/summary(lambda_model)$dispersion
    
    ### Impute Missing Values
    df$cbl[df$R == 0] <- 
      rgamma(n = sum(df$R== 0),
             shape = lambda_shape_hat,
             rate = lambda_shape_hat/exp(predict(lambda_model,newdata = df[df$R == 0,])))
    
    if(multiple_lp) {
      if(model_types$rho == 'durable') {
        rho_model <- 
          glm(smoker ~ GENDER + bbl + raceeth_hispanic + A + ptwc + cbl + bbl:ptwc,
              data = df,
              family = binomial(link = 'logit'))
      } else if(model_types$rho == 'durable_interact') {
        rho_model <- 
          glm(smoker ~ GENDER + bbl + raceeth_hispanic + A + ptwc + cbl + bbl:ptwc + A:ptwc,
              data = df,
              family = binomial(link = 'logit'))
      } else if(model_types$rho == 'aylp1_interact') {
        rho_model <- 
          glm(smoker ~ A + ptwc + cbl + A:ptwc + A:cbl + cbl:ptwc,
              data = df,
              family = binomial(link = 'logit'))
      } else if(model_types$rho == 'cbl_only') {
        rho_model <- 
          glm(smoker ~ cbl,
              data = df,
              family = binomial(link = 'logit'))
        
      } else if(model_types$rho == 'cbl_ptwc') {
        rho_model <- 
          glm(smoker ~ cbl + ptwc,
              data = df,
              family = binomial(link = 'logit'))
      } else if(model_types$rho == 'alternative' | 
                model_types$rho == 'levis_flexible') {
        rho_model <- 
          glm(smoker ~ cbl + GENDER + bbl + raceeth_hispanic,
              data = df,
              family = binomial(link = 'logit'))
        
      }
      
      ### Impute Missing Values 
      df$smoker[df$R == 0] <- 
        rbinom(n = sum(df$R == 0),
               size = 1,
               p = expit(predict(rho_model, newdata = df[df$R == 0,])))
    }
  } else if(method == 'lm') {
    
    lambda_model <- lm(cbl ~ GENDER + bbl + raceeth_hispanic + A + ptwc, data = df)
    sigma <- sqrt(sum(resid( lambda_model)^2 /  lambda_model$df.residual))
    
    
    ### Impute Missing Values
    df$cbl[df$R == 0] <- rnorm(n = sum(df$R== 0),
                               mean = predict(lambda_model, newdata = df[df$R == 0,]),
                               sd = sigma)
    
    ### 2nd Lp, if applicable
    if(multiple_lp) {
      rho_model <- 
        glm(smoker ~ GENDER + bbl + raceeth_hispanic + A + ptwc + cbl,
            data = df,
            family = binomial(link = 'logit'))
      
      df$smoker[df$R == 0] <- 
        rbinom(n = sum(df$R == 0),
               size = 1,
               p = expit(predict(rho_model, newdata = df[df$R == 0,])))
    }
    
    
  } else if(method == 'lm_interact') {
    ### Specify Pairwise Interactions
    df_tmp <- 
      df %>% 
      mutate('A_bbl' = A*bbl,
             'A_ptwc' = A*ptwc,
             'A_GENDER' = A*GENDER,
             'A_race' = A*raceeth_hispanic,
             'GENDER_bbl' = GENDER*bbl,
             'GENDER_ptwc' = GENDER*ptwc,
             'GENDER_race' = GENDER*raceeth_hispanic,
             'race_bbl' = raceeth_hispanic*bbl,
             'race_ptwc' = raceeth_hispanic*ptwc,
             'bbl_ptwc' = bbl*ptwc)
    
    
    ### Fit LASSO Model
    Y <- as.matrix(df_tmp$cbl[df$R == 1])
    X <- as.matrix(df_tmp %>% filter(R == 1) %>% select(-cbl, -R, -any_of("smoker")))
    glmnet_cv <- cv.glmnet(y = Y, x = X, alpha = 1)
    lambda_model <- 
      glmnet(y = Y, 
             x = X, 
             lambda = glmnet_cv$lambda.min, 
             alpha = 1)
    
    X0 <- as.matrix(df_tmp %>% filter(R == 0) %>% select(-cbl, -R, -any_of("smoker")))
    lambda_preds <- as.vector(predict(lambda_model, newx = X))
    lambda_preds0 <- as.vector(predict(lambda_model, newx = X0))
    
    sigma <- 
      sqrt(sum((lambda_preds - Y)^2)/(sum(df$R == 1) - (lambda_model$df + 1)))
    
    
    ### Impute Missing Values
    df$cbl[df$R == 0] <- rnorm(n = sum(df$R== 0),
                               mean = lambda_preds0,
                               sd = sigma)
    
    if(multiple_lp) {
      ### Add in all interactions w/ CBL
      df_tmp <- 
        df %>% 
        mutate('A_bbl' = A*bbl,
               'A_ptwc' = A*ptwc,
               'A_GENDER' = A*GENDER,
               'A_race' = A*raceeth_hispanic,
               'GENDER_bbl' = GENDER*bbl,
               'GENDER_ptwc' = GENDER*ptwc,
               'GENDER_race' = GENDER*raceeth_hispanic,
               'race_bbl' = raceeth_hispanic*bbl,
               'race_ptwc' = raceeth_hispanic*ptwc,
               'bbl_ptwc' = bbl*ptwc) %>% 
        mutate('A_cbl' = A*cbl,
               'ptwc_cbl' = ptwc*cbl,
               'bbl_cbl' = bbl*cbl,
               'race_cbl' = raceeth_hispanic*cbl,
               'GENDER_cbl' = GENDER*cbl)
      
      ### Fit LASSO Model (Logistic Regression)
      Y <- as.matrix(df_tmp$smoker[df$R == 1])
      X <- as.matrix(df_tmp %>% filter(R == 1) %>% select(-R, -smoker))
      glmnet_cv <- cv.glmnet(y = Y, x = X, alpha = 1, family = 'binomial')
      rho_model <- 
        glmnet(y = Y, 
               x = X, 
               lambda = glmnet_cv$lambda.min, 
               alpha = 1, 
               family = 'binomial')
      
      X0 <- as.matrix(df_tmp %>% filter(R == 0) %>% select(-R, -smoker))
      rho_preds0 <- as.vector(predict(rho_model, newx = X0, type = 'response'))
      
      
      
      ### Impute Missing Values
      df$smoker[df$R == 0] <- rbinom(n = sum(df$R== 0),
                                     size = 1,
                                     p = rho_preds0) 
    }
  }
  
  return(df)
}
