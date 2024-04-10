library(tidyverse)
library(here)
source(here('helpers.R'))

################################################################################
############################### Data Prep ######################################
################################################################################
### Load in data and create variables
### treatment/missingness indicators
### Gender (1) for (M) and (0) for (F)
### Center baseline
###
### If simulation parameters are loaded in memory, uses those specified covariate 
### shifts. Overwise, defaults to what Alex used in his paper
###
### In terms of the truth model type levis_flexible == alternative
### The difference is for what to do when estimating the Levis IF/IWOR models
load_models <- function(params = NULL) {
  df_durable <<- 
    read_csv(here('durable_missing_cbl_only.csv'), col_types = cols()) %>% 
    mutate('A' = ifelse(final_decision_merge == 'RYGB', 0, 1),
           'R' = 1 - is.na(cbl),
           'GENDER' = as.numeric(GENDER == 'M'),
           'bbl' = bbl-30, 
           'cbl' = cbl + 2)
  
  ### If we want to deal w/ multiple confounders we'll add in a smoking status
  ### variable that is not in the original data set but rather informed by 
  ### DURABLE
  if(params$multiple_lp) {
    ### [TO-DO] Base on Real Data
    ### Coefficients for Smoking Variable we will add
    if(params$models$rho == 'durable') {
      rho_star <<- list(
        'prevalence' = logit(0.2),
        'beta_GENDER' = 0.07,
        'beta_bbl' = 0.04,
        'beta_raceeth_hispanic' = 0.05,
        'beta_ptwc' = 0.25,
        'beta_cbl' = 0.05,
        'beta_bbl_ptwc' = 0.05,
        'beta_A'  = 0.03,
        'beta_A_ptwc' = 0,
        'beta_A_cbl' = 0,
        'beta_cbl_ptwc' = 0
      )
    } else if(params$models$rho == 'cbl_only') {
      rho_star <<- list(
        'prevalence' = logit(0.2),
        'beta_GENDER' = 0,
        'beta_bbl' = 0,
        'beta_raceeth_hispanic' = 0,
        'beta_ptwc' = 0,
        'beta_cbl' = 0.075,
        'beta_bbl_ptwc' = 0,
        'beta_A'  = 0,
        'beta_A_ptwc' = 0,
        'beta_A_cbl' = 0,
        'beta_cbl_ptwc' = 0
      )
    } else if(params$models$rho == 'cbl_ptwc') {
      rho_star <<- list(
        'prevalence' = logit(0.2),
        'beta_GENDER' = 0,
        'beta_bbl' = 0,
        'beta_raceeth_hispanic' = 0,
        'beta_ptwc' = -1,
        'beta_cbl' = 0.075,
        'beta_bbl_ptwc' = 0,
        'beta_A'  = 0,
        'beta_A_ptwc' = 0,
        'beta_A_cbl' = 0,
        'beta_cbl_ptwc' = 0
      )
    } else if(params$models$rho == 'durable_interact') {
      rho_star <<- list(
        'prevalence' = logit(0.2),
        'beta_GENDER' = 0.07,
        'beta_bbl' = 0.04,
        'beta_raceeth_hispanic' = 0.05,
        'beta_ptwc' = -4,
        'beta_cbl' = 0.05,
        'beta_bbl_ptwc' = 0.05,
        'beta_A'  = 0.03,
        'beta_A_ptwc' = 2,
        'beta_A_cbl' = 0,
        'beta_cbl_ptwc' = 0
      )
    } else if(params$models$rho == 'aylp1_interact') {
      rho_star <<- list(
        'prevalence' = logit(0.2),
        'beta_GENDER' = 0,
        'beta_bbl' = 0,
        'beta_raceeth_hispanic' = 0,
        'beta_ptwc' = -4,
        'beta_cbl' = 0.075,
        'beta_bbl_ptwc' = 0,
        'beta_A'  = 0,
        'beta_A_ptwc' = 2,
        'beta_A_cbl' = -0.01,
        'beta_cbl_ptwc' = 0.01) 
    } else if(params$models$rho == 'alternative' | 
              params$models$rho == 'levis_flexible') {
      ### alternative factorization Lp2 | Lc, Lp1
      rho_star <<- list(
        'prevalence' = logit(0.2),
        'beta_GENDER' = 0.03,
        'beta_bbl' = 0.05,
        'beta_raceeth_hispanic' = 0.025,
        'beta_cbl' = 0.075,
        
        'beta_ptwc' = 0,
        'beta_bbl_ptwc' = 0,
        'beta_A'  = 0,
        'beta_A_ptwc' = 0,
        'beta_A_cbl' = 0,
        'beta_cbl_ptwc' = 0)  
    }
  }
  
  ### Filter to only complete cases, we will generate the missingness
  df_cc  <<- 
    df_durable %>% 
    filter(R == 1) %>% 
    select(GENDER, bbl, raceeth_hispanic, A, cbl, R, ptwc)
  
  ################################################################################
  ########################### Model Specifications ###############################
  ################################################################################
  ### Treatment Model Specification
  if(params$models$eta == 'levis_paper') {
    eta_star <<- 
      glm(A ~ GENDER + bbl + raceeth_hispanic + I(bbl^2) + GENDER:raceeth_hispanic,
          data = df_cc, 
          family = 'binomial')
  } else if(params$models$eta == 'alternative' | 
            params$models$eta == 'levis_flexible') {
    ### alternative factorization A| Lc, Lp
    ### express lp2 effects w/ covariate shifts since smoker doesn't exist in the data
    eta_star <<- 
      glm(A ~ GENDER + bbl + raceeth_hispanic + I(bbl^2) + GENDER:raceeth_hispanic + cbl,
          data = df_cc, 
          family = 'binomial')
  }
  
  # Shifts to amplify certain effects
  if(is.null(params$shifts$eta)) {
    eta_shift <<- 
      list('GENDER' = 0.15, 
           'raceeth_hispanic' = 0.2,
           'GENDER:raceeth_hispanic' = 0.1)
  } else {
    eta_shift <<- params$shifts$eta 
  }
  
  
  ### Partial Outcome Model Specification
  if(params$models$mu == 'levis_paper') {
    mu_star <<- 
      lm(ptwc ~ A + GENDER + bbl + raceeth_hispanic +
           GENDER:raceeth_hispanic +
           A:GENDER + A:bbl + A:raceeth_hispanic,
         data = df_cc)
  } else if(params$models$mu == 'alternative' | 
            params$models$mu == 'levis_flexible') {
    ### alternative factorization Y | Lc, Lp, A
    ### express lp2 effects w/ covariate shifts since smoker doesn't exist in the data
    ### (possibly lp2:A effect as well)
    mu_star <<- 
      lm(ptwc ~ A + GENDER + bbl + raceeth_hispanic + cbl +
           GENDER:raceeth_hispanic +
           A:GENDER + A:bbl + A:raceeth_hispanic + A:cbl,
         data = df_cc)
  }
  
  # Estimated variance
  sigma_star <<- sqrt(sum(resid(mu_star)^2 / mu_star$df.residual))
  
  # Shifts to amplify certain effects
  if(is.null(params$shifts$mu)) {
    mu_shift <<- 
      list('GENDER:raceeth_hispanic' = -0.3,
           'A:GENDER' = 0.325,
           'A:raceeth_hispanic' = 0.09,
           'raceeth_hispanic' = 0.02,
           'A:bbl' = -0.003)
  } else {
    mu_shift <<- params$shifts$mu
  }
  
  ### Missingness Model Specification (uses full data, not just complete case data)
  if(params$models$pi == 'levis_paper') {
    pi_star <<- 
      glm(R ~ GENDER + bbl + raceeth_hispanic + A + ptwc +
            A:ptwc + GENDER:ptwc + bbl:raceeth_hispanic + A:bbl,
          data = df_durable,
          family="binomial")
  } else if(params$models$pi == 'alternative' | 
            params$models$pi == 'levis_flexible') {
    ### alternative factorization R | Lc, A, Y
    pi_star <<- 
      glm(R ~ GENDER + bbl + raceeth_hispanic + A + ptwc + 
            A:ptwc + GENDER:ptwc + bbl:raceeth_hispanic + A:bbl,
          data = df_durable,
          family="binomial")
  }
  
  
  # Shifts to amplify certain effects
  if(is.null(params$shifts$pi)) {
    pi_shift <<- 
      list('(Intercept)' = 1,
           'raceeth_hispanic' = -3,
           'ptwc' = 2.7,
           'A' = 2.4,
           'A:ptwc' = 2.29,
           'GENDER:ptwc' = 1.99,
           'bbl:raceeth_hispanic' = -0.14,
           'A:bbl' = 0.18)
  } else {
    pi_shift <<- params$shifts$pi 
  }
  
  ### "Imputation" Model Specification
  if(params$models$lambda == 'levis_paper') {
    lambda_star <<- 
      glm(cbl ~ GENDER + bbl + raceeth_hispanic + A + ptwc,
          data = df_cc,
          family = Gamma(link = "log"))
  } else if(params$models$lambda == 'interact_nonlinear') {
    lambda_star <<- 
      glm(cbl ~ GENDER + bbl + raceeth_hispanic + A + ptwc + 
            I(bbl^2) + bbl:ptwc ,
          data = df_cc,
          family = Gamma(link = "log"))
  } else if(params$models$lambda == 'alternative' | 
            params$models$lambda == 'levis_flexible' | 
            params$model$lambda == 'levis_flexible_linear') {
    ### alternative factorization Lp1 | Lc
    lambda_star <<- 
      glm(cbl ~ GENDER + bbl + raceeth_hispanic,
          data = df_cc,
          family = Gamma(link = "log"))
  }
  
  # Shifts to amplify certain effects
  if(is.null(params$shifts$lambda)) {
    lambda_shift <<- 
      list('A' = -0.6,
           'ptwc' = -1.3,
           'A:ptwc' = -0.4)
  } else {
    lambda_shift <<- params$shifts$lambda 
  }
  
  if(params$models$lambda_shape == 'levis_paper' | 
     params$models$lambda_shape == 'alternative' |
     params$models$lambda_shape == 'levis_flexible') {
    # Shape parameter (alpha) for gamma distribution
    lambda_shape <<- 1/summary(lambda_star)$dispersion       
  } else if(is.numeric(params$models$lambda_shape)) {
    lambda_shape <<- params$models$lambda_shape
  }
  
  if(params$multiple_lp) {
    if(is.null(params$shifts$rho)) {
      rho_shift <<- list('bbl' = 0) ### Placeholder of no shifts
    } else {
      rho_shift <<- params$shifts$rho
    }
  }
}

################################################################################
########################### Data Generation ####################################
################################################################################
### Function to generate data from the models we've fit
###
### Note in the multiple confounder setting R == 1 is the same as S == 1 (R^q == 1^q)
### as we make the simplifying assumption that subjects are either missing or have 
### missing data for each of the confounders
generate_data <- function(n_patients, data_generating_process, multiple_lp, truth = F) {
  if(data_generating_process == 'Levis') {
    ### Start from only the 3 complete covariates (sampled from empirical distribution w/ replacement)
    df_sim <- 
      df_cc %>% 
      dplyr::slice(sample(1:nrow(.), size = n_patients, replace = T)) %>% 
      select(GENDER, bbl, raceeth_hispanic)
    
    ### Generate Treatment Assignments
    df_sim <- 
      df_sim %>% 
      mutate('eta_hat' = shift_preds(preds = predict(eta_star, newdata = df_sim),
                                     shift = eta_shift,
                                     df = df_sim)) %>% 
      mutate('A' = rbinom(n = n_patients, size = 1, p = expit(eta_hat)))
    
    ### Generate Outcome 
    df_sim <- 
      df_sim %>% 
      mutate('mu_hat' = shift_preds(preds = predict(mu_star, newdata = df_sim),
                                    shift = mu_shift,
                                    df = df_sim)) %>% 
      mutate('ptwc' = rnorm(n = n_patients, mean = mu_hat, sd = sigma_star))
    
    
    ### Generate Missingness Indicator
    df_sim <- 
      df_sim %>% 
      mutate('pi_hat' = shift_preds(preds = predict(pi_star, newdata = df_sim),
                                    shift = pi_shift,
                                    df = df_sim)) %>% 
      mutate('R' = rbinom(n = n_patients, size = 1, p = expit(pi_hat)))
    
    ### Generate Lp for subject for whom not missing
    df_sim <- 
      df_sim %>% 
      mutate('lambda_hat' = shift_preds(preds = predict(lambda_star, newdata = df_sim),
                                        shift = lambda_shift,
                                        df = df_sim)) %>% 
      mutate('lambda_hat' = ifelse(R == 1, lambda_hat, NA)) %>%
      mutate('cbl' = NA) %>% 
      mutate('cbl' = replace(cbl, R == 1, rgamma(n = sum(R == 1), 
                                                 shape = lambda_shape, 
                                                 rate = lambda_shape/exp(lambda_hat[R == 1]))))
    
    if(multiple_lp) {
      df_sim <- 
        df_sim %>% 
        mutate('rho_hat' = predict_rho_star(rho_star, newdata = df_sim)) %>% 
        mutate('rho_hat' = shift_preds(preds = rho_hat,
                                       shift = rho_shift,
                                       df = df_sim)) %>% 
        mutate('smoker' = NA) %>% 
        mutate('smoker' = replace(smoker, R == 1, rbinom(n = sum(R == 1), 
                                                         size = 1,
                                                         p = expit(rho_hat[R == 1]))))
      
    }
  } else if(data_generating_process == 'Alternative') {
    ### Start from only the 3 complete covariates (sampled from empirical distribution w/ replacement)
    df_sim <- 
      df_cc %>% 
      dplyr::slice(sample(1:nrow(.), size = n_patients, replace = T)) %>% 
      select(GENDER, bbl, raceeth_hispanic)
    
    ### Generate Lp for subject for whom not missing
    df_sim <- 
      df_sim %>% 
      mutate('lambda_hat' = shift_preds(preds = predict(lambda_star, newdata = df_sim),
                                        shift = lambda_shift,
                                        df = df_sim)) %>% 
      mutate('cbl' = rgamma(n = nrow(df_sim), 
                            shape = lambda_shape, 
                            rate = lambda_shape/exp(lambda_hat)))
    
    if(multiple_lp) {
      df_sim <- 
        df_sim %>% 
        mutate('rho_hat' = predict_rho_star(rho_star, newdata = df_sim)) %>% 
        mutate('rho_hat' = shift_preds(preds = rho_hat,
                                       shift = rho_shift,
                                       df = df_sim)) %>% 
        mutate('smoker' = rbinom(n = nrow(df_sim), size = 1, p = expit(rho_hat)))
    }
    
    ### Generate Treatment Assignments
    df_sim <- 
      df_sim %>% 
      mutate('eta_hat' = shift_preds(preds = predict(eta_star, newdata = df_sim),
                                     shift = eta_shift,
                                     df = df_sim)) %>% 
      mutate('A' = rbinom(n = n_patients, size = 1, p = expit(eta_hat)))
    
    ### Generate Outcome 
    df_sim <- 
      df_sim %>% 
      mutate('mu_hat' = shift_preds(preds = predict(mu_star, newdata = df_sim),
                                    shift = mu_shift,
                                    df = df_sim)) %>% 
      mutate('ptwc' = rnorm(n = n_patients, mean = mu_hat, sd = sigma_star))
    
    ### Generate Missingness Indicator
    df_sim <- 
      df_sim %>% 
      mutate('pi_hat' = shift_preds(preds = predict(pi_star, newdata = df_sim),
                                    shift = pi_shift,
                                    df = df_sim)) %>% 
      mutate('R' = rbinom(n = n_patients, size = 1, p = expit(pi_hat))) 
    
    ### Null out Lp for Missing subjects
    if(!truth) {
      df_sim <- 
        df_sim %>% 
        mutate_at(vars(any_of(c('cbl', 'smoker'))), ~replace(.x, R == 0, NA))
    }
    
    ### Remove our intermediate prediction steps
    df_sim <- 
      df_sim %>% 
      select(-contains('hat'))
    
  }
  
  ### Remove our intermediate prediction steps
  df_sim <- 
    df_sim %>% 
    select(-contains('hat'))
  
  return(df_sim)
}
