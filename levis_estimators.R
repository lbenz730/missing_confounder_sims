library(statmod)
library(tidyverse)
library(here)
library(ranger)
library(mgcv)
source(here('helpers.R'))
source(here('levis_helpers.R'))
source(here('generate_data.R'))

### Add code for model misspecification [to-do]
compute_levis_estimators <- function(df, model_types, multiple_lp = F) {
  # Gauss-Hermite prep (for integral over gaussian outcome distribution)
  gh <- gauss.quad(n = 3, kind = "hermite")
  gh_nodes <- gh$nodes
  gh_weights <- gh$weights
  
  ### Estimate models
  # Outcome Model (mu)
  if(model_types$mu == 'levis_paper' | model_types$mu == 'alternative') {
    mu_model <- 
      lm(ptwc ~ A + GENDER + bbl + raceeth_hispanic +  A:GENDER + A:bbl + 
           A:raceeth_hispanic + GENDER:raceeth_hispanic,
         data = df)
    
  } else if(model_types$mu == 'levis_flexible') {
    mu_model <- 
      gam(ptwc ~ 
            ### Discrete terms + Interactions
            A + GENDER + raceeth_hispanic + 
            A:GENDER + A:raceeth_hispanic + GENDER:raceeth_hispanic +
            
            ### BBL Splines (by discrete interaction terms)
            s(bbl, bs = 'cs') + 
            s(bbl, bs = 'cs', by = GENDER) + 
            s(bbl, bs = 'cs', by = raceeth_hispanic) + 
            s(bbl, bs = 'cs', by = A),
          
          family = 'gaussian',
          data = df)
  }
  sigma_hat <- sqrt(sum(resid(mu_model)^2 / mu_model$df.residual)) 
  
  # Treatment Model (eta)
  if(model_types$eta == 'levis_paper' | model_types$eta == 'alternative') {
    eta_model <- 
      glm(A ~ GENDER + bbl + raceeth_hispanic + I(bbl^2) + GENDER:raceeth_hispanic,
          data = df, 
          family = 'binomial')
  } else if(model_types$eta == 'levis_flexible') {
    eta_model <- 
      gam(A ~ 
            ### Discrete terms + Interactions
            GENDER + raceeth_hispanic + GENDER:raceeth_hispanic +
            
            ### BBL Splines (by discrete interaction terms)
            s(bbl, bs = 'cs') + 
            s(bbl, bs = 'cs', by = GENDER) + 
            s(bbl, bs = 'cs', by = raceeth_hispanic), 
          
          family = 'binomial',
          data = df)
  }
  
  # Missingness mechanism (pi)
  if(model_types$pi == 'levis_paper' | model_types$pi == 'alternative') {
    pi_model <- 
      glm(R ~ A + GENDER + bbl + raceeth_hispanic + ptwc + bbl:raceeth_hispanic +
            A:ptwc + ptwc:GENDER + bbl:A, 
          data = df,
          family ="binomial")
  } else if(model_types$pi == 'levis_flexible') {
    pi_model <- 
      gam(R ~ 
            ### Discrete terms + Interactions
            A + GENDER + raceeth_hispanic +
            A:GENDER + A:raceeth_hispanic + GENDER:raceeth_hispanic +
            
            ### BBL Splines (by discrete interaction terms)
            s(bbl, bs = 'cs') + 
            s(bbl, bs = 'cs', by = GENDER) + 
            s(bbl, bs = 'cs', by = raceeth_hispanic) + 
            s(bbl, bs = 'cs', by = A) + 
            
            ### ptwc Splines (by discrete interaction terms)
            s(ptwc, bs = 'cs') + 
            s(ptwc, bs = 'cs', by = GENDER) + 
            s(ptwc, bs = 'cs', by = raceeth_hispanic) + 
            s(ptwc, bs = 'cs', by = A) + 
            
            ### Multivariate Spline
            te(ptwc, bbl, bs = 'cs'),
          
          family = 'binomial',
          data = df)
  }
  
  # Imputation model (lambda)
  if(model_types$lambda == 'levis_paper' | model_types$lambda == 'alternative') {
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
  } else if(model_types$lambda == 'levis_flexible') {
    lambda_model <- 
      gam(cbl ~ 
            ### Discrete terms + Interactions
            A + GENDER + raceeth_hispanic +
            A:GENDER + A:raceeth_hispanic + GENDER:raceeth_hispanic +
            
            ### BBL Splines (by discrete interaction terms)
            s(bbl, bs = 'cs') + 
            s(bbl, bs = 'cs', by = GENDER) + 
            s(bbl, bs = 'cs', by = raceeth_hispanic) + 
            s(bbl, bs = 'cs', by = A) + 
            
            ### ptwc Splines (by discrete interaction terms)
            s(ptwc, bs = 'cs') + 
            s(ptwc, bs = 'cs', by = GENDER) + 
            s(ptwc, bs = 'cs', by = raceeth_hispanic) + 
            s(ptwc, bs = 'cs', by = A) + 
            
            ### Multivariate Spline
            te(ptwc, bbl, bs = 'cs'),
          
          family = Gamma(link = 'log'),
          data = df[df$R == 1,])
  } else if(model_types$lambda == 'levis_flexible_linear') {
      lambda_model <- 
        gam(cbl ~ 
              ### Discrete terms + Interactions
              A + GENDER + raceeth_hispanic +
              A:GENDER + A:raceeth_hispanic + GENDER:raceeth_hispanic +
              
              ### BBL Splines (by discrete interaction terms)
              s(bbl, bs = 'cs') + 
              s(bbl, bs = 'cs', by = GENDER) + 
              s(bbl, bs = 'cs', by = raceeth_hispanic) + 
              s(bbl, bs = 'cs', by = A) + 
              
              ### ptwc Splines (by discrete interaction terms)
              s(ptwc, bs = 'cs') + 
              s(ptwc, bs = 'cs', by = GENDER) + 
              s(ptwc, bs = 'cs', by = raceeth_hispanic) + 
              s(ptwc, bs = 'cs', by = A) + 
              
              ### Multivariate Spline
              te(ptwc, bbl, bs = 'cs'),
            
            family = 'gaussian',
            data = df[df$R == 1,])
  }
  
  # Lambda Shape
  if(model_types$lambda != 'levis_flexible_linear') {
    lambda_shape_hat <- 1/summary(lambda_model)$dispersion
    lambda_sigma_hat <- lambda_shape_hat ### Alias so no errors
  } else {
    lambda_sigma_hat <- sqrt(sum(resid(lambda_model)^2 / lambda_model$df.residual)) 
    lambda_shape_hat <- lambda_sigma_hat ### Alias so no errors
  }
  
  # 2nd Imputation model (smoking, rho) if applicable
  if(multiple_lp) {
    if(model_types$rho == 'durable') {
      rho_model <- 
        glm(smoker ~ GENDER + bbl + raceeth_hispanic + A + ptwc + cbl + bbl:ptwc,
            data = df[df$R == 1,],
            family = binomial(link = 'logit'))
    } else if(model_types$rho == 'durable_interact'| model_types$rho == 'alternative') {
      rho_model <- 
        glm(smoker ~ GENDER + bbl + raceeth_hispanic + A + ptwc + cbl + bbl:ptwc + A:ptwc,
            data = df[df$R == 1,],
            family = binomial(link = 'logit'))
    } else if(model_types$rho == 'aylp1_interact') {
      rho_model <- 
        glm(smoker ~ A + ptwc + cbl + A:ptwc + A:cbl + cbl:ptwc,
            data = df[df$R == 1,],
            family = binomial(link = 'logit'))
    } else if(model_types$rho == 'cbl_only') {
      rho_model <- 
        glm(smoker ~ cbl,
            data = df[df$R == 1,],
            family = binomial(link = 'logit'))
    } else if(model_types$rho == 'cbl_ptwc') {
      rho_model <- 
        glm(smoker ~ cbl + ptwc,
            data = df[df$R == 1,],
            family = binomial(link = 'logit'))
    } else if(model_types$rho == 'levis_flexible') {
      rho_model <- 
        gam(smoker ~ 
              ### Discrete terms + Interactions
              A + GENDER + raceeth_hispanic +
              A:GENDER + A:raceeth_hispanic + GENDER:raceeth_hispanic +
              
              ### BBL Splines (by discrete interaction terms)
              s(bbl, bs = 'cs') + 
              s(bbl, bs = 'cs', by = GENDER) + 
              s(bbl, bs = 'cs', by = raceeth_hispanic) + 
              s(bbl, bs = 'cs', by = A) + 
              
              ### ptwc Splines (by discrete interaction terms)
              s(ptwc, bs = 'cs') + 
              s(ptwc, bs = 'cs', by = GENDER) + 
              s(ptwc, bs = 'cs', by = raceeth_hispanic) + 
              s(ptwc, bs = 'cs', by = A) + 
              
              ### cbl Splines (by discrete interaction terms)
              s(cbl, bs = 'cs') + 
              s(cbl, bs = 'cs', by = GENDER) + 
              s(cbl, bs = 'cs', by = raceeth_hispanic) + 
              s(cbl, bs = 'cs', by = A) + 
              
              ### Multivariate Spline
              te(ptwc, bbl, bs = 'cs') + 
              te(ptwc, cbl, bs = 'cs') + 
              te(bbl, cbl, bs = 'cs'),
            
            family = 'binomial',
            data = df[df$R == 1,])
    }
  }
  
  ### Integration over gamma imputation model
  if(model_types$lambda != 'levis_flexible_linear') {
    ### Gauss-Laguerre for Gamma Integral
    gl <- gauss.quad(n = 4, kind = "laguerre", alpha = lambda_shape_hat - 1)
    gl_nodes <- gl$nodes
    gl_weights <- gl$weights
  } else {
    ### Gauss-Herminte for Normal
    ### Keep gl to distinguish which integral (imputation) rather than (outcome)
    gl <- gauss.quad(n = 4, kind = "hermite")
    gl_nodes <- gl$nodes
    gl_weights <- gl$weights
  }
  lp1_integral <- ifelse(model_types$lambda != 'levis_flexible_linear', 'gamma', 'normal')
  
  ### Predict models (for all subjects, in one go)
  ### Note we never need Rho hat based on fully observed LP1 
  ### Rather we'll need it later on based on GL Nodes 
  df <- 
    df %>% 
    mutate('pi_hat' = predict(pi_model, newdata = df, type = 'response'),
           'eta_hat' = predict(eta_model, newdata = df, type = 'response'))
  
  df <- 
    df %>% 
    mutate('lambda_hat' = predict(lambda_model, newdata = df, type = 'response')) %>% 
    ### Predicted mu component under both counterfactuals 
    mutate('mu_a1' = predict(mu_model, newdata = mutate(df, A = 1)),
           'mu_a0' = predict(mu_model, newdata = mutate(df, A = 0))) %>% 
    ### Components needed for gamma_beta computation
    mutate('y_values_a1' = map(mu_a1, ~{sqrt(2) * sigma_hat * gh_nodes + .x}),
           'y_values_a0' = map(mu_a0, ~{sqrt(2) * sigma_hat * gh_nodes + .x})) %>% 
    mutate('lambdas_a1' = compute_lambda_preds(., lambda_model, a = 1),
           'lambdas_a0' = compute_lambda_preds(., lambda_model, a = 0)) 
  
  ### Lp2 components needed for gamma beta computation (Gauss Hermite)
  ### P(lp2 = 1 | A = a, Y = y_values_a, Lc, Lp1 = lambda_a)
  if(multiple_lp) {
    df <- 
      df %>% 
      mutate('rhos_a1' = compute_rho_preds(., rho_model, a = 1),
             'rhos_a0' = compute_rho_preds(., rho_model, a = 0),
             'rhos_gl' = compute_rho_gl(df = ., 
                                        lp1 = lambda_hat,
                                        gl_nodes = gl_nodes, 
                                        shape_hat = ifelse(lp1_integral == 'gamma', 
                                                           lambda_shape_hat, 
                                                           lambda_sigma_hat),
                                        rho_model = rho_model,
                                        integral = lp1_integral))
    
  }
  
  ### Compute Gamma/Beta for each subject 
  if(lp1_integral == 'gamma') {
    if(!multiple_lp) {
      df <- 
        df %>% 
        mutate('gammas_betas_1' = map(1:nrow(.), ~gamma_beta(y_vals = y_values_a1[[.x]],
                                                             lambdas = lambdas_a1[[.x]],
                                                             shape_hat = lambda_shape_hat,
                                                             gh_weights = gh_weights,
                                                             v = gl_nodes * lambda_hat[.x]/lambda_shape_hat)),
               'gammas_betas_0' = map(1:nrow(.), ~gamma_beta(y_vals = y_values_a0[[.x]],
                                                             lambdas = lambdas_a0[[.x]],
                                                             shape_hat = lambda_shape_hat,
                                                             gh_weights = gh_weights,
                                                             v = gl_nodes * lambda_hat[.x]/lambda_shape_hat))) %>% 
        mutate('gb1' = map(gammas_betas_1, ~{.x$beta/.x$gamma}),
               'gb0' = map(gammas_betas_0, ~{.x$beta/.x$gamma})) 
    } else {
      df <- 
        df %>% 
        mutate('gammas_betas_a1_lp2_1' = map(1:nrow(.), ~gamma_beta(y_vals = y_values_a1[[.x]],
                                                                    lambdas = lambdas_a1[[.x]],
                                                                    shape_hat = lambda_shape_hat,
                                                                    gh_weights = gh_weights,
                                                                    v = gl_nodes * lambda_hat[.x]/lambda_shape_hat,
                                                                    rhos = rhos_a1[[.x]],
                                                                    lp2 = 1)),
               'gammas_betas_a0_lp2_1' = map(1:nrow(.), ~gamma_beta(y_vals = y_values_a0[[.x]],
                                                                    lambdas = lambdas_a0[[.x]],
                                                                    shape_hat = lambda_shape_hat,
                                                                    gh_weights = gh_weights,
                                                                    v = gl_nodes * lambda_hat[.x]/lambda_shape_hat,
                                                                    rhos = rhos_a0[[.x]],
                                                                    lp2 = 1)),
               'gammas_betas_a1_lp2_0' = map(1:nrow(.), ~gamma_beta(y_vals = y_values_a1[[.x]],
                                                                    lambdas = lambdas_a1[[.x]],
                                                                    shape_hat = lambda_shape_hat,
                                                                    gh_weights = gh_weights,
                                                                    v = gl_nodes * lambda_hat[.x]/lambda_shape_hat,
                                                                    rhos = rhos_a1[[.x]],
                                                                    lp2 = 0)),
               'gammas_betas_a0_lp2_0' = map(1:nrow(.), ~gamma_beta(y_vals = y_values_a0[[.x]],
                                                                    lambdas = lambdas_a0[[.x]],
                                                                    shape_hat = lambda_shape_hat,
                                                                    gh_weights = gh_weights,
                                                                    v = gl_nodes * lambda_hat[.x]/lambda_shape_hat,
                                                                    rhos = rhos_a0[[.x]],
                                                                    lp2 = 0))) %>% 
        mutate('gb_a1_lp2_1' = map(gammas_betas_a1_lp2_1, ~{.x$beta/.x$gamma}),
               'gb_a0_lp2_1' = map(gammas_betas_a0_lp2_1, ~{.x$beta/.x$gamma}),
               'gb_a1_lp2_0' = map(gammas_betas_a1_lp2_0, ~{.x$beta/.x$gamma}),
               'gb_a0_lp2_0' = map(gammas_betas_a0_lp2_0, ~{.x$beta/.x$gamma})) 
      
      
    }
  } else if(lp1_integral == 'normal') {
    if(!multiple_lp) {
      df <- 
        df %>% 
        mutate('gammas_betas_1' = map(1:nrow(.), ~gamma_beta(y_vals = y_values_a1[[.x]],
                                                             lambdas = lambdas_a1[[.x]],
                                                             shape_hat = lambda_sigma_hat,
                                                             gh_weights = gh_weights,
                                                             v = gl_nodes * sqrt(2) * lambda_sigma_hat + lambda_hat[.x],
                                                             integral = 'normal')),
               'gammas_betas_0' = map(1:nrow(.), ~gamma_beta(y_vals = y_values_a0[[.x]],
                                                             lambdas = lambdas_a0[[.x]],
                                                             shape_hat = lambda_sigma_hat,
                                                             gh_weights = gh_weights,
                                                             v = gl_nodes * sqrt(2) * lambda_sigma_hat + lambda_hat[.x],
                                                             integral = 'normal'))) %>% 
        mutate('gb1' = map(gammas_betas_1, ~{.x$beta/.x$gamma}),
               'gb0' = map(gammas_betas_0, ~{.x$beta/.x$gamma})) 
    } else {
      df <- 
        df %>% 
        mutate('gammas_betas_a1_lp2_1' = map(1:nrow(.), ~gamma_beta(y_vals = y_values_a1[[.x]],
                                                                    lambdas = lambdas_a1[[.x]],
                                                                    shape_hat = lambda_sigma_hat,
                                                                    gh_weights = gh_weights,
                                                                    v = gl_nodes * sqrt(2) * lambda_sigma_hat + lambda_hat[.x],
                                                                    integral = 'normal',
                                                                    rhos = rhos_a1[[.x]],
                                                                    lp2 = 1)),
               'gammas_betas_a0_lp2_1' = map(1:nrow(.), ~gamma_beta(y_vals = y_values_a0[[.x]],
                                                                    lambdas = lambdas_a0[[.x]],
                                                                    shape_hat = lambda_sigma_hat,
                                                                    gh_weights = gh_weights,
                                                                    v = gl_nodes * sqrt(2) * lambda_sigma_hat + lambda_hat[.x],
                                                                    integral = 'normal',
                                                                    rhos = rhos_a0[[.x]],
                                                                    lp2 = 1)),
               'gammas_betas_a1_lp2_0' = map(1:nrow(.), ~gamma_beta(y_vals = y_values_a1[[.x]],
                                                                    lambdas = lambdas_a1[[.x]],
                                                                    shape_hat = lambda_sigma_hat,
                                                                    gh_weights = gh_weights,
                                                                    v = gl_nodes * sqrt(2) * lambda_sigma_hat + lambda_hat[.x],
                                                                    integral = 'normal',
                                                                    rhos = rhos_a1[[.x]],
                                                                    lp2 = 0)),
               'gammas_betas_a0_lp2_0' = map(1:nrow(.), ~gamma_beta(y_vals = y_values_a0[[.x]],
                                                                    lambdas = lambdas_a0[[.x]],
                                                                    shape_hat = lambda_sigma_hat,
                                                                    gh_weights = gh_weights,
                                                                    v = gl_nodes * sqrt(2) * lambda_sigma_hat + lambda_hat[.x],
                                                                    integral = 'normal',
                                                                    rhos = rhos_a0[[.x]],
                                                                    lp2 = 0))) %>% 
        mutate('gb_a1_lp2_1' = map(gammas_betas_a1_lp2_1, ~{.x$beta/.x$gamma}),
               'gb_a0_lp2_1' = map(gammas_betas_a0_lp2_1, ~{.x$beta/.x$gamma}),
               'gb_a1_lp2_0' = map(gammas_betas_a1_lp2_0, ~{.x$beta/.x$gamma}),
               'gb_a0_lp2_0' = map(gammas_betas_a0_lp2_0, ~{.x$beta/.x$gamma})) 
      
      
    }
  }
  
  ### Compute ba1 (outcome_sub_a) and ba2 (IPTW_sub_a)
  if(lp1_integral == 'gamma') {
    if(!multiple_lp) {
      df <- 
        df %>% 
        mutate('outcome_sub_1' = map_dbl(gb1, ~sum(gl_weights * .x)/gamma(lambda_shape_hat)),
               'outcome_sub_0' = map_dbl(gb0, ~sum(gl_weights * .x)/gamma(lambda_shape_hat))) %>% 
        mutate('tau_comp_1' = map2(gammas_betas_1, eta_hat, ~{.x$gamma * .y}),
               'tau_comp_0' = map2(gammas_betas_0, eta_hat, ~{.x$gamma * (1-.y)}),
               'taus' = map2(tau_comp_1, tau_comp_0, ~{.x + .y})) %>% 
        mutate('IPTW_sub1' = iptw(gl_weights, taus, ptwc, gb1, gammas_betas_1, lambda_shape_hat),
               'IPTW_sub0' = iptw(gl_weights, taus, ptwc, gb0, gammas_betas_0, lambda_shape_hat))
    } else {
      df <- 
        df %>% 
        mutate('outcome_sub_a1_lp2_1' = map2_dbl(gb_a1_lp2_1, rhos_gl, ~sum(gl_weights * .x * .y)/gamma(lambda_shape_hat)),
               'outcome_sub_a1_lp2_0' = map2_dbl(gb_a1_lp2_0, rhos_gl, ~sum(gl_weights * .x * (1-.y))/gamma(lambda_shape_hat)),
               'outcome_sub_a0_lp2_1' = map2_dbl(gb_a0_lp2_1, rhos_gl, ~sum(gl_weights * .x * .y)/gamma(lambda_shape_hat)),
               'outcome_sub_a0_lp2_0' = map2_dbl(gb_a0_lp2_0, rhos_gl, ~sum(gl_weights * .x * (1-.y))/gamma(lambda_shape_hat))) %>% 
        mutate('outcome_sub_1' = outcome_sub_a1_lp2_1 + outcome_sub_a1_lp2_0,
               'outcome_sub_0' = outcome_sub_a0_lp2_1 + outcome_sub_a0_lp2_0) %>% 
        mutate('tau_comp_a1_lp2_1' = map2(gammas_betas_a1_lp2_1, eta_hat, ~{.x$gamma * .y}),
               'tau_comp_a0_lp2_1' = map2(gammas_betas_a0_lp2_1, eta_hat, ~{.x$gamma * (1-.y)}),
               'tau_comp_a1_lp2_0' = map2(gammas_betas_a1_lp2_0, eta_hat, ~{.x$gamma * .y}),
               'tau_comp_a0_lp2_0' = map2(gammas_betas_a0_lp2_0, eta_hat, ~{.x$gamma * (1-.y)}),
               'taus_lp2_1' = map2(tau_comp_a1_lp2_1, tau_comp_a0_lp2_1, ~{.x + .y}),
               'taus_lp2_0' = map2(tau_comp_a1_lp2_1, tau_comp_a0_lp2_1, ~{.x + .y})) %>% 
        mutate('IPTW_sub_a1_lp2_1' = iptw(gl_weights = gl_weights, 
                                          taus = taus_lp2_1, 
                                          y = ptwc, 
                                          gb = gb_a1_lp2_1, 
                                          gammas_betas = gammas_betas_a1_lp2_1, 
                                          lambda_shape_hat = lambda_shape_hat,
                                          rhos = rhos_gl, 
                                          lp2 = 1),
               'IPTW_sub_a1_lp2_0' = iptw(gl_weights = gl_weights, 
                                          taus = taus_lp2_0, 
                                          y = ptwc, 
                                          gb = gb_a1_lp2_0, 
                                          gammas_betas = gammas_betas_a1_lp2_0, 
                                          lambda_shape_hat = lambda_shape_hat,
                                          rhos = rhos_gl, 
                                          lp2 = 0),
               'IPTW_sub_a0_lp2_1' = iptw(gl_weights = gl_weights, 
                                          taus = taus_lp2_1, 
                                          y = ptwc, 
                                          gb = gb_a0_lp2_1, 
                                          gammas_betas = gammas_betas_a0_lp2_1, 
                                          lambda_shape_hat = lambda_shape_hat,
                                          rhos = rhos_gl, 
                                          lp2 = 1),
               'IPTW_sub_a0_lp2_0' = iptw(gl_weights = gl_weights, 
                                          taus = taus_lp2_0, 
                                          y = ptwc, 
                                          gb = gb_a0_lp2_0, 
                                          gammas_betas = gammas_betas_a0_lp2_0, 
                                          lambda_shape_hat = lambda_shape_hat,
                                          rhos = rhos_gl, 
                                          lp2 = 0),
               'IPTW_sub0' = IPTW_sub_a0_lp2_1 + IPTW_sub_a0_lp2_0,
               'IPTW_sub1' = IPTW_sub_a1_lp2_1 + IPTW_sub_a1_lp2_0)
      
    }
  } else if(lp1_integral == 'normal') {
    if(!multiple_lp) {
      df <- 
        df %>% 
        mutate('outcome_sub_1' = map_dbl(gb1, ~sum(gl_weights * .x)/sqrt(pi)),
               'outcome_sub_0' = map_dbl(gb0, ~sum(gl_weights * .x)/sqrt(pi))) %>% 
        mutate('tau_comp_1' = map2(gammas_betas_1, eta_hat, ~{.x$gamma * .y}),
               'tau_comp_0' = map2(gammas_betas_0, eta_hat, ~{.x$gamma * (1-.y)}),
               'taus' = map2(tau_comp_1, tau_comp_0, ~{.x + .y})) %>% 
        mutate('IPTW_sub1' = iptw(gl_weights, taus, ptwc, gb1, gammas_betas_1, lambda_sigma_hat, integral = 'normal'),
               'IPTW_sub0' = iptw(gl_weights, taus, ptwc, gb0, gammas_betas_0, lambda_sigma_hat, integral = 'normal'))
    } else {
      df <- 
        df %>% 
        mutate('outcome_sub_a1_lp2_1' = map2_dbl(gb_a1_lp2_1, rhos_gl, ~sum(gl_weights * .x * .y)/sqrt(pi)),
               'outcome_sub_a1_lp2_0' = map2_dbl(gb_a1_lp2_0, rhos_gl, ~sum(gl_weights * .x * (1-.y))/sqrt(pi)),
               'outcome_sub_a0_lp2_1' = map2_dbl(gb_a0_lp2_1, rhos_gl, ~sum(gl_weights * .x * .y)/sqrt(pi)),
               'outcome_sub_a0_lp2_0' = map2_dbl(gb_a0_lp2_0, rhos_gl, ~sum(gl_weights * .x * (1-.y))/sqrt(pi))) %>% 
        mutate('outcome_sub_1' = outcome_sub_a1_lp2_1 + outcome_sub_a1_lp2_0,
               'outcome_sub_0' = outcome_sub_a0_lp2_1 + outcome_sub_a0_lp2_0) %>% 
        mutate('tau_comp_a1_lp2_1' = map2(gammas_betas_a1_lp2_1, eta_hat, ~{.x$gamma * .y}),
               'tau_comp_a0_lp2_1' = map2(gammas_betas_a0_lp2_1, eta_hat, ~{.x$gamma * (1-.y)}),
               'tau_comp_a1_lp2_0' = map2(gammas_betas_a1_lp2_0, eta_hat, ~{.x$gamma * .y}),
               'tau_comp_a0_lp2_0' = map2(gammas_betas_a0_lp2_0, eta_hat, ~{.x$gamma * (1-.y)}),
               'taus_lp2_1' = map2(tau_comp_a1_lp2_1, tau_comp_a0_lp2_1, ~{.x + .y}),
               'taus_lp2_0' = map2(tau_comp_a1_lp2_1, tau_comp_a0_lp2_1, ~{.x + .y})) %>% 
        mutate('IPTW_sub_a1_lp2_1' = iptw(gl_weights = gl_weights, 
                                          taus = taus_lp2_1, 
                                          y = ptwc, 
                                          gb = gb_a1_lp2_1, 
                                          gammas_betas = gammas_betas_a1_lp2_1, 
                                          lambda_shape_hat = lambda_sigma_hat,
                                          rhos = rhos_gl, 
                                          lp2 = 1,
                                          integral = 'normal'),
               'IPTW_sub_a1_lp2_0' = iptw(gl_weights = gl_weights, 
                                          taus = taus_lp2_0, 
                                          y = ptwc, 
                                          gb = gb_a1_lp2_0, 
                                          gammas_betas = gammas_betas_a1_lp2_0, 
                                          lambda_shape_hat = lambda_sigma_hat,
                                          rhos = rhos_gl, 
                                          lp2 = 0, 
                                          integral = 'normal'),
               'IPTW_sub_a0_lp2_1' = iptw(gl_weights = gl_weights, 
                                          taus = taus_lp2_1, 
                                          y = ptwc, 
                                          gb = gb_a0_lp2_1, 
                                          gammas_betas = gammas_betas_a0_lp2_1, 
                                          lambda_shape_hat = lambda_sigma_hat,
                                          rhos = rhos_gl, 
                                          lp2 = 1,
                                          integral = 'normal'),
               'IPTW_sub_a0_lp2_0' = iptw(gl_weights = gl_weights, 
                                          taus = taus_lp2_0, 
                                          y = ptwc, 
                                          gb = gb_a0_lp2_0, 
                                          gammas_betas = gammas_betas_a0_lp2_0, 
                                          lambda_shape_hat = lambda_sigma_hat,
                                          rhos = rhos_gl, 
                                          lp2 = 0,
                                          integral = 'normal'),
               'IPTW_sub0' = IPTW_sub_a0_lp2_1 + IPTW_sub_a0_lp2_0,
               'IPTW_sub1' = IPTW_sub_a1_lp2_1 + IPTW_sub_a1_lp2_0)
      
    }
  }
  
  if(!multiple_lp) {
    df <- 
      df %>% 
      mutate('gamma_beta_1' = map(1:nrow(.), ~gamma_beta(y_vals = y_values_a1[[.x]],
                                                         lambdas = lambdas_a1[[.x]],
                                                         shape_hat = ifelse(lp1_integral == 'gamma',
                                                                            lambda_shape_hat,
                                                                            lambda_sigma_hat),
                                                         gh_weights = gh_weights,
                                                         v = cbl[.x],
                                                         integral = lp1_integral)),
             'gamma_beta_0' = map(1:nrow(.), ~gamma_beta(y_vals = y_values_a0[[.x]],
                                                         lambdas = lambdas_a0[[.x]],
                                                         shape_hat = ifelse(lp1_integral == 'gamma',
                                                                            lambda_shape_hat,
                                                                            lambda_sigma_hat),
                                                         gh_weights = gh_weights,
                                                         v = cbl[.x],
                                                         integral = lp1_integral))) 
  } else {
    df <- 
      df %>% 
      mutate('gamma_beta_1' = map(1:nrow(.), ~gamma_beta(y_vals = y_values_a1[[.x]],
                                                         lambdas = lambdas_a1[[.x]],
                                                         shape_hat = ifelse(lp1_integral == 'gamma',
                                                                            lambda_shape_hat,
                                                                            lambda_sigma_hat),
                                                         gh_weights = gh_weights,
                                                         v = cbl[.x],
                                                         rhos = rhos_a1[[.x]],
                                                         lp2 = smoker[.x],
                                                         integral = lp1_integral)),
             'gamma_beta_0' = map(1:nrow(.), ~gamma_beta(y_vals = y_values_a0[[.x]],
                                                         lambdas = lambdas_a0[[.x]],
                                                         shape_hat = ifelse(lp1_integral == 'gamma',
                                                                            lambda_shape_hat,
                                                                            lambda_sigma_hat),
                                                         gh_weights = gh_weights,
                                                         v = cbl[.x],
                                                         rhos = rhos_a0[[.x]],
                                                         lp2 = smoker[.x],
                                                         integral = lp1_integral))) 
  }
  
  
  ### Compute IPW/IF Estimators
  df <- 
    df %>% 
    mutate('IPW_1' = pmap_dbl(.l = list(R, gamma_beta_1, pi_hat), ~ifelse(..1 == 1, ..2$beta/..2$gamma/..3, 0)),
           'IPW_0' = pmap_dbl(.l = list(R, gamma_beta_0, pi_hat), ~ifelse(..1 == 1, ..2$beta/..2$gamma/..3, 0))) %>% 
    mutate('tau' = pmap_dbl(.l = list(gamma_beta_1, gamma_beta_0, eta_hat), ~{..1$gamma*..3 + ..2$gamma*(1-..3)})) %>% 
    mutate('gb_ratio_1' = map_dbl(gamma_beta_1, ~{.x$beta/.x$gamma}),
           'gb_ratio_0' = map_dbl(gamma_beta_0, ~{.x$beta/.x$gamma}),
           'gamma_1' = map_dbl(gamma_beta_1, ~{.x$gamma}),
           'gamma_0' = map_dbl(gamma_beta_0, ~{.x$gamma})) %>% 
    mutate('IF_1' = case_when(R == 0 ~ outcome_sub_1 + (A/eta_hat) * IPTW_sub1,
                              R == 1 ~ outcome_sub_1 + (gb_ratio_1 - outcome_sub_1)/pi_hat + 
                                (A/eta_hat) * (IPTW_sub1 + ((tau/gamma_1) * (ptwc - gb_ratio_1) - IPTW_sub1)/pi_hat)),
           'IF_0' = case_when(R == 0 ~ outcome_sub_0 + ((1-A)/(1-eta_hat)) * IPTW_sub0,
                              R == 1 ~ outcome_sub_0 + (gb_ratio_0 - outcome_sub_0)/pi_hat + 
                                ((1-A)/(1-eta_hat)) * (IPTW_sub0 + ((tau/gamma_0) * (ptwc - gb_ratio_0) - IPTW_sub0)/pi_hat))) %>% 
    
    ### In some extreme cases we can't estimate a treatment effect
    filter(!is.na(IF_0), !is.na(IF_1)) %>% 
    
    ### Clean up
    mutate_at(vars(everything()), ~unname(.x)) %>% 
    select(IPW_0, IPW_1, IF_0, IF_1) %>% 
    rename('IWOR_1' = IPW_1,
           'IWOR_0' = IPW_0) %>% 
    summarise(across(everything(), mean))
  
  return(df)
}

### Function to call compute_levis_estimators over batches (helps lower 
### amount of memory being used at any one time
batch_levis_estimators <- function(datasets, batches, model_types, multiple_lp = F, seed) {
  set.seed(seed)
  seeds <- sample(1:1e7, length(batches))
  
  df <- NULL
  for(i in 1:length(batches)) {
    tmp <- 
      future_map_dfr(datasets[ batches[[i]] ], ~compute_levis_estimators(df = .x,
                                                                         model_types = model_types,
                                                                         multiple_lp = multiple_lp),
                     .options = furrr_options(seed = seeds[i])) 
    df <- bind_rows(df, tmp)
  }
  
  return(df)
}
