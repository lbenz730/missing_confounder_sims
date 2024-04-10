library(statmod)
library(tidyverse)
library(here)
source(here('helpers.R'))
source(here('levis_helpers.R'))
source(here('generate_data.R'))

### Function to compute true ATE
compute_truth <- function(params) {
  if(params$data_generating_process == 'Levis') {
    return(compute_truth_levis(params))
  } else if(params$data_generating_process == 'Alternative') {
    return(compute_truth_alternative(params))
  }
}

### Function to compute true ATE under levis decomposition
compute_truth_levis <- function(params) {
  load_models(params) 
  multiple_lp <- params$multiple_lp
  
  ### Gauss-Hermite prep (for integral over gaussian outcome distribution)
  gh <- gauss.quad(n = 15, kind = "hermite")
  gh_nodes <- gh$nodes
  gh_weights <- gh$weights
  
  ### Gauss-Laguerre (integration over gamma imputation model)
  gl <- gauss.quad(n = 15, kind = "laguerre", alpha = lambda_shape - 1)
  gl_nodes <- gl$nodes
  gl_weights <- gl$weights
  
  df_truth <- 
    df_cc %>% 
    ### Predict models (for all subjects, in one go)
    mutate('eta_hat' = expit(shift_preds(preds = predict(eta_star, newdata = df_cc),
                                         shift = eta_shift,
                                         df = df_cc))) %>% 
    ### Predicted mu component under both counterfactuals 
    mutate('mu_a1' = shift_preds(preds = predict(mu_star, newdata = mutate(df_cc, A = 1)),
                                 shift = mu_shift,
                                 df = mutate(df_cc, A = 1)),
           'mu_a0' = shift_preds(preds = predict(mu_star, newdata = mutate(df_cc, A = 0)),
                                 shift = mu_shift,
                                 df = mutate(df_cc, A = 0))) %>% 
    ### Components needed for gamma_beta computation
    mutate('y_values_a1' = map(mu_a1, ~{sqrt(2) * sigma_star * gh_nodes + .x}),
           'y_values_a0' = map(mu_a0, ~{sqrt(2) * sigma_star * gh_nodes + .x})) %>% 
    mutate('lambdas_a1' = compute_lambda_preds(., lambda_star, a = 1, pred_shift = lambda_shift),
           'lambdas_a0' = compute_lambda_preds(., lambda_star, a = 0, pred_shift = lambda_shift))
  
  ### Lp2 components needed for gamma beta computation (Gauss Hermite)
  ### P(lp2 = 1 | A = a, Y = y_values_a, Lc, Lp1 = lambda_a)
  if(multiple_lp) {
    df_truth  <- 
      df_truth  %>% 
      mutate('rhos_a1' = compute_rho_preds(., rho_star, a = 1, pred_shift = rho_shift),
             'rhos_a0' = compute_rho_preds(., rho_star, a = 0, pred_shift = rho_shift),
             'rhos_gl_a0' = compute_rho_gl(df = ., 
                                           lp1 = lambdas_a0,
                                           gl_nodes = gl_nodes, 
                                           shape_hat = lambda_shape, 
                                           rho_model = rho_star,
                                           pred_shift = rho_shift, 
                                           truth = T),
             'rhos_gl_a1' = compute_rho_gl(df = ., 
                                           lp1 = lambdas_a1,
                                           gl_nodes = gl_nodes, 
                                           shape_hat = lambda_shape, 
                                           rho_model = rho_star,
                                           pred_shift = rho_shift,
                                           truth = T))
    
  }
  
  ### Compute Gamma/Beta for each subject 
  if(!multiple_lp) {
    df_truth <- 
      df_truth %>% 
      mutate('gammas_betas_1_a1' = map(1:nrow(.), {
        function(row_id) {
          map(lambdas_a1[[row_id]],
              ~gamma_beta(y_vals = y_values_a1[[row_id]],
                          lambdas = lambdas_a1[[row_id]],
                          shape_hat = lambda_shape,
                          gh_weights = gh_weights,
                          v = gl_nodes * .x/lambda_shape))
        }}),
        'gammas_betas_0_a1' = map(1:nrow(.), {
          function(row_id) {
            map(lambdas_a0[[row_id]],
                ~gamma_beta(y_vals = y_values_a1[[row_id]],
                            lambdas = lambdas_a1[[row_id]],
                            shape_hat = lambda_shape,
                            gh_weights = gh_weights,
                            v = gl_nodes * .x/lambda_shape))
          }}),
        'gammas_betas_1_a0' = map(1:nrow(.), {
          function(row_id) {
            map(lambdas_a1[[row_id]],
                ~gamma_beta(y_vals = y_values_a0[[row_id]],
                            lambdas = lambdas_a0[[row_id]],
                            shape_hat = lambda_shape,
                            gh_weights = gh_weights,
                            v = gl_nodes * .x/lambda_shape))
          }}),
        'gammas_betas_0_a0' = map(1:nrow(.), {
          function(row_id) {
            map(lambdas_a0[[row_id]],
                ~gamma_beta(y_vals = y_values_a0[[row_id]],
                            lambdas = lambdas_a0[[row_id]],
                            shape_hat = lambda_shape,
                            gh_weights = gh_weights,
                            v = gl_nodes * .x/lambda_shape))  }})) %>% 
      mutate('gb1_a1' = map(gammas_betas_1_a1, ~{map(.x, function(l) {l$beta/l$gamma})}),
             'gb0_a1' = map(gammas_betas_0_a1, ~{map(.x, function(l) {l$beta/l$gamma})}),
             'gb1_a0' = map(gammas_betas_1_a0, ~{map(.x, function(l) {l$beta/l$gamma})}),
             'gb0_a0' = map(gammas_betas_0_a0, ~{map(.x, function(l) {l$beta/l$gamma})})) %>% 
      mutate('v_ints_1_a1' = map(gb1_a1, function(l) map_dbl(l, ~{sum(gl_weights * .x) / gamma(lambda_shape)})),
             'v_ints_0_a1' = map(gb0_a1, function(l) map_dbl(l, ~{sum(gl_weights * .x) / gamma(lambda_shape)})),
             'v_ints_1_a0' = map(gb1_a0, function(l) map_dbl(l, ~{sum(gl_weights * .x) / gamma(lambda_shape)})),
             'v_ints_0_a0' = map(gb0_a0, function(l) map_dbl(l, ~{sum(gl_weights * .x) / gamma(lambda_shape)}))) 
  } else {
    df_truth <- 
      df_truth %>% 
      mutate('gammas_betas_1_a1_lp2_1' = map(1:nrow(.), {
        function(row_id) {
          map(lambdas_a1[[row_id]],
              ~gamma_beta(y_vals = y_values_a1[[row_id]],
                          lambdas = lambdas_a1[[row_id]],
                          shape_hat = lambda_shape,
                          gh_weights = gh_weights,
                          v = gl_nodes * .x/lambda_shape,
                          rhos = rhos_a1[[row_id]],
                          lp2 = 1))
        }}),
        'gammas_betas_0_a1_lp2_1' = map(1:nrow(.), {
          function(row_id) {
            map(lambdas_a0[[row_id]],
                ~gamma_beta(y_vals = y_values_a1[[row_id]],
                            lambdas = lambdas_a1[[row_id]],
                            shape_hat = lambda_shape,
                            gh_weights = gh_weights,
                            v = gl_nodes * .x/lambda_shape,
                            rhos = rhos_a1[[row_id]],
                            lp2 = 1))
          }}),
        'gammas_betas_1_a0_lp2_1' = map(1:nrow(.), {
          function(row_id) {
            map(lambdas_a1[[row_id]],
                ~gamma_beta(y_vals = y_values_a0[[row_id]],
                            lambdas = lambdas_a0[[row_id]],
                            shape_hat = lambda_shape,
                            gh_weights = gh_weights,
                            v = gl_nodes * .x/lambda_shape,
                            rhos = rhos_a0[[row_id]],
                            lp2 = 1))
          }}),
        'gammas_betas_0_a0_lp2_1' = map(1:nrow(.), {
          function(row_id) {
            map(lambdas_a0[[row_id]],
                ~gamma_beta(y_vals = y_values_a0[[row_id]],
                            lambdas = lambdas_a0[[row_id]],
                            shape_hat = lambda_shape,
                            gh_weights = gh_weights,
                            v = gl_nodes * .x/lambda_shape,
                            rhos = rhos_a0[[row_id]],
                            lp2 = 1))  
          }}),
        'gammas_betas_1_a1_lp2_0' = map(1:nrow(.), {
          function(row_id) {
            map(lambdas_a1[[row_id]],
                ~gamma_beta(y_vals = y_values_a1[[row_id]],
                            lambdas = lambdas_a1[[row_id]],
                            shape_hat = lambda_shape,
                            gh_weights = gh_weights,
                            v = gl_nodes * .x/lambda_shape,
                            rhos = rhos_a1[[row_id]],
                            lp2 = 0))
          }}),
        'gammas_betas_0_a1_lp2_0' = map(1:nrow(.), {
          function(row_id) {
            map(lambdas_a0[[row_id]],
                ~gamma_beta(y_vals = y_values_a1[[row_id]],
                            lambdas = lambdas_a1[[row_id]],
                            shape_hat = lambda_shape,
                            gh_weights = gh_weights,
                            v = gl_nodes * .x/lambda_shape,
                            rhos = rhos_a1[[row_id]],
                            lp2 = 0))
          }}),
        'gammas_betas_1_a0_lp2_0' = map(1:nrow(.), {
          function(row_id) {
            map(lambdas_a1[[row_id]],
                ~gamma_beta(y_vals = y_values_a0[[row_id]],
                            lambdas = lambdas_a0[[row_id]],
                            shape_hat = lambda_shape,
                            gh_weights = gh_weights,
                            v = gl_nodes * .x/lambda_shape,
                            rhos = rhos_a0[[row_id]],
                            lp2 = 0))
          }}),
        'gammas_betas_0_a0_lp2_0' = map(1:nrow(.), {
          function(row_id) {
            map(lambdas_a0[[row_id]],
                ~gamma_beta(y_vals = y_values_a0[[row_id]],
                            lambdas = lambdas_a0[[row_id]],
                            shape_hat = lambda_shape,
                            gh_weights = gh_weights,
                            v = gl_nodes * .x/lambda_shape,
                            rhos = rhos_a0[[row_id]],
                            lp2 = 0))  }}))  %>% 
      mutate('gb1_a1_lp2_1' = map(gammas_betas_1_a1_lp2_1, ~{map(.x, function(l) {l$beta/l$gamma})}),
             'gb0_a1_lp2_1' = map(gammas_betas_0_a1_lp2_1, ~{map(.x, function(l) {l$beta/l$gamma})}),
             'gb1_a0_lp2_1' = map(gammas_betas_1_a0_lp2_1, ~{map(.x, function(l) {l$beta/l$gamma})}),
             'gb0_a0_lp2_1' = map(gammas_betas_0_a0_lp2_1, ~{map(.x, function(l) {l$beta/l$gamma})}),
             'gb1_a1_lp2_0' = map(gammas_betas_1_a1_lp2_0, ~{map(.x, function(l) {l$beta/l$gamma})}),
             'gb0_a1_lp2_0' = map(gammas_betas_0_a1_lp2_0, ~{map(.x, function(l) {l$beta/l$gamma})}),
             'gb1_a0_lp2_0' = map(gammas_betas_1_a0_lp2_0, ~{map(.x, function(l) {l$beta/l$gamma})}),
             'gb0_a0_lp2_0' = map(gammas_betas_0_a0_lp2_0, ~{map(.x, function(l) {l$beta/l$gamma})})) %>% 
      mutate('v_ints_1_a1_lp2_1' = map2(gb1_a1_lp2_1, rhos_gl_a1, function(l, m) map2_dbl(l, m, ~{sum(gl_weights * .x * .y) / gamma(lambda_shape)})),
             'v_ints_0_a1_lp2_1' = map2(gb0_a1_lp2_1, rhos_gl_a1, function(l, m) map2_dbl(l, m,~{sum(gl_weights * .x * .y) / gamma(lambda_shape)})),
             'v_ints_1_a0_lp2_1' = map2(gb1_a0_lp2_1, rhos_gl_a0, function(l, m) map2_dbl(l, m, ~{sum(gl_weights * .x * .y) / gamma(lambda_shape)})),
             'v_ints_0_a0_lp2_1' = map2(gb0_a0_lp2_1, rhos_gl_a0, function(l, m) map2_dbl(l, m, ~{sum(gl_weights * .x * .y) / gamma(lambda_shape)})),
             'v_ints_1_a1_lp2_0' = map2(gb1_a1_lp2_0, rhos_gl_a1, function(l, m) map2_dbl(l, m, ~{sum(gl_weights * .x * (1-.y)) / gamma(lambda_shape)})),
             'v_ints_0_a1_lp2_0' = map2(gb0_a1_lp2_0, rhos_gl_a1, function(l, m) map2_dbl(l, m, ~{sum(gl_weights * .x * (1-.y)) / gamma(lambda_shape)})),
             'v_ints_1_a0_lp2_0' = map2(gb1_a0_lp2_0, rhos_gl_a0, function(l, m) map2_dbl(l, m, ~{sum(gl_weights * .x * (1-.y)) / gamma(lambda_shape)})),
             'v_ints_0_a0_lp2_0' = map2(gb0_a0_lp2_0, rhos_gl_a0, function(l, m) map2_dbl(l, m, ~{sum(gl_weights * .x * (1-.y)) / gamma(lambda_shape)}))) %>% 
      mutate('v_ints_1_a1' = map2(v_ints_1_a1_lp2_1, v_ints_1_a1_lp2_0, ~{.x+.y}),
             'v_ints_0_a1' = map2(v_ints_0_a1_lp2_1, v_ints_0_a1_lp2_0, ~{.x+.y}),
             'v_ints_1_a0' = map2(v_ints_1_a0_lp2_1, v_ints_1_a0_lp2_0, ~{.x+.y}),
             'v_ints_0_a0' = map2(v_ints_0_a0_lp2_1, v_ints_0_a0_lp2_0, ~{.x+.y}))
    
  }
  
  df_truth <- 
    df_truth %>% 
    mutate('int_1_a1' = map_dbl(v_ints_1_a1, ~sum(.x * gh_weights)/sqrt(pi)),
           'int_0_a1' = map_dbl(v_ints_0_a1, ~sum(.x * gh_weights)/sqrt(pi)),
           'int_1_a0' = map_dbl(v_ints_1_a0, ~sum(.x * gh_weights)/sqrt(pi)),
           'int_0_a0' = map_dbl(v_ints_0_a0, ~sum(.x * gh_weights)/sqrt(pi))) %>% 
    mutate('contribution_a1' = eta_hat * int_1_a1 + (1 - eta_hat) * int_0_a1,
           'contribution_a0' = eta_hat * int_1_a0 + (1 - eta_hat) * int_0_a0) %>% 
    summarise('te_0' = mean(contribution_a0),
              'te_1' = mean(contribution_a1)) %>% 
    mutate('ate' = te_1 - te_0)
  
  return(df_truth)
}

### Function to compute true ATE under alternative data generating process
compute_truth_alternative <- function(params) {
  set.seed(13242)
  load_models(params) 
  df_mega <- 
    generate_data(n_patients = 10000000, 
                  data_generating_process = params$data_generating_process,
                  multiple_lp = params$multiple_lp, 
                  truth = T)
  
  df_truth <- 
    df_mega %>% 
    mutate('mu_hat1' = shift_preds(preds = predict(mu_star, newdata = mutate(df_mega, 'A' = 1)),
                                   shift = mu_shift,
                                   df =  mutate(df_mega, 'A' = 1))) %>% 
    mutate('mu_hat0' = shift_preds(preds = predict(mu_star, newdata = mutate(df_mega, 'A' = 0)),
                                   shift = mu_shift,
                                   df =  mutate(df_mega, 'A' = 0))) %>%
    summarise('te_1' = mean(mu_hat1),
              'te_0' = mean(mu_hat0),
              'ate' = te_1-te_0)
  
  return(df_truth)
  
}

