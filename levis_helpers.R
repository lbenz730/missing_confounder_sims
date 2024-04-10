### Helper Function to compute lambda preds for different y-values in numerical integration
compute_lambda_preds <- function(df, lambda_model, a, pred_shift = NULL) {
  if(a == 1) {
    y_vals <- df$y_values_a1
  } else {
    y_vals <- df$y_values_a0 
  }
  
  ### Turn to list of vectors, indexed by patient (1 for each gh_node)
  y_vals <- map(transpose(y_vals), unlist)
  
  ### Predictions
  if(is.null(pred_shift)) {
    preds <- map(y_vals, ~predict(lambda_model, 
                                  newdata = mutate(df, A = a, ptwc = .x), 
                                  type = 'response'))
  } else {
    preds <- map(y_vals, ~exp(shift_preds(preds = predict(lambda_model, 
                                                          newdata = mutate(df, A = a, ptwc = .x)),
                                          shift = pred_shift,
                                          df = mutate(df, A = a, ptwc = .x))))
  }
  
  ### Transform back into shape useful for data frame
  preds <- map(transpose(preds), unlist)
  
  return(preds)
}

### Helper Function to compute lp2 (rho) preds for different y-values/lp1 values in numerical integration
compute_rho_preds <- function(df, rho_model, a, pred_shift = NULL) {
  if(a == 1) {
    y_vals <- df$y_values_a1
    lambda_preds <- df$lambdas_a1
  } else {
    y_vals <- df$y_values_a0 
    lambda_preds <- df$lambdas_a0
  }
  
  ### Turn to list of vectors, indexed by patient (1 for each gh_node)
  y_vals <- map(transpose(y_vals), unlist)
  lambda_preds <- map(transpose(lambda_preds), unlist)
  
  ### Predictions
  if(is.null(pred_shift)) {
    preds <- 
      map2(y_vals, lambda_preds, ~pred_cap(predict(rho_model, 
                                                   newdata = mutate(df, A = a, ptwc = .x, cbl = .y), 
                                                   type = 'response')))
  } else {
    preds <- 
      map2(y_vals, lambda_preds, ~pred_cap(expit(shift_preds(preds = predict_rho_star(rho_model,
                                                                                      newdata = mutate(df, A = a, ptwc = .x, cbl = .y)),
                                                             shift = pred_shift,
                                                             df = mutate(df, A = a, ptwc = .x, cbl = .y)))))
  }
  
  ### Transform back into shape useful for data frame
  preds <- map(transpose(preds), unlist)
  
  return(preds)
}

### Function to predict LP2 = smoker from list of specified covariates
predict_rho_star <- function(rho_star, newdata, type = 'raw') {
  rho_log_odds <- rho_star$prevalence 
  for(i in 2:length(rho_star)) {
    beta <- rho_star[[i]]
    vars <- unlist(str_split(names(rho_star)[[i]], '_'))[-1]
    if(all(vars %in% names(newdata))) {
      if(length(vars) == 1) {
        rho_log_odds <- rho_log_odds + beta * newdata[[vars]]
      } else if(length(vars) == 2) {
        rho_log_odds <- rho_log_odds + beta * newdata[[ vars[1] ]] * newdata[[ vars[2] ]]
      }
    }
  }
  
  if(type == 'response') {
    rho_p <- expit(rho_log_odds) 
    return(rho_p)
  }
  
  return(rho_log_odds)
}

### Function to get rho values needed for collapsing gl integral
compute_rho_gl <- function(df, lp1, gl_nodes, shape_hat, rho_model, pred_shift = NULL, truth = F, integral = 'gamma') {
  ### Where we just have 1 lambda (lambda hat) to apply rho_gl preds for
  if(!is.list(lp1)) {
    if(integral == 'gamma') {
      v <- map(lp1, ~gl_nodes * .x/shape_hat)
    } else if(integral == 'normal') {
      ### If integral is normal shape_hat = sigma_hat
      v <- map(lp1, ~{gl_nodes * sqrt(2) * shape_hat + .x})
    }
    
    ### Turn to list of vectors, indexed by patient (1 for each gl_node)
    v <- map(transpose(v), unlist)
    
    ### Predictions
    if(is.null(pred_shift)) {
      preds <- map(v, ~pred_cap(predict(rho_model, newdata = mutate(df, cbl = .x), type = 'response')))
    } else {
      preds <- map(v, ~pred_cap(expit(shift_preds(preds = predict(rho_model, newdata = mutate(df, cbl = .x)),
                                                  shift = pred_shift,
                                                  df = mutate(df, cbl = .x)))))
    }
    
    ### Transform back into shape useful for data frame
    preds <- map(transpose(preds), unlist)
  } else {
    if(integral == 'gamma') {
      v <- map(lp1, function(lambda) {map(lambda, ~gl_nodes * .x/shape_hat)})
    } else if(integral == 'normal') {
      ### If integral is normal shape_hat = sigma_hat
      v <- map(lp1, function(lambda) {map(lambda, ~{gl_nodes * sqrt(2) * shape_hat + .x})})
    }
    
    v <- transpose(v)
    
    preds <- 
      map(v, function(vi) {
        vi <- map(transpose(vi), unlist)
        if(!truth) {
          p <- map(vi, ~pred_cap(expit(shift_preds(preds = predict(rho_model, newdata = mutate(df, cbl = .x)),
                                                   shift = pred_shift,
                                                   df = mutate(df, cbl = .x)))))
        } else {
          p <- map(vi, ~pred_cap(expit(shift_preds(preds = predict_rho_star(rho_model, newdata = mutate(df, cbl = .x)),
                                                   shift = pred_shift,
                                                   df = mutate(df, cbl = .x)))))
        }
        
        ### Transform back into shape useful for data frame
        map(transpose(p), unlist)
      })
    
    preds <- transpose(preds)
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

### Function to obtain gamma + beta needed to compute influence function based estimators
gamma_beta <- function(y_vals, lambdas, shape_hat, gh_weights, v, rhos = NULL, lp2 = NULL, integral = 'gamma') {
  beta <- rep(NA, length(v))
  gamma <- rep(NA, length(v))
  
  
  for(i in 1:length(v)) {
    if(integral == 'gamma') {
      lambda_preds <- dgamma(v[i], shape = shape_hat, rate = shape_hat/lambdas)
    } else if(integral == 'normal') {
      ### If integral = normal, shape_hat := sigma_hat 
      lambda_preds <- dnorm(v[i], mean = lambdas, sd = shape_hat)
    }
    
    ### Just a Single missing confounder
    if(is.null(rhos)) {
      beta[i] <- sum(y_vals * lambda_preds * gh_weights)/sqrt(pi)
      gamma[i] <- sum(lambda_preds * gh_weights)/sqrt(pi)
    } else {
      rho_preds <- dbinom(lp2, size = 1, prob = rhos)
      ### 2nd missing confounder
      beta[i] <- sum(y_vals * lambda_preds * gh_weights * rho_preds)/sqrt(pi)
      gamma[i] <- sum(lambda_preds * gh_weights * rho_preds)/sqrt(pi)
    }
  }
  
  return(list('beta' = beta,
              'gamma' = gamma))
}

iptw <- function(gl_weights, taus, y, gb, gammas_betas, lambda_shape_hat, rhos = NULL, lp2 = NULL, integral = 'gamma') {
  if(is.null(rhos)) {
    numerator1 <- map(taus, ~{gl_weights * .x}) 
  } else {
    rho_preds <- map(rhos, ~dbinom(lp2, size = 1, prob = .x))
    numerator1 <- map2(taus, rho_preds, ~{gl_weights * .x * .y}) 
  }
  numerator2 <- map2(gb, y, ~{.y - .x})
  numerator <- map2(numerator1, numerator2, ~{.x *.y})
  frac <- map2(numerator, gammas_betas, ~{.x/(.y$gamma)})
  if(integral == 'gamma') {
    weight <- map_dbl(frac, sum)/gamma(lambda_shape_hat)
  } else {
    weight <- map_dbl(frac, sum)/sqrt(pi)
  }
  return(weight)
}

