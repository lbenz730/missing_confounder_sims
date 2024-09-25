library(tidyverse)
library(knitr)
library(kableExtra)
library(glue)
options(dplyr.summarise.inform = F)
options(knitr.kable.NA = '---')



get_sim_results <- function(sim_id) {
  ### Simulation Parameters
  params <- read_rds(glue('inputs/sim_params_{sim_id}.rds'))
  
  
  df_estimates <- 
    read_rds(glue('outputs/sim_results_{sim_id}.rds')) %>% 
    filter(abs(ate_hat) < 1) %>%  ### Remove Extreme Values 
    mutate('imputation_model' = case_when(is.na(imputation_model) ~ NA_character_,
                                          imputation_model == 'True Gamma' ~ 'True DGP',
                                          imputation_model == 'OLS' ~ 'OLS/LR',
                                          imputation_model == 'OLS w/ Interactions' ~ 'OLS/LR w/ Interactions',
                                          T ~ imputation_model)
           
           
    )
  
  ### Check if any Simulations Filtered Out
  if_rm <- 5000 - sum(df_estimates$method == 'Levis IF', na.rm = T)
  print(if_rm)
  
  params <- read_rds(glue('inputs/sim_params_{sim_id}.rds'))
  
  df_eval <- 
    df_estimates %>% 
    mutate('interactions' = ifelse(interactions, 'Pairwise', 'None')) %>% 
    group_by(outcome_model, interactions, imputation_model, ipw_model, method) %>% 
    summarise('mean_rel_bias' = round(100 * mean( (ate_hat - true_ate)/true_ate)),
              'median_rel_bias' = round(100 * (median(ate_hat) - true_ate[1])/true_ate[1]),
              'std_dev' = sd(ate_hat)) %>% 
    ungroup() %>% 
    mutate('rel_uncertainty' = std_dev/std_dev[which(method == 'Levis IF')]) %>% 
    mutate('group' = case_when(grepl('Levis', method) ~ '1 \\shortstack{CCMAR-\\\\based}',
                               outcome_model == 'OLS' ~ '2 \\shortstack{Outcome\\\\Regression}',
                               outcome_model == 'GAM' ~ '3 \\shortstack{Outcome\\\\Regression}',
                               outcome_model == 'Random Forest' ~ '4 \\shortstack{Outcome\\\\Regression}',
                               ipw_model == 'Logistic Regression' ~ '5 IPW',
                               ipw_model == 'GAM' ~ '6 IPW',
                               ipw_model == 'Random Forest' ~ '7 IPW')) %>% 
    mutate('ipw_model' = case_when(ipw_model == 'GAM' ~ '\\shortstack{GAM\\\\(Logit Link)}',
                                   T ~ ipw_model)) %>% 
    mutate('method' = gsub('Levis', '', method)) %>% 
    arrange(group, method) %>% 
    mutate(group = gsub('^\\d+\\s', '', group)) %>% 
    filter(!(replace(method, is.na(method), '---') == 'Horvitz-Thompson')) %>% 
    mutate('outcome_ipw' = case_when(!is.na(outcome_model) ~ outcome_model,
                                     !is.na(ipw_model) ~ ipw_model,
                                     T ~ method)) %>% 
    mutate('outcome_ipw' = case_when(outcome_ipw == 'CCMAR-based' ~ '\\shortstack{CCMAR-\\\\based}',
                                     outcome_ipw == 'Random Forest' ~ '\\shortstack{Random\\\\Forest}',
                                     outcome_ipw == 'Logistic Regression' ~ '\\shortstack{Logistic\\\\Regression}',
                                     T ~ outcome_ipw)) %>% 
    select(group, outcome_ipw, interactions, imputation_model, mean_rel_bias, median_rel_bias, std_dev, rel_uncertainty) %>% 
    mutate('sim_id' = sim_id)
  
  return(df_eval)
  
}

### Function To Make Latex Table
get_latex <- function(paper_ids, table_ids, caption = T) {
  col_lengths <- c(4, rep(4, length(paper_ids)))
  names(col_lengths) <- c('Model', paste('Scenario', paper_ids))
  
  if(length(paper_ids) == 1) {
    col_lengths <- 8
    names(col_lengths) <- paste('Data Driven Simulation Scenario', paper_ids)
  }
  
  df_tbl <- 
    df_results %>% 
    filter(sim_id %in% table_ids) %>% 
    pivot_wider(names_from = 'sim_id',
                values_from = c('mean_rel_bias', 'median_rel_bias', 'std_dev', 'rel_uncertainty'),
                values_fn = list, 
                names_vary = 'slowest')  %>% 
    unnest(cols = everything()) 
  names(df_tbl) <- c('  ', ' ', 'Interactions', 'Imputation', 
                     rep(c('\\% Bias', '\\% M-Bias', 'SE', 'RU'), length(paper_ids)))
  
  latex <- 
    df_tbl %>% 
    kbl(align = 'c', digits = 3, escape = F, format = 'latex') %>% 
    column_spec(1, border_left = T) %>%
    column_spec(ncol(df_tbl), border_right = T) %>%
    collapse_rows(columns = c(1:3), valign = 'middle') %>%
    add_header_above(col_lengths, bold = T) %>%
    row_spec(0, bold = T) %>% 
    footnote(general = c("IF = Influence Function; IWOR = Inverse-Weighted Outcome Regression",
                         '\\\\% Bias/\\\\% M-Bias = 100 $\\\\times$ bias of mean/median ATE estimate relative to true ATE',
                         "Relative Uncertainty (RU) = Ratio of standard error to standard error of CCMAR-based IF Estimator",
                         "Imputation: Linear Regression (OLS), Logistic Regression (LR), True Data Generating Process (DGP)"),
             general_title = '', 
             escape = F) %>% 
    gsub('\\multicolumn\\{1\\}\\{c\\|\\}', '\\multicolumn\\{1\\}\\{c\\}', .) %>% 
    gsub('\\multicolumn\\{4\\}\\{c}', '\\multicolumn\\{4\\}\\{c\\}', .) %>% 
    gsub('\\multicolumn\\{8\\}\\{c}', '\\multicolumn\\{8\\}\\{c\\}', .) %>% 
    gsub('tabular\\}\\[t', 'tabular\\}\\[ht', .)
  
  
  latex <- paste0("\\begin{table}\n\\centering\\footnotesize", latex, '\n\\end{table}')
  
  if(length(paper_ids) == 1) {
    latex <- 
      gsub('\\\\textbf\\{\\s+\\}\\s+\\&\\s+\\\\textbf\\{\\s+\\}',
           '\\\\multicolumn\\{2\\}\\{c\\}\\{\\\\textbf\\{Model\\}\\}', 
           latex)
    
    latex <- gsub('\\\\end\\{table\\}', '', latex)
    latex <- paste0(latex, '\\label{table:results_', paper_ids, '}\n\\end{table}')
    
    ### Add in caption
    if(caption) {
      cap <- read_lines('paper/captions.tex')[paper_ids] 
      latex <- unlist(strsplit(latex, '\\\\end\\{tabular\\}'))
      latex <- paste0(latex[1], paste0('\\end{tabular}', cap), latex[2])
    }
  }
  
  ### Changes for v2 Submission
  latex <- gsub('cline', 'cmidrule', latex)
  latex <- gsub('tabular', 'tabularx', latex)
  latex <-
    gsub('\\[ht\\]\\{\\|>\\{\\}c\\|c\\|c\\|c\\|c\\|c\\|c\\|>\\{\\}c\\|\\}',
         '\\{\\\\textwidth\\}\\{c\\@\\{\\}c\\@\\{\\}c\\@\\{\\}c\\@\\{\\}c\\@\\{\\}cc\\@\\{\\}c}',
         latex)
  
  return(latex)
}

### All Results
df_id <- 
  tibble('paper_id' = c(1:4, 4, 4,5:19),
         'table_id' = c(4, 6, 16, 19, 20, 21,  1:3, 5, 7:15, 17:18))
sim_ids <- 1:length(dir('outputs/'))
df_results <- map_dfr(sim_ids, get_sim_results)

### Paper simulations 1-3 (sim_ids 4, 6, 16) under Levis DGP
ltx_1 <- get_latex(paper_ids = 1, table_ids = 4)
ltx_2 <- get_latex(paper_ids = 2, table_ids = 6)
ltx_3 <- get_latex(paper_ids = 3, table_ids = 16)

### Special Format for Sim 4
df_tbl <- 
  df_results %>% 
  filter(sim_id %in% c(19:21)) %>% 
  select(-sim_id, -rel_uncertainty) %>% 
  distinct(.keep_all = T) %>% 
  mutate('row_id' = row_number()) %>% 
  filter(! (row_id > 42 & !grepl('Levis', outcome_ipw)) )  %>% 
  select(-row_id) %>% 
  slice(c(1:2, 43:46, 3:42)) %>% 
  mutate('rel_uncertainty' = std_dev/std_dev[1])

df_tbl$group[1:2] <- '\\shortstack{CCMAR-\\\\based}'
df_tbl$group[3:4] <- paste0('\\shortstack{Flexible Levis\\\\(Gamma $\\lambda_1$)}')
df_tbl$group[5:6] <- paste0('\\shortstack{Flexible Levis\\\\(Guassian $\\lambda_1$)}')

names(df_tbl) <- c('  ', ' ', 'Interactions', 'Imputation', 
                   rep(c('\\% Bias', '\\% M-Bias', 'SE', 'RU'), 1))

ltx_4 <- 
  df_tbl %>% 
  kbl(align = 'c', digits = 3, escape = F, format = 'latex') %>% 
  column_spec(1, border_left = T) %>%
  column_spec(ncol(df_tbl), border_right = T) %>%
  collapse_rows(columns = c(1:3), valign = 'middle') %>%
  add_header_above(c("Data Driven Simulation Scenario 4" = 8), bold = T) %>%
  row_spec(0, bold = T) %>% 
  footnote(general = c("IF = Influence Function; IWOR = Inverse-Weighted Outcome Regression",
                       '\\\\% Bias/\\\\% M-Bias = 100 $\\\\times$ bias of mean/median ATE estimate relative to true ATE',
                       "Relative Uncertainty (RU) = Ratio of standard error to standard error of CCMAR-based IF Estimator",
                       "Imputation: Linear Regression (OLS), Logistic Regression (LR), True Data Generating Process (DGP)"),
           general_title = '', 
           escape = F) %>% 
  gsub('\\multicolumn\\{1\\}\\{c\\|\\}', '\\multicolumn\\{1\\}\\{c\\}', .) %>% 
  gsub('\\multicolumn\\{4\\}\\{c}', '\\multicolumn\\{4\\}\\{c\\}', .) %>% 
  gsub('\\multicolumn\\{8\\}\\{c}', '\\multicolumn\\{8\\}\\{c\\}', .) %>% 
  gsub('tabular\\}\\[t', 'tabular\\}\\[ht', .) %>% 
  paste0("\\begin{table}\n\\centering\n\\footnotesize", ., '\n\\end{table}') %>% 
  gsub('\\\\textbf\\{\\s+\\}\\s+\\&\\s+\\\\textbf\\{\\s+\\}',
       '\\\\multicolumn\\{2\\}\\{c\\}\\{\\\\textbf\\{Model\\}\\}', 
       .) %>% 
  paste0(., '\\label{table:results_', 4, '}') %>% 
  gsub('cline', 'cmidrule', .) %>% 
  gsub('tabular', 'tabularx', .) %>% 
  gsub('\\[ht\\]\\{\\|>\\{\\}c\\|c\\|c\\|c\\|c\\|c\\|c\\|>\\{\\}c\\|\\}',
       '\\{\\\\textwidth\\}\\{c\\@\\{\\}c\\@\\{\\}c\\@\\{\\}c\\@\\{\\}c\\@\\{\\}cc\\@\\{\\}c}',
       .)

write(paste(ltx_1, ltx_2, ltx_3, ltx_4, sep = '\n\n\n\n'), 'paper/results.tex')

### Supplement Tables
supp_id <- 
  df_id %>% 
  filter(paper_id > 4)
ltx_sup <- 
  paste(map2_chr(supp_id$paper_id, supp_id$table_id, ~get_latex(.x, .y)), collapse = '\n\n') %>% 
  gsub('\\|c\\|', 'c', .) %>% 
  gsub('\\centering\\arraybackslash CCMAR-based', '\\centering\\arraybackslash \\shortstack{CCMAR-\\based}', .)
  
write(ltx_sup, 'paper/supplement_results.tex')


### Fully Flexible Table
df_ff <- read_csv('fully_flexible/results/levis.csv')
ATE.true <-  0.2556345

df_tbl <- 
  df_ff %>% 
  group_by(method, model_type) %>% 
  summarise('mean_rel_bias' = round(100 * mean( (ate_hat - ATE.true)/ATE.true), 1),
            'median_rel_bias' = round(100 * (median(ate_hat) - ATE.true)/ATE.true, 1),
            'sd_ate' = sd(ate_hat)) %>%
  ungroup() %>% 
  mutate('rel_efficiency' = sd_ate/sd_ate[model_type == 'Parametric' & method == 'Levis IF']) %>% 
  slice(c(2,1,4,3)) %>% 
  mutate('method' = gsub('Levis', 'CCMAR', method))

names(df_tbl) <- c('Estimator', 'Nusiance Models', 
                   rep(c('\\% Bias', '\\% M-Bias', 'SE', 'Relative Uncertainty'), 1))

ff_ltx <- 
  df_tbl %>% 
  kbl(align = 'c', digits = 3, escape = F, format = 'latex') %>% 
  column_spec(1, border_left = T) %>%
  column_spec(ncol(df_tbl), border_right = T) %>%
  collapse_rows(columns = c(1:2), valign = 'middle') %>% 
  add_header_above(c('Non-Parametric Nuisance Model Simulation' = 6), bold = T) %>%
  row_spec(0, bold = T) %>% 
  footnote(general = c("IF = Influence Function; IWOR = Inverse-Weighted Outcome Regression",
                       '\\\\% Bias = 100 $\\\\times$ Relative Bias of Mean Average Treatment Effect Estimate',
                       '\\\\% M-Bias = 100 $\\\\times$ Relative Bias of Median Average Treatment Effect Estimate',
                       "Relative Uncertainty = Ratio of standard error to standard error of Parametric CCMAR-based IF Estimator"),
           general_title = '', 
           escape = F) %>% 
  gsub('\\multicolumn\\{1\\}\\{c\\|\\}', '\\multicolumn\\{1\\}\\{\\|c\\|\\}', .) %>% 
  gsub('\\multicolumn\\{4\\}\\{c}', '\\multicolumn\\{4\\}\\{\\|c\\|\\}', .) %>% 
  gsub('\\multicolumn\\{6\\}\\{c}', '\\multicolumn\\{6\\}\\{\\|c\\|\\}', .) %>% 
  gsub('tabular\\}\\[t', 'tabular\\}\\[ht', .) %>% 
  paste0("\\begin{table}\n\\centering", ., '\n\\end{table}') %>% 
  paste0(., '\\label{table:results_ff}')

write(ff_ltx, 'paper/fully_flexible_results.tex')
