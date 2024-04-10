library(tidyverse)
library(glue)
library(mgcv)
library(ranger)
library(glmnet)
library(knitr)
library(kableExtra)

source('generate_data.R')

### Sort Interactions so they are sorted
sort_term <- function(term) {
  paste(sort(unlist(strsplit(term, ':'))), collapse = ':')
}

get_model_coeff <- function(sim_id) {
  ### Load Simulation Parameters
  params <- read_rds(glue('inputs/sim_params_{sim_id}.rds'))
  
  ### Load Models
  load_models(params)
  
  ### Model Coefficients
  df_terms <- 
    map2_dfr(list(pi_star, mu_star, eta_star, lambda_star), 
             c('pi', 'mu', 'eta', 'lambda'), ~{
               tibble('term' = names(.x$coefficients),
                      'coeff' = .x$coefficients,
                      'model' = .y)
             })
  
  ### sigma_hat and lambda_hat
  df_var <- 
    tibble('term' = c('sigma_hat', 'lambda_shape_hat'),
           'coeff' = c(sigma_star, lambda_shape),
           'model' = c('mu', 'lambda'))
  df_terms <- bind_rows(df_terms, df_var)
  
  ### Lp2 model if multiple_lp = T
  if('rho' %in% names(params$models)) {
    df_rho <- 
      tibble('term' = names(unlist(rho_star)),
             'coeff' = unlist(rho_star),
             'model' = 'rho') %>% 
      mutate('term' = gsub('_', ':', gsub('beta_', '', term))) %>% 
      mutate('term' = gsub('hispanic:raceeth', 'raceeth_hispanic', term)) %>% 
      mutate('term' = gsub('raceeth:hispanic', 'raceeth_hispanic', term))
    df_terms <- bind_rows(df_terms, df_rho)
  }
  
  ### Coefficient Amplifications
  df_shifts <- 
    map2_dfr(params$shifts, names(params$shifts), ~{
      tibble('term' = names(unlist(.x)),
             'coeff' = unlist(.x),
             'model' = .y)
    })
  
  
  ### Net Coefficient will be baseline coeff + shift
  df_terms <- 
    df_terms %>% 
    bind_rows(df_shifts) %>% 
    mutate('term' = map_chr(term, sort_term)) %>% 
    group_by(model, term) %>% 
    summarise('coeff' = sum(coeff)) %>% 
    ungroup() %>% 
    mutate('sim_id' = sim_id)
  
  return(df_terms)
  
}

### Parameters for all simulations we've run
sim_ids <- 1:length(dir('outputs/'))
df_coef <- 
  map_dfr(sim_ids, get_model_coeff)

options(knitr.kable.NA = '---')
df <- 
  df_coef %>% 
  pivot_wider(names_from = 'sim_id',
              values_from = 'coeff') %>% 
  mutate('term' = case_when(term == 'prevalence' ~ '(Intercept)', 
                            term == 'GENDER' ~ '$L^{(1)}$',
                            term == 'bbl' ~ '$L^{(2)}$',
                            term == 'raceeth_hispanic' ~ '$L^{(3)}$',
                            term == 'cbl' ~ '$L^{(4)}$',
                            term == 'smoker' ~ '$L^{(5)}$',
                            term == 'A' ~ '$A$',
                            term == 'ptwc' ~ '$Y$',
                            term == 'A:ptwc' ~ '$A\\times Y$',
                            term == 'bbl:GENDER' ~ '$L^{(1)} \\times L^{(2)}$',
                            term == 'bbl:raceeth_hispanic' ~ '$L^{(2)} \\times L^{(3)}$',
                            term == 'cbl:raceeth_hispanic' ~ '$L^{(3)} \\times L^{(4)}$',
                            term == 'GENDER:raceeth_hispanic' ~ '$L^{(1)} \\times L^{(3)}$',
                            term == 'I(bbl^2)' ~ '$(L^{(2)})^2$',
                            term == 'sigma_hat' ~ '$\\sigma$',
                            term == 'lambda_shape_hat' ~'$\\alpha$',
                            term == 'bbl:ptwc' ~ '$Y \\times L^{(2)}$',
                            term == 'A:GENDER' ~ '$A \\times L^{(1)}$',
                            term == 'A:bbl' ~ '$A \\times L^{(2)}$',
                            term == 'A:raceeth_hispanic' ~ '$A \\times L^{(3)}$',
                            term == 'A:cbl' ~ '$A \\times L^{(4)}$',
                            term == 'A:smoker' ~ '$A \\times L^{(5)}$',
                            term == 'cbl:ptwc' ~ '$Y \\times L^{(4)}$',
                            term == 'GENDER:ptwc' ~ '$Y \\times L^{(1)}$',
                            term == 'cbl:smoker' ~ '$L^{(4)} \\times L^{(5)}$',
                            T ~ term)) %>% 
  mutate('term' = fct_relevel(term, '(Intercept)', '$L^{(1)}$', '$L^{(2)}$',
                              '$L^{(3)}$', '$L^{(4)}$', '$L^{(5)}$',  '$(L^{(2)})^2$',
                              '$L^{(1)} \\times L^{(2)}$', '$L^{(2)} \\times L^{(3)}$',
                              '$L^{(1)} \\times L^{(3)}$', '$L^{(4)} \\times L^{(5)}$', 
                              '$A$', '$Y$', '$A\\times Y$', '$A \\times L^{(1)}$',
                              '$A \\times L^{(2)}$', '$A \\times L^{(3)}$', 
                              '$A \\times L^{(4)}$', '$A \\times L^{(5)}$', 
                              '$Y \\times L^{(1)}$', '$Y \\times L^{(2)}$', 
                              '$Y \\times L^{(4)}$')) %>% 
  mutate('model' = fct_relevel(model, 'eta', 'mu', 'pi', 'lambda', 'rho')) %>% 
  arrange(model, term) %>% 
  mutate('model' = case_when(model == 'lambda' ~ '$\\lambda_1~\\text{or}~\\tilde\\lambda_1$',
                             model == 'rho' ~ '$\\lambda_2~\\text{or}~\\tilde\\lambda_2$',
                             model != 'pi' ~ paste0('$\\', model, '~\\text{or}~\\tilde\\', model, '$'),
                             T ~ paste0('$\\', model, '$'))) %>% 
  mutate_at(vars(everything()), ~replace(.x, .x == 0, NA)) %>% 
  rename('Model' = model,
         'Term' = term)


### For Paper (select in order)
t1_1 <- 
  df %>%
  filter(!grepl('lambda', Model)) %>% 
  select(Model, Term, '1' = '4', '2' = '6', '3' = '16', '4' = '19') %>% 
  kbl(align = 'c', digits = 3, escape = F, format = 'latex') %>% 
  column_spec(1, border_left = T) %>% 
  column_spec(6, border_right = T) %>% 
  collapse_rows(columns = 1, valign = 'middle') %>% 
  add_header_above(c(" " = 1, " " = 1, 'Data Driven Simulation Scenario' = 4), bold = T) %>% 
  row_spec(0, bold = T) %>% 
  gsub('0.000', '< 0.001', .) %>% 
  gsub('\\multicolumn\\{1\\}\\{c\\|\\}', '\\multicolumn\\{1\\}\\{\\|c\\|\\}', .) %>% 
  gsub('\\multicolumn\\{4\\}\\{c}', '\\multicolumn\\{4\\}\\{\\|c\\|\\}', .)


t1_2 <- 
  df %>%
  filter(grepl('lambda', Model)) %>% 
  select(Model, Term, '1' = '4', '2' = '6', '3' = '16', '4' = '19') %>% 
  kbl(align = 'c', digits = 3, escape = F, format = 'latex') %>% 
  column_spec(1, border_left = T) %>% 
  column_spec(6, border_right = T) %>% 
  collapse_rows(columns = 1, valign = 'middle') %>% 
  add_header_above(c(" " = 1, " " = 1, 'Data Driven Simulation Scenario' = 4), bold = T) %>% 
  row_spec(0, bold = T) %>% 
  gsub('0.000', '< 0.001', .) %>% 
  gsub('\\multicolumn\\{1\\}\\{c\\|\\}', '\\multicolumn\\{1\\}\\{\\|c\\|\\}', .) %>% 
  gsub('\\multicolumn\\{4\\}\\{c}', '\\multicolumn\\{4\\}\\{\\|c\\|\\}', .)



paper_table <- 
  paste(
    "\\begin{table}
\\begin{minipage}[t]{0.5\\textwidth}
\\centering",
    t1_1, 
    '\\end{minipage}
\\begin{minipage}[t]{0.5\\textwidth}',
    t1_2,
    '\\end{minipage}
\\end{table}',
    sep = '\n'
  )
write(paper_table, 'paper/model_param.tex')

### For Supplement
tibble('Paper Simulation' = c(1:4, 4, 4,5:19),
       'sim_id' = c(4, 6, 16, 19, 20, 21,  1:3, 5, 7:15, 17:18)) %>% 
  knitr::kable(align = 'c', format = 'markdown')

t2_1 <- 
  df %>%
  filter(!grepl('lambda', Model)) %>% 
  select(Model, Term, c('1', '2', '3', '5')) %>% 
  set_names(c('Model', 'Term', as.character(5:8))) %>% 
  kbl(align = 'c', digits = 3, escape = F, format = 'latex') %>% 
  column_spec(1, border_left = T) %>% 
  column_spec(6, border_right = T) %>% 
  collapse_rows(columns = 1, valign = 'middle') %>% 
  add_header_above(c(" " = 1, " " = 1, 'Data Driven Simulation Scenario' = 4), bold = T) %>% 
  row_spec(0, bold = T) %>% 
  gsub('0.000', '< 0.001', .) %>% 
  gsub('\\multicolumn\\{1\\}\\{c\\|\\}', '\\multicolumn\\{1\\}\\{\\|c\\|\\}', .) %>% 
  gsub('\\multicolumn\\{4\\}\\{c}', '\\multicolumn\\{4\\}\\{\\|c\\|\\}', .)


t2_2 <- 
  df %>%
  filter(grepl('lambda', Model)) %>% 
  select(Model, Term, c('1', '2', '3', '5')) %>% 
  set_names(c('Model', 'Term', as.character(5:8))) %>% 
  kbl(align = 'c', digits = 3, escape = F, format = 'latex') %>% 
  column_spec(1, border_left = T) %>% 
  column_spec(6, border_right = T) %>% 
  collapse_rows(columns = 1, valign = 'middle') %>% 
  add_header_above(c(" " = 1, " " = 1, 'Data Driven Simulation Scenario' = 4), bold = T) %>% 
  row_spec(0, bold = T) %>% 
  gsub('0.000', '< 0.001', .) %>% 
  gsub('\\multicolumn\\{1\\}\\{c\\|\\}', '\\multicolumn\\{1\\}\\{\\|c\\|\\}', .) %>% 
  gsub('\\multicolumn\\{4\\}\\{c}', '\\multicolumn\\{4\\}\\{\\|c\\|\\}', .)



sup_table_1 <- 
  paste(
    "\\begin{table}
\\begin{minipage}[t]{0.5\\textwidth}
\\centering",
    t2_1, 
    '\\end{minipage}
\\begin{minipage}[t]{0.5\\textwidth}',
    t2_2,
    '\\end{minipage}
\\end{table}',
    sep = '\n'
  )


t3_1 <- 
  df %>%
  filter(!grepl('lambda', Model)) %>% 
  select(Model, Term, c('7', '8', '9', '10')) %>% 
  set_names(c('Model', 'Term', as.character(9:12))) %>% 
  kbl(align = 'c', digits = 3, escape = F, format = 'latex') %>% 
  column_spec(1, border_left = T) %>% 
  column_spec(6, border_right = T) %>% 
  collapse_rows(columns = 1, valign = 'middle') %>% 
  add_header_above(c(" " = 1, " " = 1, 'Data Driven Simulation Scenario' = 4), bold = T) %>% 
  row_spec(0, bold = T) %>% 
  gsub('0.000', '< 0.001', .) %>% 
  gsub('\\multicolumn\\{1\\}\\{c\\|\\}', '\\multicolumn\\{1\\}\\{\\|c\\|\\}', .) %>% 
  gsub('\\multicolumn\\{4\\}\\{c}', '\\multicolumn\\{4\\}\\{\\|c\\|\\}', .)


t3_2 <- 
  df %>%
  filter(grepl('lambda', Model)) %>% 
  select(Model, Term, c('7', '8', '9', '10')) %>% 
  set_names(c('Model', 'Term', as.character(9:12))) %>% 
  kbl(align = 'c', digits = 3, escape = F, format = 'latex') %>% 
  column_spec(1, border_left = T) %>% 
  column_spec(6, border_right = T) %>% 
  collapse_rows(columns = 1, valign = 'middle') %>% 
  add_header_above(c(" " = 1, " " = 1, 'Data Driven Simulation Scenario' = 4), bold = T) %>% 
  row_spec(0, bold = T) %>% 
  gsub('0.000', '< 0.001', .) %>% 
  gsub('\\multicolumn\\{1\\}\\{c\\|\\}', '\\multicolumn\\{1\\}\\{\\|c\\|\\}', .) %>% 
  gsub('\\multicolumn\\{4\\}\\{c}', '\\multicolumn\\{4\\}\\{\\|c\\|\\}', .)



sup_table_2 <- 
  paste(
    "\\begin{table}
\\begin{minipage}[t]{0.5\\textwidth}
\\centering",
    t3_1, 
    '\\end{minipage}
\\begin{minipage}[t]{0.5\\textwidth}',
    t3_2,
    '\\end{minipage}
\\end{table}',
    sep = '\n'
  )




t4_1 <- 
  df %>%
  filter(!grepl('lambda', Model)) %>% 
  select(Model, Term, c('11', '12', '13', '14')) %>% 
  set_names(c('Model', 'Term', as.character(13:16))) %>% 
  kbl(align = 'c', digits = 3, escape = F, format = 'latex') %>% 
  column_spec(1, border_left = T) %>% 
  column_spec(6, border_right = T) %>% 
  collapse_rows(columns = 1, valign = 'middle') %>% 
  add_header_above(c(" " = 1, " " = 1, 'Data Driven Simulation Scenario' = 4), bold = T) %>% 
  row_spec(0, bold = T) %>% 
  gsub('0.000', '< 0.001', .) %>% 
  gsub('\\multicolumn\\{1\\}\\{c\\|\\}', '\\multicolumn\\{1\\}\\{\\|c\\|\\}', .) %>% 
  gsub('\\multicolumn\\{4\\}\\{c}', '\\multicolumn\\{4\\}\\{\\|c\\|\\}', .)


t4_2 <- 
  df %>%
  filter(grepl('lambda', Model)) %>% 
  select(Model, Term, c('11', '12', '13', '14')) %>% 
  set_names(c('Model', 'Term', as.character(13:16))) %>% 
  kbl(align = 'c', digits = 3, escape = F, format = 'latex') %>% 
  column_spec(1, border_left = T) %>% 
  column_spec(6, border_right = T) %>% 
  collapse_rows(columns = 1, valign = 'middle') %>% 
  add_header_above(c(" " = 1, " " = 1, 'Data Driven Simulation Scenario' = 4), bold = T) %>% 
  row_spec(0, bold = T) %>% 
  gsub('0.000', '< 0.001', .) %>% 
  gsub('\\multicolumn\\{1\\}\\{c\\|\\}', '\\multicolumn\\{1\\}\\{\\|c\\|\\}', .) %>% 
  gsub('\\multicolumn\\{4\\}\\{c}', '\\multicolumn\\{4\\}\\{\\|c\\|\\}', .)



sup_table_3 <- 
  paste(
    "\\begin{table}
\\begin{minipage}[t]{0.5\\textwidth}
\\centering",
    t4_1, 
    '\\end{minipage}
\\begin{minipage}[t]{0.5\\textwidth}',
    t4_2,
    '\\end{minipage}
\\end{table}',
    sep = '\n'
  )



t5_1 <- 
  df %>%
  filter(!grepl('lambda', Model)) %>% 
  select(Model, Term, c('15', '17', '18')) %>% 
  set_names(c('Model', 'Term', as.character(17:19))) %>% 
  kbl(align = 'c', digits = 3, escape = F, format = 'latex') %>% 
  column_spec(1, border_left = T) %>% 
  column_spec(5, border_right = T) %>% 
  collapse_rows(columns = 1, valign = 'middle') %>% 
  add_header_above(c(" " = 1, " " = 1, 'Data Driven Simulation Scenario' = 3), bold = T) %>% 
  row_spec(0, bold = T) %>% 
  gsub('0.000', '< 0.001', .) %>% 
  gsub('\\multicolumn\\{1\\}\\{c\\|\\}', '\\multicolumn\\{1\\}\\{\\|c\\|\\}', .) %>% 
  gsub('\\multicolumn\\{3\\}\\{c}', '\\multicolumn\\{3\\}\\{\\|c\\|\\}', .)


t5_2 <- 
  df %>%
  filter(grepl('lambda', Model)) %>% 
  select(Model, Term, c('15', '17', '18')) %>% 
  set_names(c('Model', 'Term', as.character(17:19))) %>% 
  kbl(align = 'c', digits = 3, escape = F, format = 'latex') %>% 
  column_spec(1, border_left = T) %>% 
  column_spec(5, border_right = T) %>% 
  collapse_rows(columns = 1, valign = 'middle') %>% 
  add_header_above(c(" " = 1, " " = 1, 'Data Driven Simulation Scenario' = 3), bold = T) %>% 
  row_spec(0, bold = T) %>% 
  gsub('0.000', '< 0.001', .) %>% 
  gsub('\\multicolumn\\{1\\}\\{c\\|\\}', '\\multicolumn\\{1\\}\\{\\|c\\|\\}', .) %>% 
  gsub('\\multicolumn\\{3\\}\\{c}', '\\multicolumn\\{3\\}\\{\\|c\\|\\}', .)



sup_table_4 <- 
  paste(
    "\\begin{table}
\\begin{minipage}[t]{0.5\\textwidth}
\\centering",
    t5_1, 
    '\\end{minipage}
\\begin{minipage}[t]{0.5\\textwidth}',
    t5_2,
    '\\end{minipage}
\\end{table}',
    sep = '\n'
  )


write(paste(sup_table_1, 
            sup_table_2,
            sup_table_3,
            sup_table_4,
            sep = '\n\n\n\n'),
      'paper/supplement_model_param.tex')
