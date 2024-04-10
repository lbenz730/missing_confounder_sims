library(tidyverse)
library(glue)
library(gt)
library(magick)

### Custom ggplot theme
theme_set(theme_bw() +
            theme(plot.title = element_text(hjust = 0.5, size = 24),
                  plot.subtitle = element_text(hjust = 0.5, size = 18),
                  axis.title = element_text(size = 20),
                  strip.text = element_text(size = 14),
                  plot.caption = element_text(size = 10),
                  legend.position = "none"))


sim_ids <- 1:length(dir('outputs/'))
for(sim_id in sim_ids) {
  if(!dir.exists(glue('figures/sim_{sim_id}/'))) {
    dir.create(glue('figures/sim_{sim_id}/'))
  }
  
  df_estimates <- 
    read_rds(glue('outputs/sim_results_{sim_id}.rds')) %>% 
    filter(abs(ate_hat) < 1) %>%  ### Remove Extreme Values 
    mutate('imputation_model' = case_when(is.na(imputation_model) ~ NA_character_,
                                          imputation_model == 'True Gamma' ~ 'True DGP',
                                          imputation_model == 'OLS' ~ 'OLS/LR',
                                          imputation_model == 'OLS w/ Interactions' ~ 'OLS/LR w/ Interactions',
                                          T ~ imputation_model)
           
           
    )
  params <- read_rds(glue('inputs/sim_params_{sim_id}.rds'))
  
  df_eval <- 
    df_estimates %>% 
    mutate('interactions' = ifelse(interactions, 'Pairwise', 'None')) %>% 
    group_by(outcome_model, interactions, imputation_model, ipw_model, method) %>% 
    summarise('mean_ate' = mean(ate_hat),
              'median_ate' = median(ate_hat),
              'mean_rel_bias' = mean( (ate_hat - true_ate)/true_ate),
              'mean_abs_bias' = mean(abs(ate_hat - true_ate)),
              'variance' = var(ate_hat),
              'std_dev' = sd(ate_hat),
              'mad' = mad(ate_hat)) %>% 
    ungroup() %>% 
    mutate('rel_uncertainty' = std_dev/std_dev[which(method == 'Levis IF')]) %>% 
    mutate('group' = case_when(grepl('Levis', method) ~ '1 Levis Estimators',
                               outcome_model == 'OLS' ~ '2 Outcome Regression (OLS)',
                               outcome_model == 'GAM' ~ '3 Outcome Regression (GAM)',
                               outcome_model == 'Random Forest' ~ '4 Outcome Regression (Random Forest)',
                               ipw_model == 'Logistic Regression' ~ '5 IPW (Simple Logistic Regression)',
                               ipw_model == 'GAM' ~ '6 IPW (GAM Logsitic Regression)',
                               ipw_model == 'Random Forest' ~ '7 IPW (Random Forest Classifier)')) %>% 
    arrange(group, method) %>% 
    mutate(group = gsub('^\\d+\\s', '', group)) %>% 
    select(method, interactions, imputation_model, everything(), -ipw_model, -outcome_model) %>% 
    filter(!(replace(method, is.na(method), '---') == 'Horvitz-Thompson')) %>%
    group_by(group)
  
  
  ### Table
  gt_analysis <-
    df_eval %>%
    gt() %>% 
    
    ### Align Columns
    cols_align(align = "center", columns = everything()) %>%
    
    ### Round Numbers
    fmt_number(columns = c('mean_ate', 'median_ate', 'rel_uncertainty'), decimals = 3, sep_mark = '') %>% 
    fmt_percent(columns = c('mean_rel_bias'), decimals = 1, sep_mark = '') %>% 
    fmt_number(columns = c('mean_abs_bias'), decimals = 3, sep_mark = '') %>% 
    fmt_scientific(columns = c('variance', 'std_dev', 'mad'), sep_mark = '') %>% 
    sub_missing(columns = c('imputation_model'), missing_text = "---") %>%
    sub_missing(columns = c('method', 'interactions'), missing_text = "") %>%
    
    ### Column Headers
    tab_spanner(label = 'ATE Estimates', columns = c('mean_ate', 'median_ate')) %>%
    tab_spanner(label = 'Bias', columns = contains('bias')) %>%
    tab_spanner(label = 'Efficiency', columns = c('variance', 'std_dev', 'mad', 'rel_uncertainty')) %>% 
    
    ### Color Columns 
    data_color(columns = c(mean_ate, median_ate),
               colors = scales::col_numeric(palette = ggsci::rgb_material('amber', n = 100),
                                            domain = range(c(df_eval$mean_ate, 
                                                             df_eval$median_ate)))) %>%
    data_color(columns = c(mean_abs_bias),
               colors = scales::col_numeric(palette = ggsci::rgb_material('amber', n = 100),
                                            domain = c(0, max(df_eval$mean_abs_bias)))) %>%
    data_color(columns = c(mean_rel_bias),
               colors = scales::col_numeric(palette = ggsci::rgb_gsea(n = 100),
                                            domain = c(-1, 1) * max(abs(df_eval$mean_rel_bias)))) %>%
    data_color(columns = c(rel_uncertainty),
               colors = scales::col_bin(palette = ggsci::rgb_gsea(n = 100),
                                        bins = unique(quantile(df_eval$rel_uncertainty, c(0, seq(0.05, 0.95, 0.01), 1))),
                                        domain = c(0, max(abs(df_eval$rel_uncertainty))))) %>%
    data_color(columns = c(variance),
               colors = scales::col_bin(palette = ggsci::rgb_material('amber', n = 100),
                                        bins = unique(quantile(df_eval$variance, c(0, seq(0.05, 0.95, 0.01), 1))),
                                        domain = c(0, max(df_eval$variance)))) %>%
    data_color(columns = c(std_dev),
               colors = scales::col_bin(palette = ggsci::rgb_material('amber', n = 100),
                                        bins = unique(quantile(df_eval$std_dev, c(0, seq(0.05, 0.95, 0.01), 1))),
                                        domain = c(0, max(df_eval$std_dev)))) %>%
    data_color(columns = c(mad),
               colors = scales::col_bin(palette = ggsci::rgb_material('amber', n = 100),
                                        bins = unique(quantile(df_eval$mad, c(0, seq(0.05, 0.95, 0.01), 1))),
                                        domain = c(0, max(df_eval$mad)))) %>% 
    ### Grid Lines
    tab_style(style = list(cell_borders(sides = "right", color = "black", weight = px(4))),
              locations = list(cells_body(columns = c('imputation_model', 'median_ate',
                                                      'mean_abs_bias', 'mad')))) %>%
    
    tab_style(style = list(cell_borders(sides = "top", color = "black", weight = px(2))),
              locations = list(cells_body(rows = interactions == 'Pairwise' & imputation_model == 'OLS/LR'))) %>%
    
    ### Names
    cols_label('method' = '',
               'interactions' = 'Interactions',
               'imputation_model' = 'Imputation Model',
               'mean_ate' = 'Mean',
               'median_ate' = 'Median',
               'mean_rel_bias' = 'Relative Bias',
               'mean_abs_bias' = 'MAE',
               'variance' = 'Variance',
               'std_dev' = 'SD',
               'mad' = 'MAD',
               'rel_uncertainty' = 'Relative Uncertainty') %>% 
    tab_header(title = md(paste('__Simulation #', sim_id,  'Analysis__')),
               subtitle = md(paste(c(paste0('__Description: ', params$description, '__'),
                                     paste0('__\\# of Simulations: ', df_estimates$n_sims[1], '__'),
                                     paste0('__\\# of Subjects: ', df_estimates$n_patients[1], '__'),
                                     paste0('__DGP: ', params$data_generating_process, '__'),
                                     paste0('__\\# of Missing Confounders: ', ifelse(params$multiple_lp, 2, 1), '__'),
                                     paste0('__True Average Treatment Effect: ', sprintf('%0.3f', df_estimates$true_ate[1]), '__')),
                                   collapse  = '<br>'))) %>%
    tab_options(column_labels.font.size = 16,
                heading.title.font.size = 40,
                heading.subtitle.font.size = 20,
                column_labels.font.weight = 'bold',
                row_group.font.weight = 'bold',
                row_group.font.size  = 22)
  
  ### Save Table
  gtExtras::gtsave_extra(gt_analysis , glue('figures/sim_{sim_id}/table_sim_{sim_id}.png'), vwidth = 3000, selector = 'table')
  
  df_estimates %>% 
    filter((method != 'Horvitz-Thompson') | is.na(method)) %>% 
    mutate('label' = case_when(outcome_model == 'OLS' ~ 'Outcome Regression (OLS)',
                               outcome_model == 'GAM' ~ 'Outcome Regression (GAM)',
                               outcome_model == 'Random Forest' ~ 'Outcome Regression (Random Forest)',
                               method == 'Levis IF' ~ 'Levis IF',
                               method == 'Levis IWOR' ~ 'Levis IWOR',
                               method == 'Hajek' & ipw_model == 'Logistic Regression' ~ 'IPW (Simple Logistic Regression)',
                               method == 'Hajek' & ipw_model == 'GAM' ~ 'IPW (GAM Logistic Regression)',
                               method == 'Hajek' & ipw_model == 'Random Forest' ~ 'IPW (Random Forest Classifier)' )) %>% 
    mutate('imputation_model' = ifelse(is.na(imputation_model), 'No Imputation', imputation_model)) %>% 
    mutate('interactions' = case_when(interactions == F ~ 'No interactions',
                                      interactions == T ~ 'Pairwise Interactions',
                                      method == 'Levis' ~ 'Correct Interactions Specified',
                                      T ~ '--')) %>% 
    group_by(label, imputation_model, interactions) %>% 
    filter(ate_hat >= quantile(ate_hat, 0.05), ate_hat <= quantile(ate_hat, 0.95)) %>%
    ungroup() %>% 
    ggplot(aes(x = ate_hat)) + 
    facet_wrap(label ~ interactions) + 
    geom_density(aes(fill = imputation_model), alpha = 0.2) + 
    geom_vline(xintercept = df_estimates$true_ate[1], lty = 2)  + 
    theme(legend.position = 'bottom') + 
    labs(x = 'Estimated ATE',
         y = 'Density',
         fill = 'Imputation Method', 
         title = 'Distribution of ATE\nInner 90% of Distribution ([5%, 95%] Quantiles)',
         subtitle = glue('Simulation #{sim_id}\n# of Subjects = {df_estimates$n_patients[1]}'))
  
  ### Save Figure
  ggsave(glue('figures/sim_{sim_id}/sim_{sim_id}_ate.png'), height = 9, width = 16)
  
  df_levis_graph <- 
    df_estimates %>% 
    filter(str_detect(method, 'Levis')) %>% 
    group_by(method) %>%
    filter(ate_hat >= quantile(ate_hat, 0.05), ate_hat <= quantile(ate_hat, 0.95)) %>%
    ungroup() %>%
    select(method, contains('hat'), contains('true')) %>% 
    pivot_longer(cols = contains('hat'),
                 names_to = 'contrast',
                 values_to = 'estimate_hat') %>% 
    pivot_longer(cols = contains('true'),
                 names_to = 'contrast_true',
                 values_to = 'truth') %>% 
    mutate('contrast' = gsub('_hat', '', contrast),
           'contrast_true' = gsub('true_', '', contrast_true)) %>% 
    filter(contrast == contrast_true) %>% 
    mutate(contrast = toupper(contrast))
  
  ggplot(df_levis_graph, aes(x = estimate_hat)) + 
    facet_grid(method ~ contrast) + 
    geom_vline(aes(xintercept = truth), lty = 2, lwd = 0.8, alpha = 0.8) + 
    geom_density(aes(fill = method), alpha = 0.2) + 
    
    labs(x = 'Estimated Contrast',
         y = 'Density',
         title = 'Distribution of Treatment Effects',
         subtitle = glue('Simulation #{sim_id}\n# of Subjects = {df_estimates$n_patients[1]}'))
  
  ### Save Figure
  ggsave(glue('figures/sim_{sim_id}/levis_sim_{sim_id}_ate.png'), height = 9, width = 16)
  
  ### Levis Scatterplot
  df_estimates %>% 
    filter(grepl('Levis', method)) %>% 
    ggplot(aes(x = te0_hat, y = te1_hat)) + 
    facet_wrap(~method) + 
    geom_point(alpha = 0.20) + 
    labs(x = 'TE0 Hat',
         y = 'TE1 Hat',
         title = 'Levis Estimators Scatterplot',
         subtitle = glue('Simulation #{sim_id}\n# of Subjects = {df_estimates$n_patients[1]}'))
  ggsave(glue('figures/sim_{sim_id}/levis_sim_{sim_id}_scatter.png'), height = 9, width = 16)
  
  
  ### Ridgeline Plot
  df_graph <- 
    df_estimates %>% 
    filter((method != 'Horvitz-Thompson') | is.na(method)) %>% 
    mutate('label' = case_when(outcome_model == 'OLS' ~ 'Outcome Regression',
                               outcome_model == 'GAM' ~ 'Outcome Regression',
                               outcome_model == 'Random Forest' ~ 'Outcome Regression',
                               method == 'Levis IF' ~ 'Levis',
                               method == 'Levis IWOR' ~ 'Levis',
                               method == 'Hajek' & ipw_model == 'Logistic Regression' ~ 'IPW',
                               method == 'Hajek' & ipw_model == 'GAM' ~ 'IPW',
                               method == 'Hajek' & ipw_model == 'Random Forest' ~ 'IPW' )) %>% 
    mutate('imputation_model' = ifelse(is.na(imputation_model), 'No Imputation', imputation_model)) %>% 
    mutate('interactions' = case_when(interactions == F ~ 'No interactions',
                                      interactions == T ~ 'Pairwise Interactions',
                                      method == 'Levis' ~ 'Correct Interactions Specified',
                                      T ~ '--')) %>% 
    mutate('y_lab' = case_when(label == 'Outcome Regression' & is.na(interactions) ~ outcome_model,
                               label == 'IPW' & is.na(interactions) ~ ipw_model,
                               label == 'Outcome Regression' & interactions == 'Pairwise Interactions' ~ paste(outcome_model, '(Interactions)'),
                               label == 'IPW' & interactions == 'Pairwise Interactions' ~ paste(ipw_model, '(Interactions)'),
                               label == 'Outcome Regression' & interactions == 'No interactions' ~ outcome_model,
                               label == 'IPW' & interactions == 'No interactions' ~ ipw_model,
                               label == 'Outcome Regression' & interactions == '--' ~ outcome_model,
                               label == 'IPW' & interactions == '--' ~ ipw_model,
                               T ~ method)) %>% 
    mutate('label' = fct_relevel(label, 'Levis', 'Outcome Regression', 'IPW')) %>% 
    mutate('y_lab' = fct_relevel(y_lab,  'Random Forest', 'GAM (Interactions)', 'GAM',
                                 'Logistic Regression (Interactions)', 'Logistic Regression',
                                 'OLS (Interactions)', 'OLS',
                                 'Levis IWOR', 'Levis IF'))
  
  
  mean_trunc <- function(x, probs) {
    df_list <- 
      df_graph %>% 
      group_by(label, y_lab, imputation_model) %>% 
      group_split()
    
    ix <- which(map_lgl(df_list, ~all(x %in% .x$ate_hat)))
    
    return(mean(df_list[[ix]]$ate_hat))
  }
  
  df_graph_ <- 
    df_graph %>% 
    group_by(label, imputation_model, interactions) %>%
    filter(ate_hat >= quantile(ate_hat, 0.05), ate_hat <= quantile(ate_hat, 0.95)) %>%
    ungroup()
  
  ggplot(df_graph_, aes(x = ate_hat, y = y_lab)) + 
    facet_wrap(~label, scales = 'free_y', ncol = 1) + 
    geom_density_ridges(aes(fill = imputation_model), 
                        alpha = 0.2, 
                        scale = 0.8, 
                        quantile_lines = T, 
                        quantile_fun = mean_trunc) + 
    geom_vline(xintercept = df_estimates$true_ate[1], lty = 2, linewidth = 2)  +
    theme(legend.position = 'bottom') + 
    labs(x = 'Estimated ATE',
         y = 'Method',
         fill = 'Imputation Method', 
         title = 'Distribution of ATE\nInner 90% of Distribution ([5%, 95%] Quantiles)',
         subtitle = glue('Simulation #{sim_id}\n# of Subjects = {df_estimates$n_patients[1]}'))
  ggsave(glue('figures/sim_{sim_id}/sim_{sim_id}_ridgeline.png'), height = 9, width = 16)
  
}

### Combine images into 1 PDF
img_files <-  glue('figures/sim_{sim_ids}/table_sim_{sim_ids}.png')
img_pdf <- reduce(map(img_files, image_read), c)
image_write(img_pdf , format = 'pdf', 'figures/sim_table.pdf')
