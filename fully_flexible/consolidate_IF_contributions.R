library(tidyverse)
library(glue)
library(ggridges)
library(furrr)
library(gt)
plan(multisession(workers = 8))

### Check if Any Files didn't process
n_batch <- 2000
setdiff(paste0('batch_', 1:n_batch, '.csv'), dir('subject_contributions/'))
setdiff(paste0('parametric_batch_', 1:n_batch, '.csv'), dir('subject_contributions/'))


### Read In All Contributions
df_contributions <- 
  future_map_dfr(paste0('subject_contributions/', intersect(paste0('batch_', 1:n_batch, '.csv'), dir('subject_contributions/'))), 
                 ~read_csv(.x, col_types = cols())) %>% 
  rename('pi.hat.1' = 'pi.hat')

df_contributions <- 
  df_contributions %>% 
  arrange(dataset_id, fold_id, batch_id, subject_uuid) 

df_contributions_parametric <- 
  future_map_dfr(paste0('subject_contributions/', intersect(paste0('parametric_batch_', 1:n_batch, '.csv'), dir('subject_contributions/'))), 
                 ~read_csv(.x, col_types = cols())) %>% 
  rename('pi.hat.1' = 'pi.hat')

df_contributions_parametric <- 
  df_contributions_parametric %>% 
  arrange(dataset_id, fold_id, batch_id, subject_uuid)

df_folds <- 
  future_map_dfr(dir('input_data/', full.names = T),
                 ~read_csv(.x, col_types = cols()))


### Merge in original dataset with computed components 
df_subject_contributions <- 
  bind_rows(
    df_contributions %>% mutate('model_type' = 'Non-Parametric'),
    df_contributions_parametric %>% mutate('model_type' = 'Parametric')
  )


df_final <- 
  df_folds %>% 
  inner_join(df_subject_contributions, by = c('dataset_id', 'subject_id', 'fold_id', 'batch_id')) 

### IWOR
df_iwor <- 
  df_final %>% 
  mutate('IOR.1' = ifelse(S == 1, (beta.hat.1/gamma.hat.1)/pi.hat.1, 0),
         'IOR.0' = ifelse(S == 1, (beta.hat.0/gamma.hat.0)/pi.hat.1, 0)) %>% 
  group_by(dataset_id, fold_id, model_type) %>% 
  summarise('ate_hat' = mean(IOR.1 - IOR.0, na.rm = T)) %>% 
  group_by(dataset_id, model_type) %>% 
  summarise('method' = 'Levis IWOR',
            'ate_hat' = mean(ate_hat)) 

### IF
df_if <- 
  df_final %>% 
  mutate('IF.0' = b.01 + ifelse(A == 0, b.02/(1 - eta.hat), 0) + 
           ifelse(S == 1, (beta.hat.0 / gamma.hat.0 - b.01 +
                             ifelse(A == 0, ((gamma.hat.0 * (1 - eta.hat) + gamma.hat.1 * eta.hat) * 
                                               (Y - beta.hat.0 / gamma.hat.0) / gamma.hat.0 - b.02) / (1 - eta.hat), 0)) / pi.hat.1, 0),
         'IF.1' = b.11 + ifelse(A == 1, b.12/eta.hat, 0) +
           ifelse(S == 1, (beta.hat.1 / gamma.hat.1 - b.11 +
                             ifelse(A == 1, ((gamma.hat.0 * (1 - eta.hat) + gamma.hat.1 * eta.hat) * 
                                               (Y - beta.hat.1 / gamma.hat.1) / gamma.hat.1 - b.12) / eta.hat, 0)) / pi.hat.1, 0)) %>% 
  group_by(dataset_id, fold_id, model_type) %>% 
  summarise('ate_hat' = mean(IF.1 - IF.0, na.rm = T)) %>% 
  group_by(dataset_id, model_type) %>% 
  summarise('method' = 'Levis IF',
            'ate_hat' = mean(ate_hat)) 


### True ATE
## parametric truth computation
delta.y <- 0.0025
y.grid <- seq(0, 1, by = delta.y)
delta.l <- 0.013875
l.grid <- seq(-3.3, 7.8, by = delta.l)

## true nuisance models
eta <- 0.5
y.pdf <- function(a, y) {
  ifelse(a == 0, dbeta(y, 2, 4), dbeta(y, 4, 2))
}
y.pdf <- Vectorize(y.pdf)
l1.pmf <- function(a, y, l1) {
  p <- plogis(-0.6 + 0.5 * a + 0.25 * y + 0.1 * a * y)
  ifelse(l1 == 1, p, 1 - p)
}
l1.pmf <- Vectorize(l1.pmf)
l2.pdf <- function(a, y, l1, l2) {
  dnorm(l2, 1.0 * a + y + 2.5 * l1 * y, 1.25)
}
l2.pdf <- Vectorize(l2.pdf)

# L1 = 0, A = 0
gamma.0.0 <- sapply(1:length(l.grid), function(j) {
  sum(l1.pmf(a = 0, y = y.grid, l1 = 0) *
        l2.pdf(a = 0, y = y.grid, l1 = 0, l2 = l.grid[j]) *
        y.pdf(a = 0, y = y.grid)) * delta.y
}, simplify = 0)
beta.0.0 <- sapply(1:length(l.grid), function(j) {
  sum(l1.pmf(a = 0, y = y.grid, l1 = 0) *
        l2.pdf(a = 0, y = y.grid, l1 = 0, l2 = l.grid[j]) *
        y.pdf(a = 0, y = y.grid) * y.grid) * delta.y
}, simplify = 0)
lambda.mat.0.0 <- sapply(1:length(y.grid), function(j) {
  l2.pdf(a = 0, y = y.grid[j], l1 = 0, l2 = l.grid)
})

# L1 = 0, A = 1
gamma.0.1 <- sapply(1:length(l.grid), function(j) {
  sum(l1.pmf(a = 1, y = y.grid, l1 = 0) *
        l2.pdf(a = 1, y = y.grid, l1 = 0, l2 = l.grid[j]) *
        y.pdf(a = 1, y = y.grid)) * delta.y
}, simplify = 0)
beta.0.1 <- sapply(1:length(l.grid), function(j) {
  sum(l1.pmf(a = 1, y = y.grid, l1 = 0) *
        l2.pdf(a = 1, y = y.grid, l1 = 0, l2 = l.grid[j]) *
        y.pdf(a = 1, y = y.grid) * y.grid) * delta.y
}, simplify = 0)
lambda.mat.0.1 <- sapply(1:length(y.grid), function(j) {
  l2.pdf(a = 1, y = y.grid[j], l1 = 0, l2 = l.grid)
})

# L1 = 1, A = 0
gamma.1.0 <- sapply(1:length(l.grid), function(j) {
  sum(l1.pmf(a = 0, y = y.grid, l1 = 1) *
        l2.pdf(a = 0, y = y.grid, l1 = 1, l2 = l.grid[j]) *
        y.pdf(a = 0, y = y.grid)) * delta.y
}, simplify = 0)
beta.1.0 <- sapply(1:length(l.grid), function(j) {
  sum(l1.pmf(a = 0, y = y.grid, l1 = 1) *
        l2.pdf(a = 0, y = y.grid, l1 = 1, l2 = l.grid[j]) *
        y.pdf(a = 0, y = y.grid) * y.grid) * delta.y
}, simplify = 0)
lambda.mat.1.0 <- sapply(1:length(y.grid), function(j) {
  l2.pdf(a = 0, y = y.grid[j], l1 = 1, l2 = l.grid)
})


# L1 = 1, A = 1
gamma.1.1 <- sapply(1:length(l.grid), function(j) {
  sum(l1.pmf(a = 1, y = y.grid, l1 = 1) *
        l2.pdf(a = 1, y = y.grid, l1 = 1, l2 = l.grid[j]) *
        y.pdf(a = 1, y = y.grid)) * delta.y
}, simplify = 0)
beta.1.1 <- sapply(1:length(l.grid), function(j) {
  sum(l1.pmf(a = 1, y = y.grid, l1 = 1) *
        l2.pdf(a = 1, y = y.grid, l1 = 1, l2 = l.grid[j]) *
        y.pdf(a = 1, y = y.grid) * y.grid) * delta.y
}, simplify = 0)
lambda.mat.1.1 <- sapply(1:length(y.grid), function(j) {
  l2.pdf(a = 1, y = y.grid[j], l1 = 1, l2 = l.grid)
})

ATE.true.0 <- 
  eta * sum(delta.y * y.pdf(a = 0, y = y.grid) *
              (l1.pmf(a = 0, y = y.grid, l1 = 0) *
                 sapply(1:length(y.grid), function(j) {
                   sum(beta.0.0 / gamma.0.0 * 
                         lambda.mat.0.0[,j] * delta.l)
                 }, simplify = 0) +
                 l1.pmf(a = 0, y = y.grid, l1 = 1) *
                 sapply(1:length(y.grid), function(j) {
                   sum(beta.1.0 / gamma.1.0 * 
                         lambda.mat.1.0[,j] * delta.l)
                 }, simplify = 0))) +
  eta * sum(delta.y * y.pdf(a = 1, y = y.grid) *
              (l1.pmf(a = 1, y = y.grid, l1 = 0) *
                 sapply(1:length(y.grid), function(j) {
                   sum(beta.0.0 / gamma.0.0 * 
                         lambda.mat.0.1[,j] * delta.l)
                 }, simplify = 0) +
                 l1.pmf(a = 1, y = y.grid, l1 = 1) *
                 sapply(1:length(y.grid), function(j) {
                   sum(beta.1.0 / gamma.1.0 * 
                         lambda.mat.1.1[,j] * delta.l)
                 }, simplify = 0)))

ATE.true.1 <- 
  eta * sum(delta.y * y.pdf(a = 0, y = y.grid) *
              (l1.pmf(a = 0, y = y.grid, l1 = 0) *
                 sapply(1:length(y.grid), function(j) {
                   sum(beta.0.1 / gamma.0.1 * 
                         lambda.mat.0.0[,j] * delta.l)
                 }, simplify = 0) +
                 l1.pmf(a = 0, y = y.grid, l1 = 1) *
                 sapply(1:length(y.grid), function(j) {
                   sum(beta.1.1 / gamma.1.1 * 
                         lambda.mat.1.0[,j] * delta.l)
                 }, simplify = 0))) +
  eta * sum(delta.y * y.pdf(a = 1, y = y.grid) *
              (l1.pmf(a = 1, y = y.grid, l1 = 0) *
                 sapply(1:length(y.grid), function(j) {
                   sum(beta.0.1 / gamma.0.1 * 
                         lambda.mat.0.1[,j] * delta.l)
                 }, simplify = 0) +
                 l1.pmf(a = 1, y = y.grid, l1 = 1) *
                 sapply(1:length(y.grid), function(j) {
                   sum(beta.1.1 / gamma.1.1 * 
                         lambda.mat.1.1[,j] * delta.l)
                 }, simplify = 0)))

ATE.true <- ATE.true.1 - ATE.true.0   


### Summary of results
bind_rows(df_iwor, df_if) %>% 
  ungroup() %>% 
  write_csv('results/levis.csv')

df_levis <- read_csv('results/levis.csv')

gt_levis <- 
  df_levis %>% 
  group_by(method, model_type) %>% 
  summarise('mean_ate' = mean(ate_hat),
            'median_ate' = median(ate_hat),
            'mean_bias' = mean(ate_hat - ATE.true),
            'rel_bias' = mean((ate_hat - ATE.true)/ATE.true),
            'sd_ate' = sd(ate_hat)) %>%
  ungroup() %>% 
  mutate('rel_efficiency' = sd_ate/sd_ate[model_type == 'Parametric' & method == 'Levis IF']) %>% 
  gt() %>% 
  cols_align('center') %>% 
  fmt_number(contains('ate'), decimals = 3) %>% 
  fmt_number('rel_efficiency', decimals = 3) %>% 
  fmt_number('mean_bias', decimals = 4) %>% 
  fmt_percent('rel_bias', decimals = 1)  %>% 
  cols_label('method' = 'Estimator',
             'model_type' = 'Model Type', 
             'mean_ate' = 'Mean',
             'median_ate' = 'Median',
             'mean_bias' = 'Bias',
             'rel_bias' = '% Bias',
             'sd_ate' = 'SD', 
             'rel_efficiency' = 'Rel. Uncertainty') %>% 
  tab_header(title = md('**Fully Flexible Simulation**')) %>% 
  tab_options(column_labels.font.size = 16,
              heading.title.font.size = 40,
              heading.subtitle.font.size = 20,
              column_labels.font.weight = 'bold',
              row_group.font.weight = 'bold',
              row_group.font.size  = 22)

gtsave(gt_levis, 'figures/levis_results.png')


### Custom ggplot theme
theme_set(theme_bw() +
            theme(plot.title = element_text(hjust = 0.5, size = 24),
                  plot.subtitle = element_text(hjust = 0.5, size = 18),
                  axis.title = element_text(size = 20),
                  strip.text = element_text(size = 14),
                  plot.caption = element_text(size = 10),
                  legend.position = "none"))


### Plot
ggplot(df_levis, aes(x = ate_hat)) +
  geom_density(aes(fill =paste(method, model_type)), alpha = 0.5) +
  facet_grid(model_type ~ method) +
  geom_vline(xintercept = ATE.true, lty = 2) + 
  labs(x = 'Estimated ATE',
       y = 'Density',
       title = 'Distribution of Levis Estimators',
       subtitle = 'Fully Flexible Simulation')
ggsave('figures/levis_ff.png', height = 9/1.2, width = 16/1.2)
ggsave('figures/levis_ff.pdf', height = 9/1.2, width = 16/1.2)
