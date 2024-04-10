library(tidyverse)
library(here)

### Custom ggplot theme
theme_set(theme_bw() +
            theme(plot.title = element_text(hjust = 0.5, size = 24),
                  plot.subtitle = element_text(hjust = 0.5, size = 18),
                  axis.title = element_text(size = 20),
                  strip.text = element_text(size = 14),
                  plot.caption = element_text(size = 10),
                  legend.position = "bottom"))


df_estimates <- read_csv(here('eda/levis_skew.csv'))

df_levis_graph <-
  df_estimates %>%
  group_by(method, n_patients) %>%
  filter(ate_hat >= quantile(ate_hat, 0.05), ate_hat <= quantile(ate_hat, 0.95)) %>%
  ungroup() %>%
  select(method, n_patients, contains('hat'), contains('true')) %>%
  pivot_longer(cols = contains('hat'),
               names_to = 'contrast',
               values_to = 'estimate_hat') %>%
  pivot_longer(cols = contains('true'),
               names_to = 'contrast_true',
               values_to = 'truth') %>%
  mutate('contrast' = gsub('_hat', '', contrast),
         'contrast_true' = gsub('true_', '', contrast_true)) %>%
  filter(contrast == contrast_true) %>%
  mutate(contrast = toupper(contrast)) %>% 
  mutate('contrast' = case_when(contrast == 'ATE' ~ 'Average Treatment Effect (ATE)',
                                contrast == 'TE0' ~ 'E[Y(0)]',
                                contrast == 'TE1' ~ 'E[Y(1)]')) %>% 
  mutate(method = gsub('Levis', 'CCMAR', method))

ggplot(df_levis_graph, aes(x = estimate_hat)) +
  facet_grid(method ~ contrast, scales = 'free') +
  geom_vline(aes(xintercept = truth), lty = 2, lwd = 0.8, alpha = 0.8) +
  geom_density(aes(fill = as.factor(n_patients)), alpha = 0.35, ) +
  labs(x = 'Estimated Contrast',
       y = 'Density',
       col = '# of Patients',
       fill = '# of Patients',
       title = 'Distribution of Treatment Effect Estimates',
       subtitle = 'Data Driven Simulation Scenario 1') +
  scale_fill_discrete(labels = function(x) scales::number(as.numeric(x), big.mark = ','))
ggsave(here('paper/levis_skew_by_contrast.pdf'), height = 9, width = 16)

