# Chignik weir and sonar full hour comparisons
# Sarah Power
# 2019-09

# load ----
source('code/helper.r')
source('code/functions.r')

# data ----
read_csv('data/chigWeirDidson201618sjp.csv') %>%
  dplyr::select(-starts_with("prop")) %>%
  gather(species, abundance, sockeye_10:total_60) %>%
  separate(species, c("species", "period"), sep ="_") %>%
  mutate(date = as.Date(date, "%m/%d/%Y"),
         year = base::factor(format(date, "%Y")),
         abundance = str_replace(abundance, "-", "0"),
         abundance = as.numeric(abundance),
         period = str_replace(period, "10", "ten_minute"),
         period = str_replace(period, "60", "sixty_minute")) %>% 
  rename(count_method = method) %>% 
  filter(species %in% c("sockeye", "coho", "total")) -> data


data %>%
  filter_all(any_vars(is.na(.)))

# model ----

data %>% 
  spread(count_method, abundance) %>% 
  drop_na %>% 
  group_by(species, year) %>%
  nest() %>% 
  mutate(fit = purrr::map(data, ~ lm(sonar ~ weir, data = .)),
         slopes = purrr::map(fit, ~ slope_eq_1(.), data = .),
         shapiro = purrr::map(fit, ~shapiro.test(.$residuals)),
         tidy = purrr::map(fit, ~ tidy(.x)),
         glance = purrr::map(fit, ~glance(.x)),
         shapiro = purrr::map(shapiro, ~tidy(.x))) -> fit_sw

data %>% 
  filter(period == "sixty_minute") %>% 
  group_by(species, year) %>% 
  nest() %>% 
  mutate(wilcox = map(data, ~ wilcox.test(abundance ~ count_method, 
                                          data = .,
                                          paired = TRUE, 
                                          alternative = "two.sided")),
         wilcox = map(wilcox, tidy)) %>% 
  unnest(wilcox) %>% 
  dplyr::select(species, year, wilcox = p.value) %>% 
  ungroup %>% 
  mutate(reject_wilcox = ifelse(wilcox >(0.05 / n()), "no", "yes")) -> sw_wilcox

# output ----
fit_sw %>% 
  unnest(glance) %>% 
  dplyr::select(-c(data, fit, slopes, shapiro, tidy)) -> sw_summary

fit_sw %>% 
  unnest(tidy) %>% 
  dplyr::select(-c(data, fit, slopes, shapiro, glance)) -> ten_sixty_params

fit_sw %>% 
  unnest(slopes) %>% 
  dplyr::select(species, year, slopes) -> sw_slope

fit_sw %>%
  unnest(shapiro) %>% 
  dplyr::select(species, year, shapiro = p.value) %>% 
  ungroup %>% 
  mutate(normal = ifelse(shapiro > 0.05, "yes", "no"),
         bonf_adj = 0.05 / sum(normal=="yes")) -> sw_shapiro

#sw_wilcox %>% 
#left_join(sw_shapiro) %>% 
sw_shapiro %>% 
  left_join(sw_slope) %>% 
  mutate(slope_eq_1 = case_when(normal=="yes" & slopes>= bonf_adj ~ "yes", 
                                normal=="yes" & slopes< bonf_adj ~ "no")) -> slope_test
# figures ----
data %>% 
  spread(count_method, abundance) %>% 
  ggplot(aes(weir, sonar)) +
  geom_point() +
  stat_smooth(method = "lm") +
  facet_grid(year ~ species, scales = "free")

data %>% 
  spread(count_method, abundance) %>% 
  ggplot(aes(weir, sonar)) +
  geom_point() +
  stat_smooth(method = "lm") +
  facet_wrap(year ~ species, scales = "free")


# using facet_grid -----------------------------------------------------------

# retain r2 values for labelling
sw_summary %>% 
  ungroup %>% 
  mutate(Species = case_when(species=="coho" ~ "Coho",
                             species=="sockeye" ~ "Sockeye",
                             TRUE ~ "Total"),
         r.squared = format(r.squared, digits = 3)) %>% 
  dplyr::select(Species, year, r.squared) -> rsq

# retain non-normal and significant
# cheater method since wilcoxin test isn't operating currently

data.frame(Species = rep(c("Coho", "Sockeye", "Total"),each = 3),
           year = rep(2016:2018, 3),
           sig = c(rep("*", 3), rep(c("", "", "*"), 2))) -> sigs

data %>% 
  spread(count_method, abundance) %>% 
  mutate(Species = case_when(species=="coho" ~ "Coho",
                             species=="sockeye" ~ "Sockeye",
                             TRUE ~ "Total")) %>% 
  drop_na %>% 
  group_split(Species, year) %>%
  map_df(~{fit = lm(sonar ~ weir, data = .)
  
  data.frame(., ci = predict(fit, ., interval = 'confidence'),
             pi = predict(fit, ., interval = 'prediction'))
  }) %>% 
  mutate(pi.lwr = ifelse(pi.lwr<0, 0, pi.lwr)) %>% 
  ggplot(aes(weir, ci.fit)) + 
  geom_point(aes(y = sonar), alpha = 0.2) +
  geom_line() +
  geom_ribbon(aes(ymin = (ci.lwr), ymax = (ci.upr)), alpha = 0.2, color = NA) +
  geom_ribbon(aes(ymin = (pi.lwr), ymax = (pi.upr)), alpha = 0.1, color = NA) +
  geom_abline(slope = 1, lty = 3, alpha = 0.5) +
  facet_wrap(year~Species, 
             scales = "free",
             labeller = label_wrap_gen(multi_line=FALSE)) +
  scale_y_continuous(labels = scales::comma, name = "Sonar counts") +
  scale_x_continuous(labels = scales::comma, name = "Weir counts") +
  expand_limits(y=0, x=0) +
  theme(strip.text = element_text(hjust = 0)) #+
  # geom_text(data = rsq, aes(x = Inf, y = -Inf, label = paste("r^2==", r.squared)), 
  #           parse = TRUE, size = 3, vjust = -1, hjust = 1) +
  # geom_text(data = sigs, aes(x = Inf, y = Inf, label =  sig), 
  #           parse = FALSE, size = 5, vjust = 2, hjust = 2) +
  # theme(strip.background = element_blank(),
  #   strip.text.x = element_blank())

ggsave("figs/sonar_weir_compare.png", width = 9, height = 6.5, units = "in")


