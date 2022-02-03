library(rjson)
library(ggplot2)
library(tibble)
library(dplyr)
library(tidyr)
library(stringr)
library(ggpubr)

####read data####

# dir = dirname(rstudioapi::getActiveDocumentContext()$path)
# dir = paste(dir,"/data_sim2",sep = "")
dir = "C:/Users/Clem/build-arc_evo-Desktop_Qt_6_1_0_MinGW_64_bit-Release/release"
setwd(dir)
all_simple_res = data.frame()
pattern = "*json$"
for (i in  list.files(path = '.', pattern = pattern))
{
results <- fromJSON(file = i)
simple_res = rowid_to_column(as_tibble(results[c("m_avg_fitnesses",
                                                 "m_env_functions",
                                                 "m_var_fitnesses")]),
                             var = "gen")

i = str_replace(i, "weights_and_activation", "weightsandactivation")
ID = data.frame(i) %>% 
 separate(i, c("mut_type","architecture","mut_rate_act","mut_rate_dup","change_freq", "selection_strength", "max_arc","seed"), sep = '_')%>% 
  separate(seed, c("seed",NA))

ID$architecture = as.factor(ID$architecture)
ID$seed = as.factor(ID$seed)
ID$max_arc = as.factor(ID$max_arc)
ID$change_freq = as.factor(ID$change_freq)
ID$mut_rate_act = as.factor(ID$mut_rate_act)
ID$mut_rate_dup = as.factor(ID$mut_rate_dup)
ID$selection_strength = as.factor(ID$selection_strength)


simple_res = cbind(simple_res, ID)
all_simple_res = rbind(all_simple_res, simple_res)
}


####save load####
save(all_simple_res, file = "all_simple_res.R")
load("all_simple_res.R")
#### Plot ####

ggplot(data = all_simple_res 
       %>%filter(mut_type == "NRduplication")  
       #%>% slice_min(gen,n = 1000)
       ) +
  geom_rect(aes(xmin = gen - 1, xmax = gen,
                ymin = 0, ymax = 1.5,
                fill = as.factor(m_env_functions),
                alpha = 0.5))+
  geom_line(aes(x = gen, y = m_avg_fitnesses)) +
  facet_grid(architecture ~ seed)

####Adaptation time

d = all_simple_res %>%
  filter(change_freq == "0.005000") %>% 
  group_by(mut_type, seed) %>% 
  mutate(
    change = case_when(
      m_env_functions != lag(m_env_functions) ~ TRUE,
      TRUE ~ FALSE
    ),
    n_change = cumsum(change),
    adapted = case_when(
      m_avg_fitnesses > 0.9 ~ TRUE, 
      TRUE ~ FALSE
    )
  ) %>%
select(-change) %>%
group_by(mut_type, seed, n_change) %>%
summarise(env = as.factor(unique(m_env_functions)),
  gen = min(gen),
  time = sum(m_avg_fitnesses < 0.9),
          time_in = sum(m_avg_fitnesses > 0.9)) %>%
subset(time_in > 0) %>%
select(-time_in)

options(scipen=999)

ggplot(d,
       aes(x = gen, y = time, color = env, fill = env)) +
 geom_col() +
  geom_smooth(method='lm')+
  stat_regline_equation(aes(label = ..eq.label..)) +
  stat_regline_equation(aes(label = ..rr.label..), label.x.npc = 0.65,label.y.npc = 0.945)+
  facet_grid(mut_type ~ seed)


######### Proportion of well-adapted time

d = all_simple_res %>%
  filter(change_freq == "0.005") %>% 
  group_by(mut_type) %>% 
  mutate(
    change = case_when(
      m_env_functions != lag(m_env_functions) ~ TRUE,
      TRUE ~ FALSE
    ),
    n_change = cumsum(change),
    adapted = case_when(
      m_avg_fitnesses > 0.9 ~ TRUE, 
      TRUE ~ FALSE
    )
  ) %>%
  group_by(mut_type, n_change) %>%
    mutate (n_adapted = cumsum (adapted),
            n_gen = cumsum (change == F) + 1,
            adapt_prop = n_adapted/n_gen)%>%
  slice_tail(n=1)%>% 
  select(-change, -adapted)


ggplot(d, aes(x = n_change, y = adapt_prop, color = as.factor(m_env_functions), fill = as.factor(m_env_functions))) +
         geom_point()+ 
         facet_grid(mut_type ~ .)+
        stat_smooth(method='lm')+
          stat_regline_equation(aes(label = ..eq.label..),label.y.npc = 0.3) +
         stat_regline_equation(aes(label = ..rr.label..), label.x.npc = 0.55,label.y.npc = 0.25)


fitted_models = d %>% 
                group_by(m_env_functions,mut_type) %>% 
                do(model = summary(lm(adapt_prop~n_change, data = .))) %>% 
                unlist(fitted_models$model, recursive=F)




######### Plotting the slope and r squared of linear regressions 

fitted_models = d %>% 
  group_by(m_env_functions,mut_type) %>% 
  mutate(coeffs = summary(lm(adapt_prop~n_change, data = .))$coefficients)
