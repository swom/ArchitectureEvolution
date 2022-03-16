library(rjson)
library(ggplot2)
library(tibble)
library(dplyr)
library(tidyr)
library(stringr)
library(ggpubr)

####read data####

dir = dirname(rstudioapi::getActiveDocumentContext()$path)
# dir = paste(dir,"/data_sim2",sep = "")
# dir = "C:/Users/Clem/build-arc_evo-Desktop_Qt_6_1_0_MinGW_64_bit-Release"
setwd(dir)
all_simple_res = data.frame()
pattern = '^m.*json$'
list.files(path = '.', pattern = pattern)
for (i in  list.files(path = '.', pattern = pattern))
{
results <- fromJSON(file = i)
simple_res = rowid_to_column(as_tibble(results[c("m_avg_fitnesses",
                                                 "m_env_functions",
                                                 "m_var_fitnesses")]),
                             var = "gen")

i = str_replace(i, "weights_and_activation", "weightsandactivation")
a_p = results$m_params

as_tibble(a_p$s_p)
as_tibble(a_p$e_p)
as_tibble(a_p$p_p)
as_tibble(a_p$i_p)

ID = data.frame(NA)
ID$architecture <- toString(a_p$i_p$net_par$net_arc)
ID$seed = as.factor(a_p$s_p$seed)
ID$max_arc = toString(a_p$i_p$net_par$max_arc)
ID$change_freq_A = as.factor(a_p$s_p$change_freq_A)
ID$change_freq_B = as.factor(a_p$s_p$change_freq_B)
ID$mut_rate_act = as.factor(a_p$p_p$mut_rate_activation)
ID$mut_rate_dup = as.factor(a_p$p_p$mut_rate_duplication)
ID$selection_strength = as.factor(a_p$s_p$selection_strength)
ID$change_type = as.factor(a_p$s_p$change_freq_type)
ID = ID[-1]

simple_res_ID = cbind(simple_res, ID)
all_simple_res = rbind(all_simple_res, simple_res_ID)
}


####save load####
save(all_simple_res, file = "all_simple_res.R")
load("all_simple_res.R")
#### Plot ####

ggplot(data = all_simple_res 
       #%>%filter(mut_type == "NRduplication")  
       #%>% slice_min(gen,n = 1000)
       ) +
  geom_rect(aes(xmin = gen - 1, xmax = gen,
                ymin = 0, ymax = 1.5,
                fill = as.factor(m_env_functions),
                alpha = 0.5))+
  geom_line(aes(x = gen, y = m_avg_fitnesses)) +
  geom_smooth(method='lm',aes(x = gen, y = m_avg_fitnesses))+
  stat_regline_equation(aes(label = ..eq.label.., x = gen, y = m_avg_fitnesses),label.y.npc = 0.9) +
  stat_regline_equation(aes(label = ..rr.label.., x = gen, y = m_avg_fitnesses), label.x.npc = 0.55,label.y.npc = 0.9)+
  facet_grid(change_freq_A~.)

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
  filter(architecture == "2-8-8-8-1") %>% 
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
  group_by(mut_type, seed, n_change) %>%
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






######### Plotting the slope and r squared of linear regressions 

fitted_models = d %>% 
  group_by(m_env_functions,mut_type, seed) %>% 
  do(rsq = summary(lm(adapt_prop~n_change, data = .))$r.squared, b = summary(lm(adapt_prop~n_change, data = .))$coefficients[1], a=summary(lm(adapt_prop~n_change, data = .))$coefficients[2])

fitted_models$rsq = unlist(fitted_models$rsq)
fitted_models$a = unlist(fitted_models$a)
fitted_models$b = unlist(fitted_models$b)
fitted_models$m_env_functions = as.factor(fitted_models$m_env_functions)

ggplot(fitted_models, aes(x = mut_type, y = seed, color = a, fill = a, size = rsq)) +
  geom_point()+
  facet_grid(m_env_functions ~ .)

####Plotting what's been going on since last env change
d =   all_simple_res %>%
  filter(change_type == "regular") %>%
  mutate(time = 1/as.numeric(as.character(change_freq_A)), avg = NA, sd = NA)

for(i in 1:nrow(d)){
  if(d$gen[i] %% d$time[i] == 0){
    d$avg[i] = mean(d$m_avg_fitnesses[(i - d$time[i] + 1) : i])
    d$sd[i] = sd(d$m_avg_fitnesses[(i - d$time[i] + 1) : i])
  }
}

d = filter(d, gen%%time ==0) %>% distinct()
 
ggplot(d, aes(x = gen, y = avg, color = as.factor(m_env_functions), fill = as.factor(m_env_functions))) +
  geom_point()+ 
  facet_grid(change_freq_A ~ .)+
  stat_smooth(method='lm')+
  stat_regline_equation(aes(label = ..eq.label..),label.y.npc = 0.3) +
  stat_regline_equation(aes(label = ..rr.label..), label.x.npc = 0.55,label.y.npc = 0.25)+
  geom_errorbar(aes(ymin=avg-sd, ymax=avg+sd), width=.2,
                position=position_dodge(.9)) 
