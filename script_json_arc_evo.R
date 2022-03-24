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
dir = "C:/Users/p288427/Desktop/data_dollo_++_3_17_22/3_23_22"
setwd(dir)
all_simple_res = data.frame()
pattern = '^m.*json$'
# pattern = 'mut_type_weights_start_arc1-2-2-2-1_act_r0.001000_dup_r0.000500_ch_A0.000000_ch_B0.010000_ch_typesymmetrical_ch_typeregular_sel_str2.0_max_arc1-2-2-2-1_sel_typesporadic_sel_freq100_1'
# pattern =  'mut_t_weights_sel_t_sporadic_sym_t_symmetrical_fr_t_regular_a_p_off_arc_1-2-2-2-1_m_arc_1-2-2-2-1_act_r_0.001_dup_r_0.000_ch_A_0.000_ch_B_0.010_s_st_2.0_s_f_0_seed0.json'
  
list.files(path = '.', pattern = pattern)
for (i in  list.files(path = '.', pattern = pattern))
{
results <- fromJSON(file = i)
simple_res = rowid_to_column(as_tibble(results[c("m_avg_fitnesses",
                                                 "m_env_functions",
                                                 "m_var_fitnesses")]),
                             var = "gen")

results$m_params$i_p$net_par$max_arc = toString(results$m_params$i_p$net_par$max_arc)
results$m_params$i_p$net_par$net_arc = toString(results$m_params$i_p$net_par$net_arc)
ID = as.data.frame(results$m_params)

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
  facet_grid(s_p.selection_freq~s_p.seed)

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
