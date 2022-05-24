library(rjson)
library(ggplot2)
library(tibble)
library(dplyr)
library(tidyr)
library(stringr)
library(ggpubr)
library(patchwork)

# load data.table
library(data.table)

extract_sensibilities <- function (x) {
  as.data.frame(do.call(rbind,do.call(cbind, x$m_sensibilities))) %>%
    mutate(across(everything(), ~ as.numeric(.x)))
}

####read data####

# dir = dirname(rstudioapi::getActiveDocumentContext()$path)
# dir = paste(dir,"/data_sim2",sep = "")
dir = "C:/Users/p288427/Desktop/data_dollo_++/5_19_22_test_new/"
setwd(dir)
all_simple_res = data.frame()

pattern = '^m.*json$'
# pattern =  'mut_t_weights_sel_t_sporadic_sym_t_symmetrical_fr_t_regular_a_p_off_arc_1-2-2-2-1_m_arc_1-2-2-2-1_act_r_0.001_dup_r_0.000_ch_A_0.000_ch_B_0.010_s_st_2.0_s_f_0_seed0.json'
# pattern = "mut_t_weights_sel_t_sporadic_sym_t_symmetrical_fr_t_regular_a_p_on_arc_1-2-2-2-1_m_arc_1-2-2-2-1_act_r_0.001_dup_r_0.000_ch_A_0.000_ch_B_0.010_s_st_2.0_s_f_0_seed2"
# pattern = "mut_t_weights_sel_t_sporadic_sym_t_symmetrical_fr_t_regular_a_p_off_arc_1-2-2-2-1_m_arc_1-2-2-2-1_act_r_0.001_dup_r_0.000_ch_A_0.000_ch_B_0.010_s_st_2.0_s_f_0_seed2"
# pattern = "mut_t_weights_sel_t_sporadic_sym_t_symmetrical_fr_t_regular_a_p_off_arc_1-2-2-2-1_m_arc_1-2-2-2-1_act_r_0.001_dup_r_0.000_ch_A_0.000_ch_B_0.010_s_st_0.1_s_f_0_seed3.json"
# pattern = "mut_t_weights_sel_t_spo_sym_t_sym_fr_t_reg_a_p_on_r_t_con_arc_1-2-2-2-1_m_arc_1-2-2-2-1_act_r_0.001_dup_r_0.000_ch_A_0.000_ch_B_0.010_s_st_1.0_s_f_100_seed2.json"

####save load####
if(file.exists("all_simple_res.Rds") && file.exists("all_sensibilities.Rds")){
  all_simple_res <- readRDS("all_simple_res.Rds")
  all_sensibilities <- readRDS("all_sensibilities.Rds")
}else{
  
  filepaths = list.files(path = ".", pattern = pattern)
  m_files = length(filepaths)
  
  all_sensibilities = list()
  
  for (i in  filepaths)
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
    
    # m fit suitability
    tmp_ = results$m_fit_phen_mut_sensibility
    
    # convert to data.table with list column
    tmp_ = lapply(tmp_, as.data.table)
    
    # deal with all gens
    tmp_ = lapply(tmp_, function(df) {
      df[, m_sensibilities := lapply(m_sensibilities, function(le) {
        as.data.table(le)
      })]
      
      df = df[, rbindlist(m_sensibilities), by = "m_generation"]
      df
    })
    
    tmp_ = rbindlist(tmp_)
    tmp_[, names(ID) := ID]
    
    tmp_
    all_sensibilities[[i]] = tmp_
    
  }
  
  all_sensibilities = rbindlist(all_sensibilities)
  
  saveRDS(all_simple_res, file = "all_simple_res.Rds")
  saveRDS(all_sensibilities, file = "all_sensibilities.Rds")
}
#### Plot ####
jpeg("fitness_plots.jpg",
     width = 700,
     height = 700)

filter_gen = 1000
show_last_n_gen = 1000000
wanted_freqs = c(1)
wanted_sel_str = c(0.1, 0.5, 1)
p <- all_simple_res %>% 
  filter(gen > max(gen) - show_last_n_gen #&
         # s_p.seed %in% wanted_seed &
         # gen %% filter_gen == 0 &
         # s_p.selection_strength %in% wanted_sel_str
  ) %>% 
  ggplot() +
  ###print all environments
  geom_rect(data = . %>%
              group_by(s_p.change_freq_B, s_p.adaptation_per) %>%
              distinct(gen, m_env_functions) %>%
              filter(gen == min(gen) |
                       gen == max(gen) |
                       m_env_functions != lag(m_env_functions)) %>%
              mutate(gen_max = lead(gen)) %>%
              mutate(gen_max = ifelse(is.na(gen_max),
                                      max(gen),
                                      gen_max)),
            aes(xmin = gen, xmax = gen_max,
                ymin = 0, ymax = 1,
                fill = as.factor(m_env_functions),
            )
  ) +
  geom_line(data = . %>% filter(gen %% filter_gen == 0),
            aes(x = gen, y = m_avg_fitnesses)
  ) +
  # geom_smooth(method='lm',aes(x = gen, y = m_avg_fitnesses))+
  # stat_regline_equation(aes(label = ..eq.label.., x = gen, y = m_avg_fitnesses),label.y.npc = 0.9) +
  # stat_regline_equation(aes(label = ..rr.label.., x = gen, y = m_avg_fitnesses), label.x.npc = 0.55,label.y.npc = 0.9)+
  # facet_grid(s_p.change_freq_B + s_p.selection_strength ~
  #              s_p.seed + s_p.adaptation_per + i_p.net_par.resp_type)
  facet_grid(s_p.selection_freq + s_p.selection_strength ~
               s_p.seed + s_p.adaptation_per)

print(p)
dev.off()


########Sensibilities


all_sensibilities$m_generation = as.factor(all_sensibilities$m_generation)
for(gen in levels(all_sensibilities$m_generation))
{
  phen_x_lim = c(0,0.3)
  fit_x_lim = c(-0.1,0.1)
  y_lim = c(0,50)
  
  n_bins = 100
  
 gen = "251999"

  gen_sens = all_sensibilities %>%
    filter(m_generation == gen) %>% 
    filter(s_p.adaptation_per == 0) %>% 
    filter(s_p.seed == 1) %>% 
    filter(s_p.selection_strength == 1) %>% 
    filter(s_p.selection_freq == 100)
  
  p1 = ggplot(data = gen_sens) +
    geom_histogram(aes(m_fitness), bins = n_bins) +
    xlim(fit_x_lim) +
    ylim(y_lim)
  
  
  p2 = ggplot(data = gen_sens) +
    geom_histogram(aes(m_phenotype), bins = n_bins) +
    xlim(phen_x_lim) +
    ylim(y_lim) +
    coord_flip()
  
  p3 = 
    ggplot(data = gen_sens) +
    geom_point(aes(x = m_fitness, y = m_phenotype, color = m_rank, alpha = 1), alpha = 0.005) +
    xlim(fit_x_lim) +
    ylim(phen_x_lim)
  
  layout <- "
AA##
AA##
CCDD
CCDD" 

p1 + p3 + p2  + 
    plot_layout(design = layout,guides = 'collect', widths = 1) 

  ggsave(paste(paste("phen_fit_sens_plot",gen,sep = "_"),".png"), device = "png")
}




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
