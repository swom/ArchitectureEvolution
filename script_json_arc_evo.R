library(rjson)
library(ggplot2)
library(tibble)
library(dplyr)
library(tidyr)
library(stringr)



####Function declaration and definition ####
Timetochangecalc = function (threshold, simple_res){
  curenv=0
  curchange=0
  gen_change = 0
  timing=FALSE
  
  stor=data.frame(matrix(ncol = 3, nrow = 0))
  x <- c("time", "change_n", "gen_change")
  colnames(stor) <- x
  
  for (i in 1: nrow(simple_res)){
   
     if (curenv!=simple_res$m_env_values[i]) {
      curenv=simple_res$m_env_values[i]
      curchange = curchange + 1
      gen_change = i
      timing=TRUE
      elapsed=0
    }
    
    if (timing == TRUE){
      elapsed=elapsed + 1
      if (simple_res$m_avg_fitnesses[i] > threshold){
        timing=FALSE
        stor= rbind(stor, data.frame(time = elapsed, 
                                     change_n = curchange,
                                     gen_change = gen_change))
      }
    }  
    
  }
  
  return (data.frame(stor))
}

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

ggplot(data = all_simple_res #%>%
         #filter(change_freq == "0.005000")  
       # %>% slice_min(gen,n = 1000)
       ) +
  geom_rect(aes(xmin = gen - 1, xmax = gen,
                ymin = 0, ymax = 1.5,
                fill = as.factor(m_env_functions),
                alpha = 0.5))+
  geom_line(aes(x = gen, y = m_avg_fitnesses)) #+
  #facet_grid(change_freq ~ architecture)



d = all_simple_res %>%
  group_by(architecture) %>% 
  mutate(
    change = case_when(
      m_env_functions != lag(m_env_functions) ~ TRUE,
      TRUE ~ FALSE
    ),
    n_change = cumsum(change)
  ) %>%
  select(-change) %>% 
  group_by(architecture, n_change) %>% 
  summarise(env = as.factor(unique(m_env_functions)),
    gen = min(gen),
    time = sum(m_avg_fitnesses < 0.9),
            time_in = sum(m_avg_fitnesses > 0.9)) %>% 
  subset(time_in > 0) %>% 
  select(-time_in)


ggplot(d %>% slice_min(gen, n = 10000),
       aes(x = gen, y = time, color = env, fill = env)) +
 geom_col() +
  facet_grid(architecture ~ . )

#Timo's code
output = Timetochangecalc(0.9,simple_res)

barplot(output$time)
