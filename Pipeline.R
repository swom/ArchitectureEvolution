library(rjson)
library(tidyr)
library(tidyverse)
library(ggplot2)
library(stringi)
library(rlist)
library(igraph)
library(ggraph)
library(ggnetwork)
library(networkD3)
library(magick)
library(tibble)
library(stringr)
library(ggpubr)
library(foreach)
library(doParallel)
library(viridis)
options(scipen=999)
registerDoParallel(detectCores()-2)

##Recovering all json files in working directory 
pattern = "*json$"
files = list.files(path = '.', pattern = pattern)

###Function to create .Rda network file from json 
create_net_file=function(i){
  if(!file.exists(paste("top_inds_net", str_replace(i, ".json", ".Rda"), sep = "_"))){
    
    results <- fromJSON(file = i)
    names(results$m_top_inds) = seq(from=0, by=1000, length.out=length(results$m_top_inds))
    results_df = results$m_top_inds %>% 
      lapply(function(l) l[[2]])%>% 
      lapply(function(l) l[[1]]) 
    results_df = as_tibble(results_df)%>%
      add_column(var = c("fitness", "input_values", "network"))%>%
      gather(key = gen, value = value, 1:length(results$m_top_inds)) %>% 
      pivot_wider(names_from = var, values_from = value)
    
    i = str_replace(i, "weights_and_activation", "weightsandactivation")
    ID = data.frame(i) %>% 
      separate(i, c("mut_type","architecture","mut_rate_act","mut_rate_dup","change_freq_A","change_freq_B","change_type", "selection_strength", "max_arc","seed"), sep = '_')%>% 
      separate(seed, c("seed",NA))
    
    ID$architecture = as.factor(ID$architecture)
    ID$seed = as.factor(ID$seed)
    ID$max_arc = as.factor(ID$max_arc)
    ID$change_freq_A = as.factor(ID$change_freq_A)
    ID$change_freq_B = as.factor(ID$change_freq_B)
    ID$mut_rate_act = as.factor(ID$mut_rate_act)
    ID$mut_rate_dup = as.factor(ID$mut_rate_dup)
    ID$selection_strength = as.factor(ID$selection_strength)
    ID$change_type = as.factor(ID$change_type)
    
    name1 = paste("top_inds", str_replace(i, ".json", ""), sep = "_")
    assign(name1, cbind(results_df, ID))
    
    ###Making a number vector out of the architecture
    if(levels(get(name1)$max_arc)[1] == ""){
      architecture = strsplit(levels(get(name1)$architecture)[1], "-")[[1]]
    } else {
      architecture = strsplit(levels(get(name1)$max_arc)[1], "-")[[1]]
    }
    architecture = as.integer(architecture)
    
    ###Keeping only network architecture, expanding to have each connection as a row
    
    top_inds_net = unnest_wider(get(name1), col = "network")%>%
      unnest_wider(col = "m_network_weights", names_sep = "_layer_")%>%
      mutate(m_input_size = NULL, input_values = NULL, fitness = NULL)%>%
      pivot_longer(cols = sprintf("m_network_weights_layer_%s", seq(1:(length(architecture)-1))), names_to = "layer")%>%
      unnest_wider(col = "value", names_sep = "_node_")%>%
      pivot_longer(cols = sprintf("value_node_%s", seq(1:(max(architecture)))), names_to = "node")%>%
      drop_na()%>%
      unnest_wider(col = "value", names_sep = "_node_")%>%
      unnest_wider(col = "value_node_m_weights", names_sep = "_")%>%
      pivot_longer(cols = sprintf("value_node_m_weights_%s", seq(1:(max(architecture)))), names_to = "weight")%>%
      drop_na()%>%
      unnest_wider(col = "value")%>%
      mutate(w_sign = if_else(m_weight < 0, 1, 2))
    top_inds_net$gen = as.factor(top_inds_net$gen)
    top_inds_net$w_sign = as.factor(top_inds_net$w_sign)
    
    name2 = paste("top_inds_net", str_replace(i, ".json", ""),"functionB", sep = "_")
    assign(name2, top_inds_net)
    get(name2)%>%saveRDS(file=paste(name2, ".Rda", sep=""))
    rm(list=c(name1,name2))
    print(i)
  }else print("skipped")
  
}

###Function to create .Rda simple_res file from json 
create_simple_res_file=function(i){
  if(!file.exists(paste("simple_res", str_replace(i, ".json", ".Rda"), sep = "_"))){
  results <- fromJSON(file = i)
  simple_res = rowid_to_column(as_tibble(results[c("m_avg_fitnesses",
                                                   "m_env_functions",
                                                   "m_var_fitnesses")]),
                               var = "gen")
    
    i = str_replace(i, "weights_and_activation", "weightsandactivation")
    ID = data.frame(i) %>% 
      separate(i, c("mut_type","architecture","mut_rate_act","mut_rate_dup","change_freq_A","change_freq_B","change_type", "selection_strength", "max_arc","seed"), sep = '_')%>% 
      separate(seed, c("seed",NA))
    
    ID$architecture = as.factor(ID$architecture)
    ID$seed = as.factor(ID$seed)
    ID$max_arc = as.factor(ID$max_arc)
    ID$change_freq_A = as.factor(ID$change_freq_A)
    ID$change_freq_B = as.factor(ID$change_freq_B)
    ID$mut_rate_act = as.factor(ID$mut_rate_act)
    ID$mut_rate_dup = as.factor(ID$mut_rate_dup)
    ID$selection_strength = as.factor(ID$selection_strength)
    ID$change_type = as.factor(ID$change_type)
    
    
    
    name1 = paste("simple_res", str_replace(i, ".json", ""), "functionB", sep = "_")
    assign(name1, cbind(simple_res, ID))
    get(name1)%>%saveRDS(file=paste(name1, ".Rda", sep=""))
    rm(list=c(name1))
    print(i)
  }else print("skipped")
}

###Looping through the two previous functions for all files, using parallel computing
foreach(i=1:length(files), .packages=.packages()) %dopar% create_simple_res_file(files[i])
foreach(i=1:length(files), .packages=.packages()) %dopar% create_net_file(files[i])


############


##Creating the big data frames where the data will be stored
all_simple_res = data.frame()
all_aggregate_data=data.frame()

##Aggregating the fitness data
res_files_fb=list.files(path = '.', pattern =glob2rx("simple_res*functionB*"))
res_files=list.files(path = '.', pattern =glob2rx("simple_res*"))
res_files=res_files[!res_files %in% res_files_fb]

simple_res_fb=map_dfr(res_files_fb, readRDS)%>%
mutate(func_B = T)

simple_res=map_dfr(res_files, readRDS)%>%
  mutate(func_B = F)

all_simple_res=rbind(simple_res,simple_res_fb)

##Aggregating the network data
net_files=list.files(path = '.', pattern =glob2rx("top_inds_net_NRaddition_1-2-2-2-1_0.001000_0.000100*"))

for(i in net_files){
  name = str_replace(i, ".Rda", "")
  assign(name,readRDS(i))}


for (i in str_replace(net_files,".Rda","")){
  
   ###Making a number vector out of the architecture
  if(get(i)$max_arc[1] == ""){
    architecture = strsplit(get(i)$architecture[1], "-")[[1]]
  } else {
    architecture = strsplit(get(i)$max_arc, "-")[[1]]
  }
  architecture = as.integer(architecture)
  
  #generate coordinates for positioning nodes correctly in plot
  x = vector()
  y = vector()
  for(j in 1:length(architecture)){
    x = c(x, rep(j, architecture[j]))
    y = c(y, seq(1:architecture[j]) + ((max(architecture)-architecture[j])/2))
  }
  
  l =   cbind(x, y)
  
  #Make the list of nodes
  layer = vector()
  node = vector()
  for(j in 1:length(architecture)){
    layer = c(layer, rep(j, architecture[j]))
    node = c(node, seq(1:architecture[j]))
  }
  id = seq(1:length(layer))
  
  node_tibble = as_tibble(cbind(id, layer, node))
  
  #Make the list of edges
  
  from = vector()
  to = vector()
  
  
  
  for(j in 1:length(levels(as.factor(node_tibble$layer)))){
    from = c(from, rep(filter(node_tibble, layer == j)$id, nrow(filter(node_tibble, layer == j+1))))
    to = c(to, sort(rep(filter(node_tibble, layer == j+1)$id, nrow(filter(node_tibble, layer == j))), decreasing = F))
  }
  
  #Now let's loop through generations
  
  aggregate_data = data.frame(gen=rep(0,length(levels(get(i)$gen))),active_nodes=rep(0,length(levels(get(i)$gen))), inactive_connections=rep(0,length(levels(get(i)$gen))), normalized_inactive=rep(0,length(levels(get(i)$gen))))
  count=0
  for(j in levels(get(i)$gen)){
    count = count + 1
    #adding weights, weight sign, activation to the edge list
    edge_tibble = as_tibble(cbind(from, to))
    ind = filter(get(i), gen == j)
    edge_tibble = cbind(edge_tibble, ind$m_weight, ind$m_is_active, ind$w_sign)
    
    node_tibble = as_tibble(cbind(id, layer, node))
    node_active = c(rep(TRUE,architecture[1]), subset(ind, ind$weight == "value_node_m_weights_1")$value_node_m_active)
    node_tibble = cbind(node_tibble, node_active)
    
    #Plot number of active nodes
    aggregate_data$active_nodes[count]=sum(node_tibble$node_active)
    aggregate_data$inactive_connections[count]=sum(edge_tibble$`ind$m_is_active`==F)
   
    con_to_active = 0
    
    for(k in 1:(length(architecture)-1)){
      con_to_active = con_to_active + sum(node_tibble$node_active[node_tibble$layer == k])*sum(node_tibble$node_active[node_tibble$layer == k+1])
    }
    
     aggregate_data$normalized_inactive[count] = aggregate_data$inactive_connections[count] / con_to_active
     aggregate_data$gen[count] = j
    
  }
  all_aggregate_data=rbind(all_aggregate_data,aggregate_data)
}


  print(i)
  aggregate_data = cbind (aggregate_data,ID)
  all_aggregate_data=rbind(all_aggregate_data, aggregate_data)
  all_aggregate_data$gen = as.numeric(all_aggregate_data$gen)
  rm(list=name1,name2)
}

save(all_simple_res, file = "all_simple_res_final.RData")
save(all_aggregate_data, file = "all_aggregate_data.RData")

#####PLOTTING TIME

####Plotting adaptation dynamics
d =   all_simple_res %>%
  filter(change_type == "regular") %>%
  mutate(time = 1/if_else(as.numeric(as.character(change_freq_A))==0,0.005,as.numeric(as.character(change_freq_A))) , avg = NA, sd = NA)%>%
  group_by(change_freq_A,mut_type,func_B,seed) %>%
  mutate(do = if_else(gen%%time ==0, T, F))%>%
  mutate(timeto = cumsum(do))%>%
  group_by(mut_rate_act,mut_rate_dup, mut_type, architecture, change_freq_A, seed,time,timeto,func_B)%>%
  summarise_at(vars(m_avg_fitnesses),list(mean = mean, sd = sd))%>%
  mutate(gen=(timeto+1)*time)%>%
  select(-timeto)%>%
  drop_na(sd)



# ###Plot with data points every 25000 generations
# d = filter(d, gen%%time ==0) %>% distinct() %>%
#   filter (gen%%25000 == 0)
# 
# ggplot(d, aes(ymin = 0, ymax = 1.2, x = gen, y = mean, color = seed, fill = seed, size=sd)) +
#   geom_point()+ 
#   facet_grid(mut_rate_act ~ mut_rate_dup)+
#   stat_smooth(method='lm')+
#   stat_regline_equation(aes(label = ..eq.label..),label.y.npc = 0.3) +
#   stat_regline_equation(aes(label = ..rr.label..), label.x.npc = 0.55,label.y.npc = 0.25)

###Plots at specific time points
d = filter(d,change_freq_A=="0.000000")

ggplot(d, aes(ymin = 0.7, ymax = 1, x = gen, y = mean)) +
  geom_line(aes(color=seed)) +
  facet_grid(func_B ~ mut_type)+
  scale_color_viridis(discrete=TRUE, option="viridis")

ggplot(d, aes(ymin = 0.25, ymax = 0, x = gen, y = sd, color=seed)) +
  geom_line() +
  facet_grid(func_B ~ mut_type)

ggplot(d, aes(ymin = 0, ymax = 1.2, x = gen, y = avg, color = as.factor(m_env_functions), fill = as.factor(m_env_functions), size=sd)) +
  geom_point(color="darkslategray3")+ 
  facet_grid(seed ~ .)+
  stat_smooth(method='lm')+
  stat_regline_equation(aes(label = ..eq.label..),label.y.npc = 0.3) +
  stat_regline_equation(aes(label = ..rr.label..), label.x.npc = 0.55,label.y.npc = 0.25)

#Plot the regression slope and last sd in circle plot
fitted_models = d %>% 
  group_by(m_env_functions, seed) %>% 
  do(rsq = summary(lm(adapt_prop~n_change, data = .))$r.squared, b = summary(lm(adapt_prop~n_change, data = .))$coefficients[1], a=summary(lm(adapt_prop~n_change, data = .))$coefficients[2])

fitted_models$rsq = unlist(fitted_models$rsq)
fitted_models$a = unlist(fitted_models$a)
fitted_models$b = unlist(fitted_models$b)
fitted_models$m_env_functions = as.factor(fitted_models$m_env_functions)

ggplot(fitted_models, aes(x = m_env_functions, y = seed, color = a, fill = a, size = sd)) + ##Size should be last sd 
  geom_point()#+
# facet_grid(m_env_functions)

#####Plotting adaptation time

d = all_simple_res %>%
  filter(as.numeric(as.character(change_freq_A))!=0) %>%
  group_by(change_freq_A, seed) %>% 
  mutate(
    change = case_when(
      m_env_functions != lag(m_env_functions) ~ TRUE,
      TRUE ~ FALSE
    ),
    n_change = cumsum(change),
    adapted = case_when(
      m_avg_fitnesses > 0.95 ~ TRUE, 
      TRUE ~ FALSE
    )
  ) %>%
  select(-change) %>%
  group_by(change_freq_A, seed, n_change) %>%
  mutate(first_adapted=which(adapted)[1])%>%
  group_by(change_freq_A, seed, n_change, first_adapted) %>%
  summarise(gen = min(gen))


ggplot(d,
       aes(x = gen, y = first_adapted, color = seed, fill = seed)) +
  geom_point() +
  # geom_smooth(method='lm')+
  # stat_regline_equation(aes(label = ..eq.label..)) +
  # stat_regline_equation(aes(label = ..rr.label..), label.x.npc = 0.65,label.y.npc = 0.945)+
  facet_grid(change_freq_A~.)

















#####Plotting badly-adapted time after a change in changing environment

d = all_simple_res %>%
  filter(as.numeric(as.character(change_freq_A))!=0) %>%
  group_by(change_freq_A, seed) %>% 
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
  group_by(change_freq_A, seed, n_change) %>%
  summarise(env = as.factor(unique(m_env_functions)),
            gen = min(gen),
            time = sum(m_avg_fitnesses < 0.9),
            time_in = sum(m_avg_fitnesses > 0.9)) #%>%
  subset(time_in > 0) %>%
  select(-time_in)


ggplot(d,
       aes(x = gen, y = time, color = seed, fill = seed)) +
  geom_point() +
  # geom_smooth(method='lm')+
  # stat_regline_equation(aes(label = ..eq.label..)) +
  # stat_regline_equation(aes(label = ..rr.label..), label.x.npc = 0.65,label.y.npc = 0.945)+
  facet_grid(change_freq_A~.)

#Plot the regression slope and last sd in circle plot
fitted_models = d %>% 
  group_by(m_env_functions, seed) %>% 
  do(rsq = summary(lm(adapt_prop~n_change, data = .))$r.squared, b = summary(lm(adapt_prop~n_change, data = .))$coefficients[1], a=summary(lm(adapt_prop~n_change, data = .))$coefficients[2])

fitted_models$rsq = unlist(fitted_models$rsq)
fitted_models$a = unlist(fitted_models$a)
fitted_models$b = unlist(fitted_models$b)
fitted_models$m_env_functions = as.factor(fitted_models$m_env_functions)

ggplot(fitted_models, aes(x = m_env_functions, y = seed, color = a, fill = a, size = rsq)) +
  geom_point()#+
# facet_grid(m_env_functions)





###Active nodes
ggplot(data = all_aggregate_data
       #%>%filter(mut_type == "NRduplication")  
       #%>% slice_min(gen,n = 1000)
) +
  geom_line(aes(x = gen, y = active_nodes))+
  facet_grid(seed~.)

###Inactive connections
ggplot(data = all_aggregate_data
       #%>%filter(mut_type == "NRduplication")  
       #%>% slice_min(gen,n = 1000)
) +
  geom_line(aes(x = gen, y = normalized_inactive))+
  facet_grid(seed~.)
