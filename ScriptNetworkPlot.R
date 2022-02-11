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


dir = "C:/Users/Clem/build-arc_evo-Desktop_Qt_6_1_0_MinGW_64_bit-Release/"
setwd(dir)
options(scipen=999)

results=list()
pattern = "*json$"
all_aggregate_data=data.frame()
for (i in  list.files(path = '.', pattern = pattern)){
  
  ###Making a data tibble with all top individuals' data 
  
  results <- fromJSON(file = i)
  
  simple_res = rowid_to_column(as_tibble(results[c("m_avg_fitnesses",
                                                   "m_env_functions",
                                                   "m_var_fitnesses")]),
                               var = "gen")
  
  
  names(results$m_top_inds) = seq(from=0, by=1000, length.out=length(results$m_top_inds))
  results_df = results$m_top_inds %>% 
  unlist(recursive=FALSE)#%>%
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
  simple_res = cbind(simple_res, ID)
  all_simple_res = rbind(all_simple_res, simple_res)
  
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
  
  name2 = paste("top_inds_net", str_replace(i, ".json", ""), sep = "_")
  assign(name2, top_inds_net)
  
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
  
  aggregate_data = data.frame(gen=rep(0,length(levels(get(name2)$gen))),active_nodes=rep(0,length(levels(get(name2)$gen))), inactive_connections=rep(0,length(levels(get(name2)$gen))))
  count=0
  for(j in levels(get(name2)$gen)){
  count = count + 1
  #adding weights, weight sign, activation to the edge list
  edge_tibble = as_tibble(cbind(from, to))
  ind = filter(get(name2), gen == j)
  edge_tibble = cbind(edge_tibble, ind$m_weight, ind$m_is_active, ind$w_sign)
  
  node_tibble = as_tibble(cbind(id, layer, node))
  node_active = c(rep(TRUE,architecture[1]), subset(ind, ind$weight == "value_node_m_weights_1")$value_node_m_active)
  node_tibble = cbind(node_tibble, node_active)
  
  #Plot number of active nodes
  aggregate_data$active_nodes[count]=sum(node_tibble$node_active)
  aggregate_data$inactive_connections[count]=sum(edge_tibble$`ind$m_is_active`==F)
  aggregate_data$gen[count] = j
  
  # ###plot network####
  ##create igraph or ggraph object
  network_d <- igraph::graph_from_data_frame(d = edge_tibble,
                                             vertices = node_tibble,
                                             directed = T)
  
  E(network_d)$color = as.factor(edge_tibble$`ind$w_sign`)
  E(network_d)$weight = if_else(edge_tibble$`ind$m_is_active` == T,  edge_tibble$`ind$m_weight`, 0)
  network_d = network_d - E(network_d)[E(network_d)$weight == 0]
  V(network_d)$color = factor(node_tibble$node_active, levels=c("FALSE", "TRUE"))
  
  
  jpeg(paste("Plot",ind$mut_type, "s",ind$seed,"arch",ind$architecture,
             "cycle", paste(paste(rep("0",max(nchar(levels(get(name2)$gen)))-(nchar(j))), collapse=""),j, sep=""),"changefreq", 
             ind$change_freq,"duprate",ind$mut_rate_dup,"actrate", ind$mut_rate_act, ".png", sep = "_")
       ,width = 700,
       height = 700)
  
  title = paste("Generation ", ind$gen[1])
  
  plot(network_d, layout = l,
       edge.arrow.size = 0.5,                           # Arrow size, defaults to 1
       edge.arrow.width = 0.7,                          # Arrow width, defaults to 1
       edge.arrow.height = 0.9,                          # Arrow width, defaults to 1
       edge.lty = c("solid"),
       edge.width = abs(E(network_d)$weight/max(abs(E(network_d)$weight)) * 10), 
       main = title,
       vertex.label = NA
  )
  dev.off()
  
  }  
  ####Create gif
  ## list file names and read in
  imgs = intersect(intersect(intersect(intersect(intersect(list.files(pattern = "*png$", full.names = T), list.files(pattern = levels(get(name1)$seed)[1], full.names =  T)),
                                                 list.files(pattern = paste("duprate_",levels(get(name1)$mut_rate_dup)[1], sep=""), full.names =  T)),
                                       list.files(pattern = paste("changefreq_", levels(get(name1)$change_freq)[1], sep=""), full.names =  T)),
                             list.files(pattern = levels(as.factor(get(name1)$mut_type))[1], full.names =  T)),
                   list.files(pattern = paste("actrate_", levels(get(name1)$mut_rate_act)[1], sep=""), full.names =  T))
  img_list = lapply(imgs, image_read)
  
  ## join the images together
  img_joined <- image_join(img_list)
  
  ## animate at 2 frames per second
  img_animated <- image_animate(img_joined, fps = 2)
  
  ## save to disk
  path = paste("Gif",get(name1)$seed[1], get(name1)$mut_rate_act[1], get(name1)$mut_rate_dup[1], get(name1)$change_freq[1], get(name1)$architecture[1], ".gif", sep = "_")
  image_write(image = img_animated,
              path = path)
  }
  aggregate_data = cbind (aggregate_data,ID)
  all_aggregate_data=rbind(all_aggregate_data, aggregate_data)
  all_aggregate_data$gen = as.numeric(all_aggregate_data$gen)
  all_aggregate_data=mutate(all_aggregate_data, normalized_inactive = inactive_connections/active_nodes)

  

ggplot(data = all_aggregate_data
       #%>%filter(mut_type == "NRduplication")  
       #%>% slice_min(gen,n = 1000)
) +
  geom_line(aes(x = gen, y = active_nodes))+
  facet_grid(seed~.)

ggplot(data = all_aggregate_data
       #%>%filter(mut_type == "NRduplication")  
       #%>% slice_min(gen,n = 1000)
) +
  geom_line(aes(x = gen, y = normalized_inactive))+
  facet_grid(seed~.)
   

 