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

dir = "C:/Users/Clem/build-arc_evo-Desktop_Qt_6_1_0_MinGW_64_bit-Release"
setwd(dir)

results=list()
pattern = "*json$"
for (i in  list.files(path = '.', pattern = pattern)){
  ###Making a data tibble with all top individuals' data 
  
  results <- fromJSON(file = i)
  names(results$m_top_inds) = seq(from=0, by=1000, length.out=length(results$m_top_inds))
  results_df = results$m_top_inds %>% 
  unlist(recursive=FALSE)%>%
  as_tibble()%>%
  add_column(var = c("fitness", "input_values", "network"))%>%
  gather(key = gen, value = value, 1:length(results$m_top_inds)) %>% 
  pivot_wider(names_from = var, values_from = value)
  
  ID = data.frame(i) %>% 
    separate(i, c("architecture", "seed"), sep = '_')%>% 
    separate(seed, c("seed",NA))
  
  ID$architecture = as.factor(ID$architecture)
  ID$seed = as.factor(ID$seed)
  
  name1 = paste("top_inds", str_replace(i, ".json", ""), sep = "_")
  assign(name1, cbind(results_df, ID))
  
  ###Making a number vector out of the architecture
  architecture = strsplit(levels(get(name1)$architecture)[1], "-")[[1]]%>%
    as.integer()
  
  
  ###Keeping only network architecture, expanding to have each connection as a row
  
  top_inds_net = unnest_wider(get(name1), col = "network")%>%
    unnest_wider(col = "m_network_weights", names_sep = "_layer_")%>%
    mutate(m_input_size = NULL, input_values = NULL, fitness = NULL)%>%
    pivot_longer(cols = sprintf("m_network_weights_layer_%s", seq(1:(length(architecture)-1))), names_to = "layer")%>%
    unnest_wider(col = "value", names_sep = "_node_")%>%
    pivot_longer(cols = sprintf("value_node_%s", seq(1:(max(architecture)))), names_to = "node")%>%
    drop_na()%>%
    unnest_wider(col = "value", names_sep = "_weight_")%>%
    pivot_longer(cols = sprintf("value_weight_%s", seq(1:(max(architecture)))), names_to = "weight")%>%
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
  for(i in 1:length(architecture)){
    x = c(x, rep(i, architecture[i]))
    y = c(y, seq(1:architecture[i]) + ((max(architecture)-architecture[i])/2))
  }
  
  l =   cbind(x, y)
  
  #Make the list of nodes
  layer = vector()
  node = vector()
  for(i in 1:length(architecture)){
    layer = c(layer, rep(i, architecture[i]))
    node = c(node, seq(1:architecture[i]))
  }
  id = seq(1:length(layer))
  
  node_tibble = as_tibble(cbind(id, layer, node))
  
  #Make the list of edges
  
  from = vector()
  to = vector()
  
  for(i in 1:nrow(node_tibble)){
    from = c(from, rep(node_tibble$id[i], nrow(filter(node_tibble, layer == (node_tibble$layer[i]+1)))))
    to = c(to, filter(node_tibble, layer == (node_tibble$layer[i]+1))$id)
  }
  
  #Now let's loop through generations
  
  for(j in levels(get(name2)$gen)){
  
  #adding weights, weight sign, activation to the edge list
  edge_tibble = as_tibble(cbind(from, to))
  ind = filter(get(name2), gen == j)
  edge_tibble = cbind(edge_tibble, ind$m_weight, ind$m_is_active, ind$w_sign)
  
  # ###plot network####
  ##create igraph or ggraph object
  network_d <- igraph::graph_from_data_frame(d = edge_tibble,
                                             vertices = node_tibble,
                                             directed = T)
  
  E(network_d)$color = as.factor(edge_tibble$`ind$w_sign`)
  E(network_d)$weight = if_else(edge_tibble$`ind$m_is_active` == T,  edge_tibble$`ind$m_weight`, 0)
  network_d = network_d - E(network_d)[E(network_d)$weight == 0]
  

  jpeg(paste("Plot","s",ind$seed,"arch",ind$architecture,
             "cycle", paste(rep("0",max(nchar(levels(ind$gen)))-(nchar(j))),j, sep=""),".png", sep = "_")
       ,width = 700,
       height = 700)
  
  title = paste("Generation ", ind$gen[1])
  
  plot(network_d, layout = l,
       edge.arrow.size = 0.5,                           # Arrow size, defaults to 1
       edge.arrow.width = 0.7,                          # Arrow width, defaults to 1
       edge.arrow.height = 0.9,                          # Arrow width, defaults to 1
       edge.lty = c("solid"),
       edge.width = abs(E(network_d)$weight/max(E(network_d)$weight) * 10), 
       main = title)
  dev.off()
  
  }  

  ####Create gif
  ## list file names and read in
  imgs = intersect(list.files(pattern = "*png$", full.names = T), list.files(pattern = levels(get(name1)$architecture)[1], full.names =  T)) 
  img_list = lapply(imgs, image_read)
  
  ## join the images together
  img_joined <- image_join(img_list)
  
  ## animate at 2 frames per second
  img_animated <- image_animate(img_joined, fps = 2)
  
  ## save to disk
  path = paste("Gif",get(name1)$seed[1], get(name1)$architecture[1], ".gif", sep = "_")
  image_write(image = img_animated,
              path = path)
  

}
 