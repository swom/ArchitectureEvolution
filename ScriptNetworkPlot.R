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
library(patchwork)

dir = "C:/Users/p288427/Desktop/data_dollo_++/28+29/"
setwd(dir)

results=list()
# pattern = "*json$"
# pattern = 'mut_type_weights_start_arc1-2-2-2-1_act_r0.001000_dup_r0.000500_ch_A0.000000_ch_B0.010000_ch_typesymmetrical_ch_typeregular_sel_str2.0_max_arc1-2-2-2-1_sel_typesporadic_sel_freq1000_1'
# pattern = "mut_t_weights_sel_t_sporadic_sym_t_symmetrical_fr_t_regular_a_p_on_arc_1-2-2-2-1_m_arc_1-2-2-2-1_act_r_0.001_dup_r_0.000_ch_A_0.000_ch_B_0.010_s_st_2.0_s_f_100_seed1.json"
pattern = "mut_t_weights_sel_t_sporadic_sym_t_symmetrical_fr_t_regular_a_p_off_arc_1-2-2-2-1_m_arc_1-2-2-2-1_act_r_0.001_dup_r_0.000_ch_A_0.000_ch_B_0.010_s_st_0.5_s_f_1000_seed10"
list.files(path = '.', pattern = pattern)

for (i in  list.files(path = '.', pattern = pattern)){
  
  ###Making a data tibble with all top individuals' data 
  results <- fromJSON(file = i)
  
  results_unnest = as.data.frame(do.call(rbind,do.call(rbind, results$m_top_inds)))
  results_unnest$generation = do.call(rbind,results_unnest$generation)
 
  reac_norms = results_unnest %>% 
    select(c(generation, m_reac_norm))

  m_ind = as.data.frame(do.call(rbind, results_unnest$m_ind)) 
  results_df = results_unnest %>% 
    select(-c(m_ind, m_reac_norm)) %>%
    cbind(m_ind) %>%
    rename(network = "m_network")
  
  i = str_replace(i, "weights_and_activation", "weightsandactivation")
  
  results$m_params$i_p$net_par$max_arc = toString(results$m_params$i_p$net_par$max_arc)
  results$m_params$i_p$net_par$net_arc = toString(results$m_params$i_p$net_par$net_arc)
  ID = as.data.frame(results$m_params)
  
  name1 = paste("top_inds", str_replace(i, ".json", ""), sep = "_")
  assign(name1, cbind(results_df, ID))
  
  architecture = as.integer(strsplit(get(name1)$i_p.net_par.max_arc, ",")[[1]])
  
  
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
  top_inds_net$generation = as.factor(top_inds_net$generation)
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
  edge_tibble = as_tibble(cbind(from, to))
  
  ####Create empty data frames and global variables to be used in the plotting loop below
  mismatch = data.frame(generation = 0, mismatch = 0)
  sum_of_weights = data.frame(generation = 0, sum_of_weights = 0)
  highest_weights = data.frame(generation = 0, max_weigth = 0)
  top_ind_fit = data.frame(generation = 0, top_ind_fit = 0)
  o_rn = reac_norms %>% slice_max(gen, n = 1) %>% select(m_reac_norm)
  optimal_rn = as.data.frame(do.call(rbind,do.call(rbind,o_rn$m_reac_norm))) %>% mutate(m_y = as.numeric(m_x) * as.numeric(m_x))
  worst_rn = o_rn %>% mutate(m_y = -1)
  max_dist = as.numeric(dist(rbind(as.array(worst_rn$m_y), as.array(optimal_rn$m_y))))
  
  #Now let's loop through generations
  for(gen in levels(get(name2)$generation)){
    #adding weights, weight sign, activation to the edge list
    ind = filter(get(name2), generation == gen)
    edge_tibble_ind = cbind(edge_tibble, ind$m_weight, ind$m_is_active, ind$w_sign)
    
    node_active = c(TRUE, subset(ind, ind$weight == "value_node_m_weights_1")$value_node_m_active)
    node_tibble_ind = cbind(node_tibble, node_active)
    
    # ###plot network####
    ##create igraph or ggraph object
    network_d <- igraph::graph_from_data_frame(d = edge_tibble_ind,
                                               vertices = node_tibble_ind,
                                               directed = T)
    
    E(network_d)$color = as.factor(edge_tibble_ind$`ind$w_sign`)
    E(network_d)$weight = if_else(edge_tibble_ind$`ind$m_is_active` == T, 
                                  edge_tibble_ind$`ind$m_weight`, 
                                  0)
    network_d = network_d - E(network_d)[E(network_d)$weight == 0]
    V(network_d)$color = factor(node_tibble_ind$node_active, levels=c("FALSE", "TRUE"))
    
    ###create data fora plot for the highest value of the weights
    highest_weights = rbind(highest_weights, 
                            data.frame(generation = gen,
                                       max_weigth = max(abs(edge_tibble_ind$`ind$m_weight`))
                            )
    ) 
    ###create a plot for the highest value of the weights
    top_ind_fit = rbind(top_ind_fit, 
                            data.frame(generation = gen,
                                       top_ind_fit = ind$m_fitness[[1]]
                                       )
                            )
    ###plot reaction norm
    rn = reac_norms %>% filter(generation == as.numeric(gen)) %>% select(m_reac_norm)
    rn_d = as.data.frame(do.call(rbind,do.call(rbind,rn$m_reac_norm)))
    
    ###create data for plot that shows mismatch level between reaction norm and optimal funciton 
    mismatch = rbind(mismatch,
                     data.frame(generation = gen,
                                mismatch = 1 - as.numeric(dist(rbind(as.array(rn_d$m_y), as.array(optimal_rn$m_y))))/max_dist))
    ###create data for plot that shows the sum of the weights of the network
    sum_of_weights = rbind(sum_of_weights,
                           data.frame(generation = gen,
                                      sum_of_weights = sum(edge_tibble_ind$`ind$m_weight`)))

    jpeg(paste("Plot",
               ind$i_p.m_mutation_type,
               "s",ind$s_p.seed,
               "arch",ind$i_p.net_par.max_arc,
               "cycle", as.numeric(gen),
               "changefreq", ind$s_p.change_freq_A,
               ".png", sep = "_")
         ,width = 700,
         height = 700)
    
    title = paste("Generation ", ind$generation[1])
    
    par(mfrow=c(3,2))
  
      plot(highest_weights$generation, highest_weights$max_weigth, type = 'l',
         main = "max_weight")  
   
     plot(sum_of_weights$generation, sum_of_weights$sum_of_weights, type = 'l',
         main = "sum of weights")
  
      plot(top_ind_fit$generation, top_ind_fit$top_ind_fit, type = 'l',
         main = "top_ind_fit")
   
     plot(mismatch$generation, mismatch$mismatch, type = 'l',
         main = "Match")
  
      plot(rn_d$m_x, rn_d$m_y,
         xlim=c(-1,1),
         ylim=c(-1,1),
         type = 'l',
              main = "reac_norm") 
    lines(rn_d$m_x,optimal_rn$m_y,col="green")

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
  # 
  # ####Create gif
  # ## list file names and read in
  # imgs = intersect(intersect(intersect(intersect(intersect(list.files(pattern = "*png$", full.names = T), list.files(pattern = levels(get(name1)$architecture)[1], full.names =  T)),
  #                                                list.files(pattern = paste("duprate_",levels(get(name1)$mut_rate_dup)[1], sep=""), full.names =  T)),
  #                                      list.files(pattern = paste("changefreq_", levels(get(name1)$change_freq)[1], sep=""), full.names =  T)),
  #                            list.files(pattern = levels(as.factor(get(name1)$mut_type))[1], full.names =  T)),
  #                  list.files(pattern = paste("actrate_", levels(get(name1)$mut_rate_act)[1], sep=""), full.names =  T))
  # img_list = lapply(imgs, image_read)
  # 
  # ## join the images together
  # img_joined <- image_join(img_list)
  # 
  # ## animate at 2 frames per second
  # img_animated <- image_animate(img_joined, fps = 2)
  # 
  # ## save to disk
  # path = paste("Gif",get(name1)$seed[1], get(name1)$mut_rate_act[1], get(name1)$mut_rate_dup[1], get(name1)$change_freq[1], get(name1)$architecture[1], ".gif", sep = "_")
  # image_write(image = img_animated,
  #             path = path)
  
  
}

