library(rjson)
library(tidyr)
library(tidyverse)
library(ggplot2)
library(stringi)
library(rlist)
library(igraph)
library(tidygraph)
library(ggraph)
library(ggnetwork)
library(networkD3)
library(magick)
library(patchwork)
library(colorspace)
dir = "C:/Users/p288427/Desktop/data_dollo_++/6_16_22_sampled_short/"
setwd(dir)

results=list()
# pattern = "*json$"
pattern = "mut_t_weights_sel_t_spo_sym_t_sym_fr_t_reg_a_p_off_r_t_con_arc_1-2-2-2-1_m_arc_1-2-2-2-1_act_r_0.001_dup_r_0.000_ch_A_0.000_ch_B_0.010_s_st_1.0_s_f_100_seed1"
for (i in  list.files(path = '.', pattern = pattern)){
  
###Making a data tibble with all top individuals' data 
  results <- fromJSON(file = i)
  
  results_unnest_top_inds = as.data.frame(do.call(rbind,do.call(rbind, results$m_top_inds)))
  results_unnest_top_inds$generation = do.call(rbind,results_unnest_top_inds$generation)
  
  reac_norms = results_unnest_top_inds %>% 
    select(c(generation, m_reac_norm))
  
  m_ind = as.data.frame(do.call(rbind, results_unnest_top_inds$m_ind)) 
  results_df_top_inds = results_unnest_top_inds %>% 
    select(-c(m_ind, m_reac_norm)) %>%
    cbind(m_ind) %>%
    rename(network = "m_network")
  
  i = str_replace(i, "weights_and_activation", "weightsandactivation")
  
  results$m_params$i_p$net_par$max_arc = toString(results$m_params$i_p$net_par$max_arc)
  results$m_params$i_p$net_par$net_arc = toString(results$m_params$i_p$net_par$net_arc)
  ID = as.data.frame(results$m_params)
  
  name1 = paste("top_inds", str_replace(i, ".json", ""), sep = "_")
  assign(name1, cbind(results_df_top_inds, ID))
  
  architecture = as.integer(strsplit(get(name1)$i_p.net_par.max_arc, ",")[[1]])
  
  #if response type is plastic add one extra node in the input layer to architecture that takes the environmental function as input
  if(length(ID$i_p.net_par.resp_type)){
    if(ID$i_p.net_par.resp_type == 0){
      architecture[1] = architecture[1] + 1
    }
  }
  
  
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
  
  layout = as.data.frame(cbind(x, y)) 
  
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
  match_a = data.frame(generation = 0, match = 0)
  match_b = data.frame(generation = 0, match = 0)
  
  
  sum_of_weights = data.frame(generation = 0, sum_of_weights = 0)
  
  
  
  highest_weights = data.frame(generation = 0, max_weigth = 0)
  
  phenotype_robustness = data.frame(generation = 0, phen_robust = 0)
  
  fitness_robustness = data.frame(generation = 0, fit_robust = 0)
  
  top_ind_fit = data.frame(generation = 0, top_ind_fit = 0)
  
  o_rn = reac_norms %>% slice_max(generation, n = 1) %>% select(m_reac_norm)
  optimal_rn_func_a = as.data.frame(do.call(rbind,do.call(rbind,o_rn$m_reac_norm))) %>% mutate(m_y = as.numeric(m_x) * as.numeric(m_x))
  optimal_rn_func_b = as.data.frame(do.call(rbind,do.call(rbind,o_rn$m_reac_norm))) %>% mutate(m_y = as.numeric(m_x) * as.numeric(m_x) * as.numeric(m_x))
  worst_rn_a = optimal_rn_func_a %>% mutate(m_y = -1)
  worst_rn_b = optimal_rn_func_b %>% mutate(m_y = ifelse(m_x <= 0,1,-1))
  
  inputs = data.frame(generation = c(), input= c())
  outputs = data.frame(generation = c(), outputs= c())
  
  #create directory where to save images
  subdir = paste(
    "s",results$m_params$s_p$seed,
    "change_freq", results$m_params$s_p$change_freq_A,
    "s_f", results$m_params$s_p$selection_freq,
    "s_s", results$m_params$s_p$selection_strength,
    "a_p", results$m_params$s_p$adaptation_per,
    sep = "_")
  dir.create(file.path(dir, subdir), showWarnings = FALSE)
  
  #Now let's loop through generations
  for(gen in levels(get(name2)$generation)){
    #adding weights, weight sign, activation to the edge list
    ind = filter(get(name2), generation == gen)
    edge_tibble_ind = cbind(edge_tibble, ind$m_weight, ind$m_is_active, ind$w_sign)
    
    node_active = c(TRUE, subset(ind, ind$weight == "value_node_m_weights_1")$value_node_m_active)
    
    #if response type is plastic add one T value for the extra node in the input layer (always active) that takes the environmental function as input
    if(length(ID$i_p.net_par.resp_type)){
      if(ID$i_p.net_par.resp_type == 0)
      {
        node_active = append(node_active, 1, TRUE)
      }
    }
    
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
    
    ###you can apply layout in ggraph!!!!
    # network_d_t = tidygraph::as_tbl_graph(network_d) %>%
    #   mutate(x = layout$x) %>%
    #   mutate(y = layout$y)
    # 
    #   n_colors = 100
    #   rbg<- colorRampPalette(c("red", "blue", "green"))
    #      
    #   ggraph(network_d_t, x = x, y = y)+
    #   geom_node_point(
    #                   aes(size = 16),
    #                   show.legend = F) +
    #   geom_edge_link(arrow = grid::arrow(),
    #                  aes(edge_width = abs(`ind$m_weight`/max(abs(`ind$m_weight`))),
    #                      edge_colour =`ind$w_sign`,
    #                      # label = as.character(weight),
    #                      start_cap = circle(),
    #                      end_cap = circle()),
    #                  show.legend = F) +  
    #     scale_edge_colour_manual(values =  rbg(2)) +
    #   ggtitle(as.character(gen))
    
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
    
    ###create data for plot that shows match level between reaction norm and optimal funciton A 
    max_dist_a = as.numeric(dist(rbind(as.array(worst_rn_a$m_y), as.array(optimal_rn_func_a$m_y))))
    match_a = rbind(match_a,
                    data.frame(generation = gen,
                               match = 1 - as.numeric(dist(rbind(as.array(rn_d$m_y), as.array(optimal_rn_func_a$m_y))))/max_dist_a))
    
    ###create data for plot that shows match level between reaction norm and optimal funciton B
    max_dist_b = as.numeric(dist(rbind(as.array(worst_rn_b$m_y), as.array(optimal_rn_func_b$m_y))))
    match_b = rbind(match_b,
                    data.frame(generation = gen,
                               match = 1 - as.numeric(dist(rbind(as.array(rn_d$m_y), as.array(optimal_rn_func_b$m_y))))/max_dist_b))
    ###create data for plot that shows the sum of the weights of the network
    sum_of_weights = rbind(sum_of_weights,
                           data.frame(generation = gen,
                                      sum_of_weights = sum(edge_tibble_ind$`ind$m_weight`)))
    
    phenotype_robustness = rbind(phenotype_robustness,
                                 data.frame(generation = gen,
                                            phen_robust = ind$m_sensibilities[[1]]$m_phenotype))
    fitness_robustness = rbind(fitness_robustness,
                                 data.frame(generation = gen,
                                            fit_robust = ind$m_sensibilities[[1]]$m_fitness))
    #### create data to plot inputs on reaction norm plot
    # inputs = do.call(rbind,do.call(cbind,results$m_input[as.numeric(gen)]))[,1]
    # outputs = do.call(cbind,results$m_optimal[as.numeric(gen)])
    
    
    jpeg(paste(subdir,paste("s",ind$s_p.seed,
                            "arch",ind$i_p.net_par.max_arc,
                            "cycle", as.numeric(gen),
                            "change_freq", ind$s_p.change_freq_A,
                            "s_f", ind$s_p.selection_freq,
                            "s_s", ind$s_p.selection_strength,
                            "a_p", ind$s_p.adaptation_per,
                            ".png", sep = "_"),sep = "/")
         ,width = 2000,
         height = 800)
    
    title = paste("Generation ", ind$generation[1])
    
    par(mfrow=c(2,5))
    
    plot(highest_weights$generation, highest_weights$max_weigth, type = 'l',
         main = "max_weight")  
    
    plot(sum_of_weights$generation, sum_of_weights$sum_of_weights, type = 'l',
         main = "sum of weights")
    
    plot(match_a$generation, match_a$match, type = 'l',
         main = "Match_a") 
    
    plot(match_b$generation, match_b$match, type = 'l',
         main = "Match_b")
      
      #generation starts from 0 but arrays in R start from 1 so add 1 to gen
    if(results$m_env_functions[as.numeric(gen) + 1] == 65){
      plot(rn_d$m_x, rn_d$m_y,
           xlim=c(-1,1),
           ylim=c(-1,1),
           type = 'l',
           main = "reac_norm") 
      clip(-1, 1, -1, 1)
      lines(rn_d$m_x,optimal_rn_func_a$m_y,col="green")
      if(nrow(inputs) > 0){
        abline(v=inputs, col = 'blue', ylim=c(-1,1))
        points(inputs, outputs)
      }
    } else if(results$m_env_functions[as.numeric(gen) + 1] == 66){    #generation starts from 0 but arrays in R start from 1 so add 1 to gen
      
      plot(rn_d$m_x, rn_d$m_y,
           xlim=c(-1,1),
           ylim=c(-1,1),
           type = 'l',
           main = "reac_norm") 
      clip(-1, 1, -1, 1)
      lines(rn_d$m_x,optimal_rn_func_b$m_y,col="red")
      if(nrow(inputs) > 0){
        abline(v=inputs, col = 'blue', ylim=c(-1,1))
        points(inputs, outputs)
      }
    }     
    
    plot(phenotype_robustness$generation, phenotype_robustness$phen_robust,
         type = 'l',
         main = "Phenotypic robustness")
   
    plot(fitness_robustness$generation, fitness_robustness$fit_robust,
         type = 'l',
         main = "Fitness robustness")
 
    
    plot(top_ind_fit$generation, top_ind_fit$top_ind_fit, type = 'l',
         main = "top_ind_fit")
    
    plot(network_d, layout = as.matrix(layout),
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

min_gen = 800000
max_gen = 900000
mismatch_sub = mismatch %>% filter(as.numeric(generation) > min_gen) %>% filter(as.numeric(generation) < max_gen)
top_ind_fit_sub = top_ind_fit %>% filter(as.numeric(generation) > min_gen) %>% filter(as.numeric(generation) < max_gen)
par(mfrow=c(1,1))
plot(top_ind_fit_sub$generation,
     top_ind_fit_sub$top_ind_fit,
     type = 'l',
     col = 'red',
     main = "match(black) and fitness(red)")
lines(mismatch_sub$generation, mismatch_sub$mismatch)
