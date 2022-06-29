library(rjson)
library(tidyr)
library(tidyverse)
library(ggplot2)
library(stringi)
library(data.table)
library(igraph)
library(tidygraph)
library(ggraph)
library(patchwork)
library(RColorBrewer)

dir = "C:/Users/p288427/Desktop/data_dollo_++/6_28_22_sampled_short_good/"
setwd(dir)

results=list()
pattern = '*json$'
# pattern = "mut_t_weights_sel_t_spo_sym_t_sym_fr_t_reg_a_p_off_r_t_con_arc_1-2-2-2-1_m_arc_1-2-2-2-1_act_r_0.001_dup_r_0.000_ch_A_0.000_ch_B_0.010_s_st_1.0_s_f_100_seed1"
for (i in  list.files(path = '.', pattern = pattern)){
  
  ###Making a data tibble with all top individuals' data 
  results <- fromJSON(file = i)
  
  #extract the sampled individuals
  #individuals in the same generations 
  #are ranked by row number 
  #the lower the row number -> the higher robustness -> the higher the rank
  
  results_unnest_sampled_inds = as.data.frame(do.call(rbind,do.call(rbind, results$m_sampled_inds))) %>% 
    mutate(generation = as.numeric(generation)) %>%
    group_by(generation) %>% 
    mutate(rank = rank(row_number())) %>% 
    ungroup()
  
  #extract the reaction norms of the sampled individuals
  reac_norms = results_unnest_sampled_inds %>% 
    select(c(generation, m_reac_norm, rank))
  
  #extract the sensibilities of the sampled individuals
  sensibilities = results_unnest_sampled_inds %>% 
    select(c(generation, m_sensibilities, rank))
  
  # extract the fintess rank from sensibilities
  ranks = sensibilities %>% 
    mutate(m_sensibilities = as.data.frame(rbindlist(m_sensibilities))) %>% 
    mutate(m_sensibilities = m_sensibilities$m_rank) %>% 
    rename(fit_rank = "m_sensibilities")
  
  #extract the network
  m_ind = as.data.frame(do.call(rbind, results_unnest_sampled_inds$m_ind))
  
  results_df_top_inds = results_unnest_sampled_inds %>% 
    select(-c(m_ind, m_reac_norm, m_sensibilities)) %>%
    cbind(m_ind) %>%
    rename(network = "m_network") %>% 
    left_join(ranks)
  
  #extract simulation parameters into data frame with arcs as strings
  results$m_params$i_p$net_par$max_arc = toString(results$m_params$i_p$net_par$max_arc)
  results$m_params$i_p$net_par$net_arc = toString(results$m_params$i_p$net_par$net_arc)
  ID = as.data.frame(results$m_params)
  
  #assign name to this specific result
  #result is combined data frame of parameters and networks
  i = str_replace(i, "weights_and_activation", "weightsandactivation")
  name1 = paste("top_inds", str_replace(i, ".json", ""), sep = "_")
  assign(name1, cbind(results_df_top_inds, ID))
  
  #create data frame that represents max architecture of network
  #each column is the number of nodes per layer
  architecture = as.integer(strsplit(get(name1)$i_p.net_par.max_arc, ",")[[1]])
  
  #if response type is plastic add one extra node in the input layer to architecture that takes the environmental function as input
  if(length(ID$i_p.net_par.resp_type)){
    if(ID$i_p.net_par.resp_type == 0){
      architecture[1] = architecture[1] + 1
    }
  }
  
  
  ###Keeping only networks architecture, expanding to have each connection as a row
  
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
  
  #create directory where to save images
  subdir = paste(
    "s",results$m_params$s_p$seed,
    "change_freq", results$m_params$s_p$change_freq_A,
    "s_f", results$m_params$s_p$selection_freq,
    "s_s", results$m_params$s_p$selection_strength,
    "a_p", results$m_params$s_p$adaptation_per,
    sep = "_")
  dir.create(file.path(dir, subdir), showWarnings = FALSE)
  
  #make plot of avg_fitness over time
  avg_fitnesses = data.frame(avg_fitness = results$m_avg_fitnesses,
                             var_fitness = results$m_var_fitnesses) %>%
    rownames_to_column() %>% 
    rename(generation = rowname) %>% 
    mutate(generation = as.numeric(generation)) %>% 
    filter(generation %in% get(name2)$generation)
  
  avg_fit_plot = ggplot(data = avg_fitnesses, 
                        aes(x = generation,
                            y = avg_fitness)) +
    geom_line() +
    geom_ribbon(aes(x = generation,
                     y = avg_fitness,
                     ymax = avg_fitness + var_fitness,
                     ymin = avg_fitness - var_fitness),
                 alpha = 0.5)+
    xlab("Generations")
 
  all_sens_summary =  results$m_fit_phen_mut_sensibility %>%
    rbindlist() %>% 
    cbind(as.data.frame(rbindlist(.$m_sensibilities),
                        col.names = names(m_sensibilities))) %>% 
    select(-c(m_sensibilities)) %>% 
    group_by(m_generation) %>% 
    summarise(mean_phen_sens = mean(m_phenotype_sens),
              mean_fit_sens = mean(m_fitness_sens),
              var_phen_sens = sd(m_phenotype_sens),
              var_fit_sens = sd(m_fitness_sens))
  
  phen_sens_plot = ggplot(data = all_sens_summary,
                          aes(x = m_generation,
                              y = mean_phen_sens)) +
    geom_line() +
    geom_ribbon(aes(x = m_generation,
                    y = mean_phen_sens,
                    ymax = mean_phen_sens + var_phen_sens,
                    ymin = mean_phen_sens - var_phen_sens),
                alpha = 0.5) +
    xlab("Generations")
  
  fit_sens_plot = ggplot(data = all_sens_summary,
                         aes(x = m_generation,
                             y = mean_fit_sens)) +
    geom_line() +
    geom_ribbon(aes(x = m_generation,
                    y = mean_fit_sens,
                    ymax = mean_fit_sens + var_fit_sens,
                    ymin = mean_fit_sens - var_fit_sens),
                alpha = 0.5)+
    xlab("Generations")
  
  
  #Now let's loop through generations
  for(gen in levels(get(name2)$generation)){
    
    #adding weights, weight sign, activation to the edge list
    ind = filter(get(name2), generation == gen)
    edge_tibble_ind = cbind(edge_tibble, ind$m_weight, ind$m_is_active, ind$w_sign, ind$rank)
    
    #finding the highest weight value among the three networks
    max_weight = max(abs(edge_tibble_ind$`ind$m_weight`))
    
    #since each weight is a row if we need to extract node 
    #it suffices to look at the first weight of a node "value_node_m_weights_1"
    #to extract that property
    #extracting if node is active and its bias
    nodes_properties = ind %>% 
      filter(ind$weight == "value_node_m_weights_1") %>% 
      select(c(rank, value_node_m_active, value_node_m_bias)) %>% 
      #add node for input 
      group_by(rank) %>%
      group_modify(~ add_row(.x,
                             .before=0)) %>% 
      mutate(value_node_m_active =  replace(value_node_m_active, 1, TRUE)) %>% 
      mutate(value_node_m_bias =  replace(value_node_m_bias, 1, 0))  
    
    
    
    #if response type is plastic add one T value for the extra node in the input layer (always active) that takes the environmental function as input
    if(length(ID$i_p.net_par.resp_type)){
      if(ID$i_p.net_par.resp_type == 0)
      {
        nodes_properties = nodes_properties %>% 
          group_by(rank) %>% 
          group_modify(~ add_row(.x, .before=0)) %>% 
          mutate(value_node_m_active =  replace(value_node_m_active, 1, TRUE)) %>% 
          mutate(value_node_m_bias =  replace(value_node_m_bias, 1, 0))  
      }
    }
    
    node_tibble_ind = cbind(node_tibble, nodes_properties)
    
    # ###plot networks####
    # create color palette for networks
    n_colors = 100
    rbg<- colorRampPalette(c("red", "blue", "green"))
    
    #create limits for plot of sensibilities
    fit_x_lim = c(-0.1,0.1)
    phen_x_lim = c(0,0.5)
    
    #create color palette for sensibilities
    myPalette <- colorRampPalette(rev(brewer.pal(11, "Spectral")))
    
    # create list where to store plots
    plot_list <- list()
    for(ind_rank in unique(edge_tibble_ind$`ind$rank`))
    {
      ##create igraph object for each
      network_d <- igraph::graph_from_data_frame(d = edge_tibble_ind %>%
                                                   filter(ind$rank == ind_rank),
                                                 vertices = node_tibble_ind %>% 
                                                   filter(rank == ind_rank),
                                                 directed = T)
      
      E(network_d)$color = as.factor(edge_tibble_ind$`ind$w_sign`)
      E(network_d)$weight = if_else(edge_tibble_ind$`ind$m_is_active` == T, 
                                    edge_tibble_ind$`ind$m_weight`, 
                                    0)
      network_d = network_d - E(network_d)[E(network_d)$weight == 0]
      V(network_d)$color = factor(node_tibble_ind$node_active, levels=c("FALSE", "TRUE"))
      
      ###you can apply layout in ggraph!!!!
      network_d_t = tidygraph::as_tbl_graph(network_d) %>%
        mutate(x = layout$x) %>%
        mutate(y = layout$y)
      
      p <- ggraph(network_d_t, x = x, y = y) +
        geom_edge_link(arrow = grid::arrow(),
                       aes(edge_width = abs(`ind$m_weight`/ max_weight),
                           edge_colour =`ind$w_sign`
                       ),
                       show.legend = F) +
        geom_node_point(shape = 21,
                        aes(size = value_node_m_bias,
                            fill = value_node_m_bias > 0),
                        show.legend = F) +
        scale_edge_colour_manual(values =  rbg(2)) +
        ggtitle(paste("gen",as.character(gen),
                      "sens_rank",as.character(ind_rank),
                      "fit_rank", as.character(ind$fit_rank[ind$rank == ind_rank])))
      
      #add plot to plot_list
      plot_list <- c(plot_list, list(p))    
    }
    # dispose on a row plot_list with patchwork
    nets <- patchwork::wrap_plots(plot_list, nrow=1)
    
    #find sensibilities of all individuals in that generation
    all_sens_gen = results$m_fit_phen_mut_sensibility %>%
      subset(m_generation == gen) %>%
      rbindlist() %>% 
      cbind(as.data.frame(rbindlist(.$m_sensibilities),
                          col.names = names(m_sensibilities))) %>% 
      select(-c(m_sensibilities))
    
    #plot sensibilities of sampled individuals
    plot_sens <- ggplot(sensibilities %>% 
                          filter(generation == as.numeric(gen)) %>%
                          mutate(sensibilities = rbindlist(m_sensibilities))) +
      geom_point(data = all_sens_gen,
                 shape = 21, 
                 aes(x = m_fitness_sens,
                     y = m_phenotype_sens,
                     # fill = m_rank,
                     colour = m_fitness
                 ), alpha = 0.5) +
      geom_point(shape = 21, size = 8,
                 aes(x = sensibilities$m_fitness_sens,
                     y = sensibilities$m_phenotype_sens,
                     fill = sensibilities$m_fitness
                 ), show.legend = F) +
      geom_text(aes(x = sensibilities$m_fitness_sens,
                    y = sensibilities$m_phenotype_sens,
                    label = rank), 
                hjust = 0.1, 
                nudge_x = -0.008) +   
      geom_text(aes(x = sensibilities$m_fitness_sens,
                    y = sensibilities$m_phenotype_sens,
                    label = paste(sensibilities$m_rank,":",round(sensibilities$m_fitness,4))), 
                hjust = -0.05, 
                nudge_x = 0.003) +
      xlim(fit_x_lim) +
      ylim(phen_x_lim) +
      scale_fill_gradientn(colours = myPalette(1000), limits=c(0,1)) +
      scale_colour_gradientn(colours = myPalette(1000), limits=c(0,1))
    
    avg_time_marked = avg_fit_plot + geom_point(aes(x = as.numeric(gen), y = avg_fitness[generation == gen]), colour ="red")
phen_time_marked = phen_sens_plot + geom_point(aes(x = as.numeric(gen), y = mean_phen_sens[m_generation == gen]), colour ="red")
fit_time_marked = fit_sens_plot + geom_point(aes(x = as.numeric(gen), y = mean_fit_sens[m_generation == gen]), colour ="red")

    nets / (plot_sens + (avg_time_marked / phen_time_marked / fit_time_marked)
              )
    
    ggsave(paste(subdir,paste(paste("sampled_nets_fit_sens_plot", gen, sep = "_"),".png"), sep = '/'),
           device = "png", 
           width = 30,
           height = 15)
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

# min_gen = 800000
# max_gen = 900000
# mismatch_sub = mismatch %>% filter(as.numeric(generation) > min_gen) %>% filter(as.numeric(generation) < max_gen)
# top_ind_fit_sub = top_ind_fit %>% filter(as.numeric(generation) > min_gen) %>% filter(as.numeric(generation) < max_gen)
# par(mfrow=c(1,1))
# plot(top_ind_fit_sub$generation,
#      top_ind_fit_sub$top_ind_fit,
#      type = 'l',
#      col = 'red',
#      main = "match(black) and fitness(red)")
# lines(mismatch_sub$generation, mismatch_sub$mismatch)
