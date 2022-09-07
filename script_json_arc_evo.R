library(rjson)
library(ggplot2)
library(tibble)
library(dplyr)
library(tidyr)
library(stringr)
library(ggpubr)
library(patchwork)
library(data.table)
library(RColorBrewer)
library(tidygraph)
library(ggraph)

#remove scientific notation 
options(scipen=999)
#declare the 4 different optimal functions
func_1 <- function(x){return(x^2)}
func_2 <- function(x){return(x^3)}
func_3 <- function(x){return(3 * x - 4 * x^3)}
func_4 <- function(x){return(-0.8 + 9.32 * x^2- 15.84 * x^4 + 6.33 * x^6)}

#create optimal reaction norm based on name of function from data
#and range of the reaction norm of individuals
produce_current_optimal_func <- function(func_name, reac_norm){
  optimal_rn = reac_norm
  if(func_name == "1"){
    for(x in optimal_rn$x){
      optimal_rn$y[x == optimal_rn$x] = func_1(x)
    }
  } else if(func_name == "2"){
    for(x in optimal_rn$x){
      optimal_rn$y[x == optimal_rn$x] = func_2(x)
    }
  } else if(func_name == "3"){
    for(x in optimal_rn$x){
      optimal_rn$y[x == optimal_rn$x] = func_3(x)
    }
  } else if(func_name == "4"){
    for(x in optimal_rn$x){
      optimal_rn$y[x == optimal_rn$x] = func_4(x)
    }
  }
  return(optimal_rn)
}



####read data####
dir = "C:/Users/p288427/Desktop/data_dollo_++/9_6_22/1-2-1"

setwd(dir)

pattern = '^m.*json$'

####save load####
if(file.exists("all_simple_res.Rds") && 
   file.exists("all_sensibilities.Rds") &&
   file.exists("obs_params.Rds") &&
   file.exists("all_inds_rns.Rds")){
  obs_params <- readRDS("obs_params.Rds")
  all_simple_res <- readRDS("all_simple_res.Rds")
  all_sensibilities <- readRDS("all_sensibilities.Rds")
  all_inds_rns <- readRDS("all_inds_rns.Rds")
}else{
  
  filepaths = list.files(pattern = pattern)
  
  all_sensibilities = list()
  all_simple_res = data.frame()
  all_inds_rns = list()
  for (i in  filepaths)
  {
    
    results <- fromJSON(file = i)
    
    ###simple results
    simple_res = rowid_to_column(as_tibble(results[c("m_avg_fitnesses",
                                                     "m_env_functions",
                                                     "m_var_fitnesses")]),
                                 var = "gen")
    results$m_params$i_p$net_par$max_arc = toString(results$m_params$i_p$net_par$max_arc)
    results$m_params$i_p$net_par$net_arc = toString(results$m_params$i_p$net_par$net_arc)
    ID = as.data.frame(results$m_params)
    simple_res_ID = cbind(simple_res, ID)
    all_simple_res = rbind(all_simple_res, simple_res_ID)
    
    ###sensibilities results
    # m fit suitability
    tmp_ = results$m_fit_phen_mut_sensibility
    # convert to data.table with list column
    tmp_ = lapply(tmp_, as.data.table)
    # deal with all gens
    tmp_ = lapply(tmp_, function(df) {
      df = df[, rbindlist(m_sensibilities), by = "m_generation"]
    })
    tmp_ = rbindlist(tmp_)
    tmp_[, names(ID) := ID]
    all_sensibilities[[i]] = tmp_
    
    ###all_reaction norms
    # Keep only what is necessary
    data <- results$m_all_inds_rn
    
    # Extract the generations
    gens <- map_dbl(data, ~ .x$generation)
    names(data) <- gens
    
    # Extract x-values and individuals
    xvalues <- map_dbl(data[[1]]$m_reac_norm[[1]], ~ .x$m_x)
    inds <- seq(data[[1]]$m_reac_norm)
    
    # Make combinations of generation and individual...
    newdata <- expand_grid(gen = seq(gens), ind = inds) %>%
      mutate(
        
        # ... and for each combination ...
        newdata = map2(gen, ind, function(gen, ind) {
          
          # Extract y-values and assemble them with x-values
          tibble(
            x = xvalues,
            y = map_dbl(data[[gen]]$m_reac_norm[[ind]], last)
          )
          
        })
      ) %>%
      unnest(newdata)
    
    # Prepare the renaming of generations
    gens_labs <- unname(gens)
    gens <- as.character(seq(gens))
    names(gens) <- gens_labs
    
    # Rename generations
    newdata <- newdata %>%
      mutate(
        gen = fct_recode(as.character(gen), !!!gens),
        gen = as.numeric(as.character(gen))
      )
    
    # Pivot wider if needed
    newdata %>%
      pivot_wider(names_from = "ind", values_from = "y")
   
    
    # ###compact version
    #  tmp_ = results$m_all_inds_rn
    # tmp_ = lapply(tmp_, as.data.table)
    # tmp_ = rbindlist(
    #   lapply(tmp_, function(df){
    #   df = df[, rbindlist(lapply(df$m_reac_norm, rbindlist)), by = "generation"] 
    #   })
    #   )
    # tmp_[, names(ID) := ID]
    # all_inds_rns[[i]] = tmp_
    
    gc()
  }
  
  all_sensibilities = rbindlist(all_sensibilities, fill =T) %>% 
    #add 1 to all generations to sync with all_simple_res
    mutate(m_generation = m_generation + 1)
  # all_inds_rns = rbindlist(all_inds_rns, fill = T) %>% 
  #   rename(m_generation = generation) %>% 
  #   mutate(m_generation = m_generation + 1)
  
  obs_params = results$m_obs_param
  # all_sensibilities$m_generation = as.factor(all_sensibilities$m_generation)
  saveRDS(obs_params, file = "obs_params.Rds")
  saveRDS(all_simple_res, file = "all_simple_res.Rds")
  saveRDS(all_sensibilities, file = "all_sensibilities.Rds")
  saveRDS(all_inds_rns, file = "all_inds_rns.Rds")
  gc()
}
### Plot ####
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
  # geom_rect(data = . %>%
  #             group_by(s_p.change_freq_B, s_p.adaptation_per) %>%
  #             distinct(gen, m_env_functions) %>%
  #             filter(gen == min(gen) |
  #                      gen == max(gen) |
  #                      m_env_functions != lag(m_env_functions)) %>%
  #             mutate(gen_max = lead(gen)) %>%
  #             mutate(gen_max = ifelse(is.na(gen_max),
  #                                     max(gen),
  #                                     gen_max)),
#           aes(xmin = gen, xmax = gen_max,
#               ymin = 0, ymax = 1,
#               fill = as.factor(m_env_functions),
#           )
# ) +
geom_line(data = . %>% filter(gen %% filter_gen == 0),
          aes(x = gen, y = m_avg_fitnesses)
) +
  # geom_smooth(method='lm',aes(x = gen, y = m_avg_fitnesses))+
  # stat_regline_equation(aes(label = ..eq.label.., x = gen, y = m_avg_fitnesses),label.y.npc = 0.9) +
  # stat_regline_equation(aes(label = ..rr.label.., x = gen, y = m_avg_fitnesses), label.x.npc = 0.55,label.y.npc = 0.9)+
  # facet_grid(s_p.change_freq_B + s_p.selection_strength ~
  #              s_p.seed + s_p.adaptation_per + i_p.net_par.resp_type)
  facet_grid(s_p.selection_freq + 
               s_p.selection_strength + 
               e_p.name_func_A +
               i_p.net_par.net_arc ~
               s_p.seed + s_p.adaptation_per)

print(p)
dev.off()

########Sensibilities

#creating palette for plots
fitness_palette <- colorRampPalette(rev(brewer.pal(11, "Spectral")))
fitness_gradient <- scale_colour_gradientn(colours = fitness_palette(1000), limits=c(0,1))
fitness_sens_palette <- colorRampPalette(rev(brewer.pal(9, "Blues")))
fitness_sensitivity_palette <- scale_colour_gradientn(colours = fitness_sens_palette(10000),
                                                      limits=c(-0.2, 0.2))
outdegree_gradient <- scale_fill_gradientn(colours = fitness_palette(30))

#setting global variables outside the loop
phen_x_lim = c(0,0.3)
fit_x_lim = c(-0.1,0.1)
y_lim = c(0,50)
n_bins = 1000


sens_summary_all = all_sensibilities %>%
  group_by(m_generation,
           s_p.selection_freq,  
           s_p.selection_strength,  
           e_p.name_func_A, 
           i_p.m_mutation_type,
           i_p.net_par.max_arc,
           s_p.seed,
           s_p.adaptation_per) %>% 
  summarise(mean_phen_sens = mean(m_phenotype_sens),
            mean_fit_sens = mean(m_fitness_sens),
            mean_fit = mean(m_fitness),
            var_phen_sens = sd(m_phenotype_sens),
            var_fit_sens = sd(m_fitness_sens),
            var_fit = sd(m_fitness))

jpeg("phen_sens_fit_plots.jpg",
     width = 700,
     height = 700)

phen_sens_plus_fit_plot = ggplot(data = sens_summary_all) +
  geom_line(aes(x = m_generation,
                 y = mean_phen_sens)) +
  geom_line(aes(x = m_generation,
                y = mean_fit), color = "red") +
  geom_line(aes(x = m_generation,
                y = mean_fit_sens), color = "blue") +
  geom_ribbon(aes(x = m_generation,
                  y = mean_phen_sens,
                  ymax = mean_phen_sens + var_phen_sens,
                  ymin = mean_phen_sens - var_phen_sens),
              alpha = 0.5) +
  ylim(c(-0.1,1)) +
  xlab("Generations") +
  facet_grid(s_p.selection_freq + 
               s_p.selection_strength + 
               e_p.name_func_A +
               i_p.m_mutation_type + i_p.net_par.max_arc ~
               s_p.seed + s_p.adaptation_per)

print(phen_sens_plus_fit_plot)
dev.off()

fit_sens_plot = ggplot(data = sens_summary_all,
                       aes(x = m_generation,
                           y = mean_fit_sens)) +
  geom_line() +
  geom_ribbon(aes(x = m_generation,
                  y = mean_fit_sens,
                  ymax = mean_fit_sens + var_fit_sens,
                  ymin = mean_fit_sens - var_fit_sens),
              alpha = 0.5)+
  xlab("Generations") +
  facet_grid(s_p.selection_freq + 
               s_p.selection_strength + 
               e_p.name_func_A +
               i_p.m_mutation_type + i_p.net_par.max_arc ~
               s_p.seed + s_p.adaptation_per)

adapt_levels = levels(as.factor(all_simple_res$s_p.adaptation_per))
seed_levels = levels(as.factor(all_simple_res$s_p.seed))
sel_str_levels = levels(as.factor(all_simple_res$s_p.selection_strength))
sel_freq_levels = levels(as.factor(all_simple_res$s_p.selection_freq))
func_name_levels = levels(as.factor(all_simple_res$e_p.name_func_A))
mut_types_levels = levels(as.factor(all_simple_res$i_p.m_mutation_type))
max_arc_levels = levels(as.factor(all_simple_res$i_p.net_par.max_arc))
for(adapt_per in adapt_levels){
  for(seed in seed_levels){
    for(sel_str in sel_str_levels){
      for(sel_freq in sel_freq_levels){
        for(func_name in func_name_levels){
          for(mut_type in mut_types_levels){
            for(max_arc in max_arc_levels){

                            ###subset to a specific simulation for now
              sim_sens = all_sensibilities %>%
                filter(s_p.seed == seed) %>% 
                filter(s_p.adaptation_per == adapt_per) %>% 
                filter(s_p.selection_strength == sel_str) %>% 
                filter(s_p.selection_freq == sel_freq) %>% 
                filter(e_p.name_func_A == func_name) %>% 
                filter(i_p.m_mutation_type == mut_type) %>% 
                filter(i_p.net_par.max_arc == max_arc)
              
              ###get and plot fitnesses of the specific simulation 
              sim_fitness = all_simple_res %>%
                filter(s_p.adaptation_per == adapt_per) %>% 
                filter(s_p.seed == seed) %>% 
                filter(s_p.selection_strength == sel_str) %>% 
                filter(s_p.selection_freq == sel_freq) %>% 
                filter(e_p.name_func_A == func_name) %>% 
                filter(i_p.m_mutation_type == mut_type) %>% 
                filter(i_p.net_par.max_arc == max_arc)
              
              ###subset rns for all inds
              sim_all_inds_rns = all_inds_rns %>% 
                filter(s_p.adaptation_per == adapt_per) %>% 
                filter(s_p.seed == seed) %>% 
                filter(s_p.selection_strength == sel_str) %>% 
                filter(s_p.selection_freq == sel_freq) %>% 
                filter(e_p.name_func_A == func_name) %>% 
                filter(i_p.m_mutation_type == mut_type) %>% 
                filter(i_p.net_par.max_arc == max_arc)
              
              optimal_rn = produce_current_optimal_func(func_name = func_name, 
                                                        data.frame(x = unique(sim_all_inds_rns$m_x),
                                                                   y = length(unique(sim_all_inds_rns$m_x))))

              rn_cloud = ggplot(sim_all_inds_rns, aes(x = m_x,
                                   y = m_y)) +
                stat_density2d(geom="tile", aes(fill = ..count..), contour = FALSE) +
                scale_fill_viridis_c() +
                geom_line(data = optimal_rn, aes(x = x, y = y), color = "white") +
                facet_grid( . ~ m_generation)
              
              sens_summary = sim_sens %>%
                group_by(m_generation) %>% 
                summarise(mean_phen_sens = mean(m_phenotype_sens),
                          mean_fit_sens = mean(m_fitness_sens),
                          var_phen_sens = sd(m_phenotype_sens),
                          var_fit_sens = sd(m_fitness_sens))
              
              sens_dens_plot = ggplot(sim_sens, aes(x = m_generation,
                                                    y = m_phenotype_sens)) +
                stat_density2d(geom="tile", aes(fill = ..count..), contour = FALSE)
              
              fit_sens_dens_plot = ggplot(sim_sens, aes(x = m_generation,
                                                        y = m_fitness_sens)) +
                stat_density2d(geom="tile", aes(fill = ..count..), contour = FALSE)
              
              # sens_dens_plot / fit_sens_dens_plot
              
              fit_plot = ggplot(data = sim_fitness %>% 
                                  filter(gen %in% sens_summary$m_generation)) +
                geom_line(aes(x = gen, y = m_avg_fitnesses)) +
                geom_ribbon( aes(x = gen, y = m_avg_fitnesses,
                                 ymax = m_avg_fitnesses + m_var_fitnesses,
                                 ymin = m_avg_fitnesses - m_var_fitnesses), alpha = 0.5)
              
              phen_sens_plot = ggplot(data = sens_summary,
                                      aes(x = m_generation,
                                          y = mean_phen_sens)) +
                geom_line() +
                geom_ribbon(aes(x = m_generation,
                                y = mean_phen_sens,
                                ymax = mean_phen_sens + var_phen_sens,
                                ymin = mean_phen_sens - var_phen_sens),
                            alpha = 0.5) +
                xlab("Generations")
              
              fit_sens_plot = ggplot(data = sens_summary,
                                     aes(x = m_generation,
                                         y = mean_fit_sens)) +
                geom_line() +
                geom_ribbon(aes(x = m_generation,
                                y = mean_fit_sens,
                                ymax = mean_fit_sens + var_fit_sens,
                                ymin = mean_fit_sens - var_fit_sens),
                            alpha = 0.5)+
                xlab("Generations")
              
              #create directory where to save images
              subdir = paste(
                "phen_fit_sens_",
                "s", seed,
                "s_f", sel_freq,
                "s_s", sel_str,
                "a_p", adapt_per,
                "func", func_name,
                "mut_t", mut_type,
                sep = "_")
              dir.create(file.path(dir, subdir), showWarnings = FALSE)
              
              generations = unique(sim_sens$m_generation)
              plot_every_n_gen = 50000
              
              for(generation in generations[generations %% plot_every_n_gen == 0]){
                record_freq = as.numeric(obs_params$m_top_ind_reg_freq)
                selection_duration = as.numeric(as.character(unique(sim_sens$s_p.selection_duration)))
                
                #skip generations that are after the selection period
                #they will be included in the same plots
                #of the generations pre-selection period
                if((generation) %% record_freq == 0){
                  
                  #get the generations 
                  #across an entire cycle of
                  #drift - selection - drift - selection
                  post_sel_gen = generations[match(generation, generations) + 1]
                  post_drift_gen = generations[match(generation, generations) + 2]
                  post_drift_post_sel_gen = generations[match(generation, generations) + 3]
                  
                  important_subsequent_generations = c(post_sel_gen,
                                                       post_drift_gen,
                                                       post_drift_post_sel_gen)
                  gen_sens = sim_sens %>%
                    filter(m_generation %in% c(generation,
                                               post_sel_gen)) 
                  if(generation %% plot_every_n_gen == 0 &&
                     generation + record_freq < max(generations))
                  {
                    gen_sens = sim_sens %>%
                      filter(m_generation %in% c(generation,
                                                 important_subsequent_generations)) 
                    
                    ###to create network of descendants 
                    #create edge data frame(Id + ancestor ID) of latest gen
                    #and node data frame
                    #(all IDs and IDs of earlier generation present in the ancestor IDs of the latest generation)
                    edge_tibble = gen_sens %>% 
                      filter(m_generation %in% important_subsequent_generations) %>% 
                      select(c(m_ancestor_ID, m_ID)) %>% 
                      rename(from = m_ancestor_ID, to = m_ID)
                    
                    ancestry_graph = tbl_graph(nodes = gen_sens,
                                               edges = edge_tibble,
                                               node_key = "m_ID") 
                    # 
                    # p_anc <- ggraph(ancestry_graph, x = m_phenotype_sens, y = m_fitness) +
                    #   geom_edge_link(alpha = 0.01)+
                    #   # geom_edge_density()+
                    #   geom_node_point(aes(color =  as.factor(m_generation),
                    #                       size = node_is_isolated(),
                    #                       alpha = ifelse(node_is_isolated(),0.01,1))
                    #   ) +
                    #   scale_size_manual(values=c(1.5,0.5), guide = 'none') +
                    #   scale_alpha(guide = 'none') +
                    #   xlab("phenotypic_sens") +
                    #   ylab("fitness") + 
                    #   theme_light() +
                    # theme(legend.position="none")
                    
                    p_anc_inter_sel <- ggraph(ancestry_graph %>%
                                                filter(m_generation %in% c(generation,
                                                                           post_sel_gen)),
                                              x = m_phenotype_sens,
                                              y = m_fitness) +
                      geom_edge_link(alpha = 0.01)+
                      # geom_edge_density()+
                      geom_node_point(shape = 21,
                                      aes(color = as.factor(m_generation),
                                          # size = centrality_degree(mode = 'out'),
                                          fill = centrality_degree(mode = 'out'),
                                          alpha = ifelse(m_generation == post_sel_gen,0,1))
                      ) +
                      outdegree_gradient +
                      scale_shape_manual(c(21,24), guide = 'none') +
                      scale_size_continuous(c(0, 1), guide = 'none') +
                      scale_alpha(c(0, 1), guide = 'none') +
                      xlab("phenotypic_sens") +
                      ylab("fitness") +  
                      xlim(phen_x_lim) + 
                      ylim(c(0,1)) +
                      theme_light() +
                      theme(legend.position="none") +
                      ggtitle(paste(generation,"gen to ", post_sel_gen))
                    
                    p_anc_inter_drift <- ggraph(ancestry_graph %>%
                                                  filter(m_generation %in% c(post_sel_gen,
                                                                             post_drift_gen)),
                                                x = m_phenotype_sens,
                                                y = m_fitness) +
                      geom_edge_link(alpha = 0.01)+
                      # geom_edge_density()+
                      geom_node_point(shape = 21,
                                      aes(color = as.factor(m_generation),
                                          shape = as.factor(m_generation),
                                          # size = centrality_degree(mode = 'out'),
                                          fill = centrality_degree(mode = 'out'),
                                          alpha = ifelse(m_generation == post_drift_gen,0,1))
                      ) +
                      outdegree_gradient +
                      scale_shape_manual(c(21,24), guide = 'none') +
                      scale_size_continuous(c(0, 1), guide = 'none') +
                      scale_alpha(c(0.5, 1.5), guide = 'none') +
                      xlab("phenotypic_sens") +
                      ylab("fitness") +  
                      xlim(phen_x_lim) + 
                      ylim(c(0,1)) +
                      theme_light() +
                      theme(legend.position="none")+
                      ggtitle(paste(post_sel_gen,"gen to ", post_drift_gen))
                    
                    
                    p_anc_inter_sel + p_anc_inter_drift 
                    ggsave(paste(subdir,paste(paste("ancestry_plot",generation,sep = "_"),".png"), sep = '/'),
                           device = "png", 
                           dpi= "screen",
                           width = 15,
                           height = 7.5)
                  } 
                  p2 = ggplot(data = gen_sens %>% 
                                filter(m_generation %in% c(generation, 
                                                           post_sel_gen))) +
                    geom_point(shape = 21, 
                               aes(x = m_phenotype_sens,
                                   y = m_fitness,
                                   colour = m_fitness_sens
                               ), alpha = 0.5) +
                    ylim(c(0,1))+
                    xlim(phen_x_lim) +
                    fitness_sensitivity_palette +
                    facet_grid(.~as.factor(m_generation)) 
                  
                  p3 = ggplot(data = gen_sens %>% 
                                filter(m_generation %in% c(generation,
                                                           post_sel_gen))) +
                    geom_point(shape = 21, 
                               aes(x = m_fitness_sens,
                                   y = m_phenotype_sens,
                                   # fill = m_rank,
                                   colour = m_fitness
                               ), alpha = 0.5)  +
                    xlim(fit_x_lim) +
                    ylim(phen_x_lim) +
                    fitness_gradient +
                    facet_grid(.~as.factor(m_generation))
                  
                  p4 = rn_cloud /
                    (fit_plot + 
                          geom_hline(yintercept = as.numeric(sim_fitness %>%
                                                               filter(gen == generation) %>% 
                                                               select(m_avg_fitnesses)),
                                     color = "red") +
                          geom_vline(xintercept = as.numeric(generation),
                                     color = "red")) /
                    (phen_sens_plot + geom_vline(xintercept = as.numeric(generation),
                                                 color = "red")) /
                    (fit_sens_plot  + geom_vline(xintercept = as.numeric(generation),
                                                 color = "red")) 
                  
                  (p3 + p2) / p4
                  
                  ggsave(paste(subdir,paste(paste("phen_fit_sens_plot",generation,sep = "_"),".png"), sep = '/'),
                         device = "png", 
                         width = 30,
                         height = 15)
                }
              }
            }
            ####Create gif
            # imgs = list.files(path = subdir, pattern = "*")
            # ## list file names and read in
            # img_list = lapply(imgs, image_read)
            # 
            # ## join the images together
            # img_joined <- image_join(img_list)
            # 
            # ## animate at 2 frames per second
            # img_animated <- image_animate(img_joined, fps = 2)
            # 
            # ## save to disk
            # path = paste("Gif",subdir,".gif", sep = "_")
            # image_write(image = img_animated,
            #             path = path)
          }
        }
      }
    }
  }
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
