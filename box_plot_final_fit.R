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
library(purrr)
library(forcats)

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
dir ="C:/Users/p288427/Desktop/data_dollo_++/10_20_22_1_10thcomplete"
# dir ="C:/Users/p288427/Github/build-ArchitectureEvolution-Desktop_Qt_6_2-_4_MSVC2019_64bit-Release/src"
setwd(dir)

####save & plot####
{
  {
    pattern = '^*.*json$'
    filepaths = list.files(pattern = pattern)
    all_simple_res = data.frame()
    for (i in  filepaths){
      tryCatch(
        {
          results <- fromJSON(file = i)
          
          ###simple results
          simple_res = tail(rowid_to_column(as_tibble(results[c("m_avg_fitnesses")]), var = "gen"), n = 10000) 
          results$m_params$i_p$net_par$max_arc = toString(results$m_params$i_p$net_par$max_arc)
          results$m_params$i_p$net_par$net_arc = toString(results$m_params$i_p$net_par$net_arc)
          ID = as.data.frame(results$m_params)
          simple_res_ID = cbind(simple_res, ID)
          all_simple_res = rbind(all_simple_res, simple_res_ID)
          gc()   
          print(paste("loading: ",i, sep = "" ))
        },
        error=function(cond){
          message(paste("could not load file: ", i, sep = ""))
          dir.create(file.path(dir, "defective"), showWarnings = FALSE)
          #file.move(i, "defective")
          message("Moving defective file into defective directory")
          message("The original message says:")
          message(cond)
        }
      )
      # tryCatch(
      #   {   
      #     tmp_ = results$m_all_inds_rn
      #     gens <- map_dbl(tmp_, ~ .x$generation)
      #     tmp_ = lapply(tmp_, as.data.table)
      #     tmp_ = rbindlist(lapply(seq_along(tmp_),function(i){
      #       tmp_[[i]]= tmp_[[i]][, ind := .I][ , rbindlist(lapply(m_reac_norm, rbindlist)), by = "ind"][, generation := gens[i]]
      #     })
      #     )
      #     tmp_[, names(ID) := ID]
      #     tmp_ = tmp_ %>%
      #       rename(m_generation = generation) %>%
      #       mutate(m_generation = m_generation + 1)
      #     tmp_[, names(ID) := ID]
      #     # all_inds_rns[[i]] = tmp_ 
      #   },
      #   error=function(cond){
      #     message(paste("could not wrangle file: ", i, sep = ""))
      #     message("The original message says:")
      #     message(cond)
      #   }
      # )
      ###Plot RN clouds
      # tryCatch(
      #   {       
      #     optimal_rn = produce_current_optimal_func(func_name = ID$e_p.name_func_A, 
      #                                               data.frame(x = unique(tmp_$m_x),
      #                                                          y = length(unique(tmp_$m_x))))
      #     rn_cloud = ggplot(tmp_, aes(x = m_x,y = m_y)) +
      #       stat_density2d(geom="tile", aes(fill = ..count..), contour = FALSE) +
      #       scale_fill_viridis_c() + 
      #       geom_line(aes(x = m_x, y = m_y, group = ind),color = "white", alpha = 0.1, size = 0.1) +
      #       geom_line(data = optimal_rn, aes(x = x, y = y), color = "red", size = 1) +
      #       # facet_grid( . ~ m_generation)+
      #       theme_pubclean() + 
      #       theme(legend.position = "none") +
      #       xlab("cue")+
      #       ylab("phenotype")
      #     
      #     # fit_plot = ggplot(data = simple_res %>% 
      #     #                     filter(gen %% 100 == 0
      #     #                            # %in% sens_summary$m_generation
      #     #                     )) +
      #     #   geom_line(aes(x = gen, y = m_avg_fitnesses))
      #     # std_fit_plot = fit_plot + ylim(c(0,1))
      #     
      #     {
      #       jpeg(paste("rn_cloud_",str_sub(i, end = -6),".jpg", sep = ""),
      #            width = 76,
      #            height = 76,
      #            units = "mm", 
      #            res = 300)
      #       p_overall = rn_cloud# /
      #       # (fit_plot) /
      #       # std_fit_plot
      #       print(p_overall)
      #       dev.off()
      #       print(paste("plotted:",str_sub(i, end = -6), sep = " "))
      #       }
      #   },
      #   error=function(cond){
      #     message(paste("could not plot file: ", i, sep = ""))
      #     message("The original message says:")
      #     message(cond)
      #   }
      # )
      }
    gc()
  }
  saveRDS(all_simple_res, file = "all_simple_res_tail.Rds")

  sum_res = all_simple_res %>% 
    group_by(s_p.seed, s_p.selection_freq, s_p.adaptation_per, s_p.selection_duration) %>% 
    summarise(mean_fit = mean(m_avg_fitnesses)) %>% 
    mutate(sel_regime = as.factor(paste(as.character(s_p.selection_duration), as.character(s_p.selection_freq),sep="/"))) %>% 
    mutate(sel_regime = recode_factor(sel_regime, "0/0" = "No_selection")) %>% 
    mutate(sel_regime = fct_relevel(sel_regime, c("1/1", "1/5", "1/10", "10/100", "1/100", "No_selection"))) 
    
  saveRDS(sum_res, file = "summarised_results.Rds")
  ### Plot ====
  jpeg("fitnesses_box_plots.jpg",
       width = 300,
       height = 300, 
       units ="mm",
       res = 900)
  
  wanted_freqs = c(0, 1, 5, 10, 20, 100)
  adapt_per = c(0,1)
  # wanted_freqs = c(20)
  # wanted_sel_str = c(0.1, 0.5, 1)
  p <- sum_res %>% 
    ggplot(aes(x = sel_regime,
               y = mean_fit,
               fill=factor(s_p.adaptation_per))) + 
    geom_boxplot() +
    theme_classic2() +
    xlab("selection regime: n rounds of selection / N generations") +
    ylab("final average fitness") +
    ylim(c(0,1)) +
    theme(legend.position = "none") +
    theme(text =  element_text(size = 25))     
  
  print(p)
  dev.off()
}



##Plotting all recatiuon norms from saved files
readRDS("summarised_results.Rds")
{
  # pattern1 = "^sut_wei_sel_spo_sym_sym_fr_reg_ap_.*_r_con_e_ful_arc_1-2-2-2-1_marc_1-2-2-2-1_wr_0.0_ar_0.0_dup_0.0_cA_0.0_cB_0.0_st_1.0_sf_.*_fA_4_g_100000_p_1000_seed1.json" 
  # pattern2 = "^mut_wei_sel_spo_sym_sym_fr_reg_ap_.*_r_con_e_ful_arc_1-2-2-2-1_marc_1-2-2-2-1_wr_0.0_ar_0.0_dup_0.0_cA_0.0_cB_0.0_st_1.0_sf_.*_fA_4_g_100000_p_1000_seed1.json" 
  pattern1 = "^mut_wei_sel_spo_sym_sym_fr_reg_ap_off_r_con_e_ful_arc_1-2-2-2-1_marc_1-2-2-2-1_wr_0.0_ar_0.0_dup_0.0_cA_0.0_cB_0.0_st_1.0_sf_100_fA_4_g_100000_p_1000_seed11.json"
  pattern2 = "^x"
  filepaths = c(list.files(pattern = pattern1),list.files(pattern = pattern2))
  # all_inds_rns = list()
  for (i in c(list.files(pattern = pattern1),list.files(pattern = pattern2))) {
    tryCatch(
      {
        results <- fromJSON(file = i)
        print(paste("loading: ",i, sep = "" ))
        results$m_params$i_p$net_par$max_arc = toString(results$m_params$i_p$net_par$max_arc)
        results$m_params$i_p$net_par$net_arc = toString(results$m_params$i_p$net_par$net_arc)
        ID = as.data.frame(results$m_params)
      },
      error=function(cond){
        message(paste("could not load file: ", i, sep = ""))
        dir.create(file.path(dir, "defective"), showWarnings = FALSE)
        #file.move(i, "defective")
        message("Moving defective file into defective directory")
        message("The original message says:")
        message(cond)
      })  
    tryCatch(
      {   
        tmp_ = results$m_all_inds_rn
        gens <- map_dbl(tmp_, ~ .x$generation)
        tmp_ = lapply(tmp_, as.data.table)
        tmp_ = rbindlist(lapply(seq_along(tmp_),function(i){
          tmp_[[i]]= tmp_[[i]][, ind := .I][ , rbindlist(lapply(m_reac_norm, rbindlist)), by = "ind"][, generation := gens[i]]
        })
        )
        tmp_[, names(ID) := ID]
        tmp_ = tmp_ %>%
          rename(m_generation = generation) %>%
          mutate(m_generation = m_generation + 1)
        tmp_[, names(ID) := ID]
        # all_inds_rns[[i]] = tmp_ 
        },
      error=function(cond){
        message(paste("could not wrangle file: ", i, sep = ""))
        message("The original message says:")
        message(cond)
      }
    )
    tryCatch(
      {       
        optimal_rn = produce_current_optimal_func(func_name = ID$e_p.name_func_A, 
                                                  data.frame(x = unique(tmp_$m_x),
                                                             y = length(unique(tmp_$m_x))))
        rn_cloud = ggplot(tmp_, aes(x = m_x,y = m_y)) +
          stat_density2d(geom="tile", aes(fill = ..count..), contour = FALSE) +
          scale_fill_viridis_c() + 
          geom_line(aes(x = m_x, y = m_y, group = ind),color = "white", alpha = 0.1, size = 0.1) +
          geom_line(data = optimal_rn, aes(x = x, y = y), color = "red", size = 1) +
          # facet_grid( . ~ m_generation)+
          theme_pubclean() + 
          theme(legend.position = "none") +
          xlab("cue")+
          ylab("phenotype")
        
        # fit_plot = ggplot(data = simple_res %>% 
        #                     filter(gen %% 100 == 0
        #                            # %in% sens_summary$m_generation
        #                     )) +
        #   geom_line(aes(x = gen, y = m_avg_fitnesses))
        # std_fit_plot = fit_plot + ylim(c(0,1))
        
        {
          jpeg(paste("rn_cloud_",str_sub(i, end = -6),".jpg", sep = ""),
             width = 76,
             height = 76,
             units = "mm", 
             res = 300)
        p_overall = rn_cloud# /
        # (fit_plot) /
        # std_fit_plot
        print(p_overall)
        dev.off()
        print(paste("plotted:",str_sub(i, end = -6), sep = " "))
          }
      },
      error=function(cond){
        message(paste("could not plot file: ", i, sep = ""))
        message("The original message says:")
        message(cond)
      }
    )
  }
  all_inds_rns = rbindlist(all_inds_rns, fill = T) 
}
