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
dir ="C:/Users/p288427/Desktop/data_dollo_++/10_19_22"
# dir ="C:/Users/p288427/Github/build-ArchitectureEvolution-Desktop_Qt_6_2_4_MSVC2019_64bit-Release/src"
setwd(dir)
pattern = '^m.*json$'

####save load####
{
  
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
  }
  saveRDS(all_simple_res, file = "all_simple_res_tail.Rds")
  gc()
}

s1_off = "mut_wei_sel_spo_sym_sym_fr_reg_ap_off_r_con_e_ful_arc_1-2-2-2-1_marc_1-2-2-2-1_wr_0.0_ar_0.0_dup_0.0_cA_0.0_cB_0.0_st_1.0_sf_100_fA_4_g_100000_p_1000_seed1.json" 
s1_on = "mut_wei_sel_spo_sym_sym_fr_reg_ap_on_r_con_e_ful_arc_1-2-2-2-1_marc_1-2-2-2-1_wr_0.0_ar_0.0_dup_0.0_cA_0.0_cB_0.0_st_1.0_sf_100_fA_4_g_100000_p_1000_seed1.json" 
s1 = c(s1_off, s1_on)

all_inds_rns = list()
for (i in s1) {
  results <- fromJSON(file = i)
  simple_res = tail(rowid_to_column(as_tibble(results[c("m_avg_fitnesses")]), var = "gen"), n = 10000) 
  results$m_params$i_p$net_par$max_arc = toString(results$m_params$i_p$net_par$max_arc)
  results$m_params$i_p$net_par$net_arc = toString(results$m_params$i_p$net_par$net_arc)
  ID = as.data.frame(results$m_params)
  
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
  all_inds_rns[[i]] = tmp_ 
  
  optimal_rn = produce_current_optimal_func(func_name = ID$e_p.name_func_A, 
                                            data.frame(x = unique(tmp_$m_x),
                                                       y = length(unique(tmp_$m_x))))
  rn_cloud = ggplot(tmp_, aes(x = m_x,y = m_y)) +
    stat_density2d(geom="tile", aes(fill = ..count..), contour = FALSE) +
    scale_fill_viridis_c() + 
    geom_line(aes(x = m_x, y = m_y, group = ind),color = "white", alpha = 0.1) +
    geom_line(data = optimal_rn, aes(x = x, y = y), color = "red", size = 1) +
    facet_grid( . ~ m_generation)+
    theme_minimal()
  
  fit_plot = ggplot(data = simple_res %>% 
                      filter(gen %% 100 == 0
                             # %in% sens_summary$m_generation
                      )) +
    geom_line(aes(x = gen, y = m_avg_fitnesses))
  std_fit_plot = fit_plot + ylim(c(0,1))
  p_overall = rn_cloud /
    (fit_plot) /
    std_fit_plot
  print(p_overall)
  
}
all_inds_rns = rbindlist(all_inds_rns, fill = T) 

sum_res = all_simple_res %>% 
  group_by(s_p.seed, s_p.selection_freq, s_p.adaptation_per) %>% 
  summarise(mean_fit = mean(m_avg_fitnesses))
### Plot ====
jpeg("fitnesses_box_plots.jpg",
     width = 1000,
     height = 1000)

wanted_freqs = c(0, 1, 5, 10, 20, 100)
adapt_per = c(0,1)
# wanted_freqs = c(20)
# wanted_sel_str = c(0.1, 0.5, 1)
p <- sum_res %>% 
  ggplot(aes(x = as.factor(1/s_p.selection_freq),
             y = mean_fit,
             fill=factor(s_p.adaptation_per))) + 
geom_boxplot()

print(p)
dev.off()
