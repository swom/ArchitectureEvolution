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
# dir ="C:/Users/p288427/Desktop/data_dollo_++/10_17_22_sqrt_avg_sqr_dist/trial"
dir ="C:/Users/p288427/Github/ArchitectureEvolution"
setwd(dir)
pattern = '^M.*json$'
# pattern = 'mut_wei_sel_spo_sym_sym_fr_reg_ap_on_r_con_e_tri_arc_1-2-2-2-1_marc_1-2-2-2-1_wr_0.0_ar_0.0_dup_0.0_cA_0.0_cB_0.0_st_1.0_sf_5_fA_3_g_100001_p_1000_seed1.json'
filepaths = list.files(pattern = pattern)

####save load####
{
  
  filepaths = list.files(pattern = pattern)
  
  all_sensibilities = list()
  all_simple_res = data.frame()
  all_inds_rns = list()
  for (i in  filepaths){
    tryCatch(
      {
        i= filepaths[1]
        results <- fromJSON(file = i)
        
       tmp_ =  results[[1]][[1]]$net_spectrum$m_net_spectrum_weights_for_weights_mutation
        ###mut_spectrum
        t = lapply(tmp_, as.data.table)
        tt = as.data.table(t[[1]])
        ttt = as.data.table(tt[[1]])[, id := .I]
        ttt1 = lapply(ttt[[1]], rbindlist)
        ttt[[1]] = rbindlist(lapply(ttt[[1]], rbindlist))
        ttt2 = rbindlist(ttt1)
        tttt = ttt[, rbindlist(lapply(ttt[[1]],rbindlist)), by = "id"]
        tttt = ttt[, rbindlist, by = "id"]
        
        tttt= lapply(seq_along(ttt),function(i){ttt[i] = rbindlist(ttt[[i]])})
        tttttt = ttt[, rbindlist(), by = "id"]
        ttttttt = rbindlist(ttt[[1]])
        
        
        
        
        tmp__ = rbindlist(lapply(seq_along(tmp_),function(i){
          tmp_[[i]]= as.data.table(tmp_[[i]])[, layer := .I][ , rbindlist(lapply(tmp_[[i]], rbindlist))]
        })
        )
        tmp_[, names(ID) := ID]
        all_inds_rns[[i]] = tmp_
        
        gc()   
        print(paste("loading: ",i, sep = "" ))
      },
      error=function(cond){
        message(paste("could not load file: ", i, sep = ""))
        dir.create(file.path(dir, "defective"), showWarnings = FALSE)
        file.move(i, "defective")
        message("Moving defective file into defective directory")
        message("The original message says:")
        message(cond)
      }
    )
  }
  
  all_sensibilities = rbindlist(all_sensibilities, fill =T) %>%
    #add 1 to all generations to sync with all_simple_res
    mutate(m_generation = m_generation + 1)
  
  all_inds_rns = rbindlist(all_inds_rns, fill = T) %>%
    rename(m_generation = generation) %>%
    mutate(m_generation = m_generation + 1)
  
  obs_params = results$m_obs_param
  all_sensibilities$m_generation = as.factor(all_sensibilities$m_generation)
  saveRDS(obs_params, file = "obs_params.Rds")
  saveRDS(all_simple_res, file = "all_simple_res.Rds")
  saveRDS(all_sensibilities, file = "all_sensibilities.Rds")
  saveRDS(all_inds_rns, file = "all_inds_rns.Rds")
  gc()
}