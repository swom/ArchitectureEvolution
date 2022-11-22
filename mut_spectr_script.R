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
dir ="C:/Users/p288427/Github/build-ArchitectureEvolution-Desktop_Qt_6_2_4_MSVC2019_64bit-Release/src"
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
        
        # tryingto extract automatically interesting part of results
        u = as.data.table(results)
        ut = data.table()
        u$V1[1][[1]]$net_spectrum$m_net_spectrum_weights_for_weights_mutation
        uu = lapply(u,  print)
        uu = unnest(results)
        
        
        #extracting manually interesting part of results
        tmp_ =  results[[1]][[1]]$net_spectrum$m_net_spectrum_weights_for_weights_mutation
        
        #for each element in each layer give name LAYER and NODE using set_names: https://stackoverflow.com/questions/35840692/r-how-to-set-names-in-a-nested-list-from-attributes
        #works, but could probably be compacted in one call
        #naming layers
        tmp_ = set_names(tmp_, lapply(seq_along(tmp_), function(l){
          paste("layer_", l, sep = "")
        }))
        #naming nodes
        tmp_ = lapply(tmp_, function(x){
          set_names(x,lapply(seq_along(x), function(n){paste("node_", n, sep = "")}))
        })
        
        #naming weights
        tmp_ = lapply(tmp_, function(x){
          lapply(x, function(y){
            set_names(y, lapply(seq_along(y), function(w){paste("weight_", w, sep = "")}))
          })
        })  
        
        #naming rn_s
        tmp_ = lapply(tmp_, function(layer){
          lapply(layer, function(node){
            lapply(node, function(weights){
              set_names(weights, lapply(seq_along(weights), function(r){paste("reac_n_", r, sep = "")}))
            })
          })
        }) 
        
        ###Transforming lists into data.tables and adding ID to reac_norms by weight ->slowish
        stt = lapply(tmp_, function(layer){
          lapply(layer, function(node){
            # node = as.data.table(node)[, id := .I]
            lapply(node, function(weight){
              lapply(weight, function(rn){
                lapply(rn, function(x_y){
                  as.data.table(x_y)
                })
              })
            })
          })
        })
    
        #adding from:https://stackoverflow.com/questions/46595080/lapply-to-all-columns-in-a-data-frame-except-one-and-replace-the-data-in-r
        #.cols <- setdiff(colnames(days), "id")
        # days[, (.cols) := lapply(.SD, round, digits = 1), .SDcols = .cols]
        unroll_node = function(node){
          .cols <- colnames(as.data.table(node))
          node = (as.data.table(node)[, id := .I][, (.cols) := lapply(.SD, function(column){lapply(column, rbindlist)}), .SDcols = .cols] %>% #add id and then trandsform into dataframe each reaction norm
                  melt(id.vars = "id", measure.vars = .cols))[,rbindlist(value), by = c("id","variable")] #reshape resulting dataframe so that instead of 1 col x weight trhere is 1 col with reaction norm  + id col + weight_id_col
        }
        
        #unrolling all nodes
        sttt = lapply(stt, function(layer){
          lapply(layer, unroll_node)
        })

        #unrolling nodes in layers into single dataframe
        #unrolling layers in network into ssingle dataframe
        
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