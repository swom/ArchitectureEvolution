library(jsonlite)
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
library(scales)

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
  
  filepaths = list.files(pattern = pattern)[3]
  all_spectrums = list()
  all_rns = list()
  for (i in  filepaths){
    tryCatch(
      {
        results <- jsonlite::fromJSON(i, simplifyVector = FALSE)
        print(paste("loading: ",i, sep = "" ))
        
        for(generation in 1:length(results)){
          result = results[[generation]]
          for(individual in 1:length(result)){
            
            # trying to extract automatically interesting part of results
            # u = as.data.table(results)
            # ut = data.table()
            # u$V1[1][[1]]$net_spectrum$m_net_spectrum_weights_for_weights_mutation
            # uu = lapply(u,  print)
            # uu = unnest(results)
            # 
            
            #extracting manually interesting part of results
            tmp_ =  result[[individual]]$net_spectrum$m_net_spectrum_weights_for_weights_mutation
            
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
            
            ###Transforming lists into data.tables 
            tmp_ = lapply(tmp_, function(layer){
              lapply(layer, function(node){
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
              node = (as.data.table(node)[, rn := .I][, (.cols) := lapply(.SD, function(column){lapply(column, rbindlist)}), .SDcols = .cols] %>% #add id and then trandsform into dataframe each reaction norm
                        melt(id.vars = "rn", measure.vars = .cols, variable.name = "weight"))[,rbindlist(value), by = c("rn","weight")]#[, weight :=  as.integer(sub("weight_", weight, replacement = ""))] reshape resulting dataframe so that instead of 1 col x weight there is 1 col with reaction_norm  + id_col + weight_id_col, make it so weight is only a number
            }
            
            #Breaking down how unroll node works on a single node step by step
            # layer_1 = stt[[2]]
            # node_1 = layer_1[[1]]
            # .cols <- colnames(as.data.table(node_1))
            # node_a = (as.data.table(node_1)[, rn := .I]) ## add ID
            # node_a1 = node_a[, (.cols) := lapply(.SD, function(column){lapply(column, rbindlist)}), .SDcols = .cols] #transform weight column from col of lists to col of df
            # node_b = melt(node_a1, id.vars = "rn", measure.vars = .cols, variable.name = "weight")
            # node_c = node_b[,rbindlist(value), by = c("rn","weight")]
            # node_c[, weight :=  as.integer(sub("weight_", weight, replacement = ""))] #substringing so that only number is left
            
            #unrolling all nodes
            tmp_ = lapply(tmp_, function(layer){
              lapply(layer, unroll_node)
            })
            
            #unrolling nodes in layers into single dataframe
            tmp_ = lapply(tmp_, function(layer){rbindlist(layer, idcol = "node")}) 
            
            #unrolling layers in network into ssingle dataframe
            tmp_ = rbindlist(tmp_, idcol = "layer")
            
            #Add ID of simulation
            tmp_[, sim_ID := i]        
            
            ###Add to all mutational_spectra
            all_spectrums[[i]][[paste("gen_",generation, sep = "")]][[paste("ind_",individual,sep = "")]] =  tmp_
            
            # ###compact version
            tmp_ = result[[individual]]$net_spectrum$m_current_reac_norm
            tmp_ = lapply(tmp_, as.data.table)
            tmp_[[1]] = tmp_[[1]][, -"ind"]
            tmp_ = rbindlist(tmp_)
            
            # z= tmp_[[1]]
            # zz = z[, ind := .I]
            # zzz = zz[ , rbindlist(lapply(m_reac_norm, rbindlist)), by = "ind"]
            # zzzz = zzz[, generation := gens[1]]
            # zzzzz = tmp_[[1]][, ind := .I][ , rbindlist(lapply(m_reac_norm, rbindlist)), by = "ind"][, generation := gens[1]]
            
            tmp_[, sim_ID := i]
            all_rns[[i]][[paste("gen_",generation, sep = "")]][[paste("ind_",individual,sep = "")]] = tmp_
            
            print(paste("loaded ind: ", individual, sep = ""))
            gc()   
          }
          print(paste("       loaded gen: ", generation, sep = ""))
        }
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
  
  obs_params = results$m_obs_param
  saveRDS(all_rns, file = "all_rns.Rds")
  saveRDS(obs_params, file = "obs_params.Rds")
  saveRDS(all_spectrums, file = "all_spectrums.Rds")
  
  gc()
}

#Create optimal reaction norm curve
optimal_rn = produce_current_optimal_func(func_name = "3", 
                                          data.frame(x = unique(all_spectrums[[1]][[1]][[1]]$m_x),
                                                     y = length(unique(all_spectrums[[1]][[1]][[1]]$m_x))))


all_rns  = readRDS(file = "all_rns.Rds")
obs_params = readRDS(file = "obs_params.Rds")
all_spectrums = readRDS(file = "all_spectrums.Rds")


###If I want to plot also optimal rn andactual rn it is bettwer to facet everything and cover plots that are not good
# all_facets = unique(expand.grid(unique(select(all_spectrums[[1]][[generation]][[individual]], c(layer, node, weight)))))
# all_spectrums[[1]][[generation]][[individual]] %>% 
#   full_join(all_facets) %>%   
#   mutate(empty=ifelse(is.na(m_x), "X", NA_character_))

backgorund_color = viridis_pal()(20)[1]
for(generation in 1:length(all_spectrums[[1]])){
  for(individual in 1:length(all_spectrums[[1]][[generation]])){
    
    jpeg(paste("mut_spectrum_gen_", generation, "_ind_", individual, ".jpg", sep = ""),
         width = 600,
         height = 300, 
         units ="mm",
         res = 300)
    
    #commented out unnecessary layers that would make a prettier visulaization but serve no purpose
    rn_cloud = ggplot(all_spectrums[[1]][[generation]][[individual]]
                      # %>% 
                      #   full_join(all_facets) %>%   
                      #   mutate(empty=ifelse(is.na(m_x), "X", "")) %>% 
                      #   filter(node != "node_2" | weight != "weight_2")
                      ,
                      aes(x = m_x,y = m_y)) +
      geom_line(aes(x = m_x, y = m_y, group = rn), color = "white", alpha = 0.1, linewidth = 0.1) +
      geom_line(data = optimal_rn, aes(x = x, y = y), color = "red", linewidth = 1) +
      geom_line(data = all_rns[[1]][[generation]][[individual]], aes(x = m_x, y = m_y), color = "green", linewidth = 0.2) +
      # geom_rect(data= .%>% 
      #             filter(is.na(m_x)) %>% 
      #             select(m_x, m_y,empty),
      #           aes(ymin= -1, ymax= 1, xmin = -1, xmax = 1, fill= empty), ) +
      # geom_text(aes(x = 0, y = 0, label=empty), colour="red", size=200) +
      ylab("phenoytpe") +
      xlab("cue") +
      theme_classic2()+
      theme(legend.position = "none") +
      theme(panel.background = element_rect(fill = backgorund_color)) +
      facet_grid(layer ~ node + weight) +
      ggtitle(paste("gen ",as.character(generation), "; ind: ", individual,  sep = "")) 
    print(rn_cloud)
    dev.off()
    print(paste("plot ind: ", individual,"; gen: ",generation, sep = ""))
  }
}
