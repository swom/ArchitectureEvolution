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

dir = dirname(rstudioapi::getActiveDocumentContext()$path)
rand_evo_dir = paste(dir,"/vG_W_2016_data/rand_evo",sep = "")
evo_dir = paste(dir,"/vG_W_2016_data/evo",sep = "")

setwd(evo_dir)
setwd(dir = "C:/Users/p288427/Desktop/research presentation/gif_s9_dd")
funders = data.frame()

# for (i in  list.files(path = '.',pattern = "funders_success_s\\d+change_\\d+"))
# {
#   replicate = read.csv(i)
#   replicate$seed = sub( "^.*s(\\d+).*",'\\1', i)
#   replicate$change = sub( "^.*_(\\d+).*",'\\1', i)
#   colnames(replicate) = colnames(funders)
#   funders = rbind(replicate, funders)
# }

name = "X:/build-simulation_logic_only-Desktop_Qt_6_0_0_MinGW_64_bit-Release/death_rand_evo_extreme_a3.000000seq_1cond_per_seq50funders_success_s9change_0.csv"
# funders = read.csv("funders_success_s9change_0.csv")
funders = read.csv(name)
#put some useful names to column
names(funders)[1]<-"cycle" 
names(funders)[2]<-"ID" 
names(funders)[length(names(funders))]<-"success" 


#find best individual
best = funders %>% 
  subset(cycle == max(cycle)) %>% 
  subset(success == max(success))

funders$ID = str_sub(funders$ID, end = - 3)
#find line of descent of best individual
best_descent =  filter(funders, stri_startswith_fixed(best$ID, funders$ID)) %>% subset(cycle > 1
                                                                                       )

###plot network####
#assuming we know network architecture:
# 3-3-2, fully connected, with hidden layer fully connected to itself
n_nodes_l1 = 3
n_nodes_l2 = 3
n_nodes_l3 = 2

#generate coordinates for positioning nodes correctly in plot

l =   cbind(
  c(
    rep(1,n_nodes_l1)
    ,rep(2,n_nodes_l2)
    ,rep(3,n_nodes_l3)
  ),
  c(seq(1:n_nodes_l1),
    seq(1:n_nodes_l2),
    seq(1:n_nodes_l3) + 0.5 ))
l[n_nodes_l1 + n_nodes_l2 / 2 + 1] = l[n_nodes_l1 + n_nodes_l2 / 2 + 1] + 0.5
colnames(l) = c("x","y")

for(row in 1:nrow(best_descent))
{
  #get clean dataframe(no characters)
  funder= best_descent[row,]
  funder_mod = funder[,3:(length(funders) - 6)]
  funder_mod$seed = 9
  funder_mod$change = 0
  
  no_ch_funder = funder_mod  %>% select_if(is.numeric)
  
  connections = data.frame()
  
  #I2H connections
  for( i in 1:(n_nodes_l1 * n_nodes_l2))
  {
    ID_sender = 1 + ((i - 1) %/% n_nodes_l2)
    ID_receiver = n_nodes_l1 + if(i %% n_nodes_l2 == 0)  n_nodes_l2 else i %% n_nodes_l2
    weight = as.numeric(no_ch_funder[i])
    node = data.frame("ID1" = ID_sender, "ID2" = ID_receiver, "weight" = weight)
    connections = rbind(connections, node)
  }
  
  #H2H connections
  #the offset from which we will start iterating in the weights vector
  offset1 =  n_nodes_l1 * n_nodes_l2
  
  for( i in 1:(n_nodes_l2 * n_nodes_l2))
  {
    #add the number of nodes in previous layer to find ID
    ID_sender = n_nodes_l1 + 1 + ((i - 1) %/% n_nodes_l2)
    #same as before
    ID_receiver = n_nodes_l1 + if(i %% n_nodes_l2 == 0)  n_nodes_l2 else i %% n_nodes_l2
    #apply the offset to the iterator
    weight = as.numeric(no_ch_funder[offset1 + i ])
    node = data.frame("ID1" = ID_sender, "ID2" = ID_receiver, "weight" = weight)
    connections = rbind(connections, node)
  }
  
  #H2O connections
  #calculate the offset from which we will start iterating in the weights vector
  offset2 =  n_nodes_l1 * n_nodes_l2 + n_nodes_l2 * n_nodes_l2
  
  for( i in 1:(n_nodes_l2 * n_nodes_l3))
  {
    #add the number of nodes in previous layer to find ID
    ID_sender = n_nodes_l1 + 1 + ((i - 1) %/% n_nodes_l3)
    #same as before
    ID_receiver = n_nodes_l1 + n_nodes_l2 + if(i %% n_nodes_l3 == 0)  n_nodes_l3 else i %% n_nodes_l3
    #apply the offset to the iterator
    weight = as.numeric(no_ch_funder[offset2 + i ])
    node = data.frame("ID1" = ID_sender, "ID2" = ID_receiver, "weight" = weight)
    connections = rbind(connections, node)
  }
  
  #create column to show if weight is above or below 0
  
  connections = connections %>% 
    mutate(w_sign = if_else(weight < 0, 1, 2))
  
  #Node ID list 
  
  nodes = as.data.frame(unique(c(connections[,1],connections[,2])))
  colnames(nodes) = "nodes"
  
  ##create igraph or ggraph object
  
  #undirected to calculate curvature
  network <- igraph::graph_from_data_frame(d = connections,
                                           vertices = nodes,
                                           directed = F)
  d = curve_multiple(network)
  
  #adjusting curvature to show connections
  d[12] = 0.5
  d[15] = 0.5
  d[13] = -0.5
  
  #directed to plot arrows
  network_d <- igraph::graph_from_data_frame(d = connections,
                                             vertices = nodes,
                                             directed = T)
  
  E(network_d)$color = as.factor(connections$w_sign)
  
  jpeg(paste("Plot","s",funder_mod$seed,
             "change",funder_mod$change,
             "cycle",row,".png", sep = "_")
       ,width = 700,
       height = 700)
  
  plot(network_d, layout = l,
       edge.curved = d,
       edge.arrow.size = 0.5,                           # Arrow size, defaults to 1
       edge.arrow.width = 0.7,                          # Arrow width, defaults to 1
       edge.arrow.height = 0.9,                          # Arrow width, defaults to 1
       edge.lty = c("solid"),
       edge.width = abs(E(network_d)$weight/max(E(network_d)$weight) * 10),             # Edge width, defaults to 1
       main = row)
  dev.off()
  
}

#Create gif
pics = paste(evo_dir,"pics", sep = "/")
## list file names and read in
imgs <- list.files(pics, full.names = TRUE)
img_list <- lapply(imgs, image_read)


## join the images together
img_joined <- image_join(img_list)

## animate at 2 frames per second
img_animated <- image_animate(img_joined, fps = 25)

## view animated image
img_animated

## save to disk
image_write(image = img_animated,
            path = "grn-evo.gif")


#####Other network plot methods#####
###None seem to inclide self-loops, meh...###
g = ggnetwork::ggnetwork(network_d, layout = l, loop = T )
g$w_sign = as.factor(g$w_sign)
ggplot(g, aes(x = x, y = y, xend = xend, yend = yend)) +
  ggnetwork::geom_edges(aes(linetype = "solid",
                            color = w_sign,
                            size = abs(weight))) +
  geom_nodes(color = "green", size = 1)+
  geom_nodetext(aes(label = name),
                fontface = "bold") 
  
  
test <- igraph::graph_from_data_frame(d = connections,
                                           vertices = nodes,
                                           directed = T)
test$layout = l
V(test)$labels = as.factor(nodes$nodes)
MetamapsDB::ig2ggplot(test)
