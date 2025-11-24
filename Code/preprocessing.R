# This script calculates local network metrics and species properties 
# from empirical interaction networks (Jauker et al., 2019, Ecology)
# stored as incidence matrices in Data folder.
#
# It outputs two dataframes into Data folder:
# - data_networks.csv - local network metrics
# - data_species.csv - local species properties

rm(list=ls())

# import libraries
library(dplyr)
library(igraph)
library(FactoMineR)


# define interaction network names
networks = c("01","02","03","04","05","06","07","08","09","10",
             "11","12","13","14","15","16","17","18","19","20",
             "21","22","23","24","25","26","27","28","29","30",
             "31","32")


# function for computing nestedness
compute_nestedness = function(B){
  
  # Get number of rows and columns
  nrows <- nrow(B)
  ncols <- ncol(B)
  
  # Compute nestedness of rows
  nestedness_rows <- 0
  for(i in 1:(nrows-1)){
    for(j in (i+1): nrows){
      
      c_ij <- sum(B[i,] * B[j,])      # Number of interactions shared by i and j
      k_i <- sum(B[i,])               # Degree of node i
      k_j <- sum(B[j,])               # Degree of node j
      
      if (k_i == 0 || k_j==0) {next}  # Handle case if a node is disconnected
      
      o_ij <- c_ij / min(k_i, k_j)    # Overlap between i and j
      
      nestedness_rows <- nestedness_rows + o_ij
    }
  }
  
  # Compute nestedness of columns
  nestedness_cols <- 0
  for(i in 1: (ncols-1)){
    for(j in (i+1): ncols){
      
      c_ij <- sum(B[,i] * B[,j])      # Number of interactions shared by i and j
      k_i <- sum(B[,i])               # Degree of node i
      k_j <- sum(B[,j])               # Degree of node j
      if (k_i == 0 || k_j==0) {next}  # Handle case if a node is disconnected.
      
      o_ij <- c_ij / min(k_i, k_j)    # Overlap between i and j
      
      nestedness_cols <- nestedness_cols + o_ij         
    }
  }
  
  # Compute nestedness of the network
  nestedness <- (nestedness_rows + nestedness_cols) / ((nrows * (nrows - 1) / 2) + (ncols * (ncols - 1) / 2))
  
  return(nestedness)
}

# function for probabilistic cell null model
# (for standardising nestedness and modularity)
cell_model = function(d,t_max){
  
  rows = nrow(d)
  columns = ncol(d)
  
  t=1
  
  null_nest = rep(NA,t_max)
  null_mod = rep(NA,t_max)
  
  while (t <= t_max){
    
    PR = matrix(0, rows, 1)
    PC = matrix(0, columns, 1)
    B = matrix(0, rows, columns)
    
    for (i in 1:rows){
      number_ones=0
      for (j in 1:columns){
        if(d[i,j] == 1){
          number_ones=number_ones+1
        }
      }
      PR[i] = number_ones/columns
    }
    
    for (j in 1:columns){
      number_ones=0
      for (i in 1:rows){
        if(d[i,j] == 1){
          number_ones=number_ones+1
        }
      }
      PC[j] = number_ones/rows
    }
    
    for (i in 1:rows){
      for (j in 1:columns){
        p = ( PR[i]+PC[j] )/2;
        r = runif(1) 
        if(r < p){  
          B[i,j] = 1;
        }
      }
      
    }
    
    # remove unconected species if present
    B = bipartite::empty(B)
    
    # skip if network has fewer than 2 rows or columns
    if(nrow(B)<2 || ncol(B)<2) {next}
    
    # compute nestedness
    null_nest[t] = compute_nestedness(B)
    
    # convert to igraph and calculate modularity
    graph = graph_from_incidence_matrix(B)
    modules = cluster_louvain(graph)
    null_mod[t] = modularity(modules)
    
    t=t+1
    
  }
  
  # unlist results 
  null_nest = unlist(null_nest)
  null_mod = unlist(null_mod)
  
  return(list(null_nest, null_mod))
  
}


# initialise dataframe for storing local network metrics
data_networks = data.frame(network=networks,
                          n_species=NA,      # number of species
                          n_interactions=NA, # number of interactions
                          connectance=NA,    # connectance
                          nestedness_obs=NA, # observed nestedness
                          nestedness_Z=NA,   # standardised nestedness
                          modularity_obs=NA, # observed modulatiry
                          modularity_Z=NA)   # standardised modularity

# initialise dataframe for storing local species data
data_species = data.frame(network=character(),
                          guild=character(),
                          species=integer(),
                          species_name=character(),
                          degree=integer())


# compute network metrics and species degree
for(n in networks){
  
  set.seed(1)
  
  # import network incidence matrix
  net = read.csv(paste0("Data/Minc_",n,".csv"))
  Minc = as.matrix(net[,2:ncol(net)])
  rownames(Minc) = net[,1]
  
  # store number of spcecies, number of interactions, connectance
  data_networks[data_networks$network==n,"n_species"] = nrow(Minc) + ncol(Minc)
  data_networks[data_networks$network==n,"n_interactions"] = sum(Minc)
  data_networks[data_networks$network==n,"connectance"] = sum(Minc) / (nrow(Minc)*ncol(Minc))
  
  # nestedness - observed
  nestedness_obs = compute_nestedness(Minc)
  
  # modularity - observed
  graph = graph_from_incidence_matrix(Minc)
  modules = cluster_louvain(graph)
  modularity_obs = modularity(modules)
  
  # run cell null model with current network and 100 replicates
  null_output = cell_model(Minc,100)
  null_nest = null_output[[1]]
  null_mod = null_output[[2]]
  
  # compute mean of nestedness and modularity
  mean_nest = mean(null_nest, na.rm=TRUE)
  mean_mod = mean(null_mod, na.rm=TRUE)
  
  # compute standard devation of nestedness and modularity
  sd_nest = sd(null_nest, na.rm=TRUE)
  sd_mod = sd(null_mod, na.rm=TRUE)
  
  # compute zscore for nestedness and modularity
  z_score_nest = (nestedness_obs - mean_nest) / sd_nest
  z_score_mod = (modularity_obs - mean_mod) / sd_mod
  
  # store nestedness and modularity
  data_networks[data_networks$network==n,"nestedness_obs"] = nestedness_obs
  data_networks[data_networks$network==n,"nestedness_Z"] = z_score_nest
  data_networks[data_networks$network==n,"modularity_obs"] = modularity_obs
  data_networks[data_networks$network==n,"modularity_Z"] = z_score_mod
  
  # store local species data
  data_species = rbind(data_species,
                       data.frame(network=n, guild="resource", species=1:nrow(Minc),
                                  species_name=rownames(Minc), 
                                  degree=rowSums(Minc), row.names=NULL),
                       data.frame(network=n, guild="consumer", species=1:ncol(Minc),
                                  species_name=colnames(Minc), 
                                  degree=colSums(Minc), row.names=NULL))
}


# Principal Component Analysis
pca = PCA(data_networks %>% select(n_species, n_interactions, connectance), graph=TRUE)
pca$eig
pca$var$coord

# store PC1 in dataframe
data_networks = data_networks %>%
  mutate(PC1 = pca$ind$coord[,1])


# write out postprocessing results
write.csv(data_networks, paste0("Data/data_networks.csv"), row.names=FALSE)
write.csv(data_species, paste0("Data/data_species.csv"), row.names=FALSE)
