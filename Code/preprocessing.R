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
library(nnet)


# define interaction network names
networks = c("01","02","03","04","05","06","07","08","09","10",
             "11","12","13","14","15","16","17","18","19","20",
             "21","22","23","24","25","26","27","28","29","30",
             "31","32")


# function for computing nestedness (Fortuna et al., 2019, Evolution)
compute_nestedness = function(B){
  
  # Get number of rows and columns
  nrows = nrow(B)
  ncols = ncol(B)
  
  # Compute nestedness of rows
  nestedness_rows = 0
  for(i in 1:(nrows-1)){
    for(j in (i+1): nrows){
      
      c_ij = sum(B[i,] * B[j,])      # Number of interactions shared by i and j
      k_i = sum(B[i,])               # Degree of node i
      k_j = sum(B[j,])               # Degree of node j
      
      if (k_i == 0 || k_j==0) {next}  # Handle case if a node is disconnected
      
      o_ij = c_ij / min(k_i, k_j)    # Overlap between i and j
      
      nestedness_rows = nestedness_rows + o_ij
    }
  }
  
  # Compute nestedness of columns
  nestedness_cols = 0
  for(i in 1: (ncols-1)){
    for(j in (i+1): ncols){
      
      c_ij = sum(B[,i] * B[,j])      # Number of interactions shared by i and j
      k_i = sum(B[,i])               # Degree of node i
      k_j = sum(B[,j])               # Degree of node j
      if (k_i == 0 || k_j==0) {next}  # Handle case if a node is disconnected.
      
      o_ij = c_ij / min(k_i, k_j)    # Overlap between i and j
      
      nestedness_cols = nestedness_cols + o_ij         
    }
  }
  
  # Compute nestedness of the network
  nestedness = (nestedness_rows + nestedness_cols) / ((nrows * (nrows - 1) / 2) + (ncols * (ncols - 1) / 2))
  
  return(nestedness)
}

# functions for computing combined NODF nestedness (Song et al., 2017, J Anim Ecol)
nestedness_NODF = function(web){
  web[web > 0] = 1
  SA = nrow(web)
  SP = ncol(web)
  N = t(web) %*% web
  num = N
  num[lower.tri(num,diag=TRUE)]=1
  den = (matrix(1,nrow=SP,ncol=1)*diag(N))%*%matrix(1,nrow=1,ncol=SP)
  dele = den - t(den)
  dele[lower.tri(dele,diag=TRUE)] = 1
  num[dele == 0] = 0
  den = pmin(den,t(den))
  den[lower.tri(den,diag=TRUE)] = 1
  nes = num/den
  nes[lower.tri(nes,diag=TRUE)] = 0
  nes[is.na(nes)] = 0
  n1 = sum(nes)
  
  N = web %*% t(web)
  num = N
  num[lower.tri(num,diag=TRUE)]=1
  den = (matrix(1,nrow=SA,ncol=1)*diag(N))%*%matrix(1,nrow=1,ncol=SA)
  dele = den - t(den)
  dele[lower.tri(dele,diag=TRUE)] = 1
  num[dele ==0 ] = 0
  den = pmin(den,t(den))
  den[lower.tri(den,diag=TRUE)]=1
  nes = num/den
  nes[lower.tri(nes,diag=TRUE)] = 0
  nes[is.na(nes)] = 0
  n2 = sum(nes)
  out = 2*(n1 + n2) / (SA*(SA-1)+SP*(SP-1))
  return(out)
}
max_nest = function(web){
  #binarize the interaction matrix
  web_binary = web
  web_binary[web_binary > 0] = 1
  #compute the number of pollinators, plants and interactions
  SA = nrow(web_binary)
  SP = ncol(web_binary)
  SI = floor(sum(web_binary))
  #initialize the interaction matrix with minimum requirements
  web_opt = matrix(0, nrow=SA, ncol=SP)  
  web_opt[1,] = 1
  web_opt[,1] = 1
  web_opt[2,2] = 1
  #counting the number of
  SI_left = SI-SP-SA
  if(SI_left>0){
    #search the best possible location
    for(j in 1:SI_left){
      #compare all possible locations and the maximum one
      position_potential = websearch_NODF(web_opt)
      nest_poten = c()
      for(i in 1:nrow(position_potential)) {
        web_poten = web_opt
        web_poten[position_potential[i,1],position_potential[i,2]] = 1
        nest_poten[i] = nestedness_NODF(web_poten)
      } 
      position_the = which.is.max(nest_poten)
      web_opt[position_potential[position_the,1],position_potential[position_the,2]] = 1
    }
    return(nestedness_NODF(web_opt))
  }
  #this is to prevent the trivial case
  else{
    return(-1)
  }
}
websearch_NODF = function(web){
  SA = nrow(web)
  SP = ncol(web)
  domain = web
  position = which(domain == 1, arr.ind=T)
  position = subset(position, position[,2] != 1)
  position = subset(position, position[,1] != 1)
  boundary = matrix(0, nrow=2*nrow(position),ncol=2)
  j=1
  #choose boundary points
  for(i in 1:nrow(position)){
    if(position[i,1]<nrow(domain)&&position[i,2]<ncol(domain))
      if(domain[position[i,1]+1,position[i,2]]+
         domain[position[i,1]-1,position[i,2]]+
         domain[position[i,1],position[i,2]+1]+
         domain[position[i,1],position[i,2]-1]<=3){
        boundary[j,1] = position[i,1]+1
        boundary[j,2] = position[i,2]
        boundary[j+1,1] = position[i,1]
        boundary[j+1,2] = position[i,2]+1
        j = j+2
      } 
  }
  #delete those with zero entries which entered as auxiliary in the first place
  keep = c()
  for(i in 1:nrow(boundary)){
    if(boundary[i,1]+boundary[i,2]>0) keep = append(keep,i)
  }
  boundary = boundary[keep,]
  #choose true boundary points
  stay = c()
  for(i in 1:nrow(boundary)){
    if(boundary[i,1]<SA&&boundary[i,2]<SP){
      if(domain[boundary[i,1]+1,boundary[i,2]]+
         domain[boundary[i,1]-1,boundary[i,2]]+
         domain[boundary[i,1],boundary[i,2]+1]+
         domain[boundary[i,1],boundary[i,2]-1]==2
         && domain[boundary[i,1],boundary[i,2]]==0){
        stay = append(stay,i)
      }
    }
  }
  boundary = boundary[stay,]
  return(boundary)
}
comb_nest = function(web,NODF,max_NODF){
  C = sum(web)/(ncol(web)*nrow(web))
  S = sqrt(ncol(web) * nrow(web) )
  out = NODF / (max_NODF * C * log10(S))
  return(out)
}


# initialise dataframe for storing local network metrics
data_networks = data.frame(network=networks,
                           n_species=NA,       # number of species
                           n_interactions=NA,  # number of interactions
                           connectance=NA,     # connectance
                           nestedness_obs=NA,  # observed nestedness
                           nestedness_NODFc=NA,# combined NODF nestedness
                           modularity_obs=NA)  # observed modulatiry

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
  
  # nestedness - observed NODF (Song et al, 2017, J Anim Ecol)
  nestedness_obs_NODF = nestedness_NODF(Minc)
  
  # nestedness - maximum NODF (Song et al, 2017, J Anim Ecol)
  max_NODF = max_nest(Minc)
  
  # nestedness - combined NODF (Song et al, 2017, J Anim Ecol)
  nestedness_comb_NODF = comb_nest(Minc,nestedness_obs_NODF,max_NODF)
  
  # modularity - observed
  graph = graph_from_incidence_matrix(Minc)
  modules = cluster_louvain(graph)
  modularity_obs = modularity(modules)
  
  # store nestedness and modularity
  data_networks[data_networks$network==n,"nestedness_obs"] = nestedness_obs
  data_networks[data_networks$network==n,"nestedness_NODFc"] = nestedness_comb_NODF
  data_networks[data_networks$network==n,"modularity_obs"] = modularity_obs

  # store local species data
  data_species = rbind(data_species,
                       data.frame(network=n, guild="resource", species=1:nrow(Minc),
                                  species_name=rownames(Minc), 
                                  degree=rowSums(Minc), row.names=NULL),
                       data.frame(network=n, guild="consumer", species=1:ncol(Minc),
                                  species_name=colnames(Minc), 
                                  degree=colSums(Minc), row.names=NULL))
}


# write out postprocessing results
write.csv(data_networks, paste0("Data/data_networks.csv"), row.names=FALSE)
write.csv(data_species, paste0("Data/data_species.csv"), row.names=FALSE)
