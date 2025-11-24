# This script simulates mutualistic coevolution on local networks.
#
# It outputs into Output folder a dataframe with reciprocity and strength 
# for all interactions.

rm(list=ls())

# import libraries
library(dplyr)


# define interaction network names
networks = c("01","02","03","04","05","06","07","08","09","10",
             "11","12","13","14","15","16","17","18","19","20",
             "21","22","23","24","25","26","27","28","29","30",
             "31","32")


# COEVOLUTION FUNCTION ----

# Function to model trait evolution in a mutualistic network (Guimaraes et al., 2017, Nature)

# Arguments:
# network - name of interaction network
# phi - selection gradient parameter
# alpha - sensitivity of selection to trait matching
# m - proportion of selection due to interactions
# tol - tolerance
# tmax - maximum number of timesteps
# reps - number of replicas

# Returns:
# df_out - dataframe with steady state reciprocity, strength and trait matching
#          for all species pairs
# df_z - dataframe with steady state trait values and theta values for all species

coevolution_mutualism = function(network, phi, alpha, m, tol, tmax, reps){
  
  # import network incidence matrix
  net = read.csv(paste0("Data/Minc_",network,".csv"))
  Minc = as.matrix(net[,2:ncol(net)])
  colnames(Minc) = NULL
  
  # number of resources, consumers and total species
  n_r = nrow(Minc)
  n_c = ncol(Minc)
  n_sp = n_r + n_c
  
  # build adjacency matrix
  f = rbind(cbind(matrix(0,n_r,n_r), Minc),
            cbind(t(Minc), matrix(0,n_c,n_c)))
  
  # initialise dataframe for storing reciprocity, strength, trait matching
  df_out = data.frame(expand.grid(replica=1:reps, resource=1:n_r, consumer=1:n_c),
                      R=NA,
                      S=NA,
                      TM=NA) %>%
    left_join(., data.frame(expand.grid(resource=1:n_r, consumer=1:n_c),
                            Minc=as.vector(Minc)), by=join_by(resource,consumer))
  
  # initialise dataframe for storing final trait values and theta values
  df_z = data.frame(rbind(expand.grid(replica=1:reps, guild="resources", species=1:n_r),
                          expand.grid(replica=1:reps, guild="consumers", species=1:n_c)),
                    z=NA,
                    theta=NA) 
  
  # loop through replicas
  for(r in 1:reps){
    
    # set seed to replica number
    set.seed(r*n_r*n_c)
    
    # sample theta values
    theta = runif(n_sp, min=0, max=10)
    
    # set initial trait values to theta
    z = theta
    
    # iterate through time
    for(t in 1:(tmax-1)){
      
      # compute differences in trait values for all species combinations
      z_dif = t(f*z) - f*z
      
      # compute the q matrix - evolutionary effects
      q = f*(exp(-alpha*(z_dif^2)))
      
      # standardize q matrix so that rows sum to 1
      q_n = q / apply(q,1,sum)
      
      # scale q_n by coevolutionary selection
      q_m = q_n * m
      
      # calculate selection differentials
      sel_dif = q_m * z_dif
      
      # calculate trait change as a result of interactions
      r_mut = phi * apply(sel_dif, 1, sum)
      
      # calculate trait change as a result of environment
      r_env = phi * (1 - m) * (theta - z)
      
      # calculate new trait value
      z_new = z + r_mut + r_env
      
      # calculate trait change from previous timestep
      dif = abs(z - z_new)
      
      # assign new trait value
      z = z_new
      
      # if trait change is smaller than tol, calculate output and stop iteration
      if(all((dif < tol))){
        
        # calculate reciprocity & strength
        
        # find mimimum value of evolutionary effect in each interaction
        qn_min = pmin(q_n,t(q_n))

        # find maximum value of evolutionary effect in each interaction
        qn_max = pmax(q_n,t(q_n))

        # calculate reciprocity
        R = qn_min / qn_max
        R[is.nan(R)] = 0
        
        # calculate strength
        Sqn = (qn_min+qn_max)/2

        # trait matching
        z_dif_all = t(matrix(1,n_sp,n_sp)*z) - matrix(1,n_sp,n_sp)*z
        TM = exp(-alpha*(z_dif_all^2))
        
        # store values
        df_out[df_out$replica==r,"R"] = as.vector(R[1:n_r,(n_r+1):n_sp])
        df_out[df_out$replica==r,"S"] = as.vector(Sqn[1:n_r,(n_r+1):n_sp])
        df_out[df_out$replica==r,"TM"] = as.vector(TM[1:n_r,(n_r+1):n_sp])

        df_z[df_z$replica==r,"z"] = z
        df_z[df_z$replica==r,"theta"] = theta
        
        # stop iterations
        break
      }
      
    } # t loop
    
  } # r loop
  
  return(list(df_out, df_z))
}


# COEVOLUTION SIMULATIONS ----

# coevolution parameters
phi = 0.5
alpha = 0.1
m = 0.5

# simulation parameters
tol = 0.0001
tmax = 10000
reps = 100

# sumulate mutualistic coevolution for all local networks
for(n in networks){
  
  # run coevolution
  out = coevolution_mutualism(n, phi, alpha, m, tol, tmax, reps)
  
  # write out result
  write.csv(out[[1]], paste0("Output/df_RS_",n,"_m",m,"_alpha",alpha,".csv"), row.names=FALSE)
  #write.csv(out[[2]], paste0("Output/df_Z_",n,"_m",m,"_alpha",alpha,".csv"), row.names=FALSE)
  
}
