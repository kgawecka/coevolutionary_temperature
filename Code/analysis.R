# This script analyses the output of the coevolution model simulations
# and plots the manuscript figures.

rm(list=ls())

# import libraries
library(dplyr)
library(data.table)
library(ggplot2)
library(ggpubr)
library(ggExtra)
library(bipartite)
library(piecewiseSEM)
library(ggeffects)


# define palettes
col_red = "#e15759ff"
col_blue = "#4e79a7ff"
col_grey = "#808080ff"
col_green_light = "#5aa05033"
col_green_dark ="#30562aff"

# define interaction network names
networks = c("01","02","03","04","05","06","07","08","09","10",
             "11","12","13","14","15","16","17","18","19","20",
             "21","22","23","24","25","26","27","28","29","30",
             "31","32")


# POSTPROCESSING ----

# import patches, network and species data
data_patches_networks = left_join(read.csv("Data/data_patches.csv"), 
                                  read.csv("Data/data_networks.csv")) %>%
  mutate(network=ifelse(network %in% c("1","2","3","4","5","6","7","8","9"),paste0("0",network),network))
data_species = read.csv("Data/data_species.csv") %>%
  mutate(network=ifelse(network %in% c("1","2","3","4","5","6","7","8","9"),paste0("0",network),network))

# dataframe with resources only
data_species_res = data_species %>%
  filter(guild=="resource") %>%
  select(network, species, degree, species_name) %>%
  rename(resource="species", degree_res="degree", name_res="species_name")
# dataframes with consumers only
data_species_con = data_species %>%
  filter(guild=="consumer") %>%
  select(network, species, degree, species_name) %>%
  rename(consumer="species", degree_con="degree", name_con="species_name")


# define coevolution parameters
alpha = 0.1
m = 0.5

# initialise dataframe for storing all simulation output
df_RS = data.frame(network=character(),
                   replica=integer(),
                   resource=integer(),
                   consumer=integer(),
                   R=double(),
                   S=double(), 
                   TM=double(),
                   Minc=integer())

# import simulation output
for(n in networks){
  
  df_RS = rbind(df_RS, 
                fread(paste0("Output/df_RS_",n,"_m",m,"_alpha",alpha,".csv")) %>%
                  mutate(network=n))
}

# calculate degree symmetry
df_RS = df_RS %>%
  mutate(network=as.factor(network),
         Minc=as.factor(Minc)) %>%
  left_join(., data_species_res) %>%
  left_join(., data_species_con) %>%
  mutate(deg_sym = pmin(degree_res,degree_con) / pmax(degree_res,degree_con))

# average reciprocity and strength for each interaction (across all replicas)
df_RS_int = df_RS %>%
  group_by(network, resource, consumer, Minc) %>%
  summarise(R_mean=mean(R, na.rm=TRUE),
            S_mean=mean(S, na.rm=TRUE),
            TM_mean=mean(TM, na.rm=TRUE),
            R_sd=sd(R, na.rm=TRUE),
            S_sd=sd(S, na.rm=TRUE),
            TM_sd=sd(TM, na.rm=TRUE),
            deg_sym=mean(deg_sym)) %>%
  left_join(., data_patches_networks) %>%
  left_join(., data_species_res) %>%
  left_join(., data_species_con) %>%
  ungroup()

# average reciprocity and strength for each network (across all interactions and replicas)
df_RS_net = df_RS %>%
  group_by(network, Minc) %>%
  summarise(R_mean=mean(R, na.rm=TRUE),
            S_mean=mean(S, na.rm=TRUE),
            TM_mean=mean(TM, na.rm=TRUE),
            deg_sym_mean=mean(deg_sym, na.rm=TRUE)) %>%
  left_join(., data_patches_networks) %>%
  ungroup()


# PLOTS ----

# patch size and network metrics ----

# patch area vs network size and connectance

summary(lm(n_species~log(area_km2), data_patches_networks))
summary(lm(n_interactions~log(area_km2), data_patches_networks))
p1 = ggplot(data=data_patches_networks, aes(x=log(area_km2)))  +
  geom_smooth(aes(y=n_species), method=lm, col=col_green_dark) +
  geom_smooth(aes(y=n_interactions), method=lm, col=col_green_light) +
  geom_point(aes(y=n_species), col=col_green_dark) +
  geom_point(aes(y=n_interactions), col=col_green_light) +
  labs(x=expression(paste("ln(patch area) [", km^2,"]")),
       y="number of species (dark green)\nnumber interactions (light green)") +
  theme(panel.background=element_rect(fill="white", colour="grey"),
        panel.grid.major=element_blank(), panel.grid.minor=element_blank(),
        axis.title=element_text(size=8), axis.text=element_text(size=6))

summary(lm(connectance~log(area_km2), data_patches_networks))
p2 = ggplot(data=data_patches_networks, aes(x=log(area_km2), y=connectance))  +
  geom_point() +
  geom_smooth(method=lm, col="black") +
  labs(x=expression(paste("ln(patch area) [", km^2,"]")),
       y="connectance") +
  theme(panel.background=element_rect(fill="white", colour="grey"),
        panel.grid.major=element_blank(), panel.grid.minor=element_blank(),
        axis.title=element_text(size=8), axis.text=element_text(size=6))

ggarrange(p1, p2, nrow=1, ncol=2, labels=c("A","B"))


# dependence of nestedness and modularity on patch area

m_nest = lm(nestedness_obs~n_species+connectance, data_patches_networks)
m_mod = lm(modularity_obs~n_species+connectance, data_patches_networks)
summary(m_nest)
summary(m_mod)

nest_species = ggpredict(m_nest, terms="n_species", type="fixed")
nest_conn    = ggpredict(m_nest, terms="connectance", type="fixed")
mod_species  = ggpredict(m_mod, terms="n_species", type="fixed")
mod_conn     = ggpredict(m_mod, terms="connectance", type="fixed")

data_patches_networks$partial_nest_nspecies = resid(m_nest) + 
  predict(m_nest, newdata=data.frame(n_species=data_patches_networks$n_species,
                                     connectance = mean(data_patches_networks$connectance)))
data_patches_networks$partial_nest_conn = resid(m_nest) +
  predict(m_nest, newdata=data.frame(n_species=mean(data_patches_networks$n_species),
                                     connectance=data_patches_networks$connectance))
data_patches_networks$partial_mod_nspecies = resid(m_mod) + 
  predict(m_mod, newdata=data.frame(n_species=data_patches_networks$n_species,
                                    connectance = mean(data_patches_networks$connectance)))
data_patches_networks$partial_mod_conn = resid(m_mod) +
  predict(m_mod, newdata=data.frame(n_species=mean(data_patches_networks$n_species),
                                    connectance=data_patches_networks$connectance))

p1 = ggplot() +
  geom_point(data=data_patches_networks, aes(x=n_species, y=partial_nest_nspecies)) +
  geom_line(data=nest_species, aes(x=x, y=predicted)) +
  geom_ribbon(data=nest_species, aes(x=x, ymin=conf.low, ymax=conf.high), alpha = 0.2) +
  labs(x="number of species", y="nestedness") +
  theme(panel.background=element_rect(fill="white", colour="grey"),
        panel.grid.major=element_blank(), panel.grid.minor=element_blank(),
        axis.title=element_text(size=8), axis.text=element_text(size=6))

p2 = ggplot() +
  geom_point(data=data_patches_networks, aes(x=connectance, y=partial_nest_conn)) +
  geom_line(data=nest_conn, aes(x=x, y=predicted)) +
  geom_ribbon(data=nest_conn, aes(x=x, ymin=conf.low, ymax=conf.high), alpha = 0.2) +
  labs(x="connectance", y="nestedness") +
  theme(panel.background=element_rect(fill="white", colour="grey"),
        panel.grid.major=element_blank(), panel.grid.minor=element_blank(),
        axis.title=element_text(size=8), axis.text=element_text(size=6))

p3 = ggplot() +
  geom_point(data=data_patches_networks, aes(x=n_species, y=partial_mod_nspecies)) +
  geom_line(data=mod_species, aes(x=x, y=predicted)) +
  geom_ribbon(data=mod_species, aes(x=x, ymin=conf.low, ymax=conf.high), alpha = 0.2) +
  labs(x="number of species", y="modularity") +
  theme(panel.background=element_rect(fill="white", colour="grey"),
        panel.grid.major=element_blank(), panel.grid.minor=element_blank(),
        axis.title=element_text(size=8), axis.text=element_text(size=6))

p4 = ggplot() +
  geom_point(data=data_patches_networks, aes(x=connectance, y=partial_mod_conn)) +
  geom_line(data=mod_conn, aes(x=x, y=predicted)) +
  geom_ribbon(data=mod_conn, aes(x=x, ymin=conf.low, ymax=conf.high), alpha = 0.2) +
  labs(x="connectance", y="modularity") +
  theme(panel.background=element_rect(fill="white", colour="grey"),
        panel.grid.major=element_blank(), panel.grid.minor=element_blank(),
        axis.title=element_text(size=8), axis.text=element_text(size=6))

ggarrange(p1, p2, p3, p4, nrow=2, ncol=2, labels=c("A","B","C","D"))


# patch locations and network sizes

ggplot(data=data_patches_networks, 
       aes(x=E/1000, y=N/1000, col=n_interactions, size=area_km2)) +
  geom_point() +
  geom_text(aes(label=network), hjust=-0.1, vjust=-0.1, col="grey", size=3) +
  scale_color_gradient(low=col_green_light, high=col_green_dark) +
  scale_size_continuous(range=c(1,5)) +
  coord_fixed() +
  scale_x_continuous(breaks=seq(min(data_patches_networks$E/1000),
                                max(data_patches_networks$E/1000),5)) +
  scale_y_continuous(breaks=seq(min(data_patches_networks$N/1000),
                                max(data_patches_networks$N/1000),5)) +
  labs(col="number of\ninteractions", size="patch area\n(km2)") +
  theme(axis.title=element_blank(), axis.text=element_blank(),
        axis.ticks=element_blank(),
        panel.background=element_rect(fill=NA, colour="lightgrey"),
        panel.grid=element_blank(),
        plot.background=element_rect(fill=NA),
        legend.title=element_text(size=8), legend.text=element_text(size=6), 
        legend.key=element_blank(), 
        legend.position="none",
        legend.box="vertical", legend.spacing.y=unit(1,"mm"))


# coevolutionary temperature mosaic ----

# define focal interaction
res_name = "Knautia.arvensis"
con_name = "Episyrphus.balteatus"

# subset data for plotting - interaction average
data_plot_int = filter(df_RS_int, Minc==1) %>%
  mutate(focal_int=ifelse(paste0(name_res,name_con)==paste0(res_name,con_name),
                          "K.arvensis-E.balteatus","other interactions"))

# subset data for plotting - network average
data_plot_net = filter(df_RS_net, Minc==1)

ggplot(data=data_plot_net, 
       aes(x=E/1000, y=N/1000, col=R_mean, size=S_mean)) +
  geom_point() +
  geom_text(aes(label=network), hjust=-0.1, vjust=-0.1, col="grey", size=3) +
  scale_color_gradient(low=col_blue, high=col_red) +
  scale_size_continuous(range=c(1,5)) +
  coord_fixed() +
  scale_x_continuous(breaks=seq(min(data_patches_networks$E/1000),
                                max(data_patches_networks$E/1000),5)) +
  scale_y_continuous(breaks=seq(min(data_patches_networks$N/1000),
                                max(data_patches_networks$N/1000),5)) +
  labs(col="average\nreciprocity", size="average\nstrength") +
  theme(axis.title=element_blank(), axis.text=element_blank(),
        axis.ticks=element_blank(),
        panel.background=element_rect(fill=NA, colour="lightgrey"),
        panel.grid=element_blank(),
        plot.background=element_rect(fill=NA),
        legend.title=element_text(size=8), legend.text=element_text(size=6), 
        legend.key=element_blank(), legend.position="bottom", 
        legend.box="vertical", legend.spacing.y=unit(1,"mm"))


plot_net = "15" # 15, 24

p = ggMarginal(
  ggplot(data=filter(data_plot_int, network==plot_net), 
         aes(x=S_mean, y=R_mean, fill=R_mean, size=S_mean)) +
    geom_point(shape=21, alpha=0.4) +
    geom_point(data=filter(data_plot_int, network==plot_net,
                           name_res==res_name,
                           name_res==con_name), shape=24) +
    coord_fixed(ratio=1) +
    scale_fill_gradient(low=col_blue, high=col_red, limits=c(0,1)) +
    scale_size_continuous(limits=c(0,1), range=c(0.1,3)) +
    scale_x_continuous(limits=c(0,1), breaks=c(0,1,1)) +
    scale_y_continuous(limits=c(0,1), breaks=c(0,1,1)) +
    labs(x="strength", y="reciprocity") +
    theme(panel.background=element_rect(fill="white", colour="grey"),
          panel.grid=element_blank(),
          axis.text=element_text(size=6, colour="grey"), axis.ticks=element_blank(),
          axis.title.x=element_text(size=8, margin=margin(t=-4)),
          axis.title.y=element_text(size=8, margin=margin(r=-2)),
          legend.position="none"),
  type="histogram", fill=col_grey, col=col_grey)

net = read.csv(paste0("Data/Minc_",plot_net,".csv"))
Minc = as.matrix(net[,2:ncol(net)])
rownames(Minc) = net[,1]

lower_color = rep("grey", nrow(Minc))
names(lower_color) = rownames(Minc)
lower_color[res_name] = "black"
higher_color = rep("grey", ncol(Minc))
names(higher_color) = colnames(Minc)
higher_color[con_name] <- "black"

plotweb(Minc,
        lower_color = lower_color, higher_color = higher_color,
        link_color = "grey",
        higher_labels = FALSE, lower_labels = FALSE)


# effect of interaction degree ratio ----

p1 = ggplot(data=data_plot_int, 
            aes(x=deg_sym, y=R_mean, 
                shape=focal_int, col=focal_int, alpha=focal_int)) +
  geom_errorbar(aes(ymin=R_mean-R_sd, ymax=R_mean+R_sd), 
                width=0, position=position_dodge()) +
  geom_point(position=position_dodge()) +
  scale_color_manual(values=c("black","grey")) +
  scale_shape_manual(values=c(17,16)) +
  scale_alpha_manual(values=c(1,0.1)) +
  guides(alpha=guide_legend(override.aes=list(alpha=1))) +
  lims(x=c(0,1), y=c(0,1)) +
  labs(x="interaction degree symmetry", y="reciprocity", 
       shape=NULL, color=NULL, alpha=NULL) +
  theme(panel.background=element_rect(fill="white", colour="grey"),
        panel.grid=element_blank(),
        axis.title=element_text(size=8), axis.text=element_text(size=6),
        legend.text=element_text(size=8), 
        legend.key=element_blank(), legend.position="bottom")

p2 = ggplot(data=data_plot_int, 
            aes(x=deg_sym, y=S_mean,
                shape=focal_int, col=focal_int, alpha=focal_int)) +
  geom_errorbar(aes(ymin=S_mean-S_sd, ymax=S_mean+S_sd), 
                width=0, position=position_dodge()) +
  geom_point(position=position_dodge()) +
  scale_color_manual(values=c("black","grey")) +
  scale_shape_manual(values=c(17,16)) +
  scale_alpha_manual(values=c(1,0.1)) +
  guides(alpha=guide_legend(override.aes=list(alpha=1))) +
  lims(x=c(0,1), y=c(0,1)) +
  labs(x="interaction degree symmetry", y="strength", 
       shape=NULL, color=NULL, alpha=NULL) +
  theme(panel.background=element_rect(fill="white", colour="grey"),
        panel.grid=element_blank(),
        axis.title=element_text(size=8), axis.text=element_text(size=6),
        legend.text=element_text(size=8), 
        legend.key=element_blank(), legend.position="bottom")

ggarrange(p1, p2, nrow=2, labels=c("A","B"), common.legend=TRUE, legend="bottom")


# structural equation model ----

# subset data
m_data = df_RS_net %>%
  filter(., Minc==1) %>% 
  mutate(area_log=log(area_km2))

# area | n_species & connectance | nestedness & modularity | R & S
sem1 = psem(
  # causal paths:
  lm(data=m_data, n_species ~ area_log),
  lm(data=m_data, connectance ~ area_log),
  lm(data=m_data, nestedness_obs ~ n_species + connectance),
  lm(data=m_data, modularity_obs ~ n_species + connectance),
  lm(data=m_data, R_mean ~ n_species + connectance + nestedness_obs + modularity_obs),
  lm(data=m_data, S_mean ~ n_species + connectance + nestedness_obs + modularity_obs),
  # covariance terms:
  n_species %~~% connectance,
  nestedness_obs %~~% modularity_obs
)
summary(sem1)
summary(sem1)$dTable


# area | n_species & connectance | nestedness & modularity | R & S
sem2 = psem(
  # causal paths:
  lm(data=m_data, n_species ~ area_log),
  lm(data=m_data, connectance ~ area_log),
  lm(data=m_data, nestedness_obs ~ area_log + n_species + connectance),
  lm(data=m_data, modularity_obs ~ area_log + n_species + connectance),
  lm(data=m_data, R_mean ~ n_species + connectance + nestedness_obs + modularity_obs),
  lm(data=m_data, S_mean ~ n_species + connectance + nestedness_obs + modularity_obs),
  # covariance terms:
  n_species %~~% connectance,
  nestedness_obs %~~% modularity_obs
)
summary(sem2)
summary(sem2)$dTable


# area | n_species & connectance | nestedness & modularity | deg_sym | R & S
sem3 = psem(
  lm(data=m_data, n_species ~ area_log),
  lm(data=m_data, connectance ~ area_log),
  lm(data=m_data, nestedness_obs ~ n_species + connectance),
  lm(data=m_data, modularity_obs ~ n_species + connectance),
  lm(data=m_data, deg_sym_mean ~ n_species + connectance + nestedness_obs + modularity_obs),
  lm(data=m_data, R_mean ~ deg_sym_mean + n_species + connectance + nestedness_obs + modularity_obs),
  lm(data=m_data, S_mean ~ deg_sym_mean + n_species + connectance + nestedness_obs + modularity_obs),
  n_species %~~% connectance,
  nestedness_obs %~~% modularity_obs
)
summary(sem3)
summary(sem3)$dTable


# area | n_species & connectance | nestedness (NODFc) & modularity | R & S
sem4 = psem(
  # causal paths:
  lm(data=m_data, n_species ~ area_log),
  lm(data=m_data, connectance ~ area_log),
  lm(data=m_data, nestedness_NODFc ~ n_species + connectance),
  lm(data=m_data, modularity_obs ~ n_species + connectance),
  lm(data=m_data, R_mean ~ n_species + connectance + nestedness_NODFc + modularity_obs),
  lm(data=m_data, S_mean ~ n_species + connectance + nestedness_NODFc + modularity_obs),
  # covariance terms:
  n_species %~~% connectance,
  nestedness_NODFc %~~% modularity_obs
)
summary(sem4)
summary(sem4)$dTable


# principal component analysis ----

pca = PCA(m_data %>% 
            select(n_species, connectance, nestedness_NODFc, modularity_obs,
                   area_log, R_mean, S_mean) %>%
            rename("patch area"=area_log, "species richness"=n_species, 
                   "nestedness (NODFc)"=nestedness_NODFc, modularity=modularity_obs, 
                   reciprocity=R_mean, strength=S_mean), 
          graph=TRUE)
