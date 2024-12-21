library(ggpubr)

set.seed(017)

indir <- "../"

source("community_and_population_functions.R")

parameter_list <- readRDS(paste(indir,"parameters.RDS",sep=""))
read_filter <- parameter_list[["read_filter"]]
sBp <- parameter_list[["sBp"]]
min_hap <- parameter_list[["min_hap"]]
sample_filter <- parameter_list[["sample_filter"]]
within_sample_filter <- parameter_list[["within_sample_filter"]]
dm_method <- parameter_list[["dm_method"]]
beta_component <- parameter_list[["beta_component"]]
amps <- parameter_list[["amps"]]
d_amp <- parameter_list[["d_amp"]]
min_sites <- parameter_list[["min_sites"]]

variables_io <- paste("assign",sBp,"filter",read_filter,sample_filter,within_sample_filter,"method",dm_method,beta_component,sep="_")

#Load RDS files
important_tables_in <- readRDS(file=paste(indir,"RDS/important_tables_",variables_io,".RDS",sep=""))
cluster_ID <- important_tables_in[["metadata"]]

community_dm_in <- readRDS(file=paste(indir,"RDS/community_dm_",variables_io,".RDS",sep=""))


#FIGURE PLOTTING COMMUNITY
neighbors <- 15
dim <- 2
com_umap_data <- list()
com_umap_plots <- list()
com_violin_plots <- list()
com_permanova_list <- list()
com_permanova_pair_list <- list()
com_dispersion_list <- list()
com_dispersion_pair_list <- list()
com_dispersion_data_list <- list()
for (amp in names(community_dm_in)){
  for (type in names(community_dm_in[[amp]])){
    for (beta_type in c(beta_component)){
      com <- community_dm_in[[amp]][[type]][[beta_type]]
      com_dist <- as.dist(as.matrix(com))
      
      #UMAP
      com_umap_data[[amp]][[type]][["default"]][[beta_type]] <- create_umap(com_dist,neighbors,dim)
      region_clusters <- cluster_ID[,"cluster_region",drop=FALSE]
      com_umap_plots[[amp]][[type]][["default"]][[beta_type]]<- create_umap_plot(com_umap_data[[amp]][[type]][["default"]][[beta_type]],region_clusters)
      
      com_umap_data[[amp]][[type]][["max"]][[beta_type]] <- create_umap(com_dist,length(labels(com))-1,dim)
      com_umap_plots[[amp]][[type]][["max"]][[beta_type]]<- create_umap_plot(com_umap_data[[amp]][[type]][["max"]][[beta_type]],region_clusters)
      
      #com_violin_plots[[amp]][[type]][[beta_type]]<- create_violin_plot(com_dist,region_clusters,rownames(com))
      
      #Statistics
      com_permanova_list[[amp]][[type]][[beta_type]] <- adonis2(com_dist~as.factor(region_clusters[rownames(com_umap_data[[amp]][[type]][["default"]][[beta_type]]),]),data=region_clusters, permutations = 9999)
      com_permanova_pair_list[[amp]][[type]][[beta_type]] <- pairwise.adonis2(com_dist, region_clusters[rownames(com_umap_data[[amp]][[type]][["default"]][[beta_type]]),], p.method = "fdr")
      dispersion <- betadisper(com_dist,region_clusters[rownames(com_umap_data[[amp]][[type]][["default"]][[beta_type]]),])
      com_dispersion_list[[amp]][[type]][[beta_type]] <- permutest(dispersion,permutations = 9999)
      com_dispersion_pair_list[[amp]][[type]][[beta_type]] <- TukeyHSD(dispersion)
      com_dispersion_data_list[[amp]][[type]][[beta_type]] <-dispersion 
    }
  }
}


#Community tables
permanova_table <- data.frame()
for (amp in names(com_permanova_list)){
  for (data_type in names(com_permanova_list[[amp]])){
    permanova_table <- rbind(permanova_table, 
                             as.data.frame(com_permanova_list[[amp]][[data_type]][[beta_component]]))
  }
}

dispersion_table <- data.frame()
for (amp in names(com_dispersion_list)){
  for (data_type in names(com_dispersion_list[[amp]])){
    dispersion_table <- rbind(dispersion_table, 
                              as.data.frame(com_dispersion_list[[amp]][[data_type]][[beta_component]][["tab"]]))
  }
}
region_pairs <- c("Northeast","Central","Southeast","Central","West","Central",
                  "Southeast","Northeast","West","Northeast","West","Southeast")
permanova_pair_table <- c()
for (amp in names(com_permanova_pair_list)){
  for (data_type in names(com_permanova_pair_list[[amp]])){
    for (x in seq(1,length(region_pairs),2)){
      current_data <-com_permanova_pair_list[[amp]][[data_type]][[beta_component]] 
      r <- region_pairs[x]
      c <- region_pairs[x+1]
      permanova_pair_table <- rbind(permanova_pair_table,
                                    c(amp,data_type,paste(region_pairs[x+1],region_pairs[x],sep=" "),
                                      current_data[["F"]][r,c],current_data[["R2"]][r,c],
                                      current_data[["p.value"]][r,c],current_data[["p.adjust"]][r,c]))
    }
  }
}

permanova_pair_table <- as.data.frame(permanova_pair_table)
colnames(permanova_pair_table) <- c("Amplicon","Type","Region pairs","F","R2","Pr(>F)","Pr(>F) adjusted")

write.csv(permanova_table,paste(indir,"Tables/TableS2_community_permanova_",variables_io,".csv",sep=""))
write.csv(dispersion_table,paste(indir,"Tables/TableS3_community_dispersion_",variables_io,".csv",sep=""))
write.csv(permanova_pair_table,paste(indir,"Tables/TableS4_community_permanova_pairs_",variables_io,".csv",sep=""))
