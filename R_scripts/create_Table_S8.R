set.seed(017)

indir <- "Biomonitoring_2.0_Refined/"

source(paste(indir,"R_scripts/community_and_population_functions.R",sep=""))

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

#Setup metadata
metadata <- read.csv(paste(indir,"Data/metadata_with_clusters_elevation.csv",sep=""),row.names=1)
rownames(metadata) <- metadata$Sitename
metadata <- subset(metadata, clusters != 5)
cluster_sample_ID <- as.data.frame(cbind(metadata$Sitename,metadata$ID,metadata$cluster_name,as.numeric(metadata$Lat),as.numeric(metadata$Long)))
colnames(cluster_sample_ID) <- c("Sample","ID","cluster_region","latitude","longitude")
rownames(cluster_sample_ID) <- cluster_sample_ID$Sample
important_tables <- list(metadata=cluster_sample_ID)


variables_io <- paste("assign",sBp,"filter",read_filter,sample_filter,within_sample_filter,"method",dm_method,beta_component,sep="_")

target_dm_in <- readRDS(file=paste(indir,"RDS/target_dm_",variables_io,".RDS",sep=""))
community_dm_in <- readRDS(file=paste(indir,"RDS/community_dm_",variables_io,".RDS",sep=""))


community_correlations <- data.frame()
for (amp1 in names(community_dm_in)){
  for (type1 in names(community_dm_in[[amp1]])){
    dm1 <- community_dm_in[[amp1]][[type1]][["all"]]
    sites1 <- labels(dm1)
    for (amp2 in names(community_dm_in)){
      for (type2 in names(community_dm_in[[amp2]])){
        dm2 <- community_dm_in[[amp2]][[type2]][["all"]]
        sites2 <- labels(dm2)
        shared_sites <- intersect(sites1,sites2)
        dm1_sub <- dist_subset(dm1,shared_sites)
        dm2_sub <- dist_subset(dm2,shared_sites)
        corr <- mantel(dm1_sub,dm2_sub,method="spearman",permutations = 9999)
        community_correlations <- rbind(community_correlations,
                                        c(amp1,type1,amp2,type2,corr[["statistic"]],corr[["signif"]]))
      }
    }
  }
}


population_correlations <- data.frame()
for (amp1 in names(target_dm_in)){
  for (type1 in names(target_dm_in[[amp1]])){
    dm1 <- as.dist(as.matrix(target_dm_in[[amp1]][[type1]]))
    sites1 <- labels(dm1)
    for (amp2 in names(target_dm_in)){
      for (type2 in names(target_dm_in[[amp2]])){
        dm2 <- as.dist(as.matrix(target_dm_in[[amp2]][[type2]]))
        sites2 <- labels(dm2)
        shared_sites <- intersect(sites1,sites2)
        dm1_sub <- dist_subset(dm1,shared_sites)
        dm2_sub <- dist_subset(dm2,shared_sites)
        corr <- mantel(dm1_sub,dm2_sub,method="spearman",permutations = 9999)
        population_correlations <- rbind(population_correlations,
                                        c(amp1,type1,amp2,type2,corr[["statistic"]],corr[["signif"]]))
      }
    }
  }
}
write.csv(community_correlations,file=paste(indir,"Tables/community_correlations.csv",sep=""))
write.csv(population_correlations,file=paste(indir,"Tables/population_correlations.csv",sep=""))
