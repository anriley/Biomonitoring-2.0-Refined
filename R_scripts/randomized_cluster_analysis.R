library(geosphere)
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

random_clusters <- readRDS(paste(indir,"RDS/randomized_clusters_",variables_io,".RDS",sep=""))
target_dm <- readRDS(paste(indir,"RDS/target_dm_",variables_io,".RDS",sep=""))
important_tables_in <- readRDS(file=paste(indir,"RDS/important_tables_",variables_io,".RDS",sep=""))

cluster_ID <- important_tables_in[["metadata"]]
lat_long <- as.data.frame(cbind(as.numeric(cluster_ID$longitude),as.numeric(cluster_ID$latitude)))
colnames(lat_long) <- c("longitude","latitude")
rownames(lat_long) <- cluster_ID$ID
lat_long_dm <- distm(lat_long, fun = distGeo)
rownames(lat_long_dm) <- rownames(lat_long)
colnames(lat_long_dm) <- rownames(lat_long)

region_clusters <- c("Northeast","Southeast","Central","West")

all_random_data <- data.frame()
for (amp in names(random_clusters)){
  site_overlap <- intersect(rownames(target_dm[[amp]][["is"]]),rownames(target_dm[[amp]][["otu"]]))
  for (data_type in names(random_clusters[[amp]])){
    for (cluster in names(random_clusters[[amp]][[data_type]])){
      random_otu_tables <- setup_cluster_tables(random_clusters[[amp]][[data_type]][[cluster]],important_tables_in[[amp]],hap=min_hap)
      random_otu_dm <- list()
      for (sp in names(random_otu_tables)){#Loop through each OTU
        current_table.dist<-create_beta_div_matrices(random_otu_tables[[sp]],method=dm_method)
        random_otu_dm[["all"]][[sp]] <- current_table.dist[["all"]]
        random_otu_dm[["nested"]][[sp]] <- current_table.dist[["nested"]]
        random_otu_dm[["turnover"]][[sp]] <- current_table.dist[["turnover"]]
      }
      random_otu_dm_selected <-list()
      for (sp in names(random_otu_dm[[beta_component]])){
        sites <- labels(random_otu_dm[[beta_component]][[sp]])
        sites_in_dm <- cluster_ID[sites,]
        region_counts <- as.data.frame(table(sites_in_dm[["cluster_region"]]))
        region_present <- subset(region_counts, Freq >= min_sites)
        if (nrow(region_present)>1){
          random_otu_dm_selected[[beta_component]][[sp]] <-random_otu_dm[[beta_component]][[sp]]
        }
      }
      
      average_random_dm_temp <- average_distance_matrices(random_otu_dm_selected[[beta_component]])
      overlap_sites2 <- intersect(site_overlap,rownames(average_random_dm_temp))
      average_random_dm <- average_random_dm_temp[overlap_sites2,overlap_sites2]
      permanova_random <-  adonis2(average_random_dm~as.factor(cluster_ID[rownames(average_random_dm),]$cluster_region),
                                   data=cluster_ID[rownames(average_random_dm),], permutations = 9999)
      all_random_data <- rbind(all_random_data,c(amp,data_type,cluster,"F",permanova_random[["F"]][1],permanova_random[["Pr(>F)"]][1]))
      average_random_dm_ordered <- average_random_dm[order(rownames(average_random_dm)),order(colnames(average_random_dm))]
      lat_long_dm_reduced <- as.matrix(dist_subset(lat_long_dm,rownames(average_random_dm)))
      lat_long_dm_ordered <- lat_long_dm_reduced[order(rownames(lat_long_dm_reduced)),order(colnames(lat_long_dm_reduced))]
      all_mantel <- mantel(as.dist(average_random_dm_ordered),as.dist(lat_long_dm_ordered),method="spearman", perm=9999)
      all_random_data <- rbind(all_random_data,c(amp,data_type,cluster,"All",all_mantel[["statistic"]],all_mantel[["signif"]]))
      for (region in region_clusters){
        r <- rownames(subset(cluster_ID,cluster_region==region))
        r2 <- intersect(r,rownames(average_random_dm))
        average_random_dm_region <- average_random_dm[r2,r2]
        average_random_dm_region_ordered <- average_random_dm_region[order(rownames(average_random_dm_region)),order(colnames(average_random_dm_region))]
        lat_long_dm_region <- as.matrix(dist_subset(lat_long_dm,rownames(average_random_dm_region_ordered)))
        lat_long_dm_region_ordered <- lat_long_dm_region[order(rownames(lat_long_dm_region)),order(colnames(lat_long_dm_region))]
        region_mantel <- mantel(as.dist(average_random_dm_region_ordered),as.dist(lat_long_dm_region_ordered),method="spearman", perm=9999)
        all_random_data <- rbind(all_random_data,c(amp,data_type,cluster,region,region_mantel[["statistic"]],region_mantel[["signif"]]))
        
      }
    }
  }
}

colnames(all_random_data) <- c("Amplicon","Type","Sample","Region","Statistic","p_value")
write.csv(all_random_data,paste(indir,"Tables/randomized_data_",variables_io,".csv",sep=""))
