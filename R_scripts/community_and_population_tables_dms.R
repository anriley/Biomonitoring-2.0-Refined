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

#Setup metadata
metadata <- read.csv(paste(indir,"Data/metadata_with_clusters_elevation.csv",sep=""),row.names=1)
rownames(metadata) <- metadata$Sitename
metadata <- subset(metadata, clusters != 5)
cluster_ID <- as.data.frame(unique(cbind(metadata$ID,metadata$cluster_name,as.numeric(metadata$Lat),as.numeric(metadata$Long))))
colnames(cluster_ID) <- c("ID","cluster_region","latitude","longitude")
rownames(cluster_ID) <- cluster_ID$ID
important_tables <- list(metadata=cluster_ID)


variables_io <- paste("assign",sBp,"filter",read_filter,sample_filter,within_sample_filter,"method",dm_method,beta_component,sep="_")

#Get community tables output
community_tables <- list()
otu_tables <- list()
isolated_species_tables <- list()

community_dm <- list()
isolated_species_dm <- list()
otu_dm <- list()

#Setup all tables and cluster dissimilarity matrices
for (amp in amps){
  #Load in files
  esv_table_in <- t(read.csv(paste(indir,"Data/ESV_tables/ESV_rcf_",read_filter,"_",amp,".csv",sep=""),row.names=1))
  swarm_clusters <- read.delim(paste(indir,"Data/Clusters/Clusters_rcf_",read_filter,"_",amp,"/clusters_rcf_",read_filter,"_d_",d_amp[[amp]],"_",amp,".txt",sep=""),
                               header=FALSE)
  #taxonomy <- read.csv(paste(indir,"Data/Taxonomy/taxonomy_rcf_",read_filter,"_",amp,".csv",sep=""))
  species_clusters <- read.csv(paste(indir,"Data/Taxonomy_clusters/taxonomy_clusters_rcf_",read_filter,"_Species_",sBp,"_",amp,"_isolated.csv",sep=""))
  
  #Get clusters
  otu_zotus <- otu_clusters(swarm_clusters)
  #taxa_zotus <- taxonomy_clusters(taxonomy,"Species",sBp)
  isolated_species_zotus <- isolate_species_clusters(species_clusters)
  
  #Manipulate ESV table (Filter, Remove sites, Combine replicates, etc.)
  esv_table <- esv_table_in[rowSums(esv_table_in[])>sample_filter,]
  esv_table[esv_table < within_sample_filter] <- 0
  esv_table <- as.data.frame(esv_table)
  esv_table_combined <- combine_samples(esv_table,metadata,"ID")
  important_tables[[amp]] <-esv_table_combined 
  
  #Combine ESVs into OTUs and Taxa for community table
  community_tables[[amp]] <- list(esv=esv_table_combined,
                                  otu=setup_OTU_table(otu_zotus,esv_table_combined),
                                  taxon=setup_OTU_table(isolated_species_zotus,esv_table_combined))
  
  #Get cluster tables
  otu_tables[[amp]] <- setup_cluster_tables(otu_zotus,esv_table_combined,hap=min_hap)
  isolated_species_tables[[amp]] <- setup_cluster_tables(isolated_species_zotus,esv_table_combined,hap=min_hap)
}

no_wetland_community_tables <- list()
wetland_community_tables<- list()
for (amp in names(community_tables)){
  for (data_type in names(community_tables[[amp]])){
    temp <- t(community_tables[[amp]][[data_type]])
    remove_unique <- temp[rowSums(temp[])>0,]
    keep_unique <- temp[rowSums(temp[])==0,]
    no_wetland_community_tables[[amp]][[data_type]] <- t(remove_unique)
    wetland_community_tables[[amp]][[data_type]] <- t(keep_unique)
  }
}
# Setup community dissimilarity matrices
community_dm <- list()
for (amp in names(community_tables)){
  for (data_type in names(community_tables[[amp]])){
    current_table.dist<-create_beta_div_matrices(community_tables[[amp]][[data_type]],method="pa")
    community_dm[[amp]][[data_type]] <- current_table.dist
  }
}

# Setup species dissimilarity matrices
for (amp in amps){#Loop through all amplicons
  for (sp in names(isolated_species_tables[[amp]])){#Loop through each isolated species
    current_table.dist<-create_beta_div_matrices(isolated_species_tables[[amp]][[sp]],method=dm_method)
    isolated_species_dm[[amp]][["all"]][[sp]] <- current_table.dist[["all"]]
    isolated_species_dm[[amp]][["nested"]][[sp]] <- current_table.dist[["nested"]]
    isolated_species_dm[[amp]][["turnover"]][[sp]] <- current_table.dist[["turnover"]]
  }
  for (sp in names(otu_tables[[amp]])){#Loop through each OTU
    current_table.dist<-create_beta_div_matrices(otu_tables[[amp]][[sp]],method=dm_method)
    otu_dm[[amp]][["all"]][[sp]] <- current_table.dist[["all"]]
    otu_dm[[amp]][["nested"]][[sp]] <- current_table.dist[["nested"]]
    otu_dm[[amp]][["turnover"]][[sp]] <- current_table.dist[["turnover"]]
  }
}


#Setup population comparison

isolated_species_dm_selected <-list()
otu_dm_selected <-list()
target_dm <- list()

for (amp in amps){
  for (sp in names(isolated_species_dm[[amp]][[beta_component]])){
    sites <- labels(isolated_species_dm[[amp]][[beta_component]][[sp]])
    sites_in_dm <- cluster_ID[sites,]
    region_counts <- as.data.frame(table(sites_in_dm[["cluster_region"]]))
    region_present <- subset(region_counts, Freq >= min_sites)
    if (nrow(region_present)>1){
      isolated_species_dm_selected[[amp]][[beta_component]][[sp]] <-isolated_species_dm[[amp]][[beta_component]][[sp]]
    }
  }
}

for (amp in amps){
  for (sp in names(otu_dm[[amp]][[beta_component]])){
    sites <- labels(otu_dm[[amp]][[beta_component]][[sp]])
    sites_in_dm <- cluster_ID[sites,]
    region_counts <- as.data.frame(table(sites_in_dm[["cluster_region"]]))
    region_present <- subset(region_counts, Freq >= min_sites)
    if (nrow(region_present)>1){
      otu_dm_selected[[amp]][[beta_component]][[sp]] <-otu_dm[[amp]][[beta_component]][[sp]]
    }
  }
}

for (amp in amps){
  target_dm[[amp]][["is"]] <- average_distance_matrices(isolated_species_dm_selected[[amp]][[beta_component]])
  target_dm[[amp]][["otu"]] <- average_distance_matrices(otu_dm_selected[[amp]][[beta_component]])
  target_dm[[amp]][["com"]] <- as.data.frame(as.matrix(community_dm[[amp]][["esv"]][[beta_component]]))
}


saveRDS(important_tables,file=paste(indir,"RDS/important_tables_",variables_io,".RDS",sep=""))
saveRDS(community_dm,file=paste(indir,"RDS/community_dm_",variables_io,".RDS",sep=""))
saveRDS(community_tables,file=paste(indir,"RDS/community_tables_",variables_io,".RDS",sep=""))
saveRDS(isolated_species_tables,file=paste(indir,"RDS/sbc_tables_",variables_io,".RDS",sep=""))
saveRDS(otu_tables,file=paste(indir,"RDS/otu_tables_",variables_io,".RDS",sep=""))

saveRDS(isolated_species_dm_selected,file=paste(indir,"RDS/sbc_dm_selected_",variables_io,".RDS",sep=""))
saveRDS(otu_dm_selected,file=paste(indir,"RDS/otu_dm_selected_",variables_io,".RDS",sep=""))
saveRDS(target_dm,file=paste(indir,"RDS/target_dm_",variables_io,".RDS",sep=""))


