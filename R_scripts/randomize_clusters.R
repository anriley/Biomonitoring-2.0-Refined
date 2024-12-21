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


randomizations <- 1000
random_clusters <- list()
for (amp in amps){
  #Randomize OTUs
  swarm_clusters <- read.delim(paste(indir,"Data/Clusters/Clusters_rcf_",read_filter,"_",amp,"/clusters_rcf_",read_filter,"_d_",d_amp[[amp]],"_",amp,".txt",sep=""),
                               header=FALSE)
  otu_zotus <- otu_clusters(swarm_clusters)
  len_otu <- c()
  all_otu_zotus <- c()
  for (otu in names(otu_zotus)){
    len_otu <- append(len_otu,length(otu_zotus[[otu]]))
    all_otu_zotus <- append(all_otu_zotus,otu_zotus[[otu]])
  }
  for (x in 1:randomizations){
    random_zotus <- sample(all_otu_zotus)
    random_otu_list <-list()
    i <- 1
    for (y in 1:length(len_otu)){
      otu_name <- paste("OTU_",y,sep="")
      j <- i + len_otu[y] -1
      random_otu_list[[otu_name]] <- random_zotus[i:j]
      i <- j +1
    }
    random_clusters[[amp]][["otu"]][[as.character(x)]] <- random_otu_list
  }
  #Randomize isolated species
  species_clusters <- read.csv(paste(indir,"Data/Taxonomy_clusters/taxonomy_clusters_rcf_",read_filter,"_Species_",sBp,"_",amp,"_isolated.csv",sep=""))
  species_zotus <- isolate_species_clusters(species_clusters)
  len_species <- c()
  all_species_zotus <- c()
  for (sp in names(species_zotus)){
    len_species <- append(len_species,length(species_zotus[[sp]]))
    all_species_zotus <- append(all_species_zotus,species_zotus[[sp]])
  }
  for (x in 1:randomizations){
    random_zotus <- sample(all_species_zotus)
    random_otu_list <-list()
    i <- 1
    for (y in 1:length(len_species)){
      otu_name <- paste("OTU_",y,sep="")
      j <- i + len_species[y] -1
      random_otu_list[[otu_name]] <- random_zotus[i:j]
      i <- j +1
    }
    random_clusters[[amp]][["is"]][[as.character(x)]] <- random_otu_list
  }
}

saveRDS(random_clusters,file=paste(indir,"RDS/randomized_clusters_",variables_io,".RDS",sep=""))
