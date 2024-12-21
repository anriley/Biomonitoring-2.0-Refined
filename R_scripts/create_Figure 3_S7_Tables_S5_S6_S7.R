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

target_dm_in <- readRDS(file=paste(indir,"RDS/target_dm_",variables_io,".RDS",sep=""))



#FIGURE PLOTTING POPULATION
neighbors <- 15
dim <- 2
umap_data <- list()
umap_plots <- list()
violin_plots <- list()
permanova_list <- list()
permanova_pair_list <- list()
dispersion_list <- list()
dispersion_pair_list <- list()
dispersion_data_list <- list()
for (amp in names(target_dm_in)){
  shared_samples <- intersect(rownames(target_dm_in[[amp]][["is"]]),rownames(target_dm_in[[amp]][["otu"]]))
  for (data_type in names(target_dm_in[[amp]])){
    current_dm <- target_dm_in[[amp]][[data_type]][shared_samples,shared_samples]
    current_dm_dist <- as.dist(as.matrix(current_dm))
    
    umap_data[[amp]][[data_type]][["default"]] <- create_umap(current_dm_dist,neighbors,dim)
    umap_permanova <- adonis2(umap_data[[amp]][[data_type]][["default"]]~as.factor(cluster_ID[shared_samples,]$cluster_region),
                              data=cluster_ID[shared_samples,], permutations = 9999, method = "euclidean")
    umap_plots[[amp]][[data_type]][["default"]] <- create_umap_plot(umap_data[[amp]][[data_type]][["default"]],
                                                                    cluster_ID,pseudo_F = umap_permanova[["F"]][1])
    
    umap_data[[amp]][[data_type]][["max"]] <- create_umap(current_dm_dist,length(shared_samples)-1,dim)
    umap_permanova <- adonis2(umap_data[[amp]][[data_type]][["max"]]~as.factor(cluster_ID[shared_samples,]$cluster_region),
                              data=cluster_ID[shared_samples,], permutations = 9999, method = "euclidean")
    umap_plots[[amp]][[data_type]][["max"]] <- create_umap_plot(umap_data[[amp]][[data_type]][["max"]],
                                                                cluster_ID,pseudo_F = umap_permanova[["F"]][1])
    
    violin_plots[[amp]][[data_type]]<- create_violin_plot(current_dm_dist,cluster_ID,shared_samples)
    
    permanova_list[[amp]][[data_type]] <- adonis2(current_dm_dist~as.factor(cluster_ID[shared_samples,]$cluster_region),data=cluster_ID[shared_samples,], permutations = 9999)
    permanova_pair_list[[amp]][[data_type]] <- pairwise.adonis2(current_dm_dist, cluster_ID[shared_samples,]$cluster_region, p.method = "fdr")
    dispersion <- betadisper(current_dm_dist,cluster_ID[shared_samples,]$cluster_region)
    dispersion_list[[amp]][[data_type]] <- permutest(dispersion,permutations = 9999)
    dispersion_pair_list[[amp]][[data_type]] <- TukeyHSD(dispersion)
    dispersion_data_list[[amp]][[data_type]] <-dispersion
  }
}

umap_plots_F230R_combined <- ggarrange(umap_plots[["F230R"]][["com"]][["default"]]+ggtitle("Community ESVs")+theme(plot.title = element_text(size = 20,hjust = 0.5),legend.title = element_text(size=16),legend.text = element_text(size=16)),
                                       umap_plots[["F230R"]][["otu"]][["default"]]+ggtitle("Intraspecific OTUs")+theme(plot.title = element_text(size = 20,hjust = 0.5)),
                                       umap_plots[["F230R"]][["is"]][["default"]]+ggtitle("Intraspecific SBCs")+theme(plot.title = element_text(size = 20,hjust = 0.5)),
                                       
                                       umap_plots[["F230R"]][["com"]][["max"]]+ggtitle(""),
                                       umap_plots[["F230R"]][["otu"]][["max"]]+ggtitle(""),
                                       umap_plots[["F230R"]][["is"]][["max"]]+ggtitle(""),
                                       
                                       violin_plots[["F230R"]][["com"]]+ggtitle(""),
                                       violin_plots[["F230R"]][["otu"]]+ggtitle(""),
                                       violin_plots[["F230R"]][["is"]]+ggtitle(""),
                                       
                                       ncol=3,nrow=3,common.legend = TRUE,legend = "right",
                                       labels=c("A","","","B","","","C",""))

umap_plots_MLJG_combined <- ggarrange(umap_plots[["MLJG"]][["com"]][["default"]]+ggtitle("Community ESVs")+theme(plot.title = element_text(size = 20,hjust = 0.5),legend.title = element_text(size=16),legend.text = element_text(size=16)),
                                      umap_plots[["MLJG"]][["otu"]][["default"]]+ggtitle("Intraspecific OTUs")+theme(plot.title = element_text(size = 20,hjust = 0.5)),
                                      umap_plots[["MLJG"]][["is"]][["default"]]+ggtitle("Intraspecific SBCs")+theme(plot.title = element_text(size = 20,hjust = 0.5)),
                                      
                                      umap_plots[["MLJG"]][["com"]][["max"]]+ggtitle(""),
                                      umap_plots[["MLJG"]][["otu"]][["max"]]+ggtitle(""),
                                      umap_plots[["MLJG"]][["is"]][["max"]]+ggtitle(""),
                                      
                                      violin_plots[["MLJG"]][["com"]]+ggtitle(""),
                                      violin_plots[["MLJG"]][["otu"]]+ggtitle(""),
                                      violin_plots[["MLJG"]][["is"]]+ggtitle(""),
                                      
                                      ncol=3,nrow=3,common.legend = TRUE,legend = "right",
                                      labels=c("A","","","B","","","C",""))


ggsave(paste(indir,"Figures/Figure3_umap_F230R_",variables_io,".png",sep=""),
       umap_plots_F230R_combined,height=24,width=32,units="cm", bg = 'white')
ggsave(paste(indir,"Figures/FigureS7_umap_MLJG_",variables_io,".png",sep=""),
       umap_plots_MLJG_combined,height=24,width=32,units="cm", bg = 'white')


#Populations tables
permanova_table_pop <- data.frame()
for (amp in names(permanova_list)){
  for (data_type in names(permanova_list[[amp]])){
    permanova_table_pop <- rbind(permanova_table_pop, 
                                 as.data.frame(permanova_list[[amp]][[data_type]]))
  }
}

dispersion_table_pop <- data.frame()
for (amp in names(dispersion_list)){
  for (data_type in names(dispersion_list[[amp]])){
    dispersion_table_pop <- rbind(dispersion_table_pop, 
                                  as.data.frame(dispersion_list[[amp]][[data_type]][["tab"]]))
  }
}

region_pairs <- c("Northeast","Central","Southeast","Central","West","Central",
                  "Southeast","Northeast","West","Northeast","West","Southeast")
permanova_pair_table_pop <- c()
for (amp in names(permanova_pair_list)){
  for (data_type in names(permanova_pair_list[[amp]])){
    for (x in seq(1,length(region_pairs),2)){
      current_data <-permanova_pair_list[[amp]][[data_type]]
      r <- region_pairs[x]
      c <- region_pairs[x+1]
      permanova_pair_table_pop <- rbind(permanova_pair_table_pop,
                                        c(amp,data_type,paste(region_pairs[x+1],region_pairs[x],sep=" "),
                                          current_data[["F"]][r,c],current_data[["R2"]][r,c],
                                          current_data[["p.value"]][r,c],current_data[["p.adjust"]][r,c]))
    }
  }
}

permanova_pair_table_pop <- as.data.frame(permanova_pair_table_pop)
colnames(permanova_pair_table_pop) <- c("Amplicon","Type","Region pairs","F","R2","Pr(>F)","Pr(>F) adjusted")

write.csv(permanova_table_pop,paste(indir,"Tables/TableS5_population_permanova_",variables_io,".csv",sep=""))
write.csv(dispersion_table_pop,paste(indir,"Tables/TableS6_population_dispersion_",variables_io,".csv",sep=""))
write.csv(permanova_pair_table_pop,paste(indir,"Tables/TableS7_population_permanova_pairs_",variables_io,".csv",sep=""))
