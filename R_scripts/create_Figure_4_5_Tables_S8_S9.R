library(ggplot2)
set.seed(017)

indir <- "../"

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

random_data <- read.csv(paste(indir,"Tables/randomized_data_assign_0.8_filter_100_0_5_method_pa_all.csv",sep=""),row.names=1)
colnames(random_data) <- c("Amplicon","Type","Sample","Region","Statistic","p_value")

#Randomized F statistics
F_pop_data_frame <-data.frame(AmpliconType = c("F230R_otu","F230R_is","MLJG_otu","MLJG_is"),
                              Statistic =c(24.6174712587324,25.8057642118536,19.6770723318228,13.3244277135726))
F_statistics <- subset(random_data, Region == "F")
F_statistics$AmpliconType <- paste(F_statistics$Amplicon,F_statistics$Type,sep="_")
F_labels <- c('F230R \nscambled OTUs', 'F230R \nscambled SBCs', 'MLJG \nscambled OTUs', 'MLJG \nscambled SBCs')

F_violin_plot <- ggplot(F_statistics,aes(x=factor(AmpliconType, level =c("F230R_otu","F230R_is","MLJG_otu","MLJG_is")),
                                                  y=as.numeric(Statistic)))+geom_violin(fill="palegreen") + theme_bw() +
  xlab("") + ylab("pseudo-F")+ ylim(12,35) +
  geom_point(F_pop_data_frame,mapping=aes(x=AmpliconType,y=as.numeric(Statistic)),pch=21,fill="white",size=3,inherit.aes = FALSE) +
  scale_x_discrete(labels= F_labels)

ggsave(paste(indir,"Figures/Figure4_scrambled_F_violin_plots.png",sep=""),
       F_violin_plot,height=12,width=16,units="cm", bg = 'white')

#Randomized Spearmen correlation
library(geosphere)
library(vegan)
library(usedist)
library(ggpubr)

important_tables_in <- readRDS(file=paste(indir,"RDS/important_tables_",variables_io,".RDS",sep=""))
target_dm_in <- readRDS(file=paste(indir,"RDS/target_dm_",variables_io,".RDS",sep=""))

cluster_ID <- important_tables_in[["metadata"]]
lat_long <- as.data.frame(cbind(as.numeric(cluster_ID$longitude),as.numeric(cluster_ID$latitude)))
colnames(lat_long) <- c("longitude","latitude")
rownames(lat_long) <- cluster_ID$ID
lat_long_dm <- distm(lat_long, fun = distGeo)
rownames(lat_long_dm) <- rownames(lat_long)
colnames(lat_long_dm) <- rownames(lat_long)
metadata <- important_tables_in[["metadata"]]
region_clusters <- c("Northeast","Southeast","Central","West")
cor_violin_plots <- list()
all_region_summaries <- data.frame()
all_region <- data.frame()
for (amp in names(target_dm_in)){
  for (type in c("is","otu")){
    current_dm <- target_dm_in[[amp]][[type]]
    region_pop_cor <- data.frame()
    for(region in region_clusters){
      current_region <- subset(metadata, cluster_region == region)
      sites <- intersect(current_region$ID,rownames(current_dm))
      bio_dm <- dist_subset(current_dm, sites)
      geo_dm <- dist_subset(lat_long_dm, sites)
      pop_cor <- mantel(geo_dm,bio_dm,method="spearman",permutations=9999)
      region_pop_cor <- as.data.frame(rbind(region_pop_cor,c(amp,type,region,pop_cor[["statistic"]],pop_cor[["signif"]])))
      current_random_region <- subset(random_data,Amplicon==amp &Type==type & Region == region)
      current_random_region_p <- subset(current_random_region, p_value < 0.05)
      outline <- c(amp,type,region,mean(as.numeric(current_random_region$Statistic)),
                   sd(as.numeric(current_random_region$Statistic)),
                   quantile(as.numeric(current_random_region$Statistic)),
                   nrow(current_random_region_p))
      all_region_summaries <-as.data.frame(rbind(all_region_summaries,outline))
    }
    colnames(region_pop_cor) <- c("Amplicon","Type","Region","Statistic")
    all_region <- as.data.frame(rbind(all_region,region_pop_cor))
    current_random <- subset(random_data,Amplicon==amp &Type==type & Region != "F" & Region != "All")
    cor_violin_plots[[amp]][[type]] <- ggplot(current_random,aes(x=Region,y=as.numeric(Statistic),fill=Region))+geom_violin() + theme_bw() +
      xlab("Regions") + ylab("Spearman's ρ")+ ylim(-0.1,1) +
      geom_point(region_pop_cor,mapping=aes(x=Region,y=as.numeric(Statistic)),pch=21,fill="white",size=3,inherit.aes = FALSE)+
      scale_fill_manual("Region Clusters", values = c("Northeast" = "#D81B60", "Southeast" = "#1E88E5",
                                                      "Central" = "#FFC107","West" = "#004D40"),guide="none")
  }
}

colnames(all_region_summaries) <- c("Amplicon","Data","Region","Mean correlation","Standard deviation",
                                    "Minimum","Q1","Median","Q3","Maximum","Significant correlations")

all_corr_violin_plots <- ggarrange(cor_violin_plots[["F230R"]][["otu"]]+ggtitle("Scrambled OTUs")+
                                     theme(plot.title = element_text(size = 16,hjust = 0.5),legend.title = element_text(size=12),legend.text = element_text(size=12)),
                                   cor_violin_plots[["F230R"]][["is"]]+ggtitle("Scrambled SBCs")+theme(plot.title = element_text(size = 16,hjust = 0.5)),
          cor_violin_plots[["MLJG"]][["otu"]]+ggtitle(""),cor_violin_plots[["MLJG"]][["is"]]+ggtitle(""),
          nrow=2,ncol=2,labels=c("A","","B"))
ggsave(paste(indir,"Figures/Figure5_scrambled_corr_violin_plots.png",sep=""),
       all_corr_violin_plots,height=16,width=16,units="cm", bg = 'white')

write.csv(all_region,paste(indir,"Tables/TableS8_randomized_data_",variables_io,".csv",sep=""))

write.csv(all_region_summaries,paste(indir,"Tables/TableS9_randomized_data_",variables_io,".csv",sep=""))
