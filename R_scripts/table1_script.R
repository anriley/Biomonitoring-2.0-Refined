library(usedist)
library(vegan)
library(ggplot2)
library(ggpubr)
library(ggtree)
library(tidytree)

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

important_tables_in <- readRDS(file=paste(indir,"RDS/important_tables_",variables_io,".RDS",sep=""))
cluster_ID <- important_tables_in[["metadata"]]

isolated_species_dm_selected <- readRDS(file=paste(indir,"RDS/sbc_dm_selected_",variables_io,".RDS",sep=""))
otu_dm_selected <- readRDS(file=paste(indir,"RDS/otu_dm_selected_",variables_io,".RDS",sep=""))
#########
#Do isolated species comparison

#Create permanova and dispersion tables
species_F230R <- data.frame()
permanova_sp_list <- list()
for (sp in names(isolated_species_dm_selected[["F230R"]][[beta_component]])){
  current_F230R_dm <- isolated_species_dm_selected[["F230R"]][[beta_component]][[sp]]
  permanova_sp_list[[sp]][["F230R"]] <- adonis2(current_F230R_dm~as.factor(cluster_ID[labels(current_F230R_dm),]$cluster_region),
                                                data=cluster_ID[labels(current_F230R_dm),], permutations = 9999)
  
  
  species_F230R <- rbind(species_F230R,c(sp,length(labels(current_F230R_dm)),
                                         permanova_sp_list[[sp]][["F230R"]][["F"]][1],
                                         permanova_sp_list[[sp]][["F230R"]][["Pr(>F)"]][1]))
  
  #}
}

colnames(species_F230R) <- c("species","Sites","F_F","F_P")
species_F230R$p_adjust <- p.adjust(as.numeric(species_F230R$F_P), method = "fdr")

check_species_F230R <- subset(species_F230R, p_adjust < 0.05)


#Create permanova and dispersion tables
species_MLJG <- data.frame()
permanova_sp_list <- list()
for (sp in names(isolated_species_dm_selected[["MLJG"]][[beta_component]])){
  current_MLJG_dm <- isolated_species_dm_selected[["MLJG"]][[beta_component]][[sp]]
  permanova_sp_list[[sp]][["MLJG"]] <- adonis2(current_MLJG_dm~as.factor(cluster_ID[labels(current_MLJG_dm),]$cluster_region),
                                                data=cluster_ID[labels(current_MLJG_dm),], permutations = 9999)
  
  
  species_MLJG <- rbind(species_MLJG,c(sp,length(labels(current_MLJG_dm)),
                                         permanova_sp_list[[sp]][["MLJG"]][["F"]][1],permanova_sp_list[[sp]][["MLJG"]][["Pr(>F)"]][1]))
  
  #}
}

colnames(species_MLJG) <- c("species","Sites","F_F","F_P")
species_MLJG$p_adjust <- p.adjust(as.numeric(species_MLJG$F_P), method = "fdr")

check_species_MLJG <- subset(species_MLJG, p_adjust < 0.05)



#Create permanova and dispersion tables
otu_F230R <- data.frame()
permanova_sp_list <- list()
for (sp in names(otu_dm_selected[["F230R"]][[beta_component]])){
  current_F230R_dm <- otu_dm_selected[["F230R"]][[beta_component]][[sp]]
  permanova_sp_list[[sp]][["F230R"]] <- adonis2(current_F230R_dm~as.factor(cluster_ID[labels(current_F230R_dm),]$cluster_region),
                                               data=cluster_ID[labels(current_F230R_dm),], permutations = 9999)
  
  
  otu_F230R <- rbind(otu_F230R,c(sp,length(labels(current_F230R_dm)),
                                         permanova_sp_list[[sp]][["F230R"]][["F"]][1],permanova_sp_list[[sp]][["F230R"]][["Pr(>F)"]][1]))
  
  #}
}

colnames(otu_F230R) <- c("species","Sites","F_F","F_P")
otu_F230R$p_adjust <- p.adjust(as.numeric(otu_F230R$F_P), method = "fdr")

check_otu_F230R <- subset(otu_F230R, p_adjust < 0.05)


#Create permanova and dispersion tables
otu_MLJG <- data.frame()
permanova_sp_list <- list()
haplotype_comparison <- data.frame()
for (sp in names(otu_dm_selected[["MLJG"]][[beta_component]])){
  current_MLJG_dm <- otu_dm_selected[["MLJG"]][[beta_component]][[sp]]
  permanova_sp_list[[sp]][["MLJG"]] <- adonis2(current_MLJG_dm~as.factor(cluster_ID[labels(current_MLJG_dm),]$cluster_region),
                                               data=cluster_ID[labels(current_MLJG_dm),], permutations = 9999)
  
  
  otu_MLJG <- rbind(otu_MLJG,c(sp,length(labels(current_MLJG_dm)),
                                         permanova_sp_list[[sp]][["MLJG"]][["F"]][1],permanova_sp_list[[sp]][["MLJG"]][["Pr(>F)"]][1]))
  
  #}
}

colnames(otu_MLJG) <- c("species","Sites","F_F","F_P")
otu_MLJG$p_adjust <- p.adjust(as.numeric(otu_MLJG$F_P), method = "fdr")

check_otu_MLJG <- subset(otu_MLJG, p_adjust < 0.05)

table_out <- data.frame(OTUs = c(nrow(check_otu_F230R),nrow(check_otu_MLJG)),
                        SBCs = c(nrow(check_species_F230R),nrow(check_species_MLJG)))
row.names(table_out) <- c("F230R","MLJG")
write.csv(table_out,paste(indir,"Tables/Table1_",variables_io,".csv",sep=""))
