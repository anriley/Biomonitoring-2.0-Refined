library(usedist)
library(vegan)
library(ggplot2)
library(ggpubr)
library(ggtree)
library(tidytree)
library(grid)
library(gridExtra)
library(seqinr)
library(ape)
set.seed(017)

#########
#For tree plotting
nodeid.tbl_tree <- utils::getFromNamespace("nodeid.tbl_tree", "tidytree")
rootnode.tbl_tree <- utils::getFromNamespace("rootnode.tbl_tree", "tidytree")
offspring.tbl_tree <- utils::getFromNamespace("offspring.tbl_tree", "tidytree")
offspring.tbl_tree_item <- utils::getFromNamespace(".offspring.tbl_tree_item", "tidytree")
child.tbl_tree <- utils::getFromNamespace("child.tbl_tree", "tidytree")
parent.tbl_tree <- utils::getFromNamespace("parent.tbl_tree", "tidytree")  

#########
#Parameter, file, and function setup

source("Region_cluster_final/map_and_haplotype_network_functions_reduced.R")

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
isolated_species_tables <- readRDS(file=paste(indir,"RDS/sbc_tables_",variables_io,".RDS",sep=""))
taxonomy <- read.csv(paste(indir,"Data/Taxonomy/taxonomy_rcf_",read_filter,"_F230R.csv",sep=""))

#########
#Do isolated species comparison

species_order_list <- list()
for (x in 1:nrow(taxonomy)){
  species_order_list[[taxonomy[x,27]]] <- taxonomy[x,18]
}
#Create permanova and dispersion tables
compare <- data.frame()
permanova_sp_list <- list()
haplotype_comparison <- data.frame()
for (sp in intersect(names(isolated_species_dm_selected[["F230R"]][[beta_component]]),names(isolated_species_dm_selected[["MLJG"]][[beta_component]]))){
  haplotype_comparison <- as.data.frame(rbind(haplotype_comparison,
                                              c(sp,ncol(isolated_species_tables[["F230R"]][[sp]]),
                                                ncol(isolated_species_tables[["MLJG"]][[sp]]),
                                                species_order_list[[sp]])))
  current_F230R_dm <- isolated_species_dm_selected[["F230R"]][[beta_component]][[sp]]
  current_MLJG_dm <- isolated_species_dm_selected[["MLJG"]][[beta_component]][[sp]]
  shared_sites <- intersect(labels(current_F230R_dm),labels(current_MLJG_dm))
  #if (length(shared_sites) >6){
  sub_current_F230R_dm <- dist_subset(current_F230R_dm,shared_sites)
  sub_current_MLJG_dm <- dist_subset(current_MLJG_dm,shared_sites)
  mantel <- mantel(sub_current_MLJG_dm,sub_current_F230R_dm,method="spearman",permutations = 9999)
  F_sites <- length(labels(current_F230R_dm)) - length(shared_sites) 
  M_sites <- length(labels(current_MLJG_dm)) - length(shared_sites) 
  permanova_sp_list[[sp]][["F230R"]] <- adonis2(sub_current_F230R_dm~as.factor(cluster_ID[labels(sub_current_F230R_dm),]$cluster_region),
                                                data=cluster_ID[labels(sub_current_F230R_dm),], permutations = 9999)
  permanova_sp_list[[sp]][["MLJG"]] <- adonis2(sub_current_MLJG_dm~as.factor(cluster_ID[labels(sub_current_MLJG_dm),]$cluster_region),
                                               data=cluster_ID[labels(sub_current_MLJG_dm),], permutations = 9999)
  
  site_sim <- (F_sites+M_sites)/(2*length(shared_sites)+F_sites+M_sites)
  
  compare <- rbind(compare,c(sp,species_order_list[[sp]],length(shared_sites),
                             site_sim,mantel[["statistic"]],mantel[["signif"]],
                             permanova_sp_list[[sp]][["F230R"]][["F"]][1],permanova_sp_list[[sp]][["F230R"]][["Pr(>F)"]][1],
                             permanova_sp_list[[sp]][["MLJG"]][["F"]][1],permanova_sp_list[[sp]][["MLJG"]][["Pr(>F)"]][1]))
  
  #}
}
colnames(compare) <- c("species","order","Sites","Site_sim","cor","P","F_F","F_P","M_F","M_P")
compare$p_adjust <- p.adjust(compare$P, method = "fdr")
compare$F_p_adjust <- p.adjust(compare$F_P, method = "fdr")
compare$M_p_adjust <- p.adjust(compare$M_P, method = "fdr")
compare$check <- ifelse(as.numeric(compare$F_p_adjust) < 0.05|as.numeric(compare$M_p_adjust) < 0.05, "One","Neither")
compare$check <- ifelse(as.numeric(compare$F_p_adjust) < 0.05&as.numeric(compare$M_p_adjust) < 0.05, "Both",compare$check)
compare$F_perm <- ifelse(as.numeric(compare$F_p_adjust) < 0.05, "Significant","Not Significant")
compare$M_perm <- ifelse(as.numeric(compare$M_p_adjust) < 0.05, "Significant","Not Significant")
compare$Correlation <- ifelse(as.numeric(compare$p_adjust) < 0.05, "Significant","Not Significant")
compare <- compare[order(as.numeric(compare$cor), decreasing = TRUE),]
compare$species <- factor(compare$species, levels=unique(compare$species))
compare$check <- factor(compare$check, levels=c('Both', 'One', 'Neither'))
compare$Correlation <- factor(compare$Correlation, levels=c("Significant","Not Significant"))

#########
#Build tree

species_match <- intersect(names(isolated_species_dm_selected[["F230R"]][[beta_component]]),
          names(isolated_species_dm_selected[["MLJG"]][[beta_component]]))

taxonomy <- read.csv(paste(indir,"Data/Taxonomy/taxonomy_rcf_",read_filter,"_F230R.csv",sep=""))
species_order_list <- list()
for (x in 1:nrow(taxonomy)){
  species_order_list[[taxonomy[x,27]]] <- taxonomy[x,18]
}

fasta_MLJG <- read.fasta(paste(indir,"Data/Sequences/nt_rcf_100_MLJG.fasta",sep=""))  
amp <- "MLJG"
sp_abundant_zotu <- list()
for (sp in species_match){
  zotu_sums <- sort(colSums(isolated_species_tables[[amp]][[sp]][,1:length(isolated_species_tables[[amp]][[sp]])]),decreasing=TRUE)
  abundant_zotu <- names(zotu_sums)[1]
  sp_abundant_zotu[[sp]] <- abundant_zotu
}

sp_seq <- list()
for (sp in species_match){
  sequence <- fasta_MLJG[[sp_abundant_zotu[[sp]]]]
  order_sp <- gsub("_"," ", sp) 
  #order_sp <- paste(species_order_list[[sp]],sp,sep=" ")
  sp_seq[[order_sp]] <- sequence
}

sp_align <- ape::as.alignment(sp_seq)
distance_alignment <- dist.alignment(sp_align, matrix = "identity")
Tree <- bionj(distance_alignment)

#########
#Plotting

trda2 <- root(Tree, node = 39, edgelabel = TRUE)
d <- data.frame(node=c(49, 56,71,45,40,44), Order=c("Diptera", "Trichoptera",
                                                   "Ephemeroptera","Ephemeroptera",
                                                   "Plecoptera","Plecoptera"))
tree_plot <- ggtree(trda2, branch.length='none')+ geom_tiplab(hjust =0.5, offset=3)+ xlim(0, 15) +
  #geom_nodelab(aes(label = node))+
  geom_hilight(data=d, aes(node=node, fill=Order))+ 
  theme(axis.text.x=element_blank(),axis.title.x=element_blank(),
        axis.ticks.x=element_blank(),legend.position = "bottom")

tree_plot

compare2 <- compare
compare2$species <- gsub("_"," ", compare2$species) 
compare2$species <- factor(compare2$species,levels=rev(ggtree::get_taxa_name(tree_plot)))
compare_tree_order <- compare2[order(compare2$species),]

cor_compare <- data.frame(compare_tree_order$species,compare_tree_order$Correlation,compare_tree_order$cor)
colnames(cor_compare) <- c('species','Significance', 'Statistic')
cor_compare$Type <- "Correlation"

F_compare <- data.frame(compare_tree_order$species,compare_tree_order$F_perm,as.numeric(compare_tree_order$F_F))
colnames(F_compare) <- c('species','Significance', 'Statistic')
F_compare$Type <- "F230R"

M_compare <- data.frame(compare_tree_order$species,compare_tree_order$M_perm,as.numeric(compare_tree_order$M_F))
colnames(M_compare) <- c('species','Significance', 'Statistic')
M_compare$Type <- "MLJG"

region_sep <- rbind(F_compare,M_compare)
region_sep$Statistic<- ifelse(region_sep$Statistic > 1000, ">1000", round(as.numeric(region_sep$Statistic),2))

cor_plot <- ggplot(cor_compare, aes(x ="ρ", y=species, fill = Significance)) +geom_tile()+ 
  geom_text(aes(label=round(as.numeric(Statistic),2))) + 
  scale_fill_manual("Significance", values = c("Significant" = "#1E88E5", "Not Significant" = "#D81B60")) + 
  theme_classic() +theme(legend.position = "bottom",axis.title.x=element_blank(),axis.title.y=element_blank())

region_sep_plot <- ggplot(region_sep, aes(x =Type, y=species, fill = Significance)) +geom_tile()+ 
  geom_text(aes(label=Statistic)) + 
  scale_fill_manual("Significance", values = c("Significant" = "#1E88E5", "Not Significant" = "#D81B60")) + 
  theme_classic() +theme(legend.position = "bottom",axis.title.x=element_blank(),axis.title.y=element_blank())

stat_combine <- ggarrange(cor_plot+ theme(axis.text.y=element_blank()),
                          region_sep_plot+ theme(axis.text.y=element_blank(),axis.line.y=element_blank(),axis.ticks.y=element_blank()),
                          common.legend = TRUE, legend = "bottom", widths = c(0.5, 1))

#https://stackoverflow.com/questions/52281227/adding-a-ggtree-object-to-already-existing-ggplot-with-shared-y-axis



final_tree <- grid.arrange(arrangeGrob(tree_plot + theme(plot.margin=margin(0,-20,0,0)),                         
                                       nullGrob(), 
                                       heights=c(0.98,0.04)),
                           arrangeGrob(stat_combine, 
                                       nullGrob(), 
                                       heights=c(0.98,0.04)), 
                           ncol=2, widths=c(2,1))

write.csv(compare,file=paste(indir,"Tables/TableS10_sbc_table_raw.csv",sep=""))
ggsave(paste(indir,"Figures/Figure6__sbc_tree_",variables_io,".png",sep=""),
       final_tree,height=18,width=26,units="cm", bg = 'white')

