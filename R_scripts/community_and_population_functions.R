library(vegan)
library(stats)
library(ggplot2)
library(ggpubr)
library(betapart)
library(reshape2)
library(uwot)
library(corrr)
library(ggpattern)
library(flexclust)
library(usedist)

#Get OTU table

#Get cluster list from SWARM file
otu_clusters <-function(clusters,cluster_column=1){
  otu_zotus <- list()
  count <- 1
  all_zotus <-strsplit(clusters[[cluster_column]], " ")
  for (zotus in all_zotus){
    otu_name <- paste("OTU_",count,sep="")
    otu_zotus[[otu_name]] <- zotus
    count <- count +1
  }
  return (otu_zotus)
}
#Get cluster list from taxonomy file
taxonomy_clusters <-function(taxonomy_file,taxonomic_level,BP_threshold){
  level_index <- which(colnames(taxonomy_file)==taxonomic_level)
  bp_index <- level_index +2
  taxonomy_assigned <- taxonomy_file[taxonomy_file[,bp_index] >= BP_threshold,]
  taxa <- unique(taxonomy_assigned[[taxonomic_level]])
  taxa_zotus <- list()
  for (tx in taxa){
    taxa_zotus[[tx]] <- taxonomy_assigned[taxonomy_assigned[,level_index] == tx,][["GlobalESV"]]
  }
  return (taxa_zotus)
}
#Get cluster list from isolated species file
isolate_species_clusters <-function(clusters,cluster_column="cluster",name_column="species"){
  is_zotus <- list()
  for (x in 1:nrow(clusters)){
    is_zotus[[clusters[x,][[name_column]]]] <- strsplit(clusters[x,][[cluster_column]], " ")[[1]]
  }
  return (is_zotus)
}

#Build table from ESV table and cluster list
setup_OTU_table <- function(cluster_list,ESV_table){
  OTU_table <- rownames(ESV_table)
  for (cluster in names(cluster_list)){
    current_zotus <- intersect(cluster_list[[cluster]],colnames(ESV_table))
    if (length(current_zotus) > 1){
      OTU_table <- cbind(OTU_table,TaxonX=rowSums(ESV_table[,current_zotus]))
    }else if (length(current_zotus) == 1){
      OTU_table <- cbind(OTU_table,TaxonX=ESV_table[,current_zotus])
    }
  }
  OTU_table <- OTU_table[,-1]
  OTU_table2 <- data.frame(apply(OTU_table, 2, function(x) as.numeric(as.character(x))))
  rownames(OTU_table2) <-rownames(ESV_table)
  colnames(OTU_table2) <-names(cluster_list)
  return (OTU_table2)
}

#Make presence absence table
setup_pa_table <- function(in_table){
  pa_table <- in_table
  pa_table[pa_table > 0] <- 1
  return (pa_table)
}

#Make relative abundance table
setup_relative_table <- function(in_table){
  relative_table <- in_table/rowSums(in_table)
  return (as.data.frame(relative_table))
}

create_beta_div_matrices <- function(in_table,method=c("absolute","relative","pa")){
  beta_dm<-list()
  if (method=="relative"){
    current_table <- setup_relative_table(in_table)
  }else if (method=="pa"){
    current_table <- setup_pa_table(in_table)
  }else{
    current_table <-in_table 
  }
  beta_dm_temp <- beta.pair.abund(current_table)
  beta_dm[["all"]] <- beta_dm_temp[["beta.bray"]]
  beta_dm[["nested"]] <- beta_dm_temp[["beta.bray.gra"]]
  beta_dm[["turnover"]] <- beta_dm_temp[["beta.bray.bal"]]
  
  return(beta_dm)
}



pull_esv_table <-function(zotus,ESV_table){
  temp <- ESV_table[,zotus]
  zotu_sites <- temp[rowSums(temp[])>0,, drop = FALSE]
  return(zotu_sites)
}
setup_cluster_tables <- function(clusters,ESV_table,hap=0){
  cluster_table_list <- list()
  for (x in names(clusters)){
    if (length(clusters[[x]])>=hap){ #Only look at species with a certain number of haplotypes
      x_ESV_table <- pull_esv_table(clusters[[x]],ESV_table)
      if (nrow(x_ESV_table)>1){ #Only collect species present in more than one site
        cluster_table_list[[x]] <-x_ESV_table 
      }
    }
  }
  return(cluster_table_list)
}

combine_samples <- function(sample_table,sample_metadata,by_variable){
  sample_table_metadata <- merge(as.data.frame(sample_table), sample_metadata[by_variable], 
                                 by = 'row.names')
  sample_table_combined <- aggregate(sample_table_metadata[,sapply(sample_table_metadata,is.numeric)],sample_table_metadata[by_variable],sum)
  rownames(sample_table_combined) <- sample_table_combined[[by_variable]]
  sample_table_combined <- sample_table_combined[,-1]
  return(sample_table_combined)
}



count_site_contribution <- function(dm_list){
  site_counts <- list()
  for (sp in names(dm_list)){
    sites <- labels(dm_list[[sp]])
    for (s in sites){
      if (s %in% names(site_counts)){
        site_counts[[s]] <- site_counts[[s]] + length(sites)-1
      }else{
        site_counts[[s]] <- length(sites)-1
      }
    }
  }
  df_out <- t(as.data.frame(site_counts))
  colnames(df_out) <- c("site_count")
  return(df_out)
}
average_distance_matrices <- function(distance_matrix_list){
  first <- TRUE
  for (sp in names(distance_matrix_list)){
    if (first){
      combine_dm_long <-melt(as.matrix(distance_matrix_list[[sp]]))
      first <- FALSE
    }else{
      combine_dm_long <- merge(combine_dm_long,melt(as.matrix(distance_matrix_list[[sp]])),by=c("Var1","Var2"), all = TRUE)
    }
  }
  combine_dm_long$mean_dist <- rowMeans(combine_dm_long[,3:ncol(combine_dm_long)], na.rm=TRUE)
  
  average_dm_long <- as.data.frame(cbind(site1 = as.character(combine_dm_long[["Var1"]]),site2=as.character(combine_dm_long[["Var2"]]),dist =combine_dm_long$mean_dist))
  
  average_dm_wide_nas <- dcast(average_dm_long, site1~site2, value.var="dist")
  rownames(average_dm_wide_nas) <-average_dm_wide_nas$site1
  #Remove NAs
  average_dm_wide_nas <- average_dm_wide_nas[,-1]
  average_dm_wide <- average_dm_wide_nas
  na_rows = rowSums(is.na(average_dm_wide))
  while(sum(na_rows) != 0){
    na_rows = rowSums(is.na(average_dm_wide))
    new_names <- rownames(average_dm_wide) != names(which.max(na_rows))
    average_dm_wide <- average_dm_wide[new_names,new_names]
  }
  return(average_dm_wide)
}



rand_index_umap <- function(umap_ord,known_clusters){
  randMeans <- c()
  randSD <- c()
  rand_indices <- c()
  for (x in 1:1000){
    cluster_1 <- kmeans(umap_ord, 4, nstart = 10 )
    cluster_2 <- known_clusters
    ct.km <- table(cluster_1$cluster, cluster_2$cluster_region)
    score <- randIndex(ct.km, correct = TRUE)
    rand_indices <- append(rand_indices, score[[1]])
  }
  output <- list(rand_mean = mean(rand_indices), rand_sd = sd(rand_indices))
  return(output)
}

create_umap <- function(dm,neigh,dim_u){
  set.seed(017)
  umap_data <- umap(dm,n_neighbors = neigh, n_components = dim_u, approx_pow = TRUE)
  return(umap_data)
}

create_umap_plot <- function(umap_data,initial_data,rand_mean="",rand_sd="",pseudo_F=""){
  umap_data_merge <- merge(umap_data,initial_data[rownames(umap_data),, drop = FALSE],by="row.names")
  if (rand_mean == "" & pseudo_F == ""){
    umap_plot <- ggplot(data=umap_data_merge,aes(x=V1,y=V2,col=as.character(cluster_region),shape=as.character(cluster_region)))  +
      geom_point(size=2)+
      xlab("UMAP1") + ylab("UMAP2")+
      scale_color_manual("Region Clusters", values = c("Wetland" = "black", "Northeast" = "#D81B60", "Southeast" = "#1E88E5",
                                                       "Central" = "#FFC107","West" = "#004D40")) +
      scale_shape_manual("Region Clusters", values = c("Wetland" = 4, "Northeast" = 15, "Southeast" = 19,
                                                       "Central" = 17,"West" = 18)) +
      theme_bw()+ theme(plot.title = element_text(hjust = 0.5))
  }else if(rand_mean == ""){
    umap_plot <- ggplot(data=umap_data_merge,aes(x=V1,y=V2,col=as.character(cluster_region),shape=as.character(cluster_region)))  +
      geom_point(size=2)+
      xlab("UMAP1") + ylab("UMAP2")+
      scale_color_manual("Region Clusters", values = c("Wetland" = "black", "Northeast" = "#D81B60", "Southeast" = "#1E88E5",
                                                       "Central" = "#FFC107","West" = "#004D40")) +
      scale_shape_manual("Region Clusters", values = c("Wetland" = 4, "Northeast" = 15, "Southeast" = 19,
                                                       "Central" = 17,"West" = 18)) +
      theme_bw()+ theme(plot.title = element_text(hjust = 0.5)) +
      labs(caption=paste("pseudo F = ",round(pseudo_F,3),sep=""))
  }else{
    umap_plot <- ggplot(data=umap_data_merge,aes(x=V1,y=V2,col=as.character(cluster_region),shape=as.character(cluster_region)))  +
      geom_point(size=2)+
      xlab("UMAP1") + ylab("UMAP2")+
      scale_color_manual("Region Clusters", values = c("Wetland" = "black", "Northeast" = "#D81B60", "Southeast" = "#1E88E5",
                                                       "Central" = "#FFC107","West" = "#004D40")) +
      scale_shape_manual("Region Clusters", values = c("Wetland" = 4, "Northeast" = 15, "Southeast" = 19,
                                                       "Central" = 17,"West" = 18)) +
      theme_bw()+ theme(plot.title = element_text(hjust = 0.5)) +
      labs(caption=paste("Adjusted Rand Index = ",round(rand_mean,3)," ± ",round(rand_sd,3),sep=""))
  }
  return(umap_plot)
}


create_violin_plot <- function(dm,initial_data,samples){
  test_melt <- melt(as.matrix(dm))
  clus <- as.data.frame(cbind(initial_data[samples,][["ID"]],initial_data[samples,][["cluster_region"]]))
  colnames(clus) <- c("Var1","cluster_name1")
  test_2 <- merge(test_melt,clus,by="Var1")
  colnames(clus) <- c("Var2","cluster_name2")
  test_3 <- merge(test_2,clus,by="Var2")
  test_3$pair <- ifelse(test_3$cluster_name1==test_3$cluster_name2,test_3$cluster_name1,
                        paste(test_3$cluster_name1,test_3$cluster_name2,sep="_"))
  test_4 <- subset(test_3, as.character(Var1) != as.character(Var2))
  test_5 <- subset(test_4, pair != "West_Northeast" & pair != "West_Southeast" &
                     pair != "West_Central" & pair != "Central_Southeast" &
                     pair != "Central_Northeast" & pair != "Southeast_Northeast")
  test_5$pair <- ifelse(test_5$pair=="Central_West","C-W",test_5$pair)
  test_5$pair <- ifelse(test_5$pair=="Southeast_West","S-W",test_5$pair)
  test_5$pair <- ifelse(test_5$pair=="Southeast_Central","S-C",test_5$pair)
  test_5$pair <- ifelse(test_5$pair=="Northeast_Southeast","N-S",test_5$pair)
  test_5$pair <- ifelse(test_5$pair=="Northeast_Central","N-C",test_5$pair)
  test_5$pair <- ifelse(test_5$pair=="Northeast_West","N-W",test_5$pair)
  
  
  test_5$pair <- factor(test_5$pair, levels=c("West","Central","Southeast","Northeast",
                                              "C-W","S-W","S-C","N-S","N-C","N-W"))
  violin_plot <- ggplot(test_5,aes(x=pair,y=as.numeric(value)))+geom_violin() +
    #geom_dotplot(binaxis='y', stackdir='center', dotsize=0.02) 
    scale_x_discrete(guide = guide_axis(angle = 90)) +
    xlab("Regions") + ylab("Dissimilarities")+ ylim(0,1)+
    geom_violin_pattern(aes(fill=cluster_name1,pattern_fill=cluster_name2),pattern_colour=NA,
                        pattern_density =0.5,pattern_angle=0) +
    scale_pattern_fill_manual("Region Clusters", values = c("Wetland" = "black", "Northeast" = "#D81B60", "Southeast" = "#1E88E5",
                                                            "Central" = "#FFC107","West" = "#004D40")) +
    scale_fill_manual("Region Clusters", values = c("Wetland" = "black", "Northeast" = "#D81B60", "Southeast" = "#1E88E5",
                                                    "Central" = "#FFC107","West" = "#004D40")) +
    theme_bw() +geom_boxplot(width=0.1)#+geom_dotplot(binaxis='y', stackdir='center', dotsize=0.015)
  return(violin_plot)
}

#pairwise.adonis2 based on
#https://uw.pressbooks.pub/appliedmultivariatestatistics/chapter/complex-models/
pairwise.adonis2 <- function(resp, fact, p.method = "none", nperm = 9999) {
  require(vegan)
  resp <- as.matrix(resp)
  fact <- factor(fact)
  fun.R <- function(i, j) {
    fact2 <- droplevels(fact[as.numeric(fact) %in% c(i, j)])
    index <- which(fact %in% levels(fact2))
    resp2 <- as.dist(resp[index, index])
    result <- adonis2(resp2 ~ fact2, permutations = nperm)
    #result$`Pr(>F)`[1]
    result$`R2`[1]
  }
  multcomp_R <- pairwise.table(fun.R, levels(fact), p.adjust.method = "none")
  fun.R <- function(i, j) {
    fact2 <- droplevels(fact[as.numeric(fact) %in% c(i, j)])
    index <- which(fact %in% levels(fact2))
    resp2 <- as.dist(resp[index, index])
    result <- adonis2(resp2 ~ fact2, permutations = nperm)
    #result$`Pr(>F)`[1]
    result$`F`[1]
  }
  multcomp_F <- pairwise.table(fun.R, levels(fact), p.adjust.method = "none")
  fun.p <- function(i, j) {
    fact2 <- droplevels(fact[as.numeric(fact) %in% c(i, j)])
    index <- which(fact %in% levels(fact2))
    resp2 <- as.dist(resp[index, index])
    result <- adonis2(resp2 ~ fact2, permutations = nperm)
    result$`Pr(>F)`[1]
    #result$`R2`[1]
  }
  multcomp_p <- pairwise.table(fun.p, levels(fact), p.adjust.method = "none")
  multcomp_p_adjust <- pairwise.table(fun.p, levels(fact), p.adjust.method = p.method)
  return(list(fact = levels(fact), F = multcomp_F, R2 = multcomp_R,
              p.value = multcomp_p, p.adjust = multcomp_p_adjust,
              p.adjust.method = p.method))
}
