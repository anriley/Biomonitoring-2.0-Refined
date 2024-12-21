library(flexclust)
library(purrr)
#library(klaR)
#library(tidyverse)

indir <- "../"

site_filename <- "../Data/metadata_with_clusters_elevation.csv"
initial_site_data <- read.csv(site_filename) #Change to directory location
site_data <- subset(initial_site_data,Sample_Type!="WETLAND")

lat_long <- cbind(site_data$Lat,site_data$Long)
{
  clust <- c(1:20)
  randMeans <- c()
  randSD <- c()
  for (c in clust){
    rand_indices <- c()
    for (x in 1:1000){
      cluster_1 <- kmeans(lat_long, c, nstart = 10 )
      cluster_2 <- kmeans(lat_long, c, nstart = 10 )
      ct.km <- table(cluster_1$cluster, cluster_2$cluster)
      score <- randIndex(ct.km)
      rand_indices <- append(rand_indices, score[[1]])
    }
    
    randMeans <- append(randMeans,mean(rand_indices))
    randSD <- append(randSD,sd(rand_indices))
  }
  
} #Adjusted Rand Index

png(file=paste(indir,"Figures/FigureS5.png",sep=""),
    width=600, height=350)
plot(clust,randMeans,ylim =c(0,1),ylab="Adjusted Rand Index",xlab="Number of clusters")
arrows(x0=clust, y0=randMeans-randSD, x1=clust, y1=randMeans+randSD, code=3, angle=90, length=0.1)
dev.off()

#Based on
#https://uc-r.github.io/kmeans_clustering
{
  wss <- function(k) {
    kmeans(lat_long, k, nstart = 10 )$tot.withinss
  }
  
  # Compute and plot wss for k = 1 to k = 20
  k.values <- 1:20
  
  # extract wss for 2-20 clusters
  wss_values <- map_dbl(k.values, wss)
  
} #within-clusters sum of squares

png(file=paste(indir,"Figures/FigureS6.png",sep=""),
    width=600, height=350)
plot(k.values, wss_values,
     type="b", pch = 19, frame = FALSE, 
     xlab="Number of clusters",
     ylab="Total within-clusters sum of squares")
dev.off()

