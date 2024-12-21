library(ggplot2)
library(raster)
library(elevatr)
library(sf)
library(tidyterra)
library(ggnewscale)
library(rnaturalearth)
library(rnaturalearthdata)
library(cowplot)

get_elevation_map <- function(long_min,long_max,lat_min,lat_max){
  target_range <- data.frame(x=c(long_min,long_max),y=c(lat_min,lat_max)) %>% st_as_sf(coords = c("x","y"), crs = 4326) 
  target_map <- get_elev_raster(target_range, z=7) #,long_min=-118,long_max=-113,lat_min=49, lat_max=52 #z=7 is ideal? z is the zoom on the map, higher z gives better resolution up to 14. A link to resolutions at different z : https://github.com/tilezen/joerd/blob/master/docs/data-sources.md#what-is-the-ground-resolution
  
  
  target_map_df <- as.data.frame(target_map, xy=TRUE)
  colnames(target_map_df) <- c("Long","Lat","Elevation")
  #plot(target_map)
  grad <- hypso.colors(10, "colombia_hypso")
  plotMap <- ggplot(data = target_map_df)+
    geom_raster(mapping=aes(x=Long, y=Lat, fill=target_map_df[,3])) +
    scale_fill_gradientn(colours= grad, name='Elevation (m)',na.value = NA) +
    theme_bw() +
    theme(axis.line = element_line(colour = "black"),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.border = element_blank(),
          panel.background = element_blank())
  return(plotMap)
}

#setup
indir <- "../"
site_filename <- paste(indir,"/Data/metadata_with_clusters_elevation.csv",sep="")
initial_site_data <- read.csv(site_filename) #Change to directory location
site_data2 <- as.data.frame(unique(cbind(ID=initial_site_data[["ID"]],Long=initial_site_data[["Long"]],Lat=initial_site_data[["Lat"]],
                                         Year=initial_site_data[["Year"]],clusters=initial_site_data[["cluster_name"]])))
rownames(site_data2) <- site_data2[["ID"]]
site_data2 <- site_data2[,-1]
site_data3 <- data.frame(apply(site_data2, 2, function(x) as.numeric(as.character(x))))
rownames(site_data3) <- rownames(site_data2)

long_min <- min(initial_site_data$Long)#http://127.0.0.1:19833/graphics/plot_zoom_png?width=1920&height=1009
long_max <- max(initial_site_data$Long)
lat_min <- min(initial_site_data$Lat)
lat_max <- max(initial_site_data$Lat)
site_map <- get_elevation_map(long_min,long_max,lat_min,lat_max)

site_map_with_sites <- site_map+new_scale("fill") +
  geom_point(aes(x=as.numeric(Long),y=as.numeric(Lat),fill=as.character(clusters),size =as.character(clusters),
                 shape = as.character(clusters),stroke =as.character(clusters)),
                                       site_data2,alpha=0.5) +
  #scale_fill_manual(values=c("red","yellow","blue","green","black"))+
  scale_fill_manual("Region clusters", values = c("Wetland" = "black", "Northeast" = "#D81B60", "Southeast" = "#1E88E5",
                                           "Central" = "#FFC107","West" = "#004D40"))+
  scale_shape_manual("Region clusters", values = c("Wetland" = 4, "Northeast" = 22, "Southeast" = 21,
                                           "Central" = 24,"West" = 23)) +
  scale_discrete_manual("Region clusters",aesthetics = "stroke",values = c("Wetland" = 2, "Northeast" = 1, "Southeast" = 1,
               "Central" = 1,"West" = 1))+
  scale_size_manual("Region clusters",values = c("Wetland" = 1, "Northeast" = 4, "Southeast" = 4,
                                                                           "Central" = 3,"West" = 4))

canada_map <- ne_states(country = "canada", returnclass = "sf")
bc_ab <- canada_map[canada_map$name == "Alberta"|canada_map$name == "British Columbia",]
map_location <- ggplot() +
  geom_sf(data = bc_ab, colour = "black", fill = "white")+
  geom_rect(aes(xmin=-118,xmax=-113,ymin=49,ymax=52),colour="red",alpha=0)+ theme_void()+
  theme(panel.background = element_rect(fill = 'lightblue', colour = 'black'))


final_site_map <- ggdraw(site_map_with_sites) + draw_plot(map_location, width = 0.28, height = 0.28, x = 0.533, y = 0.68)
ggsave(paste(indir,"Figures/Figure2_site_map.png",sep=""),
       final_site_map,height=18,width=24,units="cm", bg = 'white')

