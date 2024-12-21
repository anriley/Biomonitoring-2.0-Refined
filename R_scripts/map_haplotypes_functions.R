#Parts of this code where created by Maya Wetzl
library(raster)
library(sf)
library(elevatr)
library(ggmap)
library(ggplot2)
library(ggnewscale)
library(scatterpie)
library(pegas)
library(seqinr)
library(vegan)
library(tidyterra)

get_elevation_map <- function(long_min,long_max,lat_min,lat_max){
  target_range <- data.frame(x=c(long_min,long_max),y=c(lat_min,lat_max)) %>% st_as_sf(coords = c("x","y"), crs = 4326) 
  target_map <- get_elev_raster(target_range, z=7) #,long_min=-118,long_max=-113,lat_min=49, lat_max=52 #z=7 is ideal? z is the zoom on the map, higher z gives better resolution up to 14. A link to resolutions at different z : https://github.com/tilezen/joerd/blob/master/docs/data-sources.md#what-is-the-ground-resolution
  
  
  target_map_df <- as.data.frame(target_map, xy=TRUE)
  colnames(target_map_df) <- c("Long","Lat","Elevation")
  #plot(target_map)
  grad <- hypso.colors(10, "colombia_hypso")
  plotMap <- ggplot(data = target_map_df)+
    geom_raster(mapping=aes(x=Long, y=Lat, fill=target_map_df[,3])) +
    scale_fill_gradientn(colours= grad, name='Elevation',na.value = NA) +
    theme_bw() +
    theme(axis.line = element_line(colour = "black"),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.border = element_blank(),
          panel.background = element_blank())
  return(plotMap)
}

gg_color_hue <- function(n) {
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
} #https://stackoverflow.com/questions/8197559/emulate-ggplot2-default-color-palette


create_hap_map_pa <- function(read_data,site_data,site_map,zotu_colours,save_plots=FALSE,outfile="empty",log_base,size_coefficient){
  read_data[read_data != 0] <- 1
  data_frame_merge <- merge(read_data, site_data,by = 'row.names')
  rownames(data_frame_merge) <- data_frame_merge[,1]
  data_frame_merge <- data_frame_merge[,-1]
  a2 <- data_frame_merge
  
  a2 <- na.omit(a2)
  zotu_names <- colnames(read_data)
  hap_map_plot <- site_map +new_scale("fill")+ geom_scatterpie(aes(x=Long, y=Lat,), data=a2,cols=zotu_names,  alpha=0.7) + 
    coord_equal() + scale_fill_manual(limits = zotu_names, values=zotu_colours)
  if(save_plots){
    ggsave(outfile,hap_map_plot,height=14,width=23,units="cm")
  }
  return(hap_map_plot)
}
