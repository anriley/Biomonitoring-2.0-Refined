source("map_haplotypes_functions.R")

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
isolated_species_tables <- readRDS(file=paste(indir,"RDS/sbc_tables_",variables_io,".RDS",sep=""))
metadata <- important_tables_in[["metadata"]]

site_data2 <- as.data.frame(unique(cbind(ID=metadata[["ID"]],
                                         Long=metadata[["longitude"]],Lat=metadata[["latitude"]])))
rownames(site_data2) <- site_data2[["ID"]]
site_data2 <- site_data2[,-1]
site_data3 <- data.frame(apply(site_data2, 2, function(x) as.numeric(as.character(x))))
rownames(site_data3) <- rownames(site_data2)

#Mean region location
metadata_temp_mean <-metadata[,-1]
metadata_temp_mean$latitude <-as.numeric(as.character(metadata_temp_mean$latitude))
metadata_temp_mean$longitude <-as.numeric(as.character(metadata_temp_mean$longitude))

metadata_mean <- aggregate(.~cluster_region,metadata_temp_mean,FUN=mean)
rownames(metadata_mean) <- metadata_mean$cluster_region
metadata_mean <- metadata_mean[,-1]
colnames(metadata_mean) <- c("Lat","Long") 

#Map
lat_min <- min(as.numeric(metadata$latitude))
lat_max <- max(as.numeric(metadata$latitude))
long_min <- min(as.numeric(metadata$longitude))
long_max <- max(as.numeric(metadata$longitude))


lat_min_adj <-lat_min -(lat_max-lat_min) * 0.05
lat_max_adj <-lat_max +(lat_max-lat_min) * 0.05

long_min_adj <-long_min -(long_max-long_min) * 0.05
long_max_adj <-long_max +(long_max-long_min) * 0.05


#site_map <- get_stamen_map(long_min, long_max, lat_min, lat_max, z = 8)
#site_map_adj <- get_stamen_map(long_min_adj, long_max_adj, lat_min_adj, lat_max_adj, z = 8)

elev_map <- get_elevation_map(long_min,long_max,lat_min,lat_max)
elev_map_adj <- get_elevation_map(long_min_adj, long_max_adj, lat_min_adj, lat_max_adj)

sp_table_F230R <- isolated_species_tables[["F230R"]][["Yoraperla_brevis"]]
colnames(sp_table_F230R) <- c("YB_F230R_1","YB_F230R_2")

sp_table_MLJG <- isolated_species_tables[["MLJG"]][["Yoraperla_brevis"]]
colnames(sp_table_MLJG) <- c("YB_MLJG_1","YB_MLJG_2")

F230R_MLJG <- merge(sp_table_F230R,sp_table_MLJG,by="row.names",all=TRUE)
rownames(F230R_MLJG) <- F230R_MLJG$Row.names
F230R_MLJG <- F230R_MLJG[,-1]
F230R_MLJG[is.na(F230R_MLJG)] <- 0
{zotus <- colnames(F230R_MLJG)
  colours <- gg_color_hue(ncol(F230R_MLJG))
  colour_list <- list()
  for (zotu in zotus){
    current_colour <- colours[1]
    colours <- colours[-1]
    colour_list[[zotu]] <- current_colour
  }
  zotu_colour_list <- colour_list
}

barcode_metadata <- read.csv(file=paste(indir,"Data/yp_metadata_barcodes.csv",sep=""),row.names = 1)

read_data <- F230R_MLJG
site_data <- site_data3
site_map <- elev_map
zotu_colours <- zotu_colour_list
read_data[read_data != 0] <- 1
data_frame_merge <- merge(read_data, site_data,by = 'row.names')
rownames(data_frame_merge) <- data_frame_merge[,1]
data_frame_merge <- data_frame_merge[,-1]
a2 <- data_frame_merge

a2 <- na.omit(a2)
zotu_names <- colnames(read_data)
hap_map_plot <- site_map +geom_point(aes(x=Long,y=Lat,shape=Barcodes),data=barcode_metadata,size=3)+new_scale("fill")+ geom_scatterpie(aes(x=Long, y=Lat,), data=a2,cols=zotu_names,  alpha=0.7) + 
  coord_equal() + scale_fill_manual(limits = zotu_names, values=zotu_colours) + labs(fill='ESVs',shape="Barcodes")+ 
  guides(shape = guide_legend(order = 2),fill = guide_legend(order = 1))

ggsave(paste(indir,"Figures/Figure7_ybrevis_map.png",sep=""),
       hap_map_plot,height=14,width=24,units="cm", bg = 'white')





