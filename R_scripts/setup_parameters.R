indir <- "../"

parameter_list <- list()

#These variables can only be changed if they were run in the previous python step
parameter_list[["read_filter"]] <- 100 #Threshold for an ESV to be retained
parameter_list[["sBp"]] <- 0.8 #Threshold for species assignment, can be changed for assignment, but not isolated species

#Filters
parameter_list[["min_hap"]] <- 2 #Minimum number of haplotypes required for each species (less than 2 and no within species comparison can be made) 
parameter_list[["sample_filter"]] <- 0 #Threshold for minimum number of total reads in a sample to be retained 
parameter_list[["within_sample_filter"]] <- 5 #Threshold for minimum number of reads for an ESV to be retained in a sample

parameter_list[["dm_method"]] <- "pa" 
parameter_list[["beta_component"]] <- "all"

parameter_list[["min_sites"]] <- 3

#Amplicons
parameter_list[["amps"]] <- c("F230R","MLJG")
parameter_list[["d_amp"]] <- list("F230R"=7,"MLJG"=9) #Amplicon d threshold to use for OTUs

variables_io <- paste("assign",parameter_list[["sBp"]],"filter",parameter_list[["read_filter"]],
                      parameter_list[["sample_filter"]],parameter_list[["within_sample_filter"]],
                      "method",parameter_list[["dm_method"]],parameter_list[["beta_component"]],sep="_")

saveRDS(parameter_list,paste(indir,"parameters.RDS",sep=""))
