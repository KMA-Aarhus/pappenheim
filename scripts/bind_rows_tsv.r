library(tidyverse)

#' author: "Carl M. Kobel"

############################################################################################
# This script takes a list of tsv-files, binds them together, and spits them out of stdout #
############################################################################################



# -- Read Data --------------------------------------------

# Parse arguments:
args = commandArgs(trailingOnly=TRUE)
write("These are the args:", stderr())
write(paste("arg:", args), stderr())
write("", stderr())

# Sort the input files. 
tsv_files = args %>% sort




# Development mode:
if (FALSE) {
    rm(list = ls())
    setwd("~/GenomeDK/clinmicrocore/pappenheim-receiver/")
    
    tsv_files_string = "output/20210419.1121/20210419.1121_metadata_init.tsv output/20210423.1240/20210423.1240_metadata_init.tsv output/20210414.1219/20210414.1219_metadata_init.tsv"
    tsv_files = str_split(tsv_files_string, ' ')[[1]] # Convert input string to list}
        
}
        
  
data = tibble() # Create an empty tibble we can iteratively add tables to
    

for (file_ in tsv_files) {
    
    write(paste('Reading', file_), stderr())
    current_table = read_tsv(file = file_, col_types = cols(.default = 'c')) # Read the current file_-variable table
    write(class(current_table$batch_id), stderr())
    write(paste('Adding', dim(current_table)[1], 'records'), stderr()) 
    data = data %>% bind_rows(current_table) # Bind the current table onto the collected data table
}


    
# Write the collected data table to stdout
data %>% 
    format_tsv %>% 
    write(stdout())
