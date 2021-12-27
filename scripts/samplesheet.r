library(tidyverse)

#' author: "Carl M. Kobel"


################################################################################################################################
# This script takes variant data in, and outputs a script which helps compiling a compressed file which can be uploaded to SSI #
################################################################################################################################

args = commandArgs(trailingOnly=TRUE)

write("These are the args:", stderr())
write(paste(args), stderr())
write("", stderr())

main_batch_id = args[1]
file_metadata_pati = args[2]
dir_clean_upload = args[3]
copy_consensuses_script = args[4]
file_samplesheet_out = args[5]



# Debug variables for testing locally.
if (F) {
    rm(list = ls())
    
    setwd("~/GenomeDK/clinmicrocore/pappenheim-receiver/")
    
    #main_batch_id = "20210324.1231"
    #file_metadata_pati = "~/GenomeDK/clinmicrocore/pappenheim-receiver/output/20210324.1231/20210324.1231_metadata_pati.tsv"
    #dir_clean_upload = "../BACKUP/nanopore_sarscov2/pappenheim_clean/clean_upload_20210324.1231/"
    #copy_consensuses_script = "copy_consensuses_script.sh"
    #file_samplesheet_out = "file_samplesheet_out.tsv"
    
    
    main_batch_id = "20210324.1231"
    file_metadata_pati = "output/20210324.1231/20210324.1231_metadata_pati.tsv"
    dir_clean_upload = "../BACKUP/nanopore_sarscov2/pappenheim_clean/clean_upload_20210324.1231"
    copy_consensuses_script = "output/20210324.1231/20210324.1231_copy_consensuses.sh"
    file_samplesheet_out = "output/20210324.1231/20210324.1231_samplesheet.tsv"
    
}

date = (main_batch_id %>% str_split("\\."))[[1]][1] # Pick out the date (YYMMDD) from the main batch id.

##############################################
# Prepare files for upload to the government # 
##############################################


# Import metadata_pati
df_metadata_pati = read_tsv(file_metadata_pati, col_types = cols(.default = 'c')) %>% 
    mutate(ct = coalesce(as.character(ss_ct), as.character(ct)))

# Make the file ready for export to the government
## Filter for samples only, and select/rename the columns of interest.
# ”sample_id;cpr;sampling_date;kma_id;raw_filename;consensus_filename”.
write(paste("Writing sample sheet to", file_samplesheet_out), stderr())
df_sample_sheet = df_metadata_pati %>%
    filter(type == "sample") %>% 
    #filter(batch_control_stamp == "failed") %>%  
    rowwise() %>% 
    mutate(kma_id = "6620320",
           sample_id = full_name, #sample_id,
           #raw_full_name = paste0(batch, ".", plate, ".", moma_serial, "_", raw_sample_name),
           #raw_filename = paste0(raw_full_name, "_R", c(1, 2), ".fastq.gz", collapse = ";"),
           raw_filename = NA,
           #consensus_filename = paste0(sample_id, ".fa", collapse = " "), # hvorfor er den her collapsed?
           #consensus_filename = paste0(sample_id, ".fasta"),
           platform = "nanopore") %>% 
    ungroup() %>% 
    select(full_name, sample_id, cpr = `cprnr.`, sampling_date = afsendt, kma_id, raw_filename, consensus_filename = upload_consensus_file, platform, ct) # Consider including Ydernr/SKSnr





df_sample_sheet %>%  
    select(-full_name) %>% 
    write_tsv(file_samplesheet_out)



#######################################
# Fail if required values are missing #
######################################

# The following block has been disabled because I handle errors in the pipeline directly. See rule targz.
# # If any of the samples are missing a cpr-number, this target should fail.
# n_missing_cpr = (df_sample_sheet %>% filter(is.na(cpr) | is.na(sampling_date)) %>% dim)[1]
# if (n_missing_cpr > 0) {
#     write(paste("One or more samples are missing cpr or sampling date. Possibly because the mads-extract is out-dated."), stderr())
#     write("selected df_sample_sheet columns:", stderr())
#     write("", stderr())
#     df_sample_sheet %>% 
#         select(full_name, cpr, sampling_date) %>% 
#         format_tsv %>% 
#         write(stderr())
#     
#     stop(paste0("Update the mads-extract, or enter the missing values manually into \"", file_metadata_pati, "\"and rerun the pipeline."))   
# } else {
#     write("Sample metadata is OK", stderr())
# } 




#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# Generate a bash-script that can be used to tar the files together # 
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

# It makes sense to do it from here, because the relevant metadata is already loaded in the environment.
target_consensus = paste0("output/", main_batch_id, "/copy_consensuses/")


# Copy Consensus data
write("commanding consensuses ...", stderr())
command_consensus = df_sample_sheet %>% 
    select(full_name, consensus_filename) %>% 
    rowwise() %>% 
    mutate(source_file = paste0(dir_clean_upload, "/", consensus_filename),
           #target_dir = target_consensus,
           #target_basename = paste0(consensus_filename),
           target_file = paste0(target_consensus)) %>% 
    ungroup() %>% 
    
    #transmute(command = paste0("cp ", source_file, " ", target_dir, target_basename))
    transmute(command = paste("cp", source_file, target_file))

# Have a look at the consensus commands
#command_consensus$command[1]

command_consensus %>% select(`#!/bin/bash` = command) %>% write_tsv(copy_consensuses_script) # TODO: Make as a space-delimited file instead of tsv.





