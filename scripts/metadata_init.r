library(tidyverse)


#' author: "Carl M. Kobel"

#############################################################################################################################
# This script takes the long metadata directly from pappenheim (run on workstations) and converts the table to wide format. # 
#############################################################################################################################



#############
# Read data #
#############

# Parse arguments:
args = commandArgs(trailingOnly=TRUE)
write("These are the args:", stderr())
write(paste(args), stderr())
write("", stderr())


main_batch_id = args[1] 
file_prefix = args[2]
file_out = args[3]

    
#main_batch_id = "20211207.1456"
#file_prefix = "/media/ontseq3/ssd/sc2_sequencing/2021-12-07/no_sample/20211207_1456_MN36352_FAR57225_ad550f16/pappenheim_output/collected/20211207.1456"
#file_out = "/media/ontseq3/ssd/sc2_sequencing/2021-12-07/no_sample/20211207_1456_MN36352_FAR57225_ad550f16/pappenheim_output/20211207.1456/20211207.1456_metadata_init.tsv"
# Development mode:
if (FALSE) {
    rm(list = ls())
    setwd("~/GenomeDK/clinmicrocore/pappenheim-receiver/")
    
    
    # Test with 2021032
    #Rscript scripts/metadata_init.r 
    # main_batch_id = "20210324.1231"
    # file_prefix = "~/GenomeDK/clinmicrocore/BACKUP/nanopore_sarscov2/pappenheim_clean/clean_upload_20210324.1231/20210324.1231"
    # file_out = "testout.tsv"
    # 
    
    # Test with 20210419.1121
    #Rscript scripts/metadata_init.r
    #main_batch_id = "20210419.1121"
    #file_prefix = "../BACKUP/nanopore_sarscov2/pappenheim_clean/clean_upload_20210419.1121/20210419.1121" 
    #file_out = "output/20210419.1121/20210419.1121_metadata_init.tsv"
    
    
    # Test with 20210503.1230
    #Rscript scripts/metadata_init.r 
    main_batch_id = "20210503.1230"
    file_prefix = "../BACKUP/nanopore_sarscov2/pappenheim_clean/clean_upload_20210503.1230/20210503.1230" 
    file_out = "output/20210503.1230/20210503.1230_metadata_init.tsv"
    
    
    # Another batch I'm having problems with:
    #Rscript scripts/metadata_init.r 
    main_batch_id = "20210518.1231" 
    file_prefix = "../BACKUP/nanopore_sarscov2/pappenheim_clean/clean_upload_20210518.1231/20210518.1231" 
    file_out = "output/20210518.1231/20210518.1231_metadata_init.tsv"
    
}



###############################
# Read the long metadata file #
###############################

df_long = read_delim(paste0(file_prefix, "_all.tsv"), delim = "\t", col_types = cols(.default = col_character())) %>%
        rename(batch_id = `#batch_id`) %>% 
    mutate_at(vars(batch_id), as.character)
 write("I have used dplyr", stderr())   


# Since we can't trust that the technicians necessarily provide the needed columns, we have to force them in here.

# Define columns which must be present in the samplesheet given by the technicians.
required_columns = tibble(batch_id = NA,
                          sample_id = NA,
                          rack_id = NA,
                          position = NA,
                          #cavity_id = NA,
                          #position = NA,
                          #concentration = NA,
                          #concentrationunit = NA,
                          #volume = NA,
                          `qubit_ng/ul` = NA,
                          #fortynding = NA,
                          #`ul_prøve` = NA,
                          #korrigeret = NA,
                          #`µl_mq` = NA,
                          ct = NA,
                          note = NA,
                          barcode = NA)

unwanted_columns = tibble(cavity_id = NA,
                          concentration = NA,
                          concentrationunit = NA,
                          volume = NA,
                          fortynding = NA,
                          `ul_prøve` = NA,
                          korrigeret = NA,
                          `µl_mq` = NA,
                          userdefined1 = NA,
                          plate = NA)


# Pivot to wide format
df_wide = df_long %>%
  pivot_wider(id_cols = c(batch_id, sample_id), names_from = variable, values_from = value, values_fn = first) %>%
    # Make sure that the required columns (from the samplesheet) are present.
  add_column(!!!required_columns[!names(required_columns) %in% names(.)]) %>% 
    # Remove unwanted columns
  select(-(names(unwanted_columns[names(unwanted_columns) %in% names(.)]))) %>%
  # Remove unwanted columns. These are columns from the tecan machine that are no longer needed. 
  select(-starts_with("unnamed:_")) %>%
    mutate(ne_totalMissing_interpreted = case_when(is.na(ne_totalMissing) ~ 29903,
                                                   TRUE ~ as.numeric(ne_totalMissing)),
           sample_name_prefix = str_sub(sample_id, 1, 2), # first two characters
           sample_name_suffix = str_sub(sample_id, 3), # the rest
           full_name = paste0(batch_id, "_", sample_id),
           sample_name_prefix_converted = recode(sample_name_prefix,
                                                 "87" = paste0("L"),
                                                 "88" = paste0("R"),
                                                 "89" = paste0("I"),
                                                 "90" = paste0("V"),
                                                 "96" = paste0("P"))) %>% 
    
    rowwise() %>% 
    mutate(ya_sample_id = paste0(sample_name_prefix_converted, sample_name_suffix)) %>% 
    ungroup() %>% 
    
    
    select(-sample_name_prefix, -sample_name_suffix, -sample_name_prefix_converted) %>% 
    select(full_name, batch_id, sample_id, ya_sample_id, type, everything()) 




# Read barcode_alignment.tsv which contains the number of reads for each sample
barcode_alignment = read_delim(paste0(file_prefix, "_barcode_alignment.tsv"), delim = "\t") %>% 
    select(barcode_basename = barcode, ba_type = type, ba_reads = target_unclassified, started)

sum_reads = barcode_alignment %>% pull(ba_reads) %>% sum
unclassified_reads = barcode_alignment %>% filter(ba_type == "na") %>% pull(ba_reads)

# Join barcode_alignment onto df_wide, and add information about the number of reads.
df_wide_ba = df_wide %>% 
    left_join(barcode_alignment) %>% 
    mutate(batch_sum_reads = sum_reads,
           batch_unclassified_reads = unclassified_reads, 
           batch_unclassified_reads_prop = unclassified_reads/batch_sum_reads)



# Read final_summary.txt, which contains information about time.
terminal_summary = read_delim(paste0(file_prefix, "_final_summary.txt"), delim = "=", col_names = c("name", "value")) %>% 
    pivot_wider(names_from = name, values_from = value) %>% 

    mutate_at(vars(started, acquisition_stopped, processing_stopped), lubridate::ymd_hms) %>% 
    # Calculate run time
    #Vi
    identity




# Incorporate time information
# Parse hours sequencing
start_time = terminal_summary %>% pull(started)
end_time = terminal_summary %>% pull(acquisition_stopped) %>% head(1)
#sequencing_time = end_time-start_time
hours_sequencing = lubridate::interval(start_time, end_time) / lubridate::hours(1)

# add to the main df
df_wide_ba_hours= df_wide_ba %>% 
    mutate(batch_hours_sequencing = hours_sequencing,
           batch_million_reads_per_day = (batch_sum_reads/1000000)/(batch_hours_sequencing/24))



##################################
# Check that the controls are OK #
##################################

# Make a summary that can be used for checking the presence of controls.
type_summary = df_wide_ba %>% group_by(very_long_batch_id, type) %>%
    summarize(n = length(type), min_totalMissing = min(ne_totalMissing_interpreted), max_totalMissing = max(ne_totalMissing_interpreted)) # If there is more than one positive or negative control, we need to know the min and max-values of them.


error_messages = c() # Empty vector to collect error_messages.

# TODO: Test this thoroughly!

# Positive control presence
if(type_summary %>% filter(type == "positive_control") %>% nrow < 1) { 
    error_messages = c(error_messages, "No positive control present.")
} else if (type_summary %>% filter(type == "positive_control") %>% pull(min_totalMissing) >= 3000) { # Positive control coverage
    error_messages = c(error_messages, "The positive control totalMissing is above the critical threshold.") 
}


# Negative control presence
if(type_summary %>% filter(type == "negative_control") %>% nrow < 1) { 
    error_messages = c(error_messages, "No negative control present.")
} else if (type_summary %>% filter(type == "negative_control") %>% pull(max_totalMissing) <= 29903/2) { # Negative control coverage
    error_messages = c(error_messages, "The negative control totalMissing is below the critical threshold.")
}


if(length(error_messages) > 0) {
    error_flag = TRUE
    write("During processing of the incoming samples, the following control problems were encountered:", stderr())
    write(paste("Error:", error_messages), stderr())
} else {
    error_flag = FALSE
    error_messages = "QC OK"
    write("Info: During processing of the incoming samples, no control problems were found :)", stderr())
}

df_wide_ba_hours_stamp = df_wide_ba_hours %>% 
    mutate(batch_control_stamp = case_when(error_flag ~ "failed",
                                           TRUE ~ "passed"),
           batch_control_error_messages = paste(error_messages, collapse = ", "))


# Write the init file out.
df_wide_ba_hours_stamp %>% write_tsv(file_out)

