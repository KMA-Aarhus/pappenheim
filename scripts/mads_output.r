library(tidyverse)

#' author: "Tine S. Ebsen, Carl M. Kobel"

############################################################################
# Read the metadata file and generate a csv that can be imported into mads #
############################################################################


# Parse command line arguments:
args = commandArgs(trailingOnly = TRUE)
write("These are the args:", stderr())
write(paste(args), stderr())
write("", stderr())

main_batch_id = args[1]
metadata_file = args[2]
file_out = args[3]
obey_quality_control = args[4]



# Quick and dirty development mode solution:
if (F) {
    rm(list = ls())
    
    setwd("~/GenomeDK/clinmicrocore/pappenheim-receiver/")
    
    
    
    # Testbatch
    main_batch_id = "20210324.1231"
    metadata_file = "output/20210324.1231/20210324.1231_metadata_pati.tsv"
    file_out = "output/20210324.1231/20210324.1231_mads.csv"
    obey_quality_control = "TRUE"
    # 
    s
    
    # main_batch_id = "20210428.1218"
    # metadata_file = "output/20210428.1218/20210428.1218_metadata_pati.tsv"
    # file_out = "output/20210428.1218/20210428.1218_mads.csv"
    # obey_quality_control = "TRUE"
    # 
    
    # I'm having a periodic problem with this sample:
    #singularity run docker://rocker/tidyverse Rscript 
    #main_batch_id = scripts/mads_output.r 
    main_batch_id = "20210428.1218" 
    metadata_file = "output/20210428.1218/20210428.1218_metadata_pati.tsv"
    file_out = "output/20210428.1218/20210428.1218_mads.csv" 
    obey_quality_control = TRUE
    
}





metadata = read_tsv(metadata_file) %>% 
    select(ya_sample_id, batch_id, type, lineage = pa_lineage, clade = ne_clade, totalMissing = ne_totalMissing, totalMissing_interpreted = ne_totalMissing_interpreted, aaSubstitutions = ne_aaSubstitutions, batch_control_stamp, batch_control_error_messages) %>% 
    
    # Create a vector with the substitutions
    mutate(aaSubstitutions = str_split(aaSubstitutions, ",")) %>% 
    #unnest(aaSubstitutions) %>% 
    #filter(str_sub(aaSubstitutions, 1,2) == "S:") %>%
    
    
    
    # Remove potential duplicates
    # This step is only necessary when you're using metadata for many batches, which we are not in this setting.
    # group_by(ya_sample_id) %>% 
    # arrange(desc(batch_id)) %>% # Put the newest batch first
    # mutate(rank = row_number(ya_sample_id)) %>% 
    # filter(rank == 1) %>% # Pick the newest. Conclusion: When you later filter for batch, you will have maked sure that each sample name can only be written into mads once.
# ungroup() %>%

# After checking that there are no duplicates, we can filter by sample
filter(type == "sample") %>% 
    
    
    # Pangolin marks lineage as "None" when it is missing. I want it to be coded as NA instead
    mutate(lineage = case_when(lineage != "None" ~ lineage))



#####################################################################################
# Make a function that can be used to neatly search for tree-like pangolin lineages #
#####################################################################################

lineage_detect = function(target_column, query_string) {
    
    # An internal function that converts a dot-containing pangolin-lineage to a properly escaped string for the ICU regex engine.
    # This function follows the tree-like structure of the lineage nomenclature so that the string must proceed with another dot after the last character if the string has not already reached EOL.
    dot_escape_lineage = function(input) {
        paste0(str_replace_all(input, "\\.", "\\\\."), "($|\\.)")
    }
    
    # Test that dot_escape_lineage() behaves as we expect it to.
    stopifnot(dot_escape_lineage("B.1.1.7") == "B\\.1\\.1\\.7($|\\.)")
    
    stopifnot(str_detect("B.1.1.7", dot_escape_lineage("B.1.1.7")))
    stopifnot(! str_detect("B.1.1.77", dot_escape_lineage("B.1.1.7")))
    stopifnot(str_detect("B.1.1.7.1", dot_escape_lineage("B.1.1.7")))
    stopifnot(str_detect("B.1.1.7.2", dot_escape_lineage("B.1.1.7")))
    
    
    # Since the dot_escape_lineage-function behaves well, we can wrap it with the str_detect-function
    str_detect(target_column, dot_escape_lineage(query_string))
    
}

# Test that the wrapped function (lineage_detect()) works well.
lineage_test_table = tibble(lineage = c("B.1.1.7",
                                        "B.1.1.77",
                                        "B.1.1.7.1",
                                        "B.1.1.7.2"))
stopifnot(
    lineage_detect(lineage_test_table$lineage, "B.1.1.7") == c(T, F, T, T)
)


# Check that the quality control is OK.
# Depending on your mood, you may switch the "obey_quality_control" variable so the script doesn't fail due to inadequate controls.
unique_batch_control_stamp = metadata$batch_control_stamp %>% unique
error_message = paste("This batch does not live up to the specified requirements. The negative or positive control(s) for this batch failed with the following message(s):", metadata$batch_control_error_messages %>% unique)

if (unique_batch_control_stamp == "passed") {
    write("All samples in this batch have passed QC control.")
} else if (obey_quality_control == "TRUE") {
    warning(paste("the unique_batch_control_stamp is", unique_batch_control_stamp))
    #stop(error_message)
} else {
    warning(error_message)
} 





# Throw a warning, if NAs are present in the imported table
nalen = metadata %>% filter(is.na(ya_sample_id)) %>% 
    pull(ya_sample_id) %>% 
    length

if (nalen != 0) {
    stop(paste("Imported table contains", nalen, "missing sample names"))
} else {
    write(paste("Info: all sample names are present in imported data table"), stderr())
}


# Make the WGS-table
out = metadata %>% 
    filter(!is.na(ya_sample_id)) %>% 
    
    
    # Again, Since I'm using data available in the sample sheet (only one batch_id) I don't need to filter out other batches.
    #filter(batch_id == main_batch_id) %>% 
    
    
    # Call WGS smitsomhed
    mutate(#inconclusive = if_else((is.na(lineage) & is.na(clade)) | batch_control_stamp != "passed", T, F), # If lineage or clade is NA, the sample should answer as inconclusive
        inconclusive = case_when(
            is.na(lineage) | is.na(clade) | batch_control_stamp == "failed" | totalMissing_interpreted >= 5000 ~ T,
            TRUE ~ F
        ),
        MDSU = 32092) %>% 
    
    
    
    
    
    
    
    # This scheme is up to date with "Typningstabel_28042021"  
    rowwise() %>% # rowwise() groups the data into rows so that the %in% test works as intended.
    # Because this script fails completely if the batch does not live up to the quality thresholds, we do not need to look at the inconclusive column when triggering text.

    
    
    mutate(WGS_linje = case_when(inconclusive ~ "WGS: Inkonklusiv",
                                 TRUE ~paste0(paste0("WGS: ", sub("\\).*", "", sub(".*\\(", "", clade)), " (", sub("\\ .*", "", clade),"), ", lineage)))) %>% 
    ungroup()

# Ensure that the output is obvious to the clinic receiving it
out = out %>%
    mutate(WGS_linje = case_when(grepl("BE",WGS_linje) ~ paste0(sub("BE", "BA.5.3.1", WGS_linje), " (",lineage,")"),
                                 grepl("BG",WGS_linje) ~ paste0(sub("BG", "BA.2.12.1", WGS_linje), " (",lineage,")"),
                                 TRUE ~WGS_linje))

write("This is the results and their reasons:", stdout())
out %>% 
    select(ya_sample_id, WGS_linje) %>% 
    format_tsv() %>% write(stdout())


# Finally format the results for Mads-import.
mads_out = out %>% 
    select(`sample-id` = ya_sample_id, MDSU, WGS_linje) %>% 
    pivot_longer(c(WGS_linje)) 


# Before writing to disk, write to stdout log for easy debugging..
mads_out %>% 
    format_tsv %>% 
    write(stdout())

# Write output file to file given in variable file out.
mads_out %>% 
    write.table(file_out, quote = F, sep = ";", fileEncoding = "cp1252", row.names = F)


