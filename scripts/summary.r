library(tidyverse)

#' author: "Carl M. Kobel

#######################################################################################################################################
 
#######################################################################################################################################

args = commandArgs(trailingOnly=TRUE)

write("These are the args:", stderr())
write(paste(args), stderr())
write("", stderr())

metadata_pati_in = args[1]



if (F) {
    rm(list = ls())
    
    setwd("~/GenomeDK/clinmicrocore/pappenheim-receiver/")
    
    # Run for batch 20210422.1311
    #metadata_pati_in = "output/20210422.1311/20210422.1311_metadata_pati.tsv"
    
    # Batch med flere plader 
    metadata_pati_in ="output/20210505.0652/20210505.0652_metadata_pati.tsv"
}


write("importing metadata ...", stderr())
metadata_pati = read_tsv(metadata_pati_in)


antiparseprotection = function(input) {
    paste0("\"", input, "\"")
}

now = lubridate::today()

out = metadata_pati %>% filter(type == "sample") %>% 
    group_by(very_long_batch_id, rack_id) %>% 
    summarize(#`Pappenheim very long batch id` = very_long_batch_id,
              Maskine = str_extract(barcode_path, "ontseq.") %>% na.omit() %>% unique(),
              `Sekventeringstid (timer)` = batch_hours_sequencing %>% na.omit() %>% unique() %>% round(1),
              `Antal patientprøver` = type %>% length(),
              `Antal dækning høj` = (ne_totalMissing_interpreted <= 3000) %>% na.omit() %>% sum(),
              `Antal dækning over halv` = (ne_totalMissing_interpreted <= 29903/2) %>% na.omit() %>% sum(),
              `Antal dækning under halv` = (ne_totalMissing_interpreted > 29903/2) %>% na.omit() %>% sum(),
              `Mads uploaddato` = antiparseprotection(now),
              `Mads uploaduge` = lubridate::isoweek(now),
              `SSI uploaddato` = antiparseprotection(now), 
              `SSI uploaduge` = lubridate::isoweek(now)) %>% 

    rename(`Pappenheim very long batch id` = very_long_batch_id,
           `Registreret plade` = rack_id)


write("Processing of data is done. Printing to stdout ...", stderr())
out %>% format_tsv %>% 
    write(stdout())
