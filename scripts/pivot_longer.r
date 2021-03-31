args = commandArgs(trailingOnly=TRUE)
write("These are the args:", stderr())
write(paste(args), stderr())
write("", stderr())

input_file = args[1]
sample_id = args[2]


devel_mode = F
# devel_mode = T
if (devel_mode == T) {
    input_file = "~/repos/pappenheim/opentest/pangolin.csv"
    input_file = "~/repos/pappenheim/opentest/nextclade.tsv"
    
    sample_id = "12345"
}



if (str_detect(input_file, "\\.csv$")) {
    write("Reading .csv-format", stderr())
    delimiter = ","
} else if (str_detect(input_file, "\\.tsv$")) {
    write("Reading .tsv-format", stderr())
    delimiter = "\t"
} else {
    stop("The input file format is not supported.")
} 

input_df = read_delim(input_file, col_types = cols(.default = "c"), delim = delimiter)

input_df %>% 
    mutate(sample_id = sample_id) %>% 
    pivot_longer(-sample_id) %>% View
    
    format_tsv


