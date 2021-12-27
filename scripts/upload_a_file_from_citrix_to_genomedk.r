##################################################
# Upload thenewest 116 file (a) to the correct path on the cluster

.libPaths("H:\\Documents\\R")
library(ssh)


conn = ssh_connect("BRUGER@login.genome.au.dk")






input_file = "../Desktop/mads_116.csv"

finfo = file.info(input_file)
mtime = as.Date(finfo$mtime, format = "%y%m%d")
basename_out = paste0("mads_116_", mtime, ".csv")


out_paths = c(#"~/clinmicrocore/pipe19/batch/mads/latest/", 
              "~/clinmicrocore/pappenheim-receiver/mads/in/")
for (out_path in out_paths) {
    
    write(paste("uploading", out_path), stderr())
    file_out = paste0(out_path, basename_out)
    
    scp_upload(conn, input_file, to = file_out)
}    

ssh_disconnect(conn)
