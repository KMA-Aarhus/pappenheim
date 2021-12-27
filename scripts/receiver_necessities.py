from pathlib import Path
import os
import pandas as pd
import sys


if __name__ == "__main__":


    out_base = sys.argv[1]

    # Set the directory where clean uploads are deposited from the pappenheim workstation pipeline:
    path = sys.arg[2]

    out_path = sys.argv[3]

    bid = sys.argv[4] # batch_id


    ##################################
    # Parse clean upload directories #
    ##################################


    print("Parsing input directories from clean_dir")
    dirs = [e for e in path.iterdir() if e.is_dir()]

    df = pd.DataFrame(data = {"directory": [i for i in dirs]})
    df["directory_basename"] = df.apply(lambda row: os.path.basename(row.directory), axis=1)

    clean_upload_dir_prefix = "clean_upload_"
    df = df[df["directory_basename"].str.startswith(clean_upload_dir_prefix)]
    df["batch_id"] = df["directory_basename"].str[len(clean_upload_dir_prefix):]
    df["batch_date"] = df.apply(lambda row: str(row.batch_id).split(".")[0], axis = 1)
    df["long_metadata"] =  df.apply(lambda row: str(row.directory) + "/" + row.batch_id + "_all.tsv", axis=1)



    ##################################
    # Present the pipeline dataframe #
    ##################################

    df = df.sort_values(by=['batch_id']).reset_index()
    print()
    print("df")
    print(df[["batch_id", "directory"]])
    print("//")

    file_path = out_path + bid + "_df.csv"
    df.to_csv(path_or_buf=file_path, index=False)