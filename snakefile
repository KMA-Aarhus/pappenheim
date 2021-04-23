
__author__ = "Carl Mathias Kobel"
__version__ = "0.1"


# start_pappenheim '/run/user/1000/gvfs/smb-share:server=onerm.dk,share=nfpdata/Afdeling/AUHKLMIK/AUH/Afdelingen/Afsnit Molekyl. og Serologi/NanoPore/NanoPore Metadata/5770_seq_2021-03-23_NEB.xlsx'  '/home/ontseqa/Desktop/sc2_sequencing/COVID19-AUH-20210324-5770/' -np


import sys
import os
from os import listdir
from os.path import isfile, isdir, join
import yaml
import pandas as pd
import numpy as np
from pandas_ods_reader import read_ods
from datetime import datetime
import glob
import time
import atexit
import datetime




# When development_mode is True, the development cycle (frequency) is increased.
development_mode = True


configfile: "config.yaml"
print(config["samplesheet"])

# TODO: These variables should be set from the command line
samplesheet = "../testdata/5351_seq_2021-03-09_NEB.xlsx"
rundir = "../../GenomeDK/clinmicrocore/BACKUP/nanopore_sarscov2/COVID19-AUH-20210316-NEB/rawdata/20210316_1259_MN34697_FAP10653_02939c92/"
rundir = "../testdata"


# Actually given on the command line:
samplesheet = config["samplesheet"]
rundir = config["rundir"]


tab = "\t"
nl = "\n"



# This will make your code slower
def lag(time_ = 0.06):
    if development_mode:
        pass
    else:
        time.sleep(time_)



# Print a convincing logo, if the windows is big enough.

terminal_rows, terminal_columns = os.popen('stty size', 'r').read().split()

if int(terminal_columns) >= 126:
    print()
    print(f"          Pappenheim pipeline v{__version__}  -  Aarhus Universityhospital  -  Department of Clinical Microbiology          ")
    print()
    print("                  ██████╗  █████╗ ██████╗ ██████╗ ███████╗███╗   ██╗██╗  ██╗███████╗██╗███╗   ███╗ ")
    print("                  ██╔══██╗██╔══██╗██╔══██╗██╔══██╗██╔════╝████╗  ██║██║  ██║██╔════╝██║████╗ ████║ ")
    print("                  ██████╔╝███████║██████╔╝██████╔╝█████╗  ██╔██╗ ██║███████║█████╗  ██║██╔████╔██║ ")
    print("                  ██╔═══╝ ██╔══██║██╔═══╝ ██╔═══╝ ██╔══╝  ██║╚██╗██║██╔══██║██╔══╝  ██║██║╚██╔╝██║ ")
    print("                  ██║     ██║  ██║██║     ██║     ███████╗██║ ╚████║██║  ██║███████╗██║██║ ╚═╝ ██║ ")
    print("                  ╚═╝     ╚═╝  ╚═╝╚═╝     ╚═╝     ╚══════╝╚═╝  ╚═══╝╚═╝  ╚═╝╚══════╝╚═╝╚═╝     ╚═╝ ")
    print()
    print("                                      Press ctrl+c at any time to stop this pipeline."              )
    print()

else:
    print("Warning: Please increase your terminal window width.")





# Check that input was given.
if config["samplesheet"] == "NA":
    raise Exception("No samplesheet file was given. Please specify a samplesheet by appending --config rundir=\"path/to/samplesheet/\" to the command line call.")
if config["rundir"] == "NA":
    raise Exception("No rundir path was given. Please specify a rundir by appending --config rundir=\"path/to/rundir/\" to the command line call.")
# TODO: Implement additional input validation, like checking that the objects given are file and dir respectively.



print(f"These are the parameters given:")
print(f"  samplesheet: {samplesheet}")
print(f"  rundir: {rundir}")
print()


def check_user(prompt):
    input_OK = input(f"{prompt} [y/n]: ").strip()
    print(f"User entered \"{input_OK}\": ", end = "", flush = True)
    #sys.stdin.read(1)
    if not str(input_OK).lower()[0:1] == "y":
        print("Exiting ...")
        print(f"Hint: Press <uparrow> <enter> to rerun this pipeline with the same input files.")
        exit(1)

    print("Proceeding ...")
    print()




#########################
# Parse the samplesheet #
#########################

samplesheet_extension = samplesheet.split(".")[-1]
print(f"Reading .{samplesheet_extension}-type sample sheet \"{samplesheet}\"")
if samplesheet_extension == "ods":

    raise Exception(f"The spreadsheet file extension {samplesheet_extension} is not yet implemented.")


    # # load a sheet based on its index (1 based)
    # sheet_idx = 1
    # df = read_ods(samplesheet, sheet_idx)


    # # load a file that does not contain a header row
    # # if no columns are provided, they will be numbered
    # df = read_ods(samplesheet, 1, headers=True)

    # # Clean up the spreadsheet
    # df.columns = map(str.lower, df.columns) # Lowercase
    # df.columns = map(str.strip, df.columns) # Remove edge-spaces
    # df.columns = map(lambda x: str(x).replace(" ", "_"), df.columns) # Replace spaces with underscore
    # df = df[["barcode", "sample_id"]]
    # df = df.astype(str)


elif samplesheet_extension == "xlsx":
    # Uses openpyxl
    df = pd.read_excel(samplesheet, dtype = str)

elif samplesheet_extension == "xls":
    df = pd.read_excel(samplesheet)

else:
    raise Exception(f"The spreadsheet file extension {samplesheet_extension} is not yet implemented.")

# Clean up the spreadsheet
print("Cleaning sample sheet ...                              ", end = "", flush = True)
lag()
df.columns = map(str.lower, df.columns) # Lowercase
df.columns = map(str.strip, df.columns) # Remove edge-spaces
df.columns = map(lambda x: str(x).replace(" ", "_"), df.columns) # Replace spaces with underscore
df["barcode"] = df["barcode"].apply(np.vectorize(lambda x: str(x).strip().replace(" ", ""))) # Because we are later going to join using this column, it is necessary to strip it for spaces.
df = df.dropna(subset = ["sample_id"])# remove rows not containing a sample ID
print("✓")

#print("ee", "barcode" in list(df.columns))
# Check that the spreadsheet complies
print("Checking that the necessary columns exist ...          ", end = "", flush = True)
lag()
for i in ["barcode", "sample_id"]:
    if not i in df.columns:
        raise Exception(f"The sample sheet is missing a necessary column. The sample sheet must contain the column {i}, but it only contains {df.columns.tolist()}")
print("✓")




# df_mini is the df, with just two columns for easier handling.
print("Minimizing sample sheet ...                            ", end = "", flush = True)
lag()
df_mini = df[["barcode", "sample_id"]] # Select the only necessary columns
df_mini = df_mini.apply(np.vectorize(lambda x: str(x).strip().replace(" ", ""))) # strip whitespace and replace spaces with underscoresnothing.
df_mini = df_mini.dropna(how='all') # Drop the rows where all elements are missing.
print("✓")


acceptable_barcodes = [f"NB{i:02d}" for i in range(1,25)]
#print(acceptable_barcodes)


print("Checking that the barcodes are correctly formatted ... ", end = "", flush = True)
lag()
for i in df_mini["barcode"]:
    #print("checking", i)
    if not i in acceptable_barcodes: 
        raise Exception(f"The given barcode \'{i}\' is not an acceptable barcode. Here is a list of acceptable barcodes for inspiration:{nl} {' '.join(acceptable_barcodes)}")
print("✓")


print("Checking that the barcodes are unique ...              ", end = "", flush = True)
lag()
if not len(df_mini["barcode"]) == len(set(df_mini["barcode"])):
    counts = pd.DataFrame(df_mini['barcode'].value_counts())
    counts.columns = ["count"]
    counts = counts[counts["count"] > 1]
    #print(nl, counts)
    raise Exception(f"{nl}One or more barcodes are duplicated. Each barcode may only be used once:{nl}{counts}")
print("✓")

print("Marking sample-types following these definitions:")
print("  positive_control: the sample_id must start with \"SEQPOS\" (case insensitive).")
print("  negative_control: the sample_id must end with \"NEG\" (case insensitive).")


#print(df["sample_id"].slower())
df_mini = df_mini.assign(type = ['positive_control' if a.lower().startswith("seqpos") else ('negative_control' if a.lower().endswith("neg") else 'sample') for a in df['sample_id']])

#df_mini = df_mini.assign(type = lambda x: (x*10 if x<2 else (x**2 if x<4 else x+10)

print()
print("These are the samples from the samplesheet you have given:")
print(df_mini.to_string(index = False))
print("//")
print()

if not development_mode:
    check_user("Do you wish to proceed?")






###################
# Validate rundir #
###################

# Wait for the rundir to occur in the specified path. 
# If it doesn't occur after a specified waiting time, then stop the p

if rundir[-1] == "/":
    print("Removing trailing slash from rundir")
    rundir = rundir[0:-1]


print("Checking that the rundir exists ...                    ", end = "", flush = True)
if not os.path.isdir(rundir):
    raise Exception(f"The rundir does not exist.")
print("✓")

print(f"Looking for MinKNOW-characteristic output:") #, end = "", flush = True)

for i in range(200):
    print("  Looking ... ", end = "", flush = True)
    fastq_pass_bases = glob.glob(rundir + "/**/fastq_pass", recursive = True) # Find any occurrence of the wanted path
    if len(fastq_pass_bases) == 0:
        print("nothing found yet, waiting 10 secs ...")
        time.sleep(10) # Wait 10 seconds.
    elif(i == 10):
        print() # clean newline
        raise Exception("nothing found after 10 tries. Aborting.")
    else: 
        print(f"Found                                    ✓")
        break


if not len(fastq_pass_bases) == 1:
    raise Exception(f"There seems to be more than one fastq_pass sub-directory beneath the given rundir. These paths were found:{nl} {str(nl + ' ').join(fastq_pass_bases)}{nl}Please specify a more specific rundir.")


fastq_pass_base = fastq_pass_bases[0]
del fastq_pass_bases
print(f"Found the following fastq_pass base which will be given to rampart: {nl}  {fastq_pass_base}{nl}")





# base_dir is the place where fastq_pass, fast5_pass and the sequencing summary resides.
base_dir = os.path.dirname(fastq_pass_base) # This only works because there is NOT a trailing slash on the fastq_pass_base
print(f"This is the batch base directory:{nl}  {base_dir}")



very_long_batch_id = base_dir.split("/")[-1]
print(f"This is the very long batch id:", very_long_batch_id)

date_parse, time_parse, minion_parse, flowcell_parse, arbhash_parse = very_long_batch_id.split("_")

print("date:    ", date_parse)
print("time:    ", time_parse)
print("minion:  ", minion_parse)
print("flowcell:", flowcell_parse)
print("arbhash: ", arbhash_parse)




batch_id = ".".join(very_long_batch_id.split("_")[0:2]) # The first two words (date, time), joined by a dot.
print(f"This is the parsed batch_id:", batch_id)



out_base = os.path.join(base_dir, "pappenheim_output") # out_base is the directory where the pipeline will write its output to.



# Check that the sequence_summary.txt file exists. If it doesn't, we won't be able to polish the assemblies.
# Edit: since we're using medaka, we don't need this file up front. We can make a rule, that checks that sequencing and basecalling has stopped, and take it from there.
# print("Checking that the sequencing_summary_*.txt-file has been written to disk")
# while True:

#     sequencing_summary_file = glob.glob(base_dir + "/sequencing_summary_*.txt")
#     if len(sequencing_summary_file) == 0:
#         print("  No file yet. Waiting 15 minutes ..")
#         time.sleep(60*15)
#     else:
#         break

#     #raise Exception("sequence_summary.txt does not exist yet. Rerun the pipeline when it has been written to disk.")
# sequencing_summary_file = sequencing_summary_file[0]
# print("  The sequencing summary has been found                ✓")
# print(f"  This is the sequencing_summary_*.txt-file: \"{sequencing_summary_file.split('/')[-1]}\"")



sample_sheet_given_file = f"{fastq_pass_base}/../sample_sheet_given.tsv"
print(f"Backing up the original sample sheet ...               ", end = "", flush = True)
df.to_csv(sample_sheet_given_file, sep = "\t")
print("✓")

print()

disk_barcodes_list  = sorted(glob.glob(fastq_pass_base + "/barcode*")) # Find all fastq_pass/barcode* directories
disk_barcodes_df = pd.DataFrame({'barcode_path': disk_barcodes_list})

#df_mini = df_mini.assign(type = ['positive_control' if a.lower().startswith("seqpos") else ('negative_control' if a.lower().endswith("neg") else 'sample') for a in df['sample_id']])
disk_barcodes_df = disk_barcodes_df.assign(barcode_basename = [i.split("/")[-1] for i in disk_barcodes_df["barcode_path"]])
disk_barcodes_df = disk_barcodes_df.assign(barcode = ["NB" + i[-2:] for i in disk_barcodes_df["barcode_path"]])



print("Continuing with the following barcodes:")

# the workflow_table is the table that contains the records where the barcode could be found on the disk.
workflow_table = disk_barcodes_df.merge(df_mini, how='left', on='barcode') # left join (merge) the present barcodes onto the df_mini table.
workflow_table = workflow_table.dropna(subset = ["sample_id"])

#print(workflow_table[["barcode", "sample_id", "type"]].to_string(index = False))
print(workflow_table)
print("//")
print()





if not development_mode:
    check_user("These are the samples found on disk that match your input. Do you wish to proceed?")




####################################################################
# Finally we can run the artic protocol                            #
# https://artic.network/ncov-2019/ncov2019-bioinformatics-sop.html #
####################################################################




rampart_sub_batch_id = datetime.datetime.now().strftime("%y-%m-%dT%H%M%S")

print("this is the rampart sub batch id", rampart_sub_batch_id)




#This is the collection target, it collects all outputs from other targets. 
rule all:
   input: expand(["{out_base}/flags/{batch_id}_subbatch{rampart_sub_batch_id}_rampart.flag.ok", \
                  "{out_base}/{batch_id}_{sample_id}/consensus/{batch_id}_{sample_id}.consensus.fasta", \
                  "{out_base}/{batch_id}_{sample_id}/pangolin/{batch_id}_{sample_id}.pangolin_long.tsv", \
                  "{out_base}/{batch_id}_{sample_id}/nextclade/{batch_id}_{sample_id}.nextclade_long.tsv", \
                  "{out_base}/collected/{batch_id}_collected_nextclade_long.tsv", \
                  "{out_base}/collected/{batch_id}_collected_input_long.tsv", \
                  "{out_base}/flags/{batch_id}_clean_ready.flag.ok", \
                  "{out_base}/flags/{batch_id}_clean_uploaded.flag.ok", \
                  "{out_base}/flags/{batch_id}_raw_uploaded.flag.ok"], \
                 out_base = out_base, sample_id = workflow_table["sample_id"], batch_id = batch_id, rampart_sub_batch_id = rampart_sub_batch_id)


# # rule all for testing rampart only
# rule all:
#     input: expand("{out_base}/{batch_id}_{sample_id}/read_filtering/{batch_id}_{sample_id}.fastq", \
#                   out_base = out_base, sample_id = workflow_table["sample_id"], batch_id = batch_id)





def exit_rampart(wait_ = 6):
    print("Issuing the following command to exit Rampart:")
    command = f"""
        
        # sleep {wait_}

        # killpid=$(lsof -n -i :3000 | grep -E \"^node\" | awk '{{ print $2 }}' | grep "")
        # echo "this is the killpid $killpid"
        # # Check that any 
        # #if lsof -t -i :3000; then
        # if [ ! -z $killpid ]; then

        #     # filter for only the node-processes writing to port 3000
        #     killpid=$(lsof -n -i :3000 | grep -E \"^node\" | awk '{{ print $2 }}')

        #     #echo $killpid
        #     kill -2 $killpid
        #     echo 'Rampart has been exited'

        # else
        #     echo 'No Rampart job to kill'
        # fi

        sleep {wait_}
        kill -2 $(lsof -t -i :3000)

        """
    print(command)
    os.system(command)

exit_rampart(wait_ = 0)


# What I like about this rule is that it ensures an easy way to install rampart by using the snakemake-conda installer
# the only input needed for rampart, is the base_dir
rule start_rampart:
    output: "{out_base}/flags/{batch_id}_subbatch{rampart_sub_batch_id}_rampart.flag.ok"
    conda: "envs/rampart.yml"
    params: fastq = fastq_pass_base
    shell: """
        
        # Check that rampart actually works
        rampart --version


        # Spawn firefox with lag before calling the rampart program.
        {{ $(sleep 2; firefox localhost:3000)  & }};


        # Call rampart forked. Later this pid will be closed.
        echo "Starting rampart now."
        rampart --protocol artic-ncov2019/rampart/ --clearAnnotated --basecalledPath {params.fastq} &
 
         # Before running rampart we may touch the output such that the pipeline can finish gracefully. Of course, we then have the problem that it wont rerun when the pipeline is started again.
        touch {output}

        

    """




# This target simply waits until the sequencing and basecalling has finished.
rule wait_for_minknow:
    output: flag = "{out_base}/flags/{batch_id}_minknow_done.flag.ok",
        sequencing_summary_moved = "{out_base}/{batch_id}_sequencing_summary.txt"
    run:
        
        # heck that sequencing and basecalling has finished, by checking the existence of the sequence_summary_*.txt-file.

        minutes_wait = 10
        print("Checking that the sequencing_summary_*.txt-file has been written to disk ...")
        while True:
            sequencing_summary_file = glob.glob(base_dir + "/sequencing_summary_*.txt")
            if len(sequencing_summary_file) == 0:
                print(f"  Still sequencing/basecalling. Waiting {minutes_wait} minutes ..")
                time.sleep(60*minutes_wait)
            else:
                break

            #raise Exception("sequence_summary.txt does not exist yet. Rerun the pipeline when it has been written to disk.")
        sequencing_summary_file = sequencing_summary_file[0]
        print("  The sequencing summary has been found                ✓")
        print(f"  This is the sequencing_summary_*.txt-file: \"{sequencing_summary_file.split('/')[-1]}\"")

        os.system(f"mv {sequencing_summary_file} {output.sequencing_summary_moved}")
        os.system(f"touch {output}")



# Read filtering
# Because ARTIC protocol can generate chimeric reads, we perform length filtering.
# This step is performed for each barcode in the run.
# We first collect all the FASTQ files (typically stored in files each containing 4000 reads) into a single file.
# Because we're only using the "pass" reads we can speed up the process with skip-quality-check.
rule read_filtering:
    input:
        minknow_flag = "{out_base}/flags/{batch_id}_minknow_done.flag.ok",
        barcode_dir = directory(lambda wildcards: workflow_table[workflow_table["sample_id"] == wildcards.sample_id]["barcode_path"].values[0])
    output: "{out_base}/{batch_id}_{sample_id}/read_filtering/{batch_id}_{sample_id}.fastq" #_barcode00.fastq
    conda: "artic-ncov2019/environment.yml"
    shell: """


    artic guppyplex --skip-quality-check --min-length 400 --max-length 700 --directory {input.barcode_dir} --output {output}


    """



# Run the MinION pipeline
rule minion:
    input: "{out_base}/{batch_id}_{sample_id}/read_filtering/{batch_id}_{sample_id}.fastq"
    output: "{out_base}/{batch_id}_{sample_id}/consensus/{batch_id}_{sample_id}.consensus.fasta"
    conda: "artic-ncov2019/environment.yml"
    params:
        #workdir = "{out_base}/{batch_id}/consensus/",
        #input = "../read_filtering/{batch_id}_{sample_id}.fastq",
        #output = "../consensus/{batch_id}_{sample_id}.fasta",
        base_dir = base_dir,
        output_dir = "{out_base}/{batch_id}_{sample_id}/consensus/",
        sequencing_summary_file = "{out_base}/{batch_id}_sequencing_summary.txt"
    threads: 4
    shell: """







    # Original nanopolish
    artic minion \
        --normalise 200 \
        --threads 4 \
        --scheme-directory artic-ncov2019/primer_schemes \
        --scheme-version 3 \
        --read-file {input} \
        --fast5-directory {params.base_dir}/fast5_pass \
        --sequencing-summary {params.sequencing_summary_file} \
        nCoV-2019/V3 {wildcards.batch_id}_{wildcards.sample_id} 

    


    # # New medaka
    # artic minion \
    #     --medaka \
    #     --normalise 200 \
    #     --threads 4 \
    #     --scheme-directory artic-ncov2019/primer_schemes \
    #     --scheme-version 3 \
    #     --read-file {input} \
    #     nCoV-2019/V3 {wildcards.batch_id}_{wildcards.sample_id}





    mv {wildcards.batch_id}_{wildcards.sample_id}.* {params.output_dir}

    """





############
# Pangolin #
############

# Make sure that we have the latest pangolin environment.yml for snakemake-conda
rule pangolin_downloader:
    output: "{out_base}/flags/pangolin_downloader.flag.ok"
    shell: """


        cd pangolin
        git submodule update --remote

        touch {output}


        """

# Runs once per batch
rule pangolin_updater: 
    input: "{out_base}/flags/pangolin_downloader.flag.ok"
    output: "{out_base}/flags/pangolin_updater.flag.ok"
    conda: "pangolin/environment.yml"
    shell: """

        cd pangolin

        # Install pangolin
        python setup.py install > latest_pangolin_install_log.stdout 2> latest_pangolin_install_log.stderr

        # Check that the newest dendendencies are installed. 
        pip install git+https://github.com/cov-lineages/pangoLEARN.git --upgrade 
        pip install git+https://github.com/cov-lineages/lineages.git --upgrade 

        # Check that the install worked
        pangolin -v \
        && pangolin -pv \
        && touch "{wildcards.out_base}/flags/pangolin_install.flag.ok"
        


        touch {output}

    """

rule pangolin:
    # input: expand("{out_base}/{batch_id}_{sample_id}/consensus/{batch_id}_{sample_id}.consensus.fasta", batch_id = batch_id, sample_id = workflow_table["sample_id"]) # per batch
    input:
        pangolin_flag = "{out_base}/flags/pangolin_updater.flag.ok",
        consensuses = "{out_base}/{batch_id}_{sample_id}/consensus/{batch_id}_{sample_id}.consensus.fasta" # per sample
    output: "{out_base}/{batch_id}_{sample_id}/pangolin/{batch_id}_{sample_id}.pangolin.csv"
    conda: "pangolin/environment.yml"
    params: 
        out_dir = "{out_base}/{batch_id}_{sample_id}/pangolin"
    shell: """

        pangolin {input.consensuses} \
            --outfile {output}

        # This output should be pivoted.

    """

rule pivot_pangolin:
    input: "{out_base}/{batch_id}_{sample_id}/pangolin/{batch_id}_{sample_id}.pangolin.csv"
    output: "{out_base}/{batch_id}_{sample_id}/pangolin/{batch_id}_{sample_id}.pangolin_long.tsv"
    #conda: "envs/r-tidyverse.yml"
    run:
        wide = pd.read_csv(str(input), dtype = str)
        wide = wide.add_prefix('pa_')

        wide = wide.assign(batch_id = wildcards.batch_id, sample_id = wildcards.sample_id)

        long = pd.melt(wide, id_vars = ["batch_id", "sample_id"])
        long = long.rename(columns = {"batch_id": "#batch_id"})

        #print("This is long before assigning sample_id:", file = sys.stderr)
        #print(long)

        long.to_csv(str(output), index = False, sep = "\t")






#############
# Nextclade #
#############

# Once per batch
rule nextclade_updater:
    output: "{out_base}/flags/nextclade_updater.flag.ok"
    conda: "envs/nodejs.yml"
    shell: """


        # Install or update nextclade to the latest version.
        npm install --global @neherlab/nextclade


        touch {output}

    """

rule nextclade:
    input: 
        nextclade_flag = "{out_base}/flags/nextclade_updater.flag.ok",
        consensus = "{out_base}/{batch_id}_{sample_id}/consensus/{batch_id}_{sample_id}.consensus.fasta" # per sample
    output: "{out_base}/{batch_id}_{sample_id}/nextclade/{batch_id}_{sample_id}.nextclade.tsv"
    conda: "envs/nodejs.yml"
    shell: """

        nextclade.js \
            --input-fasta {input.consensus} \
            --output-tsv {output}

    """

# Pivot the nextclade output to long format.
rule pivot_nextclade:
    input: "{out_base}/{batch_id}_{sample_id}/nextclade/{batch_id}_{sample_id}.nextclade.tsv"
    output: "{out_base}/{batch_id}_{sample_id}/nextclade/{batch_id}_{sample_id}.nextclade_long.tsv"
    run:
        wide = pd.read_csv(str(input), dtype = str, sep = "\t")
        wide = wide.add_prefix('ne_')
        wide = wide.assign(batch_id = wildcards.batch_id, sample_id = wildcards.sample_id)

        long = pd.melt(wide, id_vars = ["batch_id", "sample_id"])
        long = long.rename(columns = {"batch_id": "#batch_id"})

        long.to_csv(str(output), index = False, sep = "\t")







###########################################
# Merge pangolin nextclade and input_data #
########################################### 

rule collect_variant_data:
    input:
        pangolin = expand("{out_base}/{batch_id}_{sample_id}/pangolin/{batch_id}_{sample_id}.pangolin_long.tsv", \
            out_base = out_base, \
            batch_id = batch_id, \
            sample_id = workflow_table["sample_id"]),
        nextclade = expand("{out_base}/{batch_id}_{sample_id}/nextclade/{batch_id}_{sample_id}.nextclade_long.tsv", \
            out_base = out_base, \
            batch_id = batch_id, \
            sample_id = workflow_table["sample_id"])
    output: 
        collected_pangolin = "{out_base}/collected/{batch_id}_collected_pangolin_long.tsv",
        collected_nextclade = "{out_base}/collected/{batch_id}_collected_nextclade_long.tsv"
    shell: """

        echo -e "#batch_id\tsample_id\tvariable\tvalue" > {output.collected_pangolin}
        cat {input.pangolin} | grep -vE "^#" >> {output.collected_pangolin}

        echo -e "#batch_id\tsample_id\tvariable\tvalue" > {output.collected_nextclade}
        cat {input.nextclade} | grep -vE "^#" >> {output.collected_nextclade}

        """

rule merge_variant_data:
    input:
        pangolin ="{out_base}/collected/{batch_id}_collected_pangolin_long.tsv",
        nextclade = "{out_base}/collected/{batch_id}_collected_nextclade_long.tsv"
    output:
        "{out_base}/collected/{batch_id}_collected_input_long.tsv"
    run:


        # df_mini contains type-information for samples which are not yes written to disk.
        merged = df.merge(df_mini, how = 'left', on = ['barcode', 'sample_id'])



        # workflow_table contains the found paths for fastq files.
        merged = df.merge(workflow_table[['barcode', 'sample_id', 'barcode_path', 'barcode_basename', 'type']], how='left', on=['barcode', 'sample_id']) # left join (merge) the present barcodes onto the df_mini table.

        
        # Couple the correct batch_id
        merged = merged.assign(batch_id = batch_id,
            very_long_batch_id = very_long_batch_id)

        merged['upload_consensus_file'] = merged.apply(lambda row: row.batch_id + "_" + row.sample_id + ".consensus.fasta", axis=1)
            


        # This data can now be long-pivoted and dumped besides pangolin and nextclade.
        merged = pd.melt(merged, id_vars = ["batch_id", "sample_id"])
        merged = merged.rename(columns = {"batch_id": "#batch_id"})

        #merged.to_csv("test_merge_out.tsv", index = False, sep = "\t")

        # Finally write the output to disk.
        merged.to_csv(str(output), index = False, sep = "\t")




rule final_merge:
    input: 
        consensuses = expand("{out_base}/{batch_id}_{sample_id}/consensus/{batch_id}_{sample_id}.consensus.fasta", out_base = out_base, batch_id = batch_id, sample_id = workflow_table["sample_id"]),
        collected_input = "{out_base}/collected/{batch_id}_collected_input_long.tsv",
        collected_pangolin = "{out_base}/collected/{batch_id}_collected_pangolin_long.tsv",
        collected_nextclade = "{out_base}/collected/{batch_id}_collected_nextclade_long.tsv"
    output:
        file = "{out_base}/collected/{batch_id}_all.tsv",
        dir = directory("{out_base}/clean_upload_{batch_id}"),
        clean_ready = "{out_base}/flags/{batch_id}_clean_ready.flag.ok"
    shell: """

        # Collect all long metadata files together in the collected-directory.
        echo -e "#batch_id\tsample_id\tvariable\tvalue" > {output.file}
        cat {input.collected_input} {input.collected_pangolin} {input.collected_nextclade} | grep -vE "^#" >> {output.file}


        # Make a copy of clean outputs
        mkdir -p {output.dir}
        rm -rf {output.dir}/* # Delete old content if any.
        

        # Copy consensus-files and metadata to the upload directory
        cp {input.consensuses} {output.dir}
        cp {output.file} {output.dir}
        cp {base_dir}/barcode_alignment*.tsv {output.dir}/{batch_id}_barcode_alignment.tsv
        cp {base_dir}/final_summary_*.txt {output.dir}/{batch_id}_final_summary.txt
        echo "machine_hostname=$(hostname)" >> {output.dir}/{batch_id}_final_summary.txt
        echo "final_merge_date=$(date --iso-8601=s)" >> {output.dir}/{batch_id}_final_summary.txt

        head -n 1 {base_dir}/throughput_*.csv > {output.dir}/{batch_id}_final_throughput.txt
        tail -n 1 {base_dir}/throughput_*.csv >> {output.dir}/{batch_id}_final_throughput.txt


        # Touch a flag to show that the clean files are OK.
        touch {output.clean_ready}


        # The output.dir now contains the essential files to backup and keep for eternity.


        """


# This is a rule that you can costumize however you want to upload your data somewhere
rule custom_upload:
    input: 
        clean_flag = "{out_base}/flags/{batch_id}_clean_ready.flag.ok"
    output:
        clean_upload_flag = "{out_base}/flags/{batch_id}_clean_uploaded.flag.ok",
        raw_upload_flag = "{out_base}/flags/{batch_id}_raw_uploaded.flag.ok"
    params:
        clean_dir = "{out_base}/clean_upload_{batch_id}"
    shell: """


        # Optionally upload the output.dir (clean data)
        touch ~/pappenheim_upload.sh
        bash ~/pappenheim_upload.sh {params.clean_dir} clinmicrocore/BACKUP/nanopore_sarscov2/pappenheim_clean/
        touch {output.clean_upload_flag}


        # Optionally upload the base_dir (raw data)
        touch ~/pappenheim_upload.sh
        bash ~/pappenheim_upload.sh {base_dir} clinmicrocore/BACKUP/nanopore_sarscov2/pappenheim_raw/
        touch {output.raw_upload_flag}


        # If all went well, we touch the final OK flag.

        """



# This function will be registered as a function to run when the pipeline is done.
# It gathers the pid for the node job that serves the rampart web site to port 3000 and closes it.



atexit.register(exit_rampart)

