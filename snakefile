
__author__ = "Tine Sneibjerg Ebsen, Carl Mathias Kobel, Benjamin L. Nichum"
__version__ = "0.2"


# start_pappenheim '/run/user/1000/gvfs/smb-share:server=onerm.dk,share=nfpdata/Afdeling/AUHKLMIK/AUH/Afdelingen/Afsnit Molekyl. og Serologi/NanoPore/NanoPore Metadata/5770_seq_2021-03-23_NEB.xlsx'  '/home/ontseqa/Desktop/sc2_sequencing/COVID19-AUH-20210324-5770/' -np
# start_pappenheim "/home/ontseq4/Desktop/9493_seq_2021-12-02_NEB_rettet_iretar.xlsx" "/media/ontseq4/ssd/sc_sequencing/2021-12-03/"

import sys
import os
from os import listdir
from os.path import isfile, isdir, join, exists, expanduser
import yaml
import pandas as pd
import numpy as np
from pandas_ods_reader import read_ods
from datetime import datetime
import glob
import time
import atexit
import datetime

# new imports
from pathlib import Path




configfile: "config.yaml"
print(config["samplesheet"])

# TODO: These variables should be set from the command line
samplesheet = "../testdata/5351_seq_2021-03-09_NEB.xlsx"
rundir = "../../GenomeDK/clinmicrocore/BACKUP/nanopore_sarscov2/COVID19-AUH-20210316-NEB/rawdata/20210316_1259_MN34697_FAP10653_02939c92/"
rundir = "../testdata"


# Actually given on the command line:
samplesheet = config["samplesheet"]
rundir = config["rundir"]
reference = config["reference"]
regions = config["regions"]
threshold = config["threshold"]
maxDepth = config["maxDepth"]


tab = "\t"
nl = "\n"



# Print a convincing logo, if the window is big enough.
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







# Check that input was given.
if config["samplesheet"] == "NA":
    raise Exception("No samplesheet file was given. Please specify a samplesheet by appending --config rundir=\"path/to/samplesheet/\" to the command line call.")
if config["rundir"] == "NA":
    raise Exception("No rundir path was given. Please specify a rundir by appending --config rundir=\"path/to/rundir/\" to the command line call.")
# TODO: Implement additional input validation, like checking that the objects given are file and dir respectively.



print(f"These are the parameters given:")
print(f"  samplesheet: {samplesheet}")
print(f"  rundir: {rundir}")
print(f"  regions: {regions}")
print(f"  threshold: {threshold}")
print(f"  maxDepth: {maxDepth}")
print()

def read_mail_list(mail_list_file):
    """ Reads a line separated file containing email addresses, 
        and returns them as a space delimited string.
    """
    mail_list = []
    with open(mail_list_file, "r") as mail_list_open:
        for line in mail_list_open:
            if line[0] in ["\n", "#"]: # Skip blank lines and comments.
                continue
            mail_list.append(line.strip())

    return ",".join(mail_list)

mail_list = read_mail_list("mail_list.txt") 
print("mail_list:", mail_list)

#########################
# Parse the samplesheet #
#########################

samplesheet_extension = samplesheet.split(".")[-1]
print(f"Reading .{samplesheet_extension}-type sample sheet \"{samplesheet}\"")
if samplesheet_extension == "ods":

    raise Exception(f"The spreadsheet file extension {samplesheet_extension} is not yet implemented.")




elif samplesheet_extension == "xlsx":
    # Uses openpyxl
    df = pd.read_excel(samplesheet, dtype = str)

elif samplesheet_extension == "xls":
    df = pd.read_excel(samplesheet)

else:
    raise Exception(f"The spreadsheet file extension {samplesheet_extension} is not yet implemented.")

# Clean up the spreadsheet
print("Cleaning sample sheet ...                              ", end = "", flush = True)
df.columns = map(str.lower, df.columns) # Lowercase
df.columns = map(str.strip, df.columns) # Remove edge-spaces
df.columns = map(lambda x: str(x).replace(" ", "_"), df.columns) # Replace spaces with underscore
df["barcode"] = df["barcode"].apply(np.vectorize(lambda x: str(x).strip().replace(" ", ""))) # Because we are later going to join using this column, it is necessary to strip it for spaces.
df = df.dropna(subset = ["sample_id"])# remove rows not containing a sample ID
print("✓")

#print("ee", "barcode" in list(df.columns))
# Check that the spreadsheet complies
print("Checking that the necessary columns exist ...          ", end = "", flush = True)
for i in ["barcode", "sample_id"]:
    if not i in df.columns:
        raise Exception(f"The sample sheet is missing a necessary column. The sample sheet must contain the column {i}, but it only contains {df.columns.tolist()}")
print("✓")




# df_mini is the df, with just two columns for easier handling.
print("Minimizing sample sheet ...                            ", end = "", flush = True)
df_mini = df[["barcode", "sample_id"]] # Select the only necessary columns
df_mini = df_mini.apply(np.vectorize(lambda x: str(x).strip().replace(" ", ""))) # strip whitespace and replace spaces with underscoresnothing.
df_mini = df_mini.dropna(how='all') # Drop the rows where all elements are missing.
print("✓")


acceptable_barcodes = [f"NB{i:02d}" for i in range(1,97)]
#print(acceptable_barcodes)


print("Checking that the barcodes are correctly formatted ... ", end = "", flush = True)
for i in df_mini["barcode"]:
    #print("checking", i)
    if not i in acceptable_barcodes: 
        raise Exception(f"The given barcode \'{i}\' is not an acceptable barcode. Here is a list of acceptable barcodes for inspiration:{nl} {' '.join(acceptable_barcodes)}")
print("✓")


print("Checking that the barcodes are unique ...              ", end = "", flush = True)
if not len(df_mini["barcode"]) == len(set(df_mini["barcode"])):
    bc_counts = pd.DataFrame(df_mini['barcode'].value_counts())
    bc_counts.columns = ["count"]
    bc_counts = bc_counts[bc_counts["count"] > 1]
    #print(nl, bc_counts)
    raise Exception(f"{nl}One or more barcodes are duplicated. Each barcode may only be used once:{nl}{bc_counts}")
print("✓")



print("Checking that the sample id's are unique ...           ", end = "", flush = True)
if not len(df_mini["sample_id"]) == len(set(df_mini["sample_id"])):
    sid_counts = pd.DataFrame(df_mini['sample_id'].value_counts())
    sid_counts.columns = ["count"]
    sid_counts = sid_counts[sid_counts["count"] > 1]
    #print(nl, sid_counts)
    raise Exception(f"{nl}One or more sample_id's are duplicated. Each sample_id may only be used once:{nl}{sid_counts}")
print("✓")

print("Checking that the sample id's are okay ...           ", end = "", flush = True)
for sample_id in df["sample_id"]:
    if "_" in sample_id:
        raise Exception(f"{nl}One or more sample_id's contain underscores (_). Sample_ids should consist only of numbers and letters. If a seperator is required, use -")
    if " " in sample_id:
        raise Exception(f"{nl}One or more sample_id's contain spaces. Sample_ids should consist only of numbers and letters. If a seperator is required, use -")
print("✓")



print("Marking sample-types following these definitions:")
print("  positive_control: the sample_id must start with \"SEQPOS\" (case insensitive).")
print("  negative_control: the sample_id must end with \"NEG\" (case insensitive).")


#print(df["sample_id"].slower())
df_mini = df_mini.assign(type = ['positive_control' if a.lower().startswith("seqpos") else ('negative_control' if a.lower().endswith("neg") else 'sample') for a in df['sample_id']])

#df_mini = df_mini.assign(type = lambda x: (x*10 if x<2 else (x**2 if x<4 else x+10)

print()
print("These are the samples from the samplesheet you have given:")
print(df_mini.to_string())
print("//")
print()





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


if config["run_monitoring"] and exists(expanduser("~/pappenheim/CoverMon.flag")):
    # Now we have all resources to start monitoring in the background
    os.system("rm ~/pappenheim/CoverMon.flag") 
    print(f"Starting monitoring in the background ... ")
    cov_mon_sh = open("start_mon.sh", "w")
    command = f"#!/bin/bash{nl}source ~/miniconda3/etc/profile.d/conda.sh{nl}cd ~/pappenheim/CoverMon{nl}conda activate covermon {nl}python seq_mon.py '{samplesheet}' {rundir} {reference} {threshold} {maxDepth} {regions}"
    cov_mon_sh.write(command)
    cov_mon_sh.close()
    # Start the sequence monitoring in a new terminal
    os.system("gnome-terminal --tab -- bash start_mon.sh")


# And here is the code from the rule wait_for_minknow
minutes_wait = 10
print("Checking that the sequencing_summary_*.txt-file has been written to disk ...")
while True:
    sequencing_summary_file = glob.glob(base_dir + "/sequencing_summary_*.txt")
    if len(sequencing_summary_file) == 0:
        print(f"  Still sequencing/basecalling; waiting {minutes_wait} minutes ...")
        time.sleep(60*minutes_wait)
    else:
        break



sequencing_summary_file = sequencing_summary_file[0]
print("  The sequencing summary has been found                ✓")
#print(f"  This is the sequencing_summary_*.txt-file: \"{sequencing_summary_file.split('/')[-1]}\"")
print(f"  This is the sequencing_summary_*.txt-file (full): \"{sequencing_summary_file}\"")








sample_sheet_given_file = f"{fastq_pass_base}/../sample_sheet_given.tsv"
print(f"Backing up the original sample sheet ...               ", end = "", flush = True)
df.to_csv(sample_sheet_given_file, sep = "\t")
print("✓")

print()

disk_barcodes_list  = sorted(glob.glob(fastq_pass_base + "/barcode*")) # Find all fastq_pass/barcode* directories
disk_barcodes_df = pd.DataFrame({'barcode_path': disk_barcodes_list})

#df_mini = df_mini.assign(type = ['positive_control' if a.lower().startswith("seqpos") else ('negative_control' if a.lower().endswith("neg") else 'sample') for a in df['sample_id']])
disk_barcodes_df = disk_barcodes_df.assign(barcode_basename = [i.split("/")[-1] for i in disk_barcodes_df["barcode_path"]])
if "RB" in df["barcode"][0]:
    disk_barcodes_df = disk_barcodes_df.assign(barcode = ["RB" + i[-2:] for i in disk_barcodes_df["barcode_path"]])
elif "NB" in df["barcode"][0]:
    disk_barcodes_df = disk_barcodes_df.assign(barcode = ["NB" + i[-2:] for i in disk_barcodes_df["barcode_path"]])
else:
    raise Exception(f"Barcodes in samplesheet are not acceptable")











print("Continuing with the following barcodes:")

# the workflow_table is the table that contains the records where the barcode could be found on the disk.
workflow_table = disk_barcodes_df.merge(df_mini, how='left', on='barcode') # left join (merge) the present barcodes onto the df_mini table.
workflow_table = workflow_table.dropna(subset = ["sample_id"])

#print(workflow_table[["barcode", "sample_id", "type"]].to_string(index = False))
print(workflow_table)
print("//")
print()




####################################################################
# Finally we can run the artic protocol                            #
# https://artic.network/ncov-2019/ncov2019-bioinformatics-sop.html #
####################################################################





#                   "{out_base}/{batch_id}/{batch_id}_metadata_init.tsv"], \
#                   "{out_base}/collected/{batch_id}_df.csv", \

#This is the collection target, it collects all outputs from other targets. 
rule all:
   input: expand(["{out_base}/{batch_id}_{sample_id}/consensus/{batch_id}_{sample_id}.consensus.fasta", \
                  "{out_base}/{batch_id}_{sample_id}/consensus/{batch_id}_{sample_id}_depths.tsv", \
                  "{out_base}/{batch_id}_{sample_id}/pangolin/{batch_id}_{sample_id}.pangolin_long.tsv", \
                  "{out_base}/{batch_id}_{sample_id}/nextclade/{batch_id}_{sample_id}.nextclade_long.tsv", \
                  "{out_base}/collected/{batch_id}_collected_nextclade_long.tsv", \
                  "{out_base}/collected/{batch_id}_collected_input_long.tsv", \
                  "{out_base}/flags/{batch_id}_clean_ready.flag.ok", \
                  "{out_base}/collected/{batch_id}_all.tsv", \
                  "{out_base}/flags/{batch_id}_clean_uploaded.flag.ok", \
                  "{out_base}/flags/{batch_id}_raw_uploaded.flag.ok"], \
                 out_base = out_base, sample_id = workflow_table["sample_id"], batch_id = batch_id)
                 
#                  "{out_base}/flags/{batch_id}_output_mail_sent.flag", \

# Read filtering
# Because ARTIC protocol can generate chimeric reads, we perform length filtering.
# This step is performed for each barcode in the run.
# We first collect all the FASTQ files (typically stored in files each containing 4000 reads) into a single file.
# Because we're only using the "pass" reads we can speed up the process with skip-quality-check.
rule read_filtering:
    input:
        #minknow_flag = "{out_base}/flags/{batch_id}_minknow_done.flag.ok",
        barcode_dir = directory(lambda wildcards: workflow_table[workflow_table["sample_id"] == wildcards.sample_id]["barcode_path"].values[0])
    output: "{out_base}/{batch_id}_{sample_id}/read_filtering/{batch_id}_{sample_id}.fastq" #_barcode00.fastq
    conda: "artic-ncov2019/environment.yml"
    shell: """


    artic guppyplex --skip-quality-check --min-length 100 --directory {input.barcode_dir} --output {output}


    """



# Run the MinION pipeline
rule minion:
    input:
        fastq ="{out_base}/{batch_id}_{sample_id}/read_filtering/{batch_id}_{sample_id}.fastq",

    output:
        consensus = "{out_base}/{batch_id}_{sample_id}/consensus/{batch_id}_{sample_id}.consensus.fasta",
        depths = ["{out_base}/{batch_id}_{sample_id}/consensus/{batch_id}_{sample_id}.coverage_mask.txt.1.depths",
                  "{out_base}/{batch_id}_{sample_id}/consensus/{batch_id}_{sample_id}.coverage_mask.txt.2.depths"]

    conda: "artic-ncov2019/environment.yml"
    params:
        base_dir = base_dir,
        output_dir = "{out_base}/{batch_id}_{sample_id}/consensus/",
        sequencing_summary_file = sequencing_summary_file
    threads: 4
    shell: """

    # Check that artic minion can run, so we are sure that we are not creating a blank output due to dependency errors.
    artic minion -h > /dev/null 2> /dev/null && echo "artic minion in itself runs fine."


    # Original nanopolish
    artic minion \
        --normalise 200 \
        --threads 4 \
        --scheme-directory primer_schemes \
        --scheme-version 1 \
        --read-file {input.fastq} \
        --fast5-directory {params.base_dir}/fast5_pass \
        --sequencing-summary {params.sequencing_summary_file} \
        ONT-Midnight-nCoV-2019 {wildcards.batch_id}_{wildcards.sample_id} \
        || echo ">{wildcards.batch_id}_{wildcards.sample_id}_notenoughdata" > {output.consensus} \
            && cat scripts/29903N.txt >> {output.consensus} \
            && touch {wildcards.batch_id}_{wildcards.sample_id}.coverage_mask.txt.1.depths \
            && touch {wildcards.batch_id}_{wildcards.sample_id}.coverage_mask.txt.2.depths

    # If there is not enough data, the job should exit gracefully and create a blank output file with 29903 N's

    # I have considered that it should only or-exit gracefully when the sample type is "positive_control" or "negative_control". But I think it is very possible that normal samples can also fail, which should not halt the complete batch in terms of the pipeline.




    mv {wildcards.batch_id}_{wildcards.sample_id}.* {params.output_dir}

    """



rule depth:
    input:
        depths = ["{out_base}/{batch_id}_{sample_id}/consensus/{batch_id}_{sample_id}.coverage_mask.txt.1.depths",
                  "{out_base}/{batch_id}_{sample_id}/consensus/{batch_id}_{sample_id}.coverage_mask.txt.2.depths"]
    output: "{out_base}/{batch_id}_{sample_id}/consensus/{batch_id}_{sample_id}_depths.tsv"
    shell: """

        # hattespase

        cat {input.depths} \
        | awk -v sample={wildcards.batch_id}_{wildcards.sample_id} '{{ print sample "\\t" $0 }}' \
        > {output}

    """




############
# Pangolin #
############

# Make sure that we have the latest pangolin environment.yml for snakemake-conda
rule pangolin_downloader:
    output: "{out_base}/flags/pangolin_downloader.flag.ok"
    shell: """

        cd pangolin
        git pull origin master

        touch {output}

        """
        

# Runs once per batch
# It seems like even though we have a new environment.yml-file for pangolin, the environment is not updated?!
# TODO: Do the git pull command (rule pangolin_downloader) before computing the job dag.
rule pangolin_updater: 
    input: "{out_base}/flags/pangolin_downloader.flag.ok"
    output: pangolin_updater_flag = "{out_base}/flags/pangolin_updater.flag.ok",
            alias_key = "{out_base}/alias_key.json"
    conda: "pangolin/environment.yml"
    shell: """

        cd pangolin

        # Update dependencies
        pip install .


        pangolin --update

        
        # Check that the install worked
        pangolin -v 
        pangolin -pv
        pangolin --alias > {output.alias_key} 
        


        touch {output.pangolin_updater_flag}

    """


rule pangolin:
    # In pangolin v4 usher is now the default method
    # input: expand("{out_base}/{batch_id}_{sample_id}/consensus/{batch_id}_{sample_id}.consensus.fasta", batch_id = batch_id, sample_id = workflow_table["sample_id"]) # per batch
    input:
        pangolin_flag = "{out_base}/flags/pangolin_updater.flag.ok",
        consensuses = "{out_base}/{batch_id}_{sample_id}/consensus/{batch_id}_{sample_id}.consensus.fasta" # per sample
    output: "{out_base}/{batch_id}_{sample_id}/pangolin/{batch_id}_{sample_id}.pangolin.csv"
    conda: "pangolin/environment.yml"
    params: 
        out_dir = "{out_base}/{batch_id}_{sample_id}/pangolin"
    shell: """

        pangolin {input.consensuses} --max-ambig 0.1672 --outfile {output}


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
    shell: """


        # Install or update nextclade to the latest version.
        cd /home/nextclade/
        wget https://github.com/nextstrain/nextclade/releases/latest/download/nextclade-x86_64-unknown-linux-gnu
        mv nextclade-x86_64-unknown-linux-gnu nextclade
        chmod +x "/home/nextclade/nextclade" 
        cd ~/pappenheim/
        nextclade dataset get --name='sars-cov-2' --output-dir='{out_base}/nextclade_files'
        touch {output}

    """

rule nextclade:
    input: 
        nextclade_flag = "{out_base}/flags/nextclade_updater.flag.ok",
        consensus = "{out_base}/{batch_id}_{sample_id}/consensus/{batch_id}_{sample_id}.consensus.fasta" # per sample
    output: "{out_base}/{batch_id}_{sample_id}/nextclade/{batch_id}_{sample_id}.nextclade.tsv"
    shell: """
    nextclade run \
    	-D  {out_base}/nextclade_files \
    	-O {out_base}/{batch_id}_{wildcards.sample_id}/nextclade/ \
    	--output-tsv {output} \
        {input.consensus} 
         
        

    if wc -l {output} < 2
    then
	nextclade run \
    	-D  {out_base}/nextclade_files \
    	-O {out_base}/{batch_id}_{wildcards.sample_id}/nextclade/ \
    	--output-tsv {output} \
        {input.consensus}
    fi

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
        depths = expand("{out_base}/{batch_id}_{sample_id}/consensus/{batch_id}_{sample_id}_depths.tsv", out_base = out_base, batch_id = batch_id, sample_id = workflow_table["sample_id"]),
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
        cat {input.depths} > {output.dir}/depths.tsv
        cp {output.file} {output.dir}
        cp {base_dir}/barcode_alignment*.tsv {output.dir}/{batch_id}_barcode_alignment.tsv
        cp {base_dir}/barcode_alignment*.tsv {out_base}/collected/{batch_id}_barcode_alignment.tsv
        cp {base_dir}/final_summary_*.txt {output.dir}/{batch_id}_final_summary.txt
        echo "machine_hostname=$(hostname)" >> {output.dir}/{batch_id}_final_summary.txt
        echo "final_merge_date=$(date --iso-8601=s)" >> {output.dir}/{batch_id}_final_summary.txt

        head -n 1 {base_dir}/throughput_*.csv > {output.dir}/{batch_id}_final_throughput.txt
        tail -n 1 {base_dir}/throughput_*.csv >> {output.dir}/{batch_id}_final_throughput.txt


        # Touch a flag to show that the clean files are OK.
        touch {output.clean_ready}


        # The output.dir now contains the essential files to backup and keep for eternity.


        """

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

        # and we can gracefully stop monitoring
        sh ~/pappenheim/scripts/stop_monitoring.sh

        # Consider also, to remove the files locally.
        # This is of course, pretty dangerous.

        """

