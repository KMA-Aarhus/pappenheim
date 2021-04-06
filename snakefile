
__author__ = "Carl Mathias Kobel"
__version__ = "0.1"
# snakemake --samplesheet path/to/samplesheet.ods --rundir path/to/minknowoutput/

import sys
from os import listdir
from os.path import isfile, isdir, join
import yaml
import pandas as pd
import numpy as np
from pandas_ods_reader import read_ods
from datetime import datetime
import glob
import time

terminal_rows, terminal_columns = os.popen('stty size', 'r').read().split()


if int(terminal_columns) < 125:
    print("Warning: Please increase your terminal window width.")


#import time
#import re
#from shutil import copyfile
#import re

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



#batch_id = datetime.now().strftime('pap%Y%m%dT%H%M')
batch_id = datetime.now().strftime('pap%Y%m%d')

#batch_id = "today"




# This will make your code slower
def lag(time_ = 0.06):
    if development_mode:
        pass
    else:
        time.sleep(time_)


# Print a convincing logo
print()
print(f"          Pappenheim pipeline v{__version__}  -  Aarhus Universityhospital  -  Department of Clinical Microbiology          ")
#print()
#print("    ▄███████▄    ▄████████    ▄███████▄    ▄███████▄    ▄████████ ███▄▄▄▄      ▄█    █▄       ▄████████  ▄█    ▄▄▄▄███▄▄▄▄   ")
#print("   ███    ███   ███    ███   ███    ███   ███    ███   ███    ███ ███▀▀▀██▄   ███    ███     ███    ███ ███  ▄██▀▀▀███▀▀▀██▄ ")
#print("   ███    ███   ███    ███   ███    ███   ███    ███   ███    █▀  ███   ███   ███    ███     ███    █▀  ███▌ ███   ███   ███ ")
#print("   ███    ███   ███    ███   ███    ███   ███    ███  ▄███▄▄▄     ███   ███  ▄███▄▄▄▄███▄▄  ▄███▄▄▄     ███▌ ███   ███   ███ ")
#print(" ▀█████████▀  ▀███████████ ▀█████████▀  ▀█████████▀  ▀▀███▀▀▀     ███   ███ ▀▀███▀▀▀▀███▀  ▀▀███▀▀▀     ███▌ ███   ███   ███ ")
#print("   ███          ███    ███   ███          ███          ███    █▄  ███   ███   ███    ███     ███    █▄  ███  ███   ███   ███ ")
#print("   ███          ███    ███   ███          ███          ███    ███ ███   ███   ███    ███     ███    ███ ███  ███   ███   ███ ")
#print("  ▄████▀        ███    █▀   ▄████▀       ▄████▀        ██████████  ▀█   █▀    ███    █▀      ██████████ █▀    ▀█   ███   █▀  ")
print()                                                                                                                            
print(f"                                      Press ctrl+c at any time to stop this pipeline.")
print()



# Check that input was given.
if config["samplesheet"] == "NA":
    raise Exception("No samplesheet file was given. Please specify a samplesheet by appending --config rundir=\"path/to/samplesheet/\" to the command line call.")
if config["rundir"] == "NA":
    raise Exception("No rundir path was given. Please specify a rundir by appending --config rundir=\"path/to/rundir/\" to the command line call.")


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




##########################
# Parse the sample sheet #
##########################



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
    raise Exception(f"There seems to be more than one fastq_pass sub-directory beneath the given rundir. These paths were found:{nl} {str(nl + ' ').join(fastq_pass_bases)}")

fastq_pass_base = fastq_pass_bases[0]
del fastq_pass_bases
#print(f"Found the following fastq_pass base: {nl}  {fastq_pass_base}{nl}  This will be regarded as the input_base directory from now on.")


# base_dir is the place where fastq_pass, fast5_pass and the sequencing summary resides.
base_dir = os.path.dirname(fastq_pass_base) # This only works because there is NOT a trailing slash on the fastq_pass_base
print(f"This is the batch base directory:{nl}  {base_dir}")



out_base = os.path.join(base_dir, "pappenheim_output") # out_base is the directory where the pipeline will write its output to.



# Check that the sequence_summary.txt file exists. If it doesn't, we won't be able to polish the assemblies.
#if not os.path.isfile(fastq_pass_base + "/../sequence_"):
print("Checking that the sequencing_summary_*.txt-file has been written to disk ... ", end = "", flush = True)
sequencing_summary_file = glob.glob(base_dir + "/sequencing_summary_*.txt")
if len(sequencing_summary_file) == 0:
    raise Exception("sequence_summary.txt does not exist yet. Rerun the pipeline when it has been written to disk.")
sequencing_summary_file = sequencing_summary_file[0]
print("✓")
print(f"  This is the sequencing_summary_*.txt-file: \"{sequencing_summary_file.split('/')[-1]}\"")




sample_sheet_given_file = f"{fastq_pass_base}/../sample_sheet_given.tsv"
print(f"Backing up the original sample sheet                   ", end = "", flush = True)
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







# This is the collection target, it collects all outputs from other targets. 
rule all:
    input: expand(["{out_base}/{batch_id}_{sample_id}/consensus/{batch_id}_{sample_id}.consensus.fasta", \
                   "{out_base}/{batch_id}_{sample_id}/pangolin/{batch_id}_{sample_id}.pangolin_long.tsv", \
                   "{out_base}/{batch_id}_{sample_id}/nextclade/{batch_id}_{sample_id}.nextclade_long.tsv", \
                   "{out_base}/collected/{batch_id}_collected_nextclade_long.tsv", \
                   "{out_base}/collected/{batch_id}_collected_input_long.tsv", \
                   "{out_base}/upload/{batch_id}_metadata.tsv"], \
                  out_base = out_base, sample_id = workflow_table["sample_id"], batch_id = batch_id)
                   #"{out_base}/{batch_id}/consensus/{batch_id}_{sample_id}.fasta"],
                  






# Read filtering
# Because ARTIC protocol can generate chimeric reads, we perform length filtering.
# This step is performed for each barcode in the run.
# We first collect all the FASTQ files (typically stored in files each containing 4000 reads) into a single file.
# Because we're only using the "pass" reads we can speed up the process with skip-quality-check.
rule read_filtering:
    input: directory(lambda wildcards: workflow_table[workflow_table["sample_id"] == wildcards.sample_id]["barcode_path"].values[0])
    output: "{out_base}/{batch_id}_{sample_id}/read_filtering/{batch_id}_{sample_id}.fastq" #_barcode00.fastq
    conda: "artic-ncov2019/environment.yml"
    shell: """


    artic guppyplex --skip-quality-check --min-length 400 --max-length 700 --directory {input} --output {output}


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
        sequencing_summary_file = sequencing_summary_file
    threads: 4
    shell: """




    artic minion \
        --normalise 200 \
        --threads 4 \
        --scheme-directory artic-ncov2019/primer_schemes \
        --read-file {input} \
        --fast5-directory {params.base_dir}/fast5_pass \
        --sequencing-summary {params.sequencing_summary_file} \
        nCoV-2019/V3 {wildcards.batch_id}_{wildcards.sample_id} 

    mv {wildcards.batch_id}_{wildcards.sample_id}.* {params.output_dir}

    """





############
# Pangolin #
############


# Make sure that we have the latest pangolin environment.yml for snakemake-conda
rule pangolin_downloader:
    output: "{out_base}/flags/pangolin_downloader.flag.ok"
    shell: """


        ping github.com -c 1 && echo "Internet OK" || echo "Warning: no internet connection"

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
        python setup.py install > plog.out 2> plog.err

        # Check that the newest dendendencies are installed. 
        pip install git+https://github.com/cov-lineages/pangoLEARN.git --upgrade 
        pip install git+https://github.com/cov-lineages/lineages.git --upgrade 

        # Check that the install worked
        pangolin -v \
        && pangolin -pv \
        && touch "{wildcards.out_base}/flags/pangolin_install.flag.ok"
        



        # Download the latest version so it can be installed next time around 
        git submodule update --remote


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

        ping github.com -c 1 && echo "Internet OK" || echo "Warning: no internet connection"


        npm install --global @neherlab/nextclade



        touch {output}

    """


rule nextclade:
    input: 
        nextclade_flag = "{out_base}/flags/nextclade_updater.flag.ok",
        consensuses = "{out_base}/{batch_id}_{sample_id}/consensus/{batch_id}_{sample_id}.consensus.fasta" # per sample
    output: "{out_base}/{batch_id}_{sample_id}/nextclade/{batch_id}_{sample_id}.nextclade.tsv"
    conda: "envs/nodejs.yml"
    shell: """


        nextclade.js \
            --input-fasta {input.consensuses} \
            --output-tsv {output}

    """


rule pivot_nextclade:
    input: "{out_base}/{batch_id}_{sample_id}/nextclade/{batch_id}_{sample_id}.nextclade.tsv"
    output: "{out_base}/{batch_id}_{sample_id}/nextclade/{batch_id}_{sample_id}.nextclade_long.tsv"
    run:
        wide = pd.read_csv(str(input), dtype = str, sep = "\t")
        wide = wide.add_prefix('ne_')
        wide = wide.assign(batch_id = wildcards.batch_id, sample_id = wildcards.sample_id)

        long = pd.melt(wide, id_vars = ["batch_id", "sample_id"])
        long = long.rename(columns = {"batch_id": "#batch_id"})

        #print("This is long before assigning sample_id:", file = sys.stderr)
        #print(long)

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
            sample_id = workflow_table["sample_id"]),


    output: 
        collected_pangolin = "{out_base}/collected/{batch_id}_collected_pangolin_long.tsv",
        collected_nextclade = "{out_base}/collected/{batch_id}_collected_nextclade_long.tsv"
    shell:
        """

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
        merged = df.merge(workflow_table[['barcode', 'sample_id', 'barcode_path']], how='left', on=['barcode', 'sample_id']) # left join (merge) the present barcodes onto the df_mini table.

        
        # Couple the correct batch_id
        merged = merged.assign(batch_id = batch_id)

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
        file = "{out_base}/upload/{batch_id}_metadata.tsv",
        
    shell: """

        echo -e "#batch_id\tsample_id\tvariable\tvalue" > {output.file}
        cat {input.collected_input} {input.collected_pangolin} {input.collected_nextclade} | grep -vE "^#" >> {output.file}
        


        # Prepare data for upload
        cp {input.consensuses} {out_base}/upload/

        # tar something


        """



        # exit()

        # # Read and pivot pangolin
        # pangolin = pd.read_csv(str(input.pangolin), sep = "\t")
        # pangolin = pangolin.pivot_table(index = ["batch_id", "sample_id"], 
        #     columns = 'variable', 
        #     values = 'value',
        #     aggfunc = 'first').reset_index() # There should only be one value per index-column combination, therefore first is fine to use as aggfunc.

        # print("This is pangolin before joining: (hidden index)")
        # print(pangolin.to_string(index = False))


        # # Read and pivot nextclade
        # nextclade = pd.read_csv(str(input.nextclade), sep = "\t")
        # nextclade = nextclade.pivot_table(index = ["batch_id", "sample_id"], 
        #     columns = 'variable', 
        #     values = 'value',
        #     aggfunc = 'first').reset_index() # There should only be one value per index-column combination, therefore first is fine to use as aggfunc.

        # print("This is nextclade before joining: (hidden index)")
        # print(nextclade.to_string(index = False))






        # # We need the batch_id so we can make sure that only the correct batched samples are integrated.

        # merged = merged.merge(pangolin, how = 'left', on = ['batch_id', 'sample_id'])
        # merged = merged.merge(nextclade, how = 'left', on = ['batch_id', 'sample_id'])
        

        # print("what")
        # # Left join nextclade onto above df.


        # # Save the df to disk and be done
        # print(merged.columns)
        # print(merged, file = sys.stderr)

        # merged.to_csv("test_merge_out.tsv", sep = "\t")





###############
# Backup data #
###############

#These jobs are optional if you want to backup your data.

#rule upload_core_data:

#rule upload_raw_data:






# # Collect all targets
# rule all:
#     input: expand(["{out_base}/metadata.tsv", \
#                    "{out_base}/samples/{sample}/{sample}.fa", \
#                    "{out_base}/samples/{sample}/prokka/{sample}.gff", \
#                    "{out_base}/roary/summary_statistics.txt", \
#                    "{out_base}/abricate/card_detail.tsv", \
#                    "{out_base}/mlst/mlst.tsv", \
#                    "{out_base}/mashtree/mashtree.newick", \
#                    "{out_base}/fasttree/fasttree.newick"], \
#                   out_base = out_base_var, sample = df["sample"]) # copy





  

# # Write the df table to the directory for later reference.
# rule metadata:
#     input: df["input_file"].tolist()
#     output: "{out_base}/metadata.tsv"
#     run: 
#         df.to_csv(str(output), index_label = "index", sep = "\t")



# # Copy the input file to its new home
# rule copy:
#     #input: "{sample}"
#     input: lambda wildcards: df[df["sample"]==wildcards.sample]["input_file"].values[0]
#     output: "{out_base}/samples/{sample}/{sample}.fa"
#     #log: "logs/{out_base}_{wildcards.sample}.out.log"

#     shell: """

#         mkdir -p logs output_asscom1

#         cp {input} {output}

#         """


