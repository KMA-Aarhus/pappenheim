
__author__ = "Carl Mathias Kobel"
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





#import time
#import re
#from shutil import copyfile
#import re

monkey_mode = False

# TODO: These variables should be set from the command line
samplesheet = "../testdata/5351_seq_2021-03-09_NEB.xlsx"
rundir = "../../GenomeDK/clinmicrocore/BACKUP/nanopore_sarscov2/COVID19-AUH-20210316-NEB/rawdata/20210316_1259_MN34697_FAP10653_02939c92/"
rundir = "../testdata"

# TODO: These variables should be given in the config






tab = "\t"
nl = "\n"
batch_date_identifier = datetime.now().strftime('pap%Y%m%dT%H%M')



# This will make your code slower
def lag(time_ = 0.06):
    time.sleep(time_)

# Print a convincing logo
print()
print("Pappenheim pipeline")
print(" ", batch_date_identifier)
print()
print("    ▄███████▄    ▄████████    ▄███████▄    ▄███████▄    ▄████████ ███▄▄▄▄      ▄█    █▄       ▄████████  ▄█    ▄▄▄▄███▄▄▄▄   ")
print("   ███    ███   ███    ███   ███    ███   ███    ███   ███    ███ ███▀▀▀██▄   ███    ███     ███    ███ ███  ▄██▀▀▀███▀▀▀██▄ ")
print("   ███    ███   ███    ███   ███    ███   ███    ███   ███    █▀  ███   ███   ███    ███     ███    █▀  ███▌ ███   ███   ███ ")
print("   ███    ███   ███    ███   ███    ███   ███    ███  ▄███▄▄▄     ███   ███  ▄███▄▄▄▄███▄▄  ▄███▄▄▄     ███▌ ███   ███   ███ ")
print(" ▀█████████▀  ▀███████████ ▀█████████▀  ▀█████████▀  ▀▀███▀▀▀     ███   ███ ▀▀███▀▀▀▀███▀  ▀▀███▀▀▀     ███▌ ███   ███   ███ ")
print("   ███          ███    ███   ███          ███          ███    █▄  ███   ███   ███    ███     ███    █▄  ███  ███   ███   ███ ")
print("   ███          ███    ███   ███          ███          ███    ███ ███   ███   ███    ███     ███    ███ ███  ███   ███   ███ ")
print("  ▄████▀        ███    █▀   ▄████▀       ▄████▀        ██████████  ▀█   █▀    ███    █▀      ██████████ █▀    ▀█   ███   █▀  ")
print()                                                                                                                            
print(f" These are the parameters given:")
print(f"   samplesheet: {samplesheet}")
print(f"        rundir: {rundir}")
print()


def check_user(prompt):
    input_OK = input(f"{prompt} [y/n]: ")
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
    #raise Exception(f"The spreadsheet file extension {samplesheet_extension} is not yet implemented.")
    # ValueError in line 61 of /Users/carkob/repos/pappenheim/Snakefile:
    # Your version of xlrd is 2.0.1. In xlrd >= 2.0, only the xls format is supported. Install openpyxl instead.
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
df = df.dropna(subset = ["sample_id"])# remove rows not containing a sample ID
print("OK")

#print("ee", "barcode" in list(df.columns))
# Check that the spreadsheet complies
print("Checking that the necessary columns exist ...          ", end = "", flush = True)
lag()
for i in ["barcode", "sample_id"]:
    if not i in df.columns:
        raise Exception(f"The sample sheet is missing a necessary column. The sample sheet must contain the column {i}, but it only contains {df.columns.tolist()}")
print("OK")





print("Minimizing sample sheet ...                            ", end = "", flush = True)
lag()
df_mini = df[["barcode", "sample_id"]] # Select the only necessary columns
df_mini = df_mini.dropna(how='all') # Drop the rows where all elements are missing.
df_mini = df_mini.apply(np.vectorize(lambda x: str(x).strip().replace(" ", "_"))) # strip whitespace and replace spaces with underscores.
print("OK")


acceptable_barcodes = [f"NB{i:02d}" for i in range(1,25)]
#print(acceptable_barcodes)
nl, tab = "\n", "\t"

print("Checking that the barcodes are correctly formatted ... ", end = "", flush = True)
lag()
for i in df_mini["barcode"]:
    #print("checking", i)
    if not i in acceptable_barcodes: 
        raise Exception(f"The given barcode \'{i}\' is not an acceptable barcode. Here is a list of acceptable barcodes for inspiration: {' '.join(acceptable_barcodes)}")
print("OK")


print("Checking that the barcodes are unique ...              ", end = "", flush = True)
lag()
if not len(df_mini["barcode"]) == len(set(df_mini["barcode"])):
    counts = pd.DataFrame(df_mini['barcode'].value_counts())
    counts.columns = ["count"]
    counts = counts[counts["count"] > 1]
    #print(nl, counts)
    raise Exception(f"{nl}One or more barcodes are duplicated. Each barcode may only be used once:{nl}{counts}")
print("OK")

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

if not monkey_mode:
    check_user("Do you wish to proceed?")






################
# Check rundir #
################

# Wait for the rundir to occur in the specified path. 
# If it doesn't occur after a specified waiting time, then stop the p

if rundir[-1] == "/":
    print("Removing trailing slash from rundir")
    rundir = rundir[0:-1]


print("Checking that the rundir exists ... ", end = "", flush = True)
if not os.path.isdir(rundir):
    raise Exception(f"The rundir does not exist.")
print("OK")


print(f"Looking for MinKNOW-characteristic output:") #, end = "", flush = True)


for i in range(200):
    print("  Looking ... ", end = "", flush = True)
    fastq_pass_bases = glob.glob(rundir + "/**/fastq_pass", recursive = True) # Find any occurrence of the wanted path
    if len(fastq_pass_bases) == 0:
        print("nothing found yet, waiting 10 secs ...")
        time.sleep(10)
    elif(i == 10):
        print() # clean newline
        raise Exception("nothing found after 10 tries. Aborting.")
    else: 
        print(f"Found the sought after output.")
        break

print("")

if not len(fastq_pass_bases) == 1:
    raise Exception(f"There seems to be more than one fastq_pass sub-directory beneath the given rundir. These paths were found: {nl.join(fastq_pass_bases)}")

fastq_pass_base = fastq_pass_bases[0]
del fastq_pass_bases
print(f"Found the following fastq_pass base: {nl}  {fastq_pass_base}{nl}  This will be regarded as the input_base directory from now on.")
print()

sample_sheet_given_file = f"{fastq_pass_base}/../sample_sheet_given.tsv"
print(f"Backing up the original sample sheet in this location: {sample_sheet_given_file}")
df.to_csv(sample_sheet_given_file, sep = "\t")


print("Listing barcodes present on the disk:")
disk_barcodes_list  = sorted(glob.glob(fastq_pass_base + "/barcode*")) # Find all fastq_pass/barcode* directories
disk_barcodes_df = pd.DataFrame({'barcode_path': disk_barcodes_list})

#df_mini = df_mini.assign(type = ['positive_control' if a.lower().startswith("seqpos") else ('negative_control' if a.lower().endswith("neg") else 'sample') for a in df['sample_id']])
disk_barcodes_df = disk_barcodes_df.assign(barcode = ["NB" + i[-2:] for i in disk_barcodes_df["barcode_path"]])



print("Continuing with these barcodes:")

workflow_table = disk_barcodes_df.merge(df_mini, how='left', on='barcode') # left join (merge) the present barcodes onto the df_mini table.


print(workflow_table.to_string(index = False))
print("//")
print()





if not monkey_mode:
    check_user("This is the samples found on disk that matches your input. Do you wish to proceed?")


#input_continue = input("continue? (y/n) ")
#if not input_continue.lower()[0] == "y":
#    exit("Quitting ...")
#


# Create the dirs


raise Exception("limit")

try:
    os.mkdir(out_base_var)
    os.mkdir("logs")
except OSError:
    print ("Creation of the directories")
else:
    print ("Successfully created the directories")


# Collect all targets
rule all:
    input: expand(["{out_base}/metadata.tsv", \
                   "{out_base}/samples/{sample}/{sample}.fa", \
                   "{out_base}/samples/{sample}/prokka/{sample}.gff", \
                   "{out_base}/roary/summary_statistics.txt", \
                   "{out_base}/abricate/card_detail.tsv", \
                   "{out_base}/mlst/mlst.tsv", \
                   "{out_base}/mashtree/mashtree.newick", \
                   "{out_base}/fasttree/fasttree.newick"], \
                  out_base = out_base_var, sample = df["sample"]) # copy


  

# Write the df table to the directory for later reference.
rule metadata:
    input: df["input_file"].tolist()
    output: "{out_base}/metadata.tsv"
    run: 
        df.to_csv(str(output), index_label = "index", sep = "\t")



# Copy the input file to its new home
rule copy:
    #input: "{sample}"
    input: lambda wildcards: df[df["sample"]==wildcards.sample]["input_file"].values[0]
    output: "{out_base}/samples/{sample}/{sample}.fa"
    #log: "logs/{out_base}_{wildcards.sample}.out.log"

    shell: """

        mkdir -p logs output_asscom1

        cp {input} {output}

        """


##################################
# Targets for each sample below: #
##################################
rule prokka:
    input: "{out_base}/samples/{sample}/{sample}.fa"
    output: "{out_base}/samples/{sample}/prokka/{sample}.gff"
    #conda: "envs/prokka.yml"
    container: "docker://staphb/prokka"
    threads: 4
    shell: """

        prokka --cpus {threads} --force --outdir {wildcards.out_base}/samples/{wildcards.sample}/prokka --prefix {wildcards.sample} {input} || echo exit 0

        """


#######################################
# Targets for the complete set below: #
#######################################
rule roary:
    input: expand("{out_base}/samples/{sample}/prokka/{sample}.gff", sample = df["sample"], out_base = out_base_var)
    output: ["{out_base}/roary/summary_statistics.txt", "{out_base}/roary/core_gene_alignment.aln", "{out_base}/roary/gene_presence_absence.csv"]
    params:
        blastp_identity = 95, # For clustering genes
        core_perc = 99  # Definition of the core genome
    #conda: "envs/roary.yml"
    threads: 8
    container: "docker://sangerpathogens/roary"
    shell: """


        # Roary is confused by the way snakemake creates directories ahead of time.
        # So I will delete it manually here before calling roary.
        rm -r {wildcards.out_base}/roary

        roary -a -r -e --mafft -p {threads} -i {params.blastp_identity} -cd {params.core_perc} -f {wildcards.out_base}/roary {input}
                
        
        """


rule abricate:
    #input: expand("{out_base}/samples/{sample}/{sample}.fa", sample = df["sample"], out_base = out_base_var)
    input: df["input_file"].tolist()
    output:
        card_detail = "{out_base}/abricate/card_detail.tsv",
        card_sum = "{out_base}/abricate/card_summary.tsv",
        plasmidfinder_detail = "{out_base}/abricate/plasmidfinder_detail.tsv",
        plasmidfinder_sum = "{out_base}/abricate/plasmidfinder_summary.tsv",
        ncbi_detail = "{out_base}/abricate/ncbi_detail.tsv",
        ncbi_sum = "{out_base}/abricate/ncbi_summary.tsv"
    container: "docker://staphb/abricate"
    shell: """


        # TODO: update these databases

        abricate --db card {input} > {output.card_detail}
        abricate --summary {output.card_detail} > {output.card_sum}
        
        abricate --db plasmidfinder {input} > {output.plasmidfinder_detail}
        abricate --summary {output.plasmidfinder_detail} > {output.plasmidfinder_sum}
        
        abricate --db ncbi {input} > {output.ncbi_detail}
        abricate --summary {output.ncbi_detail} > {output.ncbi_sum}
        


        """


rule mlst:
    #input: expand("{out_base}/samples/{sample}/{sample}.fa", sample = df["sample"], out_base = out_base_var)
    input: df["input_file"].tolist()
    output: "{out_base}/mlst/mlst.tsv"
    container: "docker://staphb/mlst"
    shell: """

        mlst {input} > {output}

        """




rule mashtree:
    #input: expand("{out_base}/samples/{sample}/{sample}.fa", sample = df["sample"], out_base = out_base_var)
    input: df["input_file"].tolist()
    output: "{out_base}/mashtree/mashtree.newick"
    container: "docker://staphb/mashtree"
    threads: 4
    shell: """

        mashtree --numcpus {threads} {input} > {output}

        """




rule fasttree:
    #input: expand("{out_base}/samples/{sample}/{sample}.fa", sample = df["sample"], out_base = out_base_var)
    input: "{out_base}/roary/core_gene_alignment.aln"
    output: "{out_base}/fasttree/fasttree.newick"
    container: "docker://staphb/fasttree"
    threads: 4
    shell: """

        OMP_NUM_THREADS={threads}

        FastTree -nt {input} > {output} 2> {output}.log

        """



rule roary_plots:
    input: genes = "{out_base}/roary/gene_presence_absence.csv",
        tree = "{out_base}/fasttree/fasttree.newick"
    output: "{out_base}/roary_plots/whatever"
    container: "docker://python" # Make our own python container with cairosvg perl etc...
    shell: """
        
        # Failing because matplotlib is missing...
        python3 scripts/roary_plots.py {input.tree} {input.genes} > hat 2> hat.err

        # TODO: add the other weird stuff from https://github.com/cmkobel/assemblycomparator/blob/61c9a891a75e2f252dc54185d74c0fbb092815e5/workflow_templates.py#L489
        """



#print(mashtree.input)







