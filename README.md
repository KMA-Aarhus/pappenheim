# Pappenheim ONT-SARS-CoV-2 (Work in progress)

Pipeline for sequencing of SARS-CoV-2 on local workstations with realtime-basecalling (GPU). This pipeline is based on snakemake and uses conda. All intermediary outputs are long-pivoted.


## What it does

Pappenheim is designed to trail a ONT-minknow sequencing run. When you have started sequencing in minknow, you can start this pipeline and it will progress through the following steps.

### On a batch level:
* Input validation of sample sheet
* Starting rampart on the rundir for realtime sequencing oversight

### On a sample level:
* Filtering of reads using **artic guppyplex**
* Consensus calling using **artic minion** via **nanopolish**
   * Using the [artic-ncov2019](https://github.com/artic-network/artic-ncov2019) nCoV-2019 V3 reference scheme
* Variant typing using **pangolin** and **nextclade**


When pappenheim runs, it automatically checks and installs the newest versions of pangolin and nextclade.



## Usage

### Sample sheet
Have a  .xlsx-sample sheet ready for input to the pipeline. This samplesheet should have at least to columns present: "barcode" and "sample id". The barcode column should contain barcode names in the EXP-NBD104 and EXP-NBD114 kits. The "sample id" column has specific requirements to the formatting, but please be aware that spaces will be truncated. Positive- and negative controls are recognized with the following specification: Positive controls must have a sample id starting with "seqpos". Negative controls must have a sample id ending with "neg". For both types of controls the recognition is case insensitive. 

### Rundir 
Because pappenheim is designed to trail the output of minknow, you must first have started the sequencing (including high accuracy basecalling) in the minknow user interface and decided where the output should be written. Pappenheim needs to know this rundir such that it can start taking the correct reads for downstream analysis.


When you have a sample sheet ready and know the rundir, you can start the pappenheim pipeline:

```
cd ~/pappenheim
snakemake --samplesheet path/to/above/file.csv --rundir path/to/minknow-output/
```

First, the pappenheim pipeline validates the sample sheet and checks that the necessary columns exists and are correctly formatted. It also checks that the barcodes are unique. It then proceeds to check that the rundir exists. If it doesn't, pappenheim waits a few minutes and tries again.

When the sequencing is done, minknow writes a specific file to the rundir: "sequencing_summary_\*.txt". This file is necessary for consensus calling and thus pappenheim can only start when minknow is done basecalling.



## Installation 


1. Install conda (if you haven't already):

    ```
    cd ~
    wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh
    chmod +x Miniconda3-latest-Linux-x86_64.sh
    ./Miniconda3-latest-Linux-x86_64.sh
    ```
    
2. We recommend to use mamba (inside conda) for installing packages onto the miniconda environment.
   Hint: mamba is simply faster

   ```
   conda install mamba -n base -c conda-forge
   ```

3. Clone this repository
    ```
    # install pappenheim in the home directory
    cd ~
    
    # Make sure you clone recursively, in order for the artic and pangolin submodules to be included in the download.
    git clone --recurse-submodules https://github.com/KMA-Aarhus/pappenheim.git 
    ```
    
4. Install snakemake and other dependencies needed for pappenheim:

    Optional: Make sure to uncomment the `conda create ...` command if you want pappenheim in its own environment.

    ```    
    # conda create -n pappenheim && conda activate pappenheim    
    
    mamba env update --file pappenheim/environment.yaml 
    ```



