# Pappenheim ONT-SARS-CoV-2

Pipeline for sequencing of SARS-CoV-2 on local workstations with realtime-basecalling (GPU). This pipeline is based on snakemake and uses conda. All intermediary outputs are long-pivoted.




## What it does

Pappenheim is designed to trail a ONT-minknow sequencing run. When you have started sequencing in minknow, you can start this pipeline and it will progress through the following steps.

### On a batch level:
* Input validation of the sample sheet
* Starting rampart on the rundir for realtime sequencing oversight

### On a sample level:

When the sequencing and basecalling is complete, the following steps are completed as well:

* Filtering of reads using **artic guppyplex**
* Consensus calling using **artic minion** via **nanopolish**
   * Using the [artic nCov-2019 V3](https://github.com/artic-network/artic-ncov2019) reference scheme
* Variant typing using **pangolin** and **nextclade**
* Optionally: Upload the variant data to a central server which can integrate patient data (pappenheim-receiver).


When pappenheim runs, it automatically checks and installs the newest versions of pangolin and nextclade. Given that you have a working internet connection at the time of starting the pipeline.





## Usage

Pappenheim needs information on the samples you are analysing and the rundir directory where the raw data is stored. Below is a description of why and how these arguments must be formatted for pappenheim.

### CLI-argument: `samplesheet`
Have a  .xlsx-sample sheet ready for input to the pipeline. This samplesheet should have at least to columns present: "barcode" and "sample id". The barcode column should contain barcode names in the EXP-NBD104 and EXP-NBD114 kits. The "sample id" column has specific requirements to the formatting, but please be aware that spaces will be truncated. Positive- and negative controls are recognized with the following specification: Positive controls must have a sample id starting with "seqpos". Negative controls must have a sample id ending with "neg". For both types of controls the recognition is case insensitive. 

![](https://github.com/KMA-Aarhus/pappenheim/blob/main/documentation/Screenshot%202021-04-16%20at%2010.03.01.png)

### CLI-argument: `rundir` 
Because pappenheim is designed to trail the output of minknow, you must first have started the sequencing (including high accuracy basecalling) in the minknow user interface and decided where the output should be written. Pappenheim needs to know this rundir such that it can start taking the correct reads for downstream analysis.


When you have a sample sheet ready and know the rundir, you can start the pappenheim pipeline:

Because pappenheim starts rampart automatically based on the given rundir, you should start pappenheim just after minknow has started sequencing.

```
cd ~/pappenheim
snakemake --config samplesheet='path/to/samplesheet.csv' rundir='path/to/minknow-output/'
```

First, the pappenheim pipeline validates the sample sheet and checks that the necessary columns exists and are correctly formatted. It also checks that the barcodes are unique. It then proceeds to check that the rundir exists. If it doesn't, pappenheim waits a few minutes and tries again.

When the sequencing is done in minknow, minknow writes a specific file to the rundir: "sequencing_summary_\*.txt". This file is necessary for consensus calling and thus pappenheim can only start when minknow is done basecalling.





## Installation 


1. Install conda (if you haven't already):

    ```
    cd ~
    wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh
    chmod +x Miniconda3-latest-Linux-x86_64.sh
    ./Miniconda3-latest-Linux-x86_64.sh
    ```
    
2. We recommend to use mamba (inside conda) for installing packages onto the miniconda environment.
   
   Hint: mamba is simply faster than the stock conda dependency solver and installer.

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

    Optionally, uncomment the `conda create ...` command if you want pappenheim and its dependencies in an enclosed environment.

    ```    
    # conda create -n pappenheim && conda activate pappenheim    
    
    mamba env update --file pappenheim/environment.yaml 
    ```
    
Please create an issue (in this repo) if you encounter any problems or unanswered questions during installation or use.





## Affiliations

Pappenheim is developed at Department of Clinical Microbiology (KMA) at Aarhus Universityhospital, Denmark, in response to the political decision of sequencing all positive SARS-CoV-2 samples.


