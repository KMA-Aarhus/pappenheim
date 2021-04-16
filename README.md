# Pappenheim ONT-SARS-CoV-2 (Work in progress)

Pipeline for sequencing of SARS-CoV-2 on local workstations with realtime-basecalling (GPU). This pipeline is based on snakemake and uses conda. All intermediary outputs are long-pivoted.



## Installation 


1. Install conda (if you haven't already):

    ```
    cd ~
    wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh
    chmod +x Miniconda3-latest-Linux-x86_64.sh
    ./Miniconda3-latest-Linux-x86_64.sh
    ```
    
2. We recommend to use mamba (inside conda) for installing packages onto the miniconda environment.

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


## Usage

Have a  .xlsx-sample sheet ready for input to the pipeline. This samplesheet should have at least to columns present: "barcode" and "sample id". The barcode column should contain barcode names in the EXP-NBD104 and EXP-NBD114 kits. The "sample id" column has specific requirements to the formatting, but please be aware that spaces will be truncated. Positive- and negative controls are recognized with the following specification: Positive controls must have a sample id starting with "seqpos". Negative controls must have a sample id ending with "neg". For both types of controls the recognition is case insensitive.




Denne fil gives til pipelinen ved kørslen ved argumentet `--samplesheet`
Der skal også specificeres en kørselsmappe fra minknow, hvor de basecallede fastq filer ligger. Denne gives ved argumentet `--rundir`.

```
snakemake --samplesheet path/to/above/file.csv --rundir path/to/minknow-output/
```



Hvis der er problemer med inputtet, gives en advarsel.
Ellers starter pipelinen og skulle gerne passe sig selv.

Efterhånden som genomerne bliver behandlet kommer outputtet i en mappe ved navn "output".

Der vil være en mappe med hvert isolat hvor der findes reads.
Der vil være en pangolin batch-kørsel.
Der vil være en nextclade batch-kørsel.

Pangolin og Nextclade kørslerne vil blive integreret med input metadata, og uploades i sidste ende til GenomeDK serveren.
På GenomeDK sker patientdata-integrationen, og derefter uploades hele baduljen til SSI.


