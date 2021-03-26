# Pappenheim ONT-SARS-CoV-2 (Work in progress)
Pipeline til sekventering af SARS-CoV-2 på lokale workstations med realtidsbasecalling. Pipelinen er baseret på snakemake og bruger i videste omfang conda. Alle mellemværende outputs er lang-pivoterede.





## Installation 


1. Install conda:

    We recommend to use mamba (instead of conda) for installing packages onto the miniconda environment system.
    ```
    cd ~
    wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh
    chmod +x Miniconda3-latest-Linux-x86_64.sh
    ./Miniconda3-latest-Linux-x86_64.sh
    conda install mamba -n base -c conda-forge
    ```

2. Clone this repository and install snakemake

    ```
    # install pappenheim in the home directory
    cd ~
    
    # Make sure you clone recursively, in order for the artic submodule to be included in the download.
    git clone --recurse-submodules https://github.com/KMA-Aarhus/pappenheim.git 

    # Install snakemake using conda
    # If you wish to install the package NOT in the base environment, create an environment, activate it, before running the following line:
    mamba env update --file pappenheim/environment.yaml 
    ```


## Usage

Hav en csv-fil klar til kørslen. Denne fil indeholder information om barcodes og samplenames.

```
barcode,sample
NB01,R012938
NB02,R382988
NB03,R328837
```

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


