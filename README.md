# ONT-SARS-CoV-2 (Work in progress)
Pipeline til sekventering af SARS-CoV-2 på lokale workstations med realtidsbasecalling. Pipelinen er baseret på snakemake og bruger i videste omfang conda. Alle mellemværende outputs er lang-pivoterede.





## Installation 


Install conda by downloading and following the [instructions](https://docs.conda.io/en/latest/miniconda.html):
```
cd
wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh
chmod +x Miniconda3-latest-Linux-x86_64.sh
./Miniconda3-latest-Linux-x86_64.sh
```

Clone this repository and install snakemake

```
cd 
# Install snakemake using conda
conda install -c conda-forge -c bioconda snakemake
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


