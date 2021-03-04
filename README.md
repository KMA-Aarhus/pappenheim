# ONT-SARS-CoV-2 (Work in progress)
Pipeline til sekventering af SARS-CoV-2 på lokale workstations med realtidsbasecalling. Pipelinen er baseret på snakemake og bruger i videste omfang conda. Alle mellemværende outputs er lang-pivoterede.





## Installation 


```
# Install conda by downloading and following the instructions:
wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh
chmod +x Miniconda3-latest-Linux-x86_64.sh
./Miniconda3-latest-Linux-x86_64.sh

# Install snakemake
mamba create -c conda-forge -c bioconda -n snakemake snakemake
```


## Brug pipelinen
```
conda activate snakemake
snakemake --help

# Hav et ark med metadata klar til kørslen:
echo "
```



 - [Installer conda](https://conda.io/projects/conda/en/latest/user-guide/install/linux.html)
 - Installer (https://snakemake.readthedocs.io/en/stable/getting_started/installation.html)
 - 
