#!/bin/bash

# Private aliases
# Should be started using the "artic-rampart" conda environment
alias start_rampart='{ sleep 5; firefox localhost:3000 & }; conda activate artic-rampart && rampart --clearAnnotated --protocol ~/repos/artic-ncov2019/rampart/ --basecalledPath'
alias survey='watch "sensors; nvidia-smi"'
alias start_rampart='{  & }; conda activate artic-rampart && rampart --protocol ~/repos/artic-ncov2019/rampart/ --clearAnnotated --basecalledPath'


#alias pappenheim='cd ~/pappenheim && snakemake --profile default --config'
start_pappenheim () {
 if [ -z "$1" ]; then
     echo "The samplesheet argument is empty. Please specify a samplesheet."
 elif [ -z "$2" ]; then
     echo "The rundir argument is empty. Please specify a rundir."
 else

     #clear
     cd ~/pappenheim && snakemake --profile default --config samplesheet="${1}" rundir="${2}" ${3} && echo "pappenheim finished successfully."
 fi
}


alias whatthe="echo it works"