#!/bin/bash


# You may source this bashrc-extension, if you wish to have some shortcuts which can be nice for development and deployment.



# Neat aliases
alias mconda='mamba'
alias survey='watch "sensors; nvidia-smi"'
alias start_rampart='{ $(sleep 5; firefox localhost:3000)  & }; conda activate artic-rampart && rampart --protocol ~/repos/artic-ncov2019/rampart/ --clearAnnotated --basecalledPath'
alias citament='git add -u && git commit -m "amend" &&  git pull && git push && echo OK'






#alias pappenheim='cd ~/pappenheim && snakemake --profile default --config'


start_pappenheim () {

	example="\n\texample:\n\t start_pappenheim path/to/my_samplesheet.xlsx path/to/my_rundir\n"

if [ -z "$1" ]; then
    echo "Input error: The samplesheet argument is empty. Please specify a samplesheet."
    echo -e $example
     
elif [ -z "$2" ]; then
    echo "Input error: The rundir argument is empty. Please specify a rundir."
    echo -e $example
else

    #clear
    cd ~/pappenheim && snakemake --profile default --config samplesheet="${1}" rundir="${2}" ${3} && echo && cowsay "pappenheim finished successfully."
fi
}



