#!/bin/bash
#SBATCH --time=0-4
#SBATCH --mem=200GB
#SBATCH --out=Integrate-Search_%j.out
#SBATCH --error=Integrate-Search_%j.err


module load r/4.1.2

source configFile.txt

Rscript RunIntegrateSearch.R ${query} ${ref} ${k} ${out_dir} ${project_name} ${cellID} ${s}
