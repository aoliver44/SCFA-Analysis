#!/bin/bash

## ML_script2.sh to write and submit scripts for taxaHFE
## Author: Andrew Oliver
## Date: Nov 1, 2023
## assumptions: working in this dir on spitfire: /share/lemaylab/aoliver/SCFA
## to run: bash ML_script2.sh

mkdir -p taxaHFE_scripts

for scfa in {acetate,propionate,new_butyrate,total_scfa,acetate_norm_ratio_dist,propionate_norm_ratio_dist,new_butyrate_norm_ratio_dist,p_acetic_acid_nmol,p_propionic_acid_nmol,p_butyric_acid_nmol,p_scfa_nmol_total}; do

echo "#!/bin/bash
#SBATCH --partition=production # partition to submit to
#SBATCH --job-name=taxahfe_${scfa} # Job name
#SBATCH --nodes=1 # single node, anything more than 1 will not run
#SBATCH --ntasks=1 # number of tasks (like array, mpi parallelization)
#SBATCH --cpus-per-task=6 # equivalent to cpus
#SBATCH --mem=24000 # in MB, memory pool all cores, default is 2GB per cpu
#SBATCH --time=0-02:00:00  # expected time of completion in days, hours, minutes, seconds, default 1-day
#SBATCH --output=taxahfe_${scfa}.out # STDOUT
#SBATCH --error=taxahfe_${scfa}.err # STDERR
##SBATCH --mail-user=aoliver # commented out so i dont get a ton of emails
##SBATCH --mail-type=ALL # commented out so i dont get a ton of emails

## go to working director
cd /share/lemaylab/aoliver/SCFA

## load module for singularity
module load singularity

## run taxaHFEv2 command
singularity run -W `pwd` --bind `pwd`:/home/docker /share/lemaylab/aoliver/software/taxaHFE_v2.sif bash -c 'taxaHFE --subject_identifier subject_id --label ${scfa} --feature_type numeric -a 0 -L 3 --nperm 80 -n 6 /home/docker/raw_data/taxahfe_${scfa}.csv /home/docker/raw_data/merged_metaphlan_v4-0-6.txt /home/docker/input_data/${scfa}_microbe_taxaHFE.csv'

" > /share/lemaylab/aoliver/SCFA/taxaHFE_scripts/taxahfe_microbe_${scfa}.sh


echo "#!/bin/bash
#SBATCH --partition=production # partition to submit to
#SBATCH --job-name=taxahfe_food_${scfa} # Job name
#SBATCH --nodes=1 # single node, anything more than 1 will not run
#SBATCH --ntasks=1 # number of tasks (like array, mpi parallelization)
#SBATCH --cpus-per-task=6 # equivalent to cpus
#SBATCH --mem=24000 # in MB, memory pool all cores, default is 2GB per cpu
#SBATCH --time=0-02:00:00  # expected time of completion in days, hours, minutes, seconds, default 1-day
#SBATCH --output=taxahfe_food_${scfa}.out # STDOUT
#SBATCH --error=taxahfe_food_${scfa}.err # STDERR
##SBATCH --mail-user=aoliver # commented out so i dont get a ton of emails
##SBATCH --mail-type=ALL # commented out so i dont get a ton of emails

## go to working director
cd /share/lemaylab/aoliver/SCFA

## load module for singularity
module load singularity

## run taxaHFEv2 command
singularity run -W `pwd` --bind `pwd`:/home/docker /share/lemaylab/aoliver/software/taxaHFE_v2.sif bash -c 'taxaHFE --subject_identifier subject_id --label ${scfa} --feature_type numeric -a 0 -L 3 --nperm 80 -n 6 /home/docker/raw_data/taxahfe_${scfa}.csv /home/docker/raw_data/fl100_otu_abundance.txt /home/docker/input_data/${scfa}_food_taxaHFE.csv'

" > /share/lemaylab/aoliver/SCFA/taxaHFE_scripts/taxahfe_food_${scfa}.sh

done