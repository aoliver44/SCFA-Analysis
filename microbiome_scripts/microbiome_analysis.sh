#!/bin/bash

## set directory locations
WORKDIR=/share/lemaylab/aoliver/FL100_metagenomic_analysis/
IFDP_DB=/share/lemaylab/aoliver/software/IFDP/ec_full.dmnd

cd ${WORKDIR}

## make some output dirs
#mkdir -p ${WORKDIR}metaphlan_outs/
#mkdir -p ${WORKDIR}metaphlan_outs/sams/
#mkdir -p ${WORKDIR}metaphlan_outs/bowtie_outs/
#mkdir -p ${WORKDIR}humann3/
#mkdir -p ${WORKDIR}humann3/logs/
#mkdir -p ${WORKDIR}IFDP/
#mkdir -p ${WORKDIR}plant_blast/

## copy IFDP scripts
cp /share/lemaylab/aoliver/software/IFDP/* ${WORKDIR}

## loop through metagenomes and run merge/humann
while read line path mc; do
echo "#!/bin/bash
#SBATCH --partition=production # partition to submit to
#SBATCH --job-name=metagenome_${line} # Job name
#SBATCH --nodes=1 # single node, anything more than 1 will not run
#SBATCH --ntasks=1 # number of tasks (like array, mpi parallelization)
#SBATCH --cpus-per-task=8 # equivalent to cpus
#SBATCH --mem=32000 # in MB, memory pool all cores, default is 2GB per cpu
#SBATCH --time=1-00:00:00  # expected time of completion in days, hours, minutes, seconds, default 1-day
#SBATCH --output=${line}_metagenome.out # STDOUT
#SBATCH --error=${line}_metagenome.err # STDERR
##SBATCH --mail-user=${USER} #
##SBATCH --mail-type=ALL #

## work in workdir
echo 'Work directory is: ${WORKDIR}'
cd ${WORKDIR}

###########################
## Subsample to test
###########################
## comment out to run full metagenomes

## load modules
#module load bbmap/38.87

## subsample to exact reads
#reformat.sh in=${path} out=${WORKDIR}${line}_subsample.fastq samplereadstarget=250000

#mv ${WORKDIR}${line}_subsample.fastq ${WORKDIR}${line}.extendedFrags.fastq

#module purge

############################
## run humann3 + metaphlan4
###########################

## run metaphlan v4.0.6 because it has
## the massively expanded OCT 22 database
## see the release notes:
## https://github.com/biobakery/MetaPhlAn/releases

module load metaphlan/4.0.6
source activate metaphlan-4.0.6

echo 'Beginning Metaphlan4 analysis...'

metaphlan ${path} \
--input_type fastq \
--bowtie2db /software/anaconda3/23.1.0/lssc0-linux/envs/metaphlan-4.0.6/db/ \
--samout ${WORKDIR}metaphlan_outs/sams/${line}.metaphlan.sam \
--bowtie2out ${WORKDIR}metaphlan_outs/bowtie_outs/${line}.bowtie2.out \
--add_viruses \
--nproc 8 \
-o ${WORKDIR}metaphlan_outs/${line}.mpa4_profile.txt

conda deactivate
module purge

## Run humann 3.7, which plays nice with
## metaphlan 4.0.6...which we use for the
## october 22 expanded taxonomic database.

## trick to get humann3.7 working on barbera
## was to install the conda env NOT in /home/.conda
## ie conda create --prefix /share/lemaylab/aoliver/software/humann_conda_env/humann3.7
## (that gets installed wherever the prefix says)
## home directory is not readable from barbera

module load anaconda3/4.9.2
source activate /share/lemaylab/aoliver/software/humann_conda_env/humann3.7

echo 'Beginning Humann analysis...'

humann -i ${path} \
-o ${WORKDIR}humann3 \
--threads 8 \
--search-mode uniref90 \
--nucleotide-database /software/humann/3.6.1/lssc0-linux/db/chocophlan/ \
--protein-database /software/humann/3.6.1/lssc0-linux/db/uniref/ \
--taxonomic-profile ${WORKDIR}metaphlan_outs/${line}.mpa4_profile.txt


echo 'Humann analysis finished for ${line} sample'

conda deactivate
module purge

mv ${WORKDIR}humann3/${line}.extendedFrags_humann_temp/${line}.extendedFrags.log ${WORKDIR}humann3/logs/

############################
## run ifdp
###########################

## load up modules needed
module load diamond/2.0.13

## Get inferred Fiber profile, run ecn_map_mc.sh
## after.
## This code is based on this work PMID: 36464700

run_sample.sh \
-p 8 \
-d ${IFDP_DB} \
-i ${path} \
-o ${WORKDIR}IFDP/${line}_output

## We normalize the fiber degradation profile
## using some of their scripts and some custom
## stuff too i made.

#ecn_map_mc.sh ${line}_output ${mc}


############################
## clean up...more
###########################

echo 'Cleaning...'

rm -r ${WORKDIR}humann3/${line}.extendedFrags_humann_temp/
rm ${WORKDIR}${line}.tmp.fasta

echo 'Success! probably moving to next sample!'

" > metagenome_${line}.sh

sbatch metagenome_${line}.sh

sleep 1

done < FL100_metageonome_paths.txt