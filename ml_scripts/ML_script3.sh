#!/bin/bash

## ML_script3.sh to write and submit scripts for dietML
## Author: Andrew Oliver
## Date: Nov 14, 2023
## assumptions: working in this dir on spitfire: /share/lemaylab/aoliver/SCFA
## to run: bash ML_script3.sh

mkdir -p dietML_results
mkdir -p dietML_scripts

while read subject label cor split engine folds score hp min shap cores input outdir; do
job_name=$(basename $input .csv)

echo "#!/bin/bash
#SBATCH --partition=production # partition to submit to
#SBATCH --job-name=${job_name} # Job name
#SBATCH --nodes=1 # single node, anything more than 1 will not run
#SBATCH --ntasks=1 # number of tasks (like array, mpi parallelization)
#SBATCH --cpus-per-task=8 # equivalent to cpus
#SBATCH --mem=32000 # in MB, memory pool all cores, default is 2GB per cpu
#SBATCH --time=0-07:00:00  # expected time of completion in days, hours, minutes, seconds, default 1-day
#SBATCH --output=${job_name}.out # STDOUT
#SBATCH --error=${job_name}.err # STDERR
##SBATCH --mail-user=aoliver # commented out so i dont get a ton of emails
##SBATCH --mail-type=ALL # commented out so i dont get a ton of emails

module load singularity

singularity run -W `pwd` --bind `pwd`:/home/docker /share/lemaylab/aoliver/software/nutrition_tools.sif \
bash -c 'dietML --subject_identifier $subject --label $label --cor_level $cor \
--train_split $split --model $engine --folds $folds --metric $score \
--type regression --seed 400 --tune_length $hp --tune_time $min \
--shap $shap --ncores $cores /home/docker/input_data/${input} /home/docker/dietML_results/${outdir}'

singularity run -W `pwd` --bind `pwd`:/home/docker /share/lemaylab/aoliver/software/nutrition_tools.sif \
bash -c 'dietML --subject_identifier $subject --label $label --cor_level $cor \
--train_split $split --model $engine --folds $folds --metric $score \
--type regression --seed 9142 --tune_length $hp --tune_time $min \
--shap $shap --ncores $cores /home/docker/input_data/${input} /home/docker/dietML_results/${outdir}'

singularity run -W `pwd` --bind `pwd`:/home/docker /share/lemaylab/aoliver/software/nutrition_tools.sif \
bash -c 'dietML --subject_identifier $subject --label $label --cor_level $cor \
--train_split $split --model $engine --folds $folds --metric $score \
--type regression --seed 607 --tune_length $hp --tune_time $min \
--shap $shap --ncores $cores /home/docker/input_data/${input} /home/docker/dietML_results/${outdir}'

singularity run -W `pwd` --bind `pwd`:/home/docker /share/lemaylab/aoliver/software/nutrition_tools.sif \
bash -c 'dietML --subject_identifier $subject --label $label --cor_level $cor \
--train_split $split --model $engine --folds $folds --metric $score \
--type regression --seed 7053 --tune_length $hp --tune_time $min \
--shap $shap --ncores $cores /home/docker/input_data/${input} /home/docker/dietML_results/${outdir}'

singularity run -W `pwd` --bind `pwd`:/home/docker /share/lemaylab/aoliver/software/nutrition_tools.sif \
bash -c 'dietML --subject_identifier $subject --label $label --cor_level $cor \
--train_split $split --model $engine --folds $folds --metric $score \
--type regression --seed 507 --tune_length $hp --tune_time $min \
--shap $shap --ncores $cores /home/docker/input_data/${input} /home/docker/dietML_results/${outdir}'

singularity run -W `pwd` --bind `pwd`:/home/docker /share/lemaylab/aoliver/software/nutrition_tools.sif \
bash -c 'dietML --subject_identifier $subject --label $label --cor_level $cor \
--train_split $split --model $engine --folds $folds --metric $score \
--type regression --seed 8636 --tune_length $hp --tune_time $min \
--shap $shap --ncores $cores /home/docker/input_data/${input} /home/docker/dietML_results/${outdir}'

singularity run -W `pwd` --bind `pwd`:/home/docker /share/lemaylab/aoliver/software/nutrition_tools.sif \
bash -c 'dietML --subject_identifier $subject --label $label --cor_level $cor \
--train_split $split --model $engine --folds $folds --metric $score \
--type regression --seed 7471 --tune_length $hp --tune_time $min \
--shap $shap --ncores $cores /home/docker/input_data/${input} /home/docker/dietML_results/${outdir}'

singularity run -W `pwd` --bind `pwd`:/home/docker /share/lemaylab/aoliver/software/nutrition_tools.sif \
bash -c 'dietML --subject_identifier $subject --label $label --cor_level $cor \
--train_split $split --model $engine --folds $folds --metric $score \
--type regression --seed 1525 --tune_length $hp --tune_time $min \
--shap $shap --ncores $cores /home/docker/input_data/${input} /home/docker/dietML_results/${outdir}'

singularity run -W `pwd` --bind `pwd`:/home/docker /share/lemaylab/aoliver/software/nutrition_tools.sif \
bash -c 'dietML --subject_identifier $subject --label $label --cor_level $cor \
--train_split $split --model $engine --folds $folds --metric $score \
--type regression --seed 7268 --tune_length $hp --tune_time $min \
--shap $shap --ncores $cores /home/docker/input_data/${input} /home/docker/dietML_results/${outdir}'

singularity run -W `pwd` --bind `pwd`:/home/docker /share/lemaylab/aoliver/software/nutrition_tools.sif \
bash -c 'dietML --subject_identifier $subject --label $label --cor_level $cor \
--train_split $split --model $engine --folds $folds --metric $score \
--type regression --seed 1378 --tune_length $hp --tune_time $min \
--shap $shap --ncores $cores /home/docker/input_data/${input} /home/docker/dietML_results/${outdir}'


" > /share/lemaylab/aoliver/SCFA/dietML_scripts/${job_name}.sh

done < /share/lemaylab/aoliver/SCFA/dietML_manifest.txt