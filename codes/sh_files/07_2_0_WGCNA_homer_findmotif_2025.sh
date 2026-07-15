#!/bin/bash
#SBATCH --export=ALL # export all environment variables to the batch job
#SBATCH -D . # set working directory to .
#SBATCH -p pq # partition to use
#SBATCH --time=50:00:00 # maximum walltime
#SBATCH -A xxxxx # research project to submit under
#SBATCH --array=1-4
#SBATCH --cpus-per-task=8
#SBATCH --exclusive
#SBATCH --job-name=homer_wgcna
#SBATCH --mem=10G
#SBATCH --error=homer_%A_%a.err
#SBATCH --output=homer_%A_%a.out

conda init
conda activate homer

ob=$(sed -n "${SLURM_ARRAY_TASK_ID}p" group_list.txt)
fd=./${ob}/wgcna_modules
md_list=./$ob/WGCNA_module_list.txt

echo "starting with modules for" $ob

while read module; do
findMotifs.pl $fd/${module}.txt \
zebrafish \
./${ob}/results/${module}_adv2kd1k \
-start -2000 -end 1000 -p 8 -mis 3 -cache 1000 \
&> ./${ob}/${module}_out2k1k
done < md_list



