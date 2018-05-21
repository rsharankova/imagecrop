#!/bin/bash
#
#SBATCH --job-name=cropss
#SBATCH --output=log_cropss.txt
#SBATCH --ntasks=400
#SBATCH --time=20:00
#SBATCH --mem-per-cpu=2000

WORKDIR=/cluster/kappa/90-days-archive/wongjiradlab/twongj01/ssnet_retrain_cropper/grid
CONTAINER=/cluster/kappa/90-days-archive/wongjiradlab/larbys/images/dllee_unified/singularity-dllee-intertool-20180216.img
JOBIDLIST=${WORKDIR}/rerunlist.txt
INPUTLISTDIR=${WORKDIR}/inputlists
OUTPUTDIR=/cluster/kappa/90-days-archive/wongjiradlab/twongj01/ssnet_retrain_cropper/output/cocktail/p03

mkdir -p ${OUTPUTDIR}
module load singularity
srun singularity exec ${CONTAINER} bash -c "cd ${WORKDIR} && source run_job.sh ${WORKDIR} ${INPUTLISTDIR} ${JOBIDLIST} ${OUTPUTDIR}"

