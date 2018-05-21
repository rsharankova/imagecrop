#!/bin/bash

CONTAINER=/cluster/kappa/90-days-archive/wongjiradlab/larbys/images/dllee_unified/singularity-dllee-unified-taggerv2alpha-20171121.img
WORKDIR=/cluster/kappa/90-days-archive/wongjiradlab/twongj01/ssnet_retrain_cropper/grid

module load singularity

singularity exec ${CONTAINER} bash -c "source /usr/local/bin/thisroot.sh && cd /usr/local/share/dllee_unified && source configure.sh && cd ${WORKDIR} && python check_jobs.py"

