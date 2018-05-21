#!/bin/sh

# REMEMBER, WE ARE IN THE CONTAINER RIGHT NOW
# This means we access the next work drives through some mounted folders

workdir=$1
inputlist_dir=$2
jobid_list=$3
output_dir=$4

echo "WORKDIR: $1"

export DLLEE_UNIFIED_BASEDIR=/usr/local/share/dllee_unified
export PATH=/cluster/kappa/90-days-archive/wongjiradlab/twongj01/ssnet_retrain_cropper/bin:${PATH}
export LD_LIBRARY_PATH=/cluster/kappa/90-days-archive/wongjiradlab/twongj01/ssnet_retrain_cropper:${LD_LIBRARY_PATH}

source /usr/local/bin/thisroot.sh
source /usr/local/share/dllee_unified/configure.sh

cd $workdir

let NUM_PROCS=`cat ${jobid_list} | wc -l`
echo "number of processes: $NUM_PROCS"
if [ "$NUM_PROCS" -lt "${SLURM_PROCID}" ]; then
    echo "No Procces ID to run."
    return
fi

let "proc_line=${SLURM_PROCID}+1"
echo "sed -n ${proc_line}p ${jobid_list}"
let jobid=`sed -n ${proc_line}p ${jobid_list}`
echo "JOBID ${jobid}"

slurm_folder=`printf slurm_genssnet_job%04d ${jobid}`
mkdir -p ${slurm_folder}
cd ${slurm_folder}/

# copy over input list
inputlist=`printf ${inputlist_dir}/input_%d.txt ${jobid}`

cp $inputlist input.txt
supera=`sed -n 1p input.txt`
mcinfo=`sed -n 2p input.txt`

logfile=`printf log_%04d.txt ${jobid}`
make_ssnet_vertexcrop $supera $mcinfo ssnetout_larcv.root >& $logfile || exit

outfile=`printf ${output_dir}/ssnet_retraining_%06d.root ${jobid}`
cp ssnetout_larcv.root $outfile

# clean up
cd ../
rm -r $slurm_folder
