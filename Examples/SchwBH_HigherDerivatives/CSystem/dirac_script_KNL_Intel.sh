#!/bin/bash
#SBATCH -p knl
#SBATCH -A DIRAC-DP128-SL4-KNL

#SBATCH --job-name=CSystemArray
#SBATCH -N 1
#SBATCH --ntasks-per-node=32
#SBATCH --cpus-per-task=2
#SBATCH --time=03:30:00
##SBATCH --mem=

INDEX=$SLURM_ARRAY_TASK_ID
TOTAL=$SLURM_ARRAY_TASK_COUNT
echo INDEX=$INDEX
echo TOTAL=$TOTAL

FOLDER=search_$INDEX
mkdir $FOLDER
cd $FOLDER
 
export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK
echo OMP_NUM_THREADS $OMP_NUM_THREADS

# Use Chombo timers
export CH_TIMER=TRUE

echo ID $SLURM_JOB_ID
echo PROCID $SLURM_PROCID
echo NAME $SLURM_JOB_NAME
echo
echo NODES $SLURM_JOB_NUM_NODES
echo TASKS_PER_NODE $SLURM_NTASKS_PER_NODE
echo CPUS_PER_TASK $SLURM_CPUS_PER_TASK
echo NTASKS $SLURM_NTASKS
echo

if [ "$SLURM_JOB_NODELIST" ]; then
    #! Create a machine file:
    export NODEFILE=`generate_pbs_nodefile`
    cat $NODEFILE | uniq > machine.file.$JOBID
    echo -e "\nNodes allocated:\n================"
    echo `cat machine.file.$JOBID | sed -e 's/\..*$//g'`
    echo
    rm -f machine.file.$JOBID
fi

module purge
module load use.own
module load grchombo_intel
module list

#python ../job_array_setup.py $INDEX $TOTAL
python ../job_array_setup.py $INDEX 900

time mpiexec ../Main_HigherDerivatives_CSystem3d_ch.Linux.64.mpiicpc.ifort.OPTHIGH.MPI.OPENMPCC.KNL.Intel.ex ../params_C.txt FILE=specific_params.txt



