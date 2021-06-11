...
#SBATCH --nodes=20
#SBATCH --ntasks-per-node=8
#SBATCH --cpus-per-task=6
...
export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK
...
 
# Fill just the following two
MY_NODES_PER_UNIT=2
MY_WORK_UNITS=128
 
# Just algebra, these need no changes
MY_WORKERS=$(( $SLURM_JOB_NUM_NODES / $MY_NODES_PER_UNIT ))
MY_TASKS_PER_WORKER=$(( $SLURM_NTASKS_PER_NODE * $MY_NODES_PER_UNIT ))
 
for (( I=1; I<=$MY_WORK_UNITS; I++ )) ; do # Count the work units that started
  while true ; do # Scan for a free worker
    SUBJOBS=`jobs -r | wc -l` # detect how many subjobs are already running
    if [ $SUBJOBS -lt $MY_WORKERS ] ; then  # submit only if at least one worker is free
      sleep 4 # wait before any submission
      # wrapper could also be an MPI program,`-c $SLURM_CPUS_PER_TASK` is only for OpenMP
      srun -N $MY_NODES_PER_UNIT -n $MY_TASKS_PER_WORKER -c $SLURM_CPUS_PER_TASK -J subjob.$I ./wrapper $I >OUT.$I &
      break # So "I" will get +1
    fi
  done
done
 
wait # for the last pool of work units
 
exit