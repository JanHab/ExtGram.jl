# SLURM

## A Collection of Commands

Submit a job

``` bash
sbatch batch_script.sh
```

Running jobs
``` bash
squeue --me
```

Cancel job
``` bash
scance --me # Cancels all of your jobs
scancel -v <JOB_ID> # Cancel single job
```

Make script usable before
``` bash
chmod 775 myscript.sh
./myscript.sh
```