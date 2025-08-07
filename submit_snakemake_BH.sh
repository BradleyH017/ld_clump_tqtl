#!/bin/bash
#BSUB -o sm_logs/snakemake_master-%J-output.log
#BSUB -e sm_logs/snakemake_master-%J-error.log 
#BSUB -q oversubscribed
#BSUB -G team152
#BSUB -n 1
#BSUB -M 10000
#BSUB -a "memlimit=True"
#BSUB -R "select[mem>10000] rusage[mem=10000] span[hosts=1]"
#BSUB -J 1

# Define some params
config_var=scripts/config.yaml
worfklow_prefix="clump_"
group="team152"
workdir=${PWD}

# Load snakemake and singulatiry
module load HGI/common/snakemake/7
module load ISG/singularity/3.11.4
which singularity


# Copy config to results
thresh=$(grep "thresh" scripts/config.yaml | awk '{print $2}')
mkdir -p results_${thresh}
cp $config_var results_${thresh}/

# Make a log dir
mkdir -p sm_logs

# Run snakemake
snakemake -j 20000 \
    --latency-wait 90 \
    --use-envmodules \
    --rerun-incomplete \
    --keep-going \
    --directory ${workdir} \
    --cluster-config ${config_var} \
    --cluster-config cluster_config.yaml \
    --use-singularity \
    --singularity-args "-B /lustre,/software,/nfs/users/nfs_b/bh18/.local/lib/python3.7/site-packages" \
    --keep-going \
    --restart-times 0 \
    --snakefile $PWD/scripts/run_ldclump.smk

# bsub -M 10000 -a "memlimit=True" -R "select[mem>10000] rusage[mem=10000] span[hosts=1]" -o sm_logs/snakemake_master-%J-output.log -e sm_logs/snakemake_master-%J-error.log -q oversubscribed -J "snakemake_master_CLUMP" < submit_snakemake_BH.sh 
