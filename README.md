### Running LD clumping to approximate eQTL aggregation

1. Prepare tests
To run LD clumping, need to generate 2 input files `input/variants_per_gene_list.txt` and `input/gene_list`, lists of unique gene x variant pairs and unique genes respectively. \
An example of how to prepare these from TensorQTL output and colocalisation is with `scripts/001_prep_tests.r`. This makes use of the master repo (https://github.com/andersonlab/IBDVerse-sc-eQTL-code) as the 'repo_dir'.

2. Adjust `scripts/config.yaml`
Specify the path to the genotype information (this must be plink converted already), the input files listed above, and set the desired LD clumping threshold.

3. Run LD clumping
Submit the workflow to an LSF cluster using:
```
bsub -M 10000 -a "memlimit=True" -R "select[mem>10000] rusage[mem=10000] span[hosts=1]" -o sm_logs/snakemake_master-%J-output.log -e sm_logs/snakemake_master-%J-error.log -q normal -J "snakemake_master_CLUMP" -G humgen-priority < submit_snakemake_BH.sh 
```