Snakemake pipeline to process Cystic Fibrosis isolates for phage discovery and evaluation

### Setup
Recommend creating a new conda environment to help avoid conda package conflicts

```shell script
conda create -n cfphage_pipe -c bioconda snakemake conda=4.8.4
conda activate cfphage_pipe
```

### Usage
```shell script
snakemake --use-conda --latency-wait 20 --cores 30 polishing --configfile config.yaml 
```
