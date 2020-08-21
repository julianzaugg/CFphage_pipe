Snakemake pipeline to process Cystic Fibrosis isolates for phage discovery and evaluation

### Setup
Recommend creating a new conda environment to help avoid conda package conflicts

```shell script
conda create -n cfphage_pipe -c bioconda snakemake conda=4.8.4
conda activate cfphage_pipe
```

### Usage
For QC
```shell script
snakemake --use-conda --latency-wait 10 --keep-going --cores 30 qc --configfile config.yaml
```
For QC + assembly
```shell script
snakemake --use-conda --latency-wait 10 --keep-going --cores 30 assemble --configfile config.yaml
```
For QC + assembly + polishing
```shell script
snakemake --use-conda --latency-wait 10 --keep-going --cores 30 polish --configfile config.yaml 
```
For reference genome coverage
```shell script
snakemake --use-conda --latency-wait 10 --keep-going --cores 30 coverage_reference_genomes_all --configfile config.yaml 
```