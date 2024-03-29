# DNAvc pipeline

a DNA variant calling pipeline from indexed `.bam` files to `.vcf` files with stats report.

(DeepVariant currently optimised for Illumina WGS)


---
### set up pipeline


before running, you have to set up the attached Docker image which includes OpenCRAVAT and 2 others (all will take quite some time and overall need >50 GB of memory, you can add/remove OpenCRAVAT modules according to your needs in the Dockerfile, currently adapted for breast cancer):
```sh
docker build -t dnavc-pipeline https://raw.githubusercontent.com/loipf/DNAvc-pipeline/master/docker/Dockerfile
docker pull google/deepvariant:1.2.0
docker pull quay.io/mlin/glnexus:v1.2.7
```

now either replace the Docker container hash (last output line from previous build command) in `nextflow.config` or run nextflow with the `-with-docker dnavc-pipeline` argument.


---
### run mapping pipeline

it can be run locally with downloaded github-repo and edited `nextflow.config` file with:
```sh
nextflow run main.nf
```

or

```sh
nextflow run loipf/DNAvc-pipeline -r v1.0 --project_dir /path/to/folder --reads_mapped_dir /path/to/samples --reference_genome /path/to/ref --num_threads 10 -with-docker dnavc-pipeline
```
for this execution to work properly, you have to be in the current project directory.


optional extendable with:
```sh
-resume
-with-report report_DNAseq-pipeline
-with-timeline timeline_DNAseq-pipeline
-w work_dir
```


by default, all output will be saved into the `data` folder of the current directory.
mapped `.bam` files must already be indexed.




