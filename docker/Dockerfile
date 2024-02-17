FROM continuumio/miniconda3@sha256:cea2f5114fdbeb2303f55ffe96cc0f4170472f7a959e564649b9adb0127f09d7

MAINTAINER loipf

#RUN apt-get update 
#RUN conda update --all

RUN apt-get update && apt install -y procps gzip

RUN conda config --add channels defaults
RUN conda config --add channels bioconda
RUN conda config --add channels conda-forge


### bioinformatic tools
RUN conda install \
		multiqc>=1.9 \
		samtools>=1.11.0 \
		bcftools>=1.11 \
		pyvcf>=0.6.8 \
		open-cravat>=2.4.2


### set up open-cravat
RUN oc module install -y install-base
RUN oc module install -y csvreporter vcfreporter 
RUN oc module install -y hg38
#RUN oc module install -y clinvar cancer_genome_interpreter go target   ### small test

RUN oc module install -y clinvar dbsnp cosmic cosmic_gene cancer_genome_interpreter go gnomad3 gwas_catalog
RUN oc module install -y target pharmgkb cancer_hotspots cgl litvar ncRNA genehancer omim chasmplus_BRCA
RUN oc module install -y segway_blood alphamissense





