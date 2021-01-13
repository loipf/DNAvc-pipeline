#!/bin/sh

### maybe set up everything with bioconda ? - dont know if they have newest versions
### export tool-pathes to .bashrc ?


TOOL_DIR="tools/"

MULTIQC_VERSION="1.9"
DEEPTOOLS_VERSION="3.5.0"


mkdir $TOOL_DIR
cd $TOOL_DIR

### nextflow v20.07.1
curl -fsSL https://get.nextflow.io | bash


### MultiQC
pip install multiqc==$MULTIQC_VERSION



### DeepTools
#wget https://github.com/deeptools/deepTools/archive/$DEEPTOOLS_VERSION.tar.gz
#tar -xzvf $DEEPTOOLS_VERSION.tar.gz
#rm -r $DEEPTOOLS_VERSION.tar.gz
#cd deepTools-$DEEPTOOLS_VERSION
#python setup.py install --prefix /User/Tools/deepTools2.0

pip install deeptools==$DEEPTOOLS_VERSION





### open-cravat
pip install open-cravat
oc module install-base
oc module install -y csvreporter vcfreporter 

# oc module ls -a -t annotator   ### check all annotators
# oc gui

### test
oc module install -y clinvar cancer_genome_interpreter go target


### reality
oc module install -y clinvar dbsnp cosmic cosmic_gene cancer_genome_interpreter go gnomad 


### wish
oc module install -y clinvar dbsnp cosmic cosmic_gene cancer_genome_interpreter go gnomad uniprot
oc module install -y target mutpanning pharmgkb cancer_hotspots cgl litvar ncRNA genehancer omim chasmplus_BRCA


oc run ./example_input -l hg38








