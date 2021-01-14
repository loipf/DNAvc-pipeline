#!/bin/bash

############################################
### specifiy directories, tools and properties

data_dir="/home/stefan/Documents/umcg/bash_dna-seq_pipeline/data"

#tool_fastqc="/home/stefan/tools/FastQC/fastqc"
#tool_cutadapt="/home/stefan/.local/bin/cutadapt"
#tool_multiqc="/home/stefan/miniconda3/bin/multiqc"
#tool_samtools="/home/stefan/tools/samtools-1.10/samtools"
#tool_bwa="/home/stefan/tools/bwa/bwa"
#tool_deeptools="/home/stefan/miniconda3/bin/deeptools"

tool_bcftools="/home/stefan/tools/bcftools-1.11/bcftools" # needs tabix and bgzip
tool_vep="/home/stefan/tools/ensembl-vep/vep"

num_threads="3"



deepvariant_version="1.0.0"
glnexus_version="v1.2.7"  

############################################

mkdir -p $data_dir/variants_vcf


### install deepvariant docker
#sudo docker pull google/deepvariant:"${deepvariant_version}"
### docker id 712e5bbcf309

### install GLnexus docker
#sudo docker pull quay.io/mlin/glnexus:"${glnexus_version}" 
### docker id 1f9c94778a08



### index reference for deepvariant
#gunzip $data_dir/Homo_sapiens.GRCh38.dna.alt.fa.gz
#$tool_samtools faidx $data_dir/Homo_sapiens.GRCh38.dna.alt.fa
#gzip -v $data_dir/Homo_sapiens.GRCh38.dna.alt.fa



############################################
### iterate over samples and map

for dir_path in $data_dir/reads_mapped/[^_]*
do
#    dir_path="/home/stefan/Documents/umcg/bash_dna-seq_pipeline/data/reads_mapped/SRR037207"

    sample_id=$(basename $dir_path)
    #bam_file=$(find $dir_path -name "*.bam")

    ### call variants
	bam_file_in_data=/reads_mapped/$sample_id/$sample_id.bam
	sample_out_dir=$data_dir/variants_vcf/$sample_id
    mkdir -p $sample_out_dir
    cd $sample_out_dir
    

### output databases version
oc module ls -t annotator > oc_databases_version.txt



### with GPU
#	docker run --gpus 1 \
###		--intermediate_results_dir=/output/deepvariant_intermediate_results \
###		--sample_name   # not implemented yet

	docker run \
		-v $data_dir:"/data" \
		-v $sample_out_dir:"/output" \
		google/deepvariant:$deepvariant_version \
		/opt/deepvariant/bin/run_deepvariant \
		--model_type=WGS \
		--ref=/data/Homo_sapiens.GRCh38.dna.alt.fa \
		--reads=/data/$bam_file_in_data \
		--output_vcf=/output/$sample_id.vcf.gz \
		--output_gvcf=/output/$sample_id.g.vcf.gz \
		--num_shards=$num_threads

#	$tool_bcftools stats --threads $num_threads $sample_out_dir/$sample_id.vcf.gz > $sample_id\_vcfstats.txt
#	$tool_vep -i $sample_out_dir/$sample_id.vcf.gz -o $sample_id\_vep.txt --database    # --offline

### VEP plugins: 
#- spliceAI/MMsplice 
#- HPO https://github.com/molgenis/vip/blob/master/plugins/vep/Hpo.pm
#- dbNSFP https://m.ensembl.org/info/docs/tools/vep/script/vep_plugins.html
#- dbscSNV
#- DisGeNET
#- GO
#- COSMIC
#- Mastermind ?
#- clinvar https://www.ensembl.org/info/docs/tools/vep/script/vep_custom.html
# use open-cravat - combination of a lot of variant callers ? https://github.com/KarchinLab/open-cravat
### driver mutation significance: https://opencravat.org/examples.html
# chasmplus
# gnomAD



done

# VEP
# annovar https://doc-openbio.readthedocs.io/projects/annovar/en/latest/
# snpeff
# (snpsift)

# GLnexus https://github.com/dnanexus-rnd/GLnexus/wiki/Getting-Started


############################################
### create common .vcf and VEP

mkdir -p $data_dir/variants_vcf/_all

#docker run --rm -i \
#	-v $data_dir/variants_vcf:"/in" \
#	quay.io/mlin/glnexus:$glnexus_version \
#    bash -c "glnexus_cli --threads ${num_threads} --config DeepVariant /in/*/*.g.vcf.gz" > $data_dir/variants_vcf/_all/glnexus_vcf_all.bcf
###--mem-gbytes 5
#$tool_bcftools view $data_dir/variants_vcf/_all/glnexus_vcf_all.bcf | bgzip -@ $num_threads -c > $data_dir/variants_vcf/_all/glnexus_vcf_all.vcf.gz




#bcftools view dv_1000G_ALDH2.bcf | bgzip -@ 4 -c > dv_1000G_ALDH2.vcf.gz




############################################
### create quality reports


#multiqc -f -o $data_dir/quality_reports/variants_vcf $data_dir/variants_vcf



















#mkdir -p $out_dir/data/quality_reports

#read_files_mapped=$(find $out_dir/data/reads_mapped -name "*.bam")


############################################
### deeptools analysis
#mkdir -p $out_dir/data/reads_mapped/_deepTools
#multiBamSummary bins -p $num_threads --smartLabels --bamfiles $read_files_mapped -o $out_dir/data/reads_mapped/_deepTools/multiBamSummary.npz
#plotCorrelation --corData $out_dir/data/reads_mapped/_deepTools/multiBamSummary.npz --corMethod spearman --whatToPlot heatmap --outFileCorMatrix $out_dir/data/reads_mapped/_deepTools/plotCorrelation_matrix.tsv
#plotPCA --corData $out_dir/data/reads_mapped/_deepTools/multiBamSummary.npz --outFileNameData $out_dir/data/reads_mapped/_deepTools/plotPCA_matrix.tsv

#plotCoverage -p $num_threads --ignoreDuplicates --bamfiles $read_files_mapped --outRawCounts $out_dir/data/reads_mapped/_deepTools/plotCoverage_rawCounts_woDuplicates.tsv > $out_dir/data/reads_mapped/_deepTools/plotCoverage_output.tsv

#bamPEFragmentSize -p $num_threads --bamfiles $read_files_mapped --table $out_dir/data/reads_mapped/_deepTools/bamPEFragment_table.tsv --outRawFragmentLengths $out_dir/data/reads_mapped/_deepTools/bamPEFragment_rawLength.tsv

#estimateReadFiltering -p $num_threads --smartLabels --bamfiles $read_files_mapped > $out_dir/data/reads_mapped/_deepTools/estimateReadFiltering_output.tsv

### TODO needs BED files and gtf ?
###plotEnrichment -p $num_threads --smartLabels --bamfiles $read_files_mapped --outRawCounts $out_dir/data/reads_mapped/_deepTools/plotEnrichment_rawCounts.tsv

#multiqc -f -o $out_dir/data/quality_reports/reads_raw $out_dir/data/reads_raw
#multiqc -f -o $out_dir/data/quality_reports/reads_prepro $out_dir/data/reads_prepro
#multiqc -f -o $out_dir/data/quality_reports/reads_mapped $out_dir/data/reads_mapped















