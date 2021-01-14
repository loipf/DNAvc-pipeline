

// need params.data_dir declared in main.nf, not implemented in DSL2 yet
params.data_dir	= "$launchDir/data"



process INDEX_REFERENCE { 
	publishDir "$params.data_dir", mode: "copy", pattern:"oc_databases_version.txt"

	input:
		path reference_genome

	output:
		tuple path("*.fa"), path("*.fa.fai"), emit: reference_genome
		path "oc_databases_version.txt"

	shell:
	'''
	reference_name=!{reference_genome}
	reference_name=${reference_name%.*}  # removes last file extension .gz

	gunzip -c !{reference_genome} > $reference_name 
	samtools faidx $reference_name -o $reference_name.fai
	
	oc module ls -t annotator > oc_databases_version.txt  ### output OpenCRAVAT database versions
	'''
}



process VARIANT_CALLING { 
	container "google/deepvariant:1.1.0"
	tag "$sample_id"
	publishDir "$params.data_dir/variants_vcf", mode: "copy", pattern:"${sample_id}.{vcf.gz,vcf.gz.tbi,visual_report.html}", overwrite: false, saveAs: { filename -> "${sample_id}/$filename" }

	input:
		tuple val(sample_id), path(reads_mapped)
		val num_threads
		path reference_genome

	output:
		tuple val("${sample_id}"), path("${sample_id}.vcf.gz"), path("${sample_id}.vcf.gz.tbi") , emit: vcf
		path "${sample_id}.g.vcf.gz" , emit: global_vcf
		path "${sample_id}.visual_report.html"

	shell:
	'''
	/opt/deepvariant/bin/run_deepvariant \
		--model_type=WGS \
		--ref=!{reference_genome[0]} \
		--reads=!{reads_mapped[0]} \
		--output_vcf=!{sample_id}.vcf.gz \
		--output_gvcf=!{sample_id}.g.vcf.gz \
		--num_shards=!{num_threads} \
		--vcf_stats_report=true
	'''
}


process VARIANT_MERGING { 
	container "quay.io/mlin/glnexus:v1.2.7"
	publishDir "$params.data_dir/variants_vcf/_all", mode: "copy", pattern:"GLnexus.DB", overwrite: false

	input:
		path global_vcf_files
		val num_threads

	output:
		path "*"
		path "pvcf_all_glnexus.bcf", emit: glnexus_bcf

	shell:
	'''
	bash -c "glnexus_cli --threads !{num_threads} --config DeepVariant !{global_vcf_files}" > pvcf_all_glnexus.bcf
	'''
}


process GLNEXUS_BCF_TO_VCF { 
	publishDir "$params.data_dir/variants_vcf/_all", mode: "copy", overwrite: false 

	input:
		path bcf_file
		val num_threads

	output:
		path "pvcf_all_glnexus.vcf.gz", emit: pvcf_file

	shell:
	'''
	bcftools view !{bcf_file} | bgzip -@ !{num_threads} -c > pvcf_all_glnexus.vcf.gz
	'''
}


process VARIANT_CALLING_STATS { 
	tag "$sample_id"
	publishDir "$params.data_dir/variants_vcf", mode: "copy", overwrite: false, saveAs: { filename -> "${sample_id}/$filename" }

	input:
		tuple val(sample_id), path(vcf_file), path(vcf_file_index) 
		val num_threads

	output:
		path "${sample_id}_vcfstats.txt", emit: vcf_stats

	shell:
	'''
	bcftools stats -f PASS --threads !{num_threads} !{vcf_file} > !{sample_id}_vcfstats.txt
	'''
}


process MULTIQC_VCF { 
	publishDir "$params.data_dir/quality_reports", mode: "copy"

	input:
		path stat_files

	output:
		path "*"

	shell:
	'''
	multiqc -f -o variants_vcf .
	'''
}



process VARIANT_ANNOTATION { 
	publishDir "$params.data_dir/variants_vcf/_all", mode: "copy", overwrite: false

	input:
		path pvcf_glnexus
		val num_threads

	output:
		path "*"

	shell:
	'''
	oc run -l hg38 -t csv -x --mp !{num_threads} !{pvcf_glnexus}
	'''
}





















