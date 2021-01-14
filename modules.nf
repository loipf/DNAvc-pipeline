

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
		tuple val("${sample_id}"), path("${sample_id}.g.vcf.gz"), path("${sample_id}.g.vcf.gz.tbi") , emit: global_vcf
		tuple val("${sample_id}"), path("${sample_id}.vcf.gz"), path("${sample_id}.vcf.gz.tbi") , emit: vcf
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








process PREPROCESS_READS { 
	tag "$sample_id"
	publishDir "$params.data_dir/reads_prepro", pattern:"*cutadapt_output.txt", mode: "copy", saveAs: { filename -> "${sample_id}/$filename" }
	stageInMode = 'copy'   // avoids permission denied error

	input:
		tuple val(sample_id), path(reads) 
		val num_threads
		path adapter_seq

	output:
		tuple val(sample_id), path("${sample_id}_prepro_1.fastq.gz"), path("${sample_id}_prepro_2.fastq.gz"), emit: reads_prepro
		path "${sample_id}_cutadapt_output.txt", emit: cutadapt

	shell:
	'''
	ADAPTER_5=$(cat !{adapter_seq} | sed -n 1p | cut -f 2)  # forward
	ADAPTER_3=$(cat !{adapter_seq} | sed -n 2p | cut -f 2)  # reverse

	cutadapt --cores=!{num_threads} --max-n 0.1 --discard-trimmed --pair-filter=any --minimum-length 10 -b $ADAPTER_5 -B $ADAPTER_3 -o !{sample_id}_prepro_1.fastq.gz -p !{sample_id}_prepro_2.fastq.gz !{reads} > !{sample_id}_cutadapt_output.txt

	'''
}



// not possible to run dynamically fastqc with same name
process FASTQC_READS_RAW { 
	tag "$sample_id"
	publishDir "$params.data_dir/reads_raw", mode: "copy", overwrite: false, saveAs: { filename -> "${sample_id}/$filename" }
	stageInMode = 'copy'   // avoids permission denied error

	input:
		tuple val(sample_id), path(reads) 
		val num_threads
		path adapter_seq

	output:
		path "*.zip", emit: reports
		path "*.html"

	shell:
	'''
	fastqc -a !{adapter_seq} -t !{num_threads} --noextract !{reads}
	'''
}



process FASTQC_READS_PREPRO { 
	tag "$sample_id"
	publishDir "$params.data_dir/reads_prepro", mode: "copy", overwrite: false, saveAs: { filename -> "${sample_id}/$filename" }

	input:
		tuple val(sample_id), path(reads_prepro) 
		val num_threads
		path adapter_seq

	output:
		path "*.zip", emit: reports
		path "*.html"

	shell:
	'''
	fastqc -a !{adapter_seq} -t !{num_threads} --noextract !{reads_prepro}
	'''
}



process CREATE_BWA_INDEX { 
	publishDir "$params.data_dir/bwa_index", mode: "copy"

	input:
		path reference_genome

	output:
		path "*.{amb,ann,bwt,pac,sa}", emit: bwa_index

	shell:
	'''
	bwa index !{reference_genome} 
	'''
}



process MAPPING_BWA { 
	tag "$sample_id"
	publishDir "$params.data_dir/reads_mapped", mode: 'copy', saveAs: { filename -> "${sample_id}/$filename" }

	input:
		tuple val(sample_id), path(reads_prepro) 
		val num_threads
		path reference_genome
		path bwa_index  // to ensure index is created


	output:
		path "${sample_id}.bam", emit: reads_mapped
		path "${sample_id}.bam.bai", emit: reads_mapped_index
		path "*", emit: all


	shell:
	'''
	bwa mem -Y -R "@RG\\tID:!{sample_id}\\tSM:!{sample_id}" -t !{num_threads} -K 100000000 !{reference_genome} !{reads_prepro} \
	| samtools view -@ !{num_threads} -h -b - \
	| samtools sort -n -@ !{num_threads} - \
	| samtools fixmate -m -@ !{num_threads} - - \
	| samtools sort -@ !{num_threads} - \
	| samtools markdup -@ !{num_threads} -f !{sample_id}_markdup_stats.txt - !{sample_id}.bam

	samtools index -b -@ !{num_threads} !{sample_id}.bam
	samtools stats -@ !{num_threads} !{sample_id}.bam > !{sample_id}_stats.txt

	'''
}




process DEEPTOOLS_ANALYSIS { 
	publishDir "$params.data_dir/reads_mapped/_deepTools", mode: 'copy'

	input:
		path reads_mapped
		path reads_mapped_index
		val num_threads

	output:
		path "*", emit:all


	shell:
	'''
	multiBamSummary bins -p !{num_threads} --smartLabels --bamfiles !{reads_mapped} -o multiBamSummary.npz

	plotCorrelation --corData multiBamSummary.npz --corMethod spearman --whatToPlot heatmap --outFileCorMatrix plotCorrelation_matrix.tsv

	plotPCA --corData multiBamSummary.npz --outFileNameData plotPCA_matrix.tsv

	plotCoverage -p !{num_threads} --ignoreDuplicates --smartLabels --bamfiles !{reads_mapped} --outRawCounts plotCoverage_rawCounts_woDuplicates.tsv > plotCoverage_output.tsv

	bamPEFragmentSize -p !{num_threads} --smartLabels --bamfiles !{reads_mapped} --table bamPEFragment_table.tsv --outRawFragmentLengths bamPEFragment_rawLength.tsv

	estimateReadFiltering -p !{num_threads} --smartLabels --bamfiles !{reads_mapped} > estimateReadFiltering_output.tsv

	'''
}






process MULTIQC_RAW { 
	publishDir "$params.data_dir/quality_reports", mode: "copy"

	input:
		path stat_files

	output:
		path "*"

	shell:
	'''
	multiqc -f -o reads_raw .
	'''
}


process MULTIQC_PREPRO { 
	publishDir "$params.data_dir/quality_reports", mode: "copy"

	input:
		path stat_files

	output:
		path "*"

	shell:
	'''
	multiqc -f -o reads_prepro .
	'''
}


process MULTIQC_MAPPED { 
	publishDir "$params.data_dir/quality_reports", mode: "copy"

	input:
		path stat_files

	output:
		path "*"

	shell:
	'''
	multiqc -f -o reads_mapped .
	'''
}


















