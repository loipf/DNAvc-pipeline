
/* 
 * DNA-VC PIPELINE 
 * for mapped reads
 */


/* 
 * import modules 
 */

nextflow.enable.dsl=2

include { 
	INDEX_REFERENCE;
	VARIANT_CALLING;
	VARIANT_MERGING;
	GLNEXUS_BCF_TO_VCF;
	VARIANT_ANNOTATION;
	VARIANT_CALLING_STATS;
	MULTIQC_VCF
} from './modules.nf' 




/*
 * default parameters
 */ 

params.dev_samples = 3

params.project_dir	= "$projectDir"
params.reads_mapped_dir	= "$params.project_dir/data/reads_mapped" 

params.reads_mapped	= "$params.reads_mapped_dir/*/*.{bam,bam.bai}"
params.data_dir		= "$params.project_dir/data"
params.scripts_dir	= "$params.project_dir/scripts"


/*
 * other parameters
 */

params.num_threads		= 3
params.reference_genome	= "$params.project_dir/data/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz"




log.info """\
DNA-VC PIPELINE
===================================================
reads_mapped		: $params.reads_mapped
data_dir		: $params.data_dir
reference_genome	: $params.reference_genome

===================================================

"""



/* 
 * main pipeline logic
 */
workflow {
	channel_reads_mapped = Channel
			.fromFilePairs( params.reads_mapped )
			.ifEmpty { error "cannot find any reads matching: ${params.reads_mapped}" }
			.take( params.dev_samples )  // only consider a few files for debugging

	INDEX_REFERENCE(params.reference_genome)

	VARIANT_CALLING(channel_reads_mapped, params.num_threads, INDEX_REFERENCE.out.reference_genome)
	VARIANT_CALLING_STATS(VARIANT_CALLING.out.vcf, params.num_threads)
	MULTIQC_VCF(VARIANT_CALLING_STATS.out.vcf_stats.collect())
	
	VARIANT_MERGING(VARIANT_CALLING.out.global_vcf.collect(), params.num_threads)
	GLNEXUS_BCF_TO_VCF(VARIANT_MERGING.out.glnexus_bcf, params.num_threads)

	VARIANT_ANNOTATION(GLNEXUS_BCF_TO_VCF.out.pvcf_file, params.num_threads)

}




workflow.onComplete { 
	println ( workflow.success ? "\ndone! check the quality reports in --> $params.data_dir/quality_reports\n" : "oops .. something went wrong" ) } 

















