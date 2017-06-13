#! /usr/bin/env nextflow

// STILL IN DEVELOPMENT !!!

// IARC - C. VOEGELE - Last update 13 Jun 2017 // Version of GATK needs to be >= 3.4

// usage : ./GATK_calling_GVCF.nf --bam_folder BAM/ --mem 32 --fasta_ref hg19.fasta
// Input files should be formatted: .bam and .bam.bai

if (params.help) {
    log.info ''
    log.info '--------------------------------------------------'
    log.info 'NEXTFLOW: GVCF VARIANT CALLING FOLLOWING GATK BEST PRACTICES'
    log.info '--------------------------------------------------'
    log.info ''
    log.info 'Usage: '
    log.info 'nextflow run GATK_calling_GVCF.nf --bam_folder BAM/ --mem 32 --fasta_ref hg19.fasta'
    log.info ''
    log.info 'Mandatory arguments:'
    log.info '    --bam_folder          FOLDER                  Folder containing BAM files to be called.'
    log.info '    --genome_ref          FILE                    Reference fasta file (with index) (excepted if in your config).'
    log.info '    --known_hapmap        FILE                    GATK_Bundle: hapmap_3.3.hg19.sites.vcf"
    log.info '    --known_omni          FILE                    GATK_Bundle: 1000G_omni2.5.hg19.sites.vcf"
    log.info '    --known_1000G         FILE                    GATK_Bundle: 1000G_phase1.snps.high_confidence.hg19.sites.vcf"
    log.info '    --known_snps          FILE                    GATK_Bundle: dbsnp_138.hg19.vcf"
    log.info '    --known_mills         FILE                    GATK_Bundle: Mills_and_1000G_gold_standard.indels.hg19.sites.vcf"
    log.info '    --java17              FOLDER                  Path to java directory.'
    log.info '    --GenomeAnalysisTK    FILE                    GenomeAnalysisTK.jar file.'
    log.info 'Optional arguments:'
    log.info '    --reserved_cpu        INTEGER                 Number of cpu reserved by nextflow (default: 8).'
    log.info '    --used_cpu            INTEGER                 Number of cpu used by bwa mem and sambamba (default: 8).'
    log.info '    --mem                 INTEGER                 Size of memory allocated by nextflow (default: 32).'
    log.info '    --out_folder          STRING                  Output folder (default: results_GVCF_calling).'
    log.info '    --bed                 FILE                    Bed file provided to GATK HaplotypeCaller.'
    log.info '    --cluster_options	    STRING		              Specific options for cluster scheduler'
    log.info ''
    exit 1
}

params.reserved_cpu = 8
params.used_cpu = 8
params.mem = 32
params.out_folder="results_GVCF_calling"
fasta_ref = file(params.fasta_ref)

// Check folder with bam files + bai files & pair them

bams = Channel.fromPath( params.bam_folder+'/*.bam' )
              .ifEmpty { error "Cannot find any bam file in: ${params.bam_folder}" }
    		.map {  path -> [ path.name.replace("bam",""), path ] }	
bais = Channel.fromPath( params.bam_folder+'/*.bam.bai' )
    		.map {  path -> [ path.name.replace("bam.bai",""), path ] }		

b_pair = bams
	.phase (bais)
	.map { pair1, pair2 -> [ pair1[1], pair2[1] ] }

//gvcf_files_list = Channel.fromPath('output_gvcf/*.vcf').toList()
//gvcf_files_idx_list = Channel.fromPath('output_gvcf/.vcf.idx').toList()

////////// STEP 01 ################### Variant calling on samples in GVCF mode - By default: --genotyping_mode = DISCOVERY --min_base_quality_score (-mbq) = 10, -stand_call_conf = 10 (vs 30 for tumors)

process creation_gvcf {
	echo "creation_gvcf"
	tag { bam_tag }
	cpus 6
	memory '12GB'
  cpus params.reserved_cpu
  clusterOptions '-R "rusage[mem=' + params.mem + '000]" -M ' + params.mem + '000'  
	publishDir 'output_gvcf'
  // publishDir params.out_folder, mode: 'move'  

input:
	file pair from b_pair
output:
	set val(bam_tag), file("${bam_tag}_raw_calls.g.vcf") into gvcf_files, count_gvcf
 	set val(bam_tag), file("${bam_tag}_raw_calls.g.vcf.idx") into gvcf_idx_files
shell:
	bam_tag = pair[0].baseName
     	'''
	!{params.java17} -jar !{params.gatk}GenomeAnalysisTK.jar -T HaplotypeCaller -nct 8 -R !{params.genome_ref} -I !{pair[0]} --emitRefConfidence GVCF -L !{params.bed} -o !{bam_tag}_raw_calls.g.vcf
       '''
}

nb_files = count_gvcf.count().val
// command_list_gvcf = ls *.gvcf

////////// STEP DATA AGGREGATION IF NB OF SAMPLES > 200-300 samples (java -jar GenomeAnalysisTK.jar -T CombineGVCFs -R reference.fasta --variant sample1.g.vcf --variant sample2.g.vcf -o cohort.g.vcf)

////////// STEP 02 ################### Joint genotyping


process GVCF {
    set val(bam_tag), file("${bam_tag}_realigned_recal.bam"), file("${bam_tag}_realigned_recal.bai") from outputs_recalibration
    file fasta_ref
    file fasta_ref_fai
    file fasta_ref_dict

    output:
    file("${bam_tag}_raw_calls.g.vcf") into output_gvcf
    file("${bam_tag}_raw_calls.g.vcf.idx") into output_gvcf_idx

    shell:
    '''
    java -jar !{params.GenomeAnalysisTK} -T HaplotypeCaller -nct !{params.used_cpu} -R !{fasta_ref} -I !{bam_tag}_realigned_recal.bam --emitRefConfidence GVCF !{intervals_gvcf} -o !{bam_tag}_raw_calls.g.vcf
    '''
}


//////// Variant quality score recalibration => VQSLOD

//process VQSR {
//	echo "creation_vqsr"
//	cpus 12
//	memory '24GB'
//	clusterOptions = params.cluster_options
//input:
	//file("All_samples.gvcf")
//output:recalibrate_SNP.recal , recalibrate_SNP.tranches, recalibrate_SNP_plots.R
//shell:

		//SNP recal (be careful: do not use -an DP for exomes)
   //	'''
		//!{params.java17} -jar !{params.gatk}GenomeAnalysisTK.jar -T VariantRecalibrator -nt 12 -R !{params.genome_ref} -input !{All_samples.gvcf} -resource:hapmap,known=false,training=true,truth=true,prior=15.0 ![params.known_hapmap} -resource:omni,known=false,training=true,truth=true,prior=12.0 ![params.known_omni} -resource:1000G,known=false,training=true,truth=false,prior=10.0 ![params.known_1000G} -resource:dbsnp,known=true,training=false,truth=false,prior=2.0 ![params.known_snps} -an QD -an FS -an SOR -an MQ -an MQRankSum -an ReadPosRankSum -an InbreedingCoeff -mode SNP -tranche 100.0 -tranche 99.9 -tranche 99.0 -tranche 90.0 -recalFile recalibrate_SNP.recal -tranchesFile recalibrate_SNP.tranches -rscriptFile recalibrate_SNP_plots.R 
		//!{params.java17} -jar !{params.gatk}GenomeAnalysisTK.jar -T ApplyRecalibration -nt 12 -R !{params.genome_ref} -input !{All_samples.gvcf} -mode SNP --ts_filter_level 99.9 -recalFile !{recalibrate_SNP.recal} -tranchesFile !{recalibrate_SNP.tranches} -o All_samples_recalibrated_snps_raw_indels.vcf

		//Indel recal
		//!{params.java17} -jar !{params.gatk}GenomeAnalysisTK.jar -T VariantRecalibrator -nt 12 -R !{params.genome_ref} -input !{All_samples_recalibrated_snps_raw_indels.vcf} -resource:mills,known=true,training=true,truth=true,prior=12.0 ![params.known_mills} -an QD -an FS -an SOR -an MQRankSum -an ReadPosRankSum -an InbreedingCoeff -mode INDEL -tranche 100.0 -tranche 99.9 -tranche 99.0 -tranche 90.0 --maxGaussians 4 -recalFile recalibrate_INDEL.recal -tranchesFile recalibrate_INDEL.tranches -rscriptFile recalibrate_INDEL_plots.R
		//!{params.java17} -jar !{params.gatk}GenomeAnalysisTK.jar -T ApplyRecalibration -nt 12 -R !{params.genome_ref} -input !{All_samples_recalibrated_snps_raw_indels.vcf} -mode INDEL --ts_filter_level 99.9 -recalFile !{recalibrate_INDEL.recal} -tranchesFile !{recalibrate_INDEL.tranches} -o All_samples_recalibrated_snps_indels.vcf
   //	'''
 //}

// process creation_joint {
   	// echo "creation_joint"
	// tag { bam_tag }
	// cpus 24
	// memory '24GB'
	// clusterOptions = params.cluster_options
// input:
    	// set val(bam_tag), file("${bam_tag}_raw_calls.g.vcf") from gvcf_files.groupTuple(size: nb_files)
    	// set val(bam_tag), file("${bam_tag}_raw_calls.g.vcf.idx") from gvcf_idx_files.groupTuple(size: nb_files)
	// file gvcf_files_list
	// file gvcf_files_idx_list


// output:
	// file("All_samples.gvcf") into uniq_gvcf_file, uniq_gvcf_file2
// shell:
	// command_list_gvcf = gvcf_files_list.collect { f -> "-V '${f}'" }.join(' ')
   	// '''
	// !{params.java17} -jar !{params.gatk}GenomeAnalysisTK.jar -T GenotypeGVCFs -nt 24 -R !{params.genome_ref} !{command_list_gvcf} -o All_samples.gvcf
   	// '''
//}
