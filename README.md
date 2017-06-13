# GATK-Calling-GVCF-nf
Performs variant calling in GVCF mode on bam files, joint genotyping and variant recalibration (SNPs and indels) following GATK best practices

The whole pipeline is made up of 6 steps: 

		1. Variant calling on samples bam files in GVCF mode

		2. Joint genotyping

		3. & 4. SNP recalibration

		5. & 6. Indel recalibration


Before using the pipeline, specify in your config file the paths to the following files and softwares:

		   genome_ref = '/appli/reference/GATK_Bundle/ucsc.hg19.fasta'
		   
		   dbsnp = '/appli/reference/GATK_Bundle/dbsnp_138.hg19.vcf'
		   
		   Mills_indels = '/appli/reference/GATK_Bundle/Mills_and_1000G_gold_standard.indels.hg19.sites.vcf'
		   
		   java17 = 'java'
		   
		   gatk = '/appli/GenomeAnalysisTK/GATK-3.4-0/'
		   
		   hapmap = '/appli57/reference/GATK_Bundle/hapmap_3.3.hg19.sites.vcf'
		   
		   omni = '/appli57/reference/GATK_Bundle/1000G_omni2.5.hg19.sites.vcf'
		   
		   ThousandG_snps = '/appli57/reference/GATK_Bundle/1000G_phase1.snps.high_confidence.hg19.sites.vcf'
		   



		   

