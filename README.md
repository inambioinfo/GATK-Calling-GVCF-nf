# GATK_Calling_GVCF
## Nextflow pipeline to run variant calling in GVCF mode on bam files, joint genotyping and variant recalibration (SNPs and indels) following GATK best practices

## Description

The whole pipeline is made up of 6 steps: 

		1. Variant calling on samples bam files in GVCF mode

		2. Joint genotyping

		3. & 4. SNP recalibration

		5. & 6. Indel recalibration

## Dependencies

1. This pipeline is based on [nextflow](https://www.nextflow.io). As we have several nextflow pipelines, we have centralized the common information in the [IARC-nf](https://github.com/IARCbioinfo/IARC-nf) repository. Please read it carefully as it contains essential information for the installation, basic usage and configuration of nextflow and our pipelines.
2. External software:
- Java
- GATK

Before using the pipeline, specify in your config file the paths to the following files and softwares:
  
		gatk = '/appli/GenomeAnalysisTK/GATK-3.4-0/'
		   
		dbsnp = '/appli/reference/GATK_Bundle/dbsnp_138.hg19.vcf'
		   
		Mills_indels = '/appli/reference/GATK_Bundle/Mills_and_1000G_gold_standard.indels.hg19.sites.vcf'
		   
		hapmap = '/appli57/reference/GATK_Bundle/hapmap_3.3.hg19.sites.vcf'
		   
		omni = '/appli57/reference/GATK_Bundle/1000G_omni2.5.hg19.sites.vcf'
		   
		ThousandG_snps = '/appli57/reference/GATK_Bundle/1000G_phase1.snps.high_confidence.hg19.sites.vcf'
		   
You can avoid installing all the external software by only installing Docker. See the [IARC-nf](https://github.com/IARCbioinfo/IARC-nf) repository for more information.

		   

