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

## Input 
  | Type      | Description     |
  |-----------|---------------|
  | bam folder    | Folder containing the bam files on which you want to run GATK variant calling (bam files and their associated index (bai files)) 

## Parameters

  * #### Mandatory
| Name      | Example value | Description     |
|-----------|---------------|-----------------| 
| --ref    |            PATH TO FILE | Path to genome of reference |
| --bed    |            PATH TO FILE | Path to file with target regions |
| --known_hapmap    |            PATH TO FILE | Path to GATK Bundle: hapmap_3.3.hg19.sites.vcf |
| --known_omni    |            PATH TO FILE | Path to GATK Bundle: 1000G_omni2.5.hg19.sites.vcf |
| --known_1000G    |            PATH TO FILE | Path to GATK Bundle: 1000G_phase1.snps.high_confidence.hg19.sites.vcf |
| --known_snps    |            PATH TO FILE | Path to GATK Bundle: dbsnp_138.hg19.vcf |
| --known_mills    |            PATH TO FILE | Path to GATK Bundle: Mills_and_1000G_gold_standard.indels.hg19.sites.vcf |

  
  * #### Optional
| Name      | Default value | Description     |
|-----------|---------------|-----------------| 
| --cpu   |            8 | Number of cpu to be allocated by the various processes |
| --mem    |            32G | Memory to be allocated to the various processes |

  * #### Flags
  
Flags are special parameters without value.

| Name      | Description     |
|-----------|-----------------| 
| --help    | Display help |
		   
## Usage 
  ```
  ...
  ```
  
## Output 
  | Type      | Description     |
  |-----------|---------------|
  | .gvcf    | GVCF File containing list of variants for each sample (input bam) |

## Detailed description
Cf. [GATK website](https://software.broadinstitute.org/gatk/best-practices/)

## Directed Acyclic Graph
To be generated

## Contributions

  | Name      | Email | Description     |
  |-----------|---------------|-----------------| 
  | VOEGELE Catherine    |            voegelec@iarc.fr | Developer|
  | DELHOMME Tiffany    |            delhommet@students.iarc.fr | Developer |
  
## References

Cf. [GATK website](https://software.broadinstitute.org/gatk/best-practices/)
