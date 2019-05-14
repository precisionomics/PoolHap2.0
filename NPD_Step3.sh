#!/bin/bash
#SBATCH --mem=40gb	 		# Job memory request

# $1 = $inter_dir/$prefix
# $2 = $input_dir/HIV_HXB2.fa

# Make the patient part of the prefix of the input BAM and output gVCF files from the job array ID. 
prefix="${1}_p$SLURM_ARRAY_TASK_ID"
inbam="${prefix}.rg.bam"
outgvcf="${prefix}.raw.g.vcf"

# echo $prefix
# echo $inbam
# echo $outgvcf
# echo "gunzip ${prefix}.bwa.read1.fastq.gz"
# echo "/gpfs/home/lmak/programs/jre1.8.0_181/bin/java -Xmx20g -XX:+UseConcMarkSweepGC -XX:ParallelGCThreads=4 -jar /gpfs/common/programs/gatk-4.0.5.1/gatk-package-4.0.5.1-local.jar HaplotypeCaller -R $2 -I ${inbam} -ERC GVCF -ploidy 150 --heterozygosity 0.1 --max-alternate-alleles 1 -O ${outgvcf}"

# Steps 2 and 3) For each pool, align the simulated reads to a reference sequence, and call variants using GATK HaplotypeCaller in gVCF mode.
gunzip ${prefix}.bwa.read1.fastq.gz
gunzip ${prefix}.bwa.read2.fastq.gz
/gpfs/common/programs/bwa-0.7.15 mem $2 ${prefix}.bwa.read1.fastq ${prefix}.bwa.read2.fastq | /gpfs/common/programs/samtools-0.1.19 view -Shub - > ${prefix}.bam
/gpfs/common/programs/samtools-0.1.19 sort ${prefix}.bam ${prefix}.srt
	
/gpfs/home/lmak/programs/jre1.8.0_181/bin/java -jar /gpfs/common/programs/gatk-4.0.5.1/gatk-package-4.0.5.1-local.jar AddOrReplaceReadGroups -I ${prefix}.srt.bam -O ${inbam} -ID ${SLURM_ARRAY_TASK_ID} -LB NPD -PL Illumina -PU NPD -SM ${SLURM_ARRAY_TASK_ID} 
/gpfs/common/programs/samtools-0.1.19 index ${inbam}

/gpfs/home/lmak/programs/jre1.8.0_181/bin/java -Xmx20g -XX:+UseConcMarkSweepGC -XX:ParallelGCThreads=4 -jar /gpfs/common/programs/gatk-4.0.5.1/gatk-package-4.0.5.1-local.jar HaplotypeCaller -R $2 -I ${inbam} -ERC GVCF -ploidy 150 --heterozygosity 0.01 --max-alternate-alleles 1 -O ${outgvcf}

# Edit the final gVCF files such that the sample name in the header line is $pool_num because otherwise, they will collide at the GenotypeGVCFs step.
# sed -i '/^#CHROM/ s/$/${job}/' ${outgvcf}
