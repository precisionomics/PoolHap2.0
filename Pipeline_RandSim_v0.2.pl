#!/usr/bin/perl

# 1) Take parameters from the first command line argument.  
my @paramlist; 
open PARAMFILE, "< $ARGV[0]" or die "Error: Couldn't open parameter file $ARGV[0].\n";
while(<PARAMFILE>) {
	chomp $_;
	my @columns = split("\t");
	if (index($columns[0], "NEXT_PROG") != -1) {next};
	my $curr_param = $columns[2];
	$curr_param =~ s/\s+$//;
	push(@paramlist,$curr_param);
}
close PARAMFILE; 

print STDOUT "Goal: Reconstruct haplotypes and frequencies from randomly simulated patient haplotypes of the HIV HXB2 pol polyprotein gene (3012 bp).\n";
print STDOUT "Reading parameters from input file...\n";
my $hdir = $paramlist[0];
my $syspdir = $paramlist[1];
my $ownpdir = $paramlist[2];
my $refseq = $paramlist[12]; 
system("$syspdir/perl/perl-5.24.0/bin/perl -pi -e 's/\r\n/\n/g' $hdir/$refseq");
system("$syspdir/samtools/samtools-0.1.19/bin/samtools faidx $hdir/$refseq"); 
my $refname = (split /\./, "$refseq")[0]; 
system("$syspdir/java/jdk1.8.0_77/bin/java -jar $syspdir/picard/picard-tools-2.6.0/picard.jar CreateSequenceDictionary R=$hdir/$refseq O=$hdir/$refname.dict"); 
my $outdir = "$hdir/$paramlist[37]"; 
system("mkdir $outdir");
my $prefix = "$outdir/simhaps"; 
print STDOUT "The output directory is $outdir.\n\n";

my $sim_haps = $paramlist[33]; 
my $sim_pts = $paramlist[34]; 

# 2a) Simulate the set of all haplotype sequences and the set of haplotypes in each patient.
print STDOUT "Simulating a total pool of $sim_haps haplotypes from $refseq using ms...\n\n";
my $theta = 2 * $paramlist[20] * $paramlist[8]; # theta = 2 * N0 * mu
my $rho = 2 * $paramlist[20] * $paramlist[8] / 2; # theta = 2 * N0 * rate_recomb
system("$ownpdir/msdir/ms $sim_haps 1 -t $theta -L -r $rho $paramlist[35] > $prefix.vc.txt"); 
system("$syspdir/java/jdk1.8.0_77/bin/java -jar $ownpdir/FullSimulator.jar $hdir/$refseq $prefix.vc.txt $paramlist[38] $paramlist[39] $outdir $sim_pts $paramlist[40]");

my $gVCFs; 

# For all patient-specific FastQ files...

# my @biggVCFs = qw( 1 2 8 9 );	# To finish final gVCFs of very large raw gVCFs.
# foreach my $p (@biggVCFs) {		#

#	my $pfile = "$outdir/p$p";	#

for (my $p = 0; $p < $sim_pts; $p++) {

	my $pfile = "$outdir/p$p";

	# 2b) Simulate patient-specific FastQs using DWGSIM.
	print STDOUT "Making patient $p FASTQ files using DWGSIM...\n";
	system("$ownpdir/DWGSIM-master/dwgsim $pfile.fa $pfile -e $paramlist[3] -E $paramlist[4] -C $paramlist[5] -1 $paramlist[6] -2 $paramlist[7] -r 0 -F $paramlist[9] -H -o $paramlist[11]");

	# 3) Align the simulated reads to a reference sequence, and sort the BAM file.
	print STDOUT "\nMaking sorted patient $p BAM files using BWA and Samtools...\n";
	system("$syspdir/bwa/bwa-0.7.13/bwa index $hdir/$refseq");
	system("gunzip $pfile.bwa.read1.fastq.gz");
	system("gunzip $pfile.bwa.read2.fastq.gz");
	system("$syspdir/bwa/bwa-0.7.13/bwa mem -R '\@RG\tID:bwa\tLB:HIVlib\tPL:ILLUMINA\tSM:p$p\tPU:HWI' $hdir/$refseq $pfile.bwa.read1.fastq $pfile.bwa.read2.fastq | $syspdir/samtools/samtools-0.1.19/bin/samtools view -Shub - > $pfile.bam");
	system("$syspdir/samtools/samtools-0.1.19/bin/samtools sort $pfile.bam $pfile.srt");		
	
	# 4) Mark all duplicated reads in the sorted BAM file.
	print STDOUT "\nMarking duplicates in sorted patient $p BAM files using Picard...\n";
	system("$syspdir/java/jdk1.8.0_77/bin/java -jar $syspdir/picard/picard-tools-2.6.0/picard.jar MarkDuplicates I=$pfile.srt.bam O=$pfile.mdup.bam  M=$outdir/p$p_mdup_metrics.txt");
	system("$syspdir/samtools/samtools-0.1.19/bin/samtools index $pfile.mdup.bam");  

	# 5a) Call variants using GATK HaplotypeCaller in gVCF mode.
	print STDOUT "\nCalling variants to raw patient $p gVCF files using GATK HaplotypeCaller...\n";
	system("$syspdir/java/jdk1.8.0_77/bin/java -Xmx$paramlist[13] -XX:+UseConcMarkSweepGC -XX:ParallelGCThreads=$paramlist[15] -jar $ownpdir/gatk-4.0.2.1/gatk-package-4.0.2.1-local.jar HaplotypeCaller -R $hdir/$refseq -I $pfile.mdup.bam -ERC GVCF -ploidy $paramlist[16] --max-alternate-alleles $paramlist[17] --minimum-mapping-quality 0 -O $pfile.raw.g.vcf");

	# 5b) Adjust base quality scores using GATK BaseRecalibrator.
	print STDOUT "\nAdjusting patient $p base quality scores using GATK BaseRecalibrator and PrintReads...\n";
	system("$syspdir/java/jdk1.8.0_77/bin/java -jar $syspdir/gatk/gatk-3.6/GenomeAnalysisTK.jar -T BaseRecalibrator -R $hdir/$refseq -I $pfile.mdup.bam -knownSites $pfile.raw.g.vcf -o $pfile.recal.table");
	system("$syspdir/java/jdk1.8.0_77/bin/java -jar $syspdir/gatk/gatk-3.6/GenomeAnalysisTK.jar -T BaseRecalibrator -R $hdir/$refseq -I $pfile.mdup.bam -knownSites $pfile.raw.g.vcf -BQSR $pfile.recal.table -o $pfile.post.table");
	system("$syspdir/java/jdk1.8.0_77/bin/java -jar $syspdir/gatk/gatk-3.6/GenomeAnalysisTK.jar -T AnalyzeCovariates -R $hdir/$refseq -before $pfile.recal.table -after $pfile.post.table -plots $pfile.plots.pdf");
	system("$syspdir/java/jdk1.8.0_77/bin/java -jar $syspdir/gatk/gatk-3.6/GenomeAnalysisTK.jar -T PrintReads -R $hdir/$refseq -I $pfile.mdup.bam -BQSR $pfile.post.table -o $pfile.bqsr.bam"); 

	# 5c) Re-call variants using the updated base quality scores in gVCF mode.
	print STDOUT "\nRe-calling variants in final patient $p gVCF files using GATK HaplotypeCaller...\n\n";
	system("$syspdir/java/jdk1.8.0_77/bin/java -Xmx$paramlist[14] -XX:+UseConcMarkSweepGC -XX:ParallelGCThreads=$paramlist[15] -jar $ownpdir/gatk-4.0.2.1/gatk-package-4.0.2.1-local.jar HaplotypeCaller -R $hdir/$refseq -I $pfile.bqsr.bam -ERC GVCF -ploidy $paramlist[16] --max-alternate-alleles $paramlist[17] --minimum-mapping-quality 0 -O $pfile.final.g.vcf");

	# 5d) Edit the final gVCF files such that the sample name in the header line is p$patient_num because otherwise, they will collide at the GenotypeGVCFs step.
	system("sed -i '/^#CHROM/ s/$/$p/' $pfile.final.g.vcf"); 

	$gVCFs .= "-V $pfile.final.g.vcf "
}

# my $gVCFs; 								# To finish the rest of the script.

# for (my $p = 0; $p < $sim_pts; $p++) {	# 
#	my $pfile = "$outdir/p$p";			#
#	$gVCFs .= "-V $pfile.final.g.vcf ";	#
# }										#

# 6a) Join all of the patient-specific gVCFs into a GenomicsDB using GenomicsDBImport.
print STDOUT "Making a multi-pool gVCF for all $sim_pts final gVCFs using GATK CombineGVCFs...\n\n";
system("$syspdir/java/jdk1.8.0_77/bin/java -Xmx10g -XX:+UseConcMarkSweepGC -XX:ParallelGCThreads=$paramlist[15] -jar $ownpdir/gatk-4.0.2.1/gatk-package-4.0.2.1-local.jar CombineGVCFs -R $hdir/$refseq $gVCFs -O $outdir/p.all.g.vcf");

# 6b) Join all of the patient-specific gVCFs in the GenomicsDB using GenotypeGVCFs.
print STDOUT "\nMerging all patient gVCF files using GATK GenotypeGVCFs...\n\n";
system("$syspdir/java/jdk1.8.0_77/bin/java -jar $ownpdir/gatk-4.0.2.1/gatk-package-4.0.2.1-local.jar GenotypeGVCFs -R $hdir/$refseq -V $outdir/p.all.g.vcf -ploidy $paramlist[16] -new-qual -O $outdir/p.all.vcf"); 

# 7) Convert each patient-specific BAM file to SAM (i.e.: text), then VEF files. 
print STDOUT "Converting all patient BAM files to VEF files and assembling patient allele counts...\n\n";
for (my $p = 0; $p < $sim_pts; $p++) { 
	my $pfile = "$outdir/p$p";
	system("$syspdir/samtools/samtools-0.1.19/bin/samtools view -ho $pfile.bqsr.sam $pfile.bqsr.bam");
 	system("$syspdir/java/jdk1.8.0_77/bin/java -jar $ownpdir/BAMFormatter.jar $pfile $outdir/p.all.vcf $paramlist[6] $paramlist[35] $outdir");
}
my $VEFFILE = "$outdir/p.all.vef.list";
open(my $vw, '>', $VEFFILE) or die "Error: Couldn't open list of VEFs file $outdir/p.all.vef.list.\n";
for (my $p = 0; $p < $sim_pts; $p++) {
	print $vw "$outdir/p$p.vef\n"; 
}
close $vw;
close $VEFFILE;
system("cat $outdir/p\*.ct > $outdir/p.all.ct"); 

# 8) Run the input VEF file through the PoolHap (GC -> rjMCMC -> EM) algorithms. 
print STDOUT "Reconstructing haplotypes and frequencies using PoolHap2.0...\n\n";
system("$syspdir/java/jdk1.8.0_77/bin/java -jar $ownpdir/PoolHap2.0.2.jar $outdir $outdir/p.all.vef.list $outdir/p.all.ct $paramlist[18] $paramlist[19] $paramlist[20] $paramlist[21] $paramlist[22] $paramlist[23] $paramlist[24] $paramlist[25] $paramlist[26] $paramlist[27] $paramlist[28] $paramlist[29] $paramlist[30] $paramlist[31] $paramlist[32] $sim_pts $outdir/p.all.pos $paramlist[36]"); 