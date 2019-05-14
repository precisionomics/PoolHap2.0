#!/usr/bin/perl
#SBATCH --mem=40gb	 		# Job memory request

my $master_dir = $ARGV[0];
my $own_pdir = $ARGV[1];
my $java_cmd = "/gpfs/home/lmak/programs/jre1.8.0_181/bin/java";
my $gatk_cmd = "/gpfs/common/programs/gatk-4.0.5.1/gatk-package-4.0.5.1-local.jar";
my $sys_pdir = "/gpfs/common/programs"; 
my $input_dir = "$master_dir/input";
my $inter_dir = "$master_dir/intermediate";
my $gs_dir = "$master_dir/gold_standard";
my $out_dir = "$master_dir/output";
my $prefix = $ARGV[2];
my $pools = $ARGV[3];

my $ref =  "$input_dir/HIV_HXB2.fa";
my $gVCFs; 
for (my $p = 0; $p < $pools; $p++) {
	$gVCFs .= "-V $inter_dir/$prefix\_p$p.raw.g.vcf "; # With BQSR, $pfile.final.g.vcf
}

# print STDOUT $gVCFs;

# Step 4) Join all of the pool-specific gVCFs into a joint gVCF file usign Combine GVCFs.
print STDOUT "Step 4a) Join all of the pool-specific gVCFs into a joint gVCF file usign Combine GVCFs.\n\n";
system("$java_cmd -jar $gatk_cmd CombineGVCFs -R $ref $gVCFs -O $inter_dir/$prefix.g.vcf");

# Step 4b) Convert the joint gVCF into a VCF file using GenotypeGVCFs.
print STDOUT "\nStep 4b) Convert the joint gVCF into a VCF file using GenotypeGVCFs.\n\n";
system("$java_cmd -jar $gatk_cmd GenotypeGVCFs -R $ref -V $inter_dir/$prefix.g.vcf -ploidy 150 -O $inter_dir/$prefix.raw.vcf"); 
system("$java_cmd -jar $gatk_cmd SelectVariants -R $ref -V $inter_dir/$prefix.raw.vcf -O $inter_dir/$prefix.vcf"); 

# Step 5) Convert each pool-specific BAM file to SAM (i.e.: text), then VEF files. 
print STDOUT "\nStep 5) Convert each pool-specific BAM file to SAM (i.e.: text), then VEF files, then paired-read VEF files.\n\n";
for (my $p = 0; $p < $pools; $p++) { 
	my $pfile = "$inter_dir/$prefix\_p$p";
	system("$sys_pdir/samtools-0.1.19 view -ho $pfile.srt.sam $pfile.srt.bam");
}
system("$java_cmd -jar $own_pdir/BF3.jar $inter_dir/ $prefix $pools");
for (my $p = 0; $p < $pools; $p++) {
    my $pfile = "$inter_dir/$prefix\_p$p";
    system("$java_cmd -jar $own_pdir/PRL.jar $pfile");
}

# Step 6) Reconstructing haplotypes from multi-pool VEF files.
print STDOUT "\nStep 6) Reconstructing haplotypes from multi-pool VEF files...\n\n";
system("$java_cmd -jar $own_pdir/PHX8.2.jar $input_dir/PHX.properties $prefix $pools");

# Step 7) Comparing reconstructed haplotypes to the gold-standard information.
print STDOUT "\nStep 7) Comparing reconstructed haplotypes to the gold-standard information...\n\n";
system("$java_cmd -jar $own_pdir/O2R3.jar $gs_dir $out_dir $prefix $pools 0.05");
system("$java_cmd -jar $own_pdir/O2R3.jar $gs_dir $out_dir $prefix $pools 0.10");
