#!/usr/bin/perl

my $c = $ARGV[0];
my $s = $ARGV[1];
my $pools = $ARGV[2];
my $c0 = $ARGV[3];
my $c1 = $ARGV[4];

print STDOUT "Start of combination $c simulation $s.\n";
my $master_dir = "/gpfs/home/lmak/PHX_NonPerfect_Data3";
my $input_dir = "$master_dir/input";
my $gs_dir = "$master_dir/gold_standard";
my $inter_dir = "$master_dir/intermediate";
my $out_dir = "$master_dir/output";
my $own_pdir = "$master_dir/programs";
my $java_cmd = "/gpfs/home/lmak/programs/jre1.8.0_181/bin/java";
# Cannot use the default system JRE because even though commandline Java is showing version 8, Slurm is accessing a different default Java which is likely <8. 
my $prefix = "$c\_$s"; 
my $ref = "$input_dir/HIV_HXB2.fa";

# Step 1) Generate coalescence-simulated haplotypes, distribute to each of the pools, and simulate reads for each pool.

print STDOUT "Step 1) Generating coalescence-simulated haplotypes, distribute to each of the pools, and simulate reads for each pool.\n\n";
# print STDOUT "$java_cmd -jar $own_pdir/FS3.jar $input_dir/FS3.properties $prefix $c0 $c1";
system("$java_cmd -jar $own_pdir/FS3.jar $input_dir/FS3.properties $prefix $c0 $c1");


print STDOUT "\nSteps 2 and 3) For each pool, align the simulated reads to a reference sequence, and call variants using GATK HaplotypeCaller in gVCF mode.\n";
my $step3 = 0;
my $array_end = $pools - 1;
my $sub3 = `sbatch --array=0-$array_end -J NPD3_$prefix -o NPD3_$prefix NPD_Step3.sh $inter_dir/$prefix $ref`;
$sub3 =~ /^Submitted batch job (\d+)/; 
$step3 = $1;
print STDOUT "The array job ID was $step3.";

system("sbatch --dependency=afterok:$step3 -J NPD4_$prefix -o NPD4_$prefix NPD_Step4.pl $master_dir $own_pdir $prefix $pools");

