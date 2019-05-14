#!/usr/bin/perl

my $c = 0; 
my $simulations = 10;
my $pools = 20;
my @seg_sites = (20, 60);
my @err_rates = (0, 0.001); 

# $1 = Master directory
# $2 = Prefix 
# $3 = Sequencing error
# $4 = Population size

for (my $p = 0; $p < scalar(@seg_sites); $p++) { # scalar(@seg_sites)
	for (my $e = 0; $e < scalar(@err_rates); $e++) { # scalar(@err_rates)
		for (my $s = 0; $s < $simulations; $s++){
			print STDOUT "Combination $c simulation $s: The number of segregating sites is $seg_sites[$p], and the Illumina per-base error rate is $err_rates[$e].\n";
			system("sbatch -J TSR4_$c\_$s -o TSR4_$c\_$s NPD_TenSQR.sh /gpfs/home/lmak/PHX_NonPerfect_Data4 $c\_$s $err_rates[$e] 1000000 $pools");
		}
		$c++; 
	}
}

