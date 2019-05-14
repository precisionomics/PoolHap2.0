#!/usr/bin/perl

my $c = 2; 
my $simulations = 5;
my $pools = 20;
my @seg_sites = (20, 60);
my @err_rates = (0, 0.001); 

for (my $p = 0; $p < scalar(@seg_sites); $p++) { # scalar(@seg_sites)
	for (my $e = 1; $e < scalar(@err_rates); $e++) { # scalar(@err_rates)
		for (my $s = 0; $s < $simulations; $s++){
			system("LD_LIBRARY_PATH=/home/lauren.mak/programs/minisat-master/lib");
			system("export LD_LIBRARY_PATH");
			print STDOUT "Combination $c simulation $s: The number of segregating sites is $seg_sites[$p], and the Illumina per-base error rate is $err_rates[$e].\n";
			system("sbatch -J NPD2_$c\_$s -o NPD2_$c\_$s NPD_Step2.pl $c $s $pools $seg_sites[$p] $err_rates[$e]");
		}
		$c++; 
	}
}
