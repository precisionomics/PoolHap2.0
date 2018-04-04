#!/usr/bin/perl

my @mutationRate = (0.00000125, 0.0000025, 0.000025);
my @hapsPerPat = (5, 10, 100);

for (my $m = 0; $m < scalar(@mutationRate); $m++) {
	for (my $h = 0; $h < scalar(@hapsPerPat); $h++) {
		system("sed 's/PARAM0/$mutationRate[$m]/g; s/PARAM1/$hapsPerPat[$h]/g;' /home/lauren.mak/PoolHap_Testing/parameters_v0.2_multi > /home/lauren.mak/PoolHap_Testing/parameters_v0.2_$m\_$h");
		system("LD_LIBRARY_PATH=/home/lauren.mak/programs/minisat-master/lib");
		system("export LD_LIBRARY_PATH");
		system("echo \"/home/lauren.mak/PoolHap_Testing/Pipeline_RandSim_v0.2.pl /home/lauren.mak/PoolHap_Testing/parameters_v0.2_$m\_$h\" | qsub -N PoolHapMultiTest -l walltime=24:00:00 -l mem=32GB -v LD_LIBRARY_PATH -k oe -m abe -M lauren.mak@ucalgary.ca");
	}
}