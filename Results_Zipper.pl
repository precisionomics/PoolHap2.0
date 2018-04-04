#!/usr/bin/perl

my @mutationRate = (0.00000125, 0.0000025, 0.000025);
my @hapsPerPat = (5, 10, 100);
my $zipNames; 

for (my $m = 0; $m < scalar(@mutationRate); $m++) {
	for (my $h = 0; $h < scalar(@hapsPerPat); $h++) {
		system("zip $mutationRate[$m]\_$hapsPerPat[$h].zip runthrough_v0.2_$mutationRate[$m]\_$hapsPerPat[$h]/simhaps.c.txt runthrough_v0.2_$mutationRate[$m]\_$hapsPerPat[$h]/simhaps.mutations.txt runthrough_v0.2_$mutationRate[$m]\_$hapsPerPat[$h]/p.all.vcf runthrough_v0.2_$mutationRate[$m]\_$hapsPerPat[$h]/p.all.results runthrough_v0.2_$mutationRate[$m]\_$hapsPerPat[$h]/p*.fa");
		$zipNames .= "$mutationRate[$m]\_$hapsPerPat[$h].zip ";
	}
}

system("zip MultiTest_results.zip $zipNames");