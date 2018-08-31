#!/usr/bin/perl


my @hapsPerPat = (5, 10);
my $replicates = 3; 
my $zipNames; 

for (my $m = 0; $m < scalar(@hapsPerPat); $m++) {
	for (my $h = 0; $h < scalar(@hapsPerPat); $h++) {
		for (my $r = 0; $r < $replicates; $r++) {
			system("zip $hapsPerPat[$m]\_$hapsPerPat[$h]\_$r.zip bmb_v0.3_$hapsPerPat[$m]\_$hapsPerPat[$h]\_$r/simhaps.vc.txt bmb_v0.3_$hapsPerPat[$m]\_$hapsPerPat[$h]\_$r/simhaps.mutations.txt bmb_v0.3_$hapsPerPat[$m]\_$hapsPerPat[$h]\_$r/p.all.vcf bmb_v0.3_$hapsPerPat[$m]\_$hapsPerPat[$h]\_$r/p.all.pos bmb_v0.3_$hapsPerPat[$m]\_$hapsPerPat[$h]\_$r/p.all.results bmb_v0.3_$hapsPerPat[$m]\_$hapsPerPat[$h]\_$r/p*.fa");
			$zipNames .= "$hapsPerPat[$m]\_$hapsPerPat[$h]\_$r.zip ";
		}
	}
}

system("zip BMB_PHX_results.zip $zipNames");