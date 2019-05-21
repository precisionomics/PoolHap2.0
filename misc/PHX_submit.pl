#!/usr/bin/perl

my @problem_0 = ("0_1", "3_4");
my @problem_1 = ("2_4", "3_1");

for (my $p = 0; $p < scalar(@problem_0); $p++){
	system("sbatch -J $problem_0[$p]_rerun -o $problem_0[$p]_rerun PHX_run.pl $problem_0[$p]");
}

for (my $p = 0; $p < scalar(@problem_1); $p++){
	system("sbatch -J $problem_1[$p]_rerun -o $problem_1[$p]_rerun PHX_run.problem_1 $problem_0[$p]");
}