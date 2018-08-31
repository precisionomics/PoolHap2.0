#!/usr/bin/perl

# my @mutationRate = (0.00000125, 0.0000025, 0.000025);
# my @hapsPerPat = (5, 10);
my @rounds = (10, 100); 
my @iterations = (1000, 10000); 
my $replicates = 50; 
my $homedir = "/gpfs/home/lmak/v0.3_burnin"; 

=begin
for (my $rp = 0; $rp < $replicates; $rp++) {
	# system("sed 's/PARAM0/$hapsPerPat[$b]/g; s/PARAM1/$hapsPerPat[$h]/g; s/PARAM2/$r/g;' TODO/parameters_v0.3_multi > TODO/parameters_v0.3_$rounds[$rd]\_$iterations[$it]\_$rp");
	my $outdir = "$homedir/0_0_$rp"; 
	system("sed 's|PARAM0|0|g; s|PARAM1|0|g; s|PARAM2|$outdir|g;' $homedir/parameters_v0.3_multi > $homedir/parameters_v0.3_0_0_$rp");
	system("LD_LIBRARY_PATH=/gpfs/home/lmak/programs/minisat-master/lib");
	# system("PARAMETERS=$homedir/parameters_v0.3_$rd\_$it\_$rp");
	system("mkdir $outdir");
	# system("export LD_LIBRARY_PATH");
	# system("echo \"/home/lauren.mak/PoolHap_Testing/Pipeline_RandSim_v0.3.pl /home/lauren.mak/PoolHap_Testing/parameters_v0.3_$rounds[$rd]\_$iterations[$it]\_$rp\" | qsub -N PHX_BMB -l walltime=24:00:00 -l mem=32GB -v LD_LIBRARY_PATH -k oe -m abe -M lauren.mak@ucalgary.ca");
	system("sbatch --mem=9000M --export=LD_LIBRARY_PATH -D $outdir -o stdout.txt $homedir/Pipeline_RandSim_v0.3.pl $homedir/parameters_v0.3_0_0_$rp");
}
=cut
for (my $rd = 0; $rd < scalar(@rounds); $rd++) {
	for (my $it = 0; $it < scalar(@iterations); $it++) {
		for (my $rp = 0; $rp < $replicates; $rp++) {
			# system("sed 's/PARAM0/$hapsPerPat[$b]/g; s/PARAM1/$hapsPerPat[$h]/g; s/PARAM2/$r/g;' TODO/parameters_v0.3_multi > TODO/parameters_v0.3_$rounds[$rd]\_$iterations[$it]\_$rp");
			my $outdir = "$homedir/$rounds[$rd]\_$iterations[$it]\_$rp"; 
			system("sed 's|PARAM0|$rounds[$rd]|g; s|PARAM1|$iterations[$it]|g; s|PARAM2|$outdir|g;' $homedir/parameters_v0.3_multi > $homedir/parameters_v0.3_$rounds[$rd]\_$iterations[$it]\_$rp");
			# system("LD_LIBRARY_PATH=/gpfs/home/lmak/programs/minisat-master/lib");	// Have to set this manually when you start the script.
			# system("PARAMETERS=$homedir/parameters_v0.3_$rd\_$it\_$rp");
			system("mkdir $outdir");
			# system("export LD_LIBRARY_PATH");
			# system("echo \"/home/lauren.mak/PoolHap_Testing/Pipeline_RandSim_v0.3.pl /home/lauren.mak/PoolHap_Testing/parameters_v0.3_$rounds[$rd]\_$iterations[$it]\_$rp\" | qsub -N PHX_BMB -l walltime=24:00:00 -l mem=32GB -v LD_LIBRARY_PATH -k oe -m abe -M lauren.mak@ucalgary.ca");
			system("sbatch --mem=9000M --export=LD_LIBRARY_PATH -D $outdir -o stdout.txt $homedir/Pipeline_RandSim_v0.3.pl $homedir/parameters_v0.3_$rounds[$rd]\_$iterations[$it]\_$rp");
		}
	}
}

# /gpfs/home/lmak/programs/jre1.8.0_181/bin/java -jar /gpfs/home/lmak/programs/Orig2Recons.jar 50 /gpfs/home/lmak/v0.3_burnin/10_1000_ 10 0 0.05
# /gpfs/home/lmak/programs/jre1.8.0_181/bin/java -jar /gpfs/home/lmak/programs/Orig2Recons.jar 50 /gpfs/home/lmak/v0.3_burnin/10_10000_ 10 0 0.05
# /gpfs/home/lmak/programs/jre1.8.0_181/bin/java -jar /gpfs/home/lmak/programs/Orig2Recons.jar 50 /gpfs/home/lmak/v0.3_burnin/100_1000_ 10 0 0.05
# /gpfs/home/lmak/programs/jre1.8.0_181/bin/java -jar /gpfs/home/lmak/programs/Orig2Recons.jar 50 /gpfs/home/lmak/v0.3_burnin/100_10000_ 10 0 0.05