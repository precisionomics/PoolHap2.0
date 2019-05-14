#!/usr/bin/perl

my $combinations = 4; 
my @simulations = ( 
       [0,2,3,4],
       [0,1,2,3,4],
       [0,1,2,3],
       [0,2,3]
);

for (my $c = 0; $c < $combinations; $c++) {
	for (my $s = 0; $s < scalar(@{$simulations[$c]}); $s++){
		$sim = $simulations[$c][$s]; 
		system("/gpfs/home/lmak/programs/jre1.8.0_181/bin/java -jar /gpfs/home/lmak/PHX_NonPerfect_Data3/programs/O2R3.jar /gpfs/home/lmak/PHX_NonPerfect_Data3/gold_standard /gpfs/home/lmak/PHX_NonPerfect_Data3/output $c\_$sim 20");
	}
}
