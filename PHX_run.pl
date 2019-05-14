#!/usr/bin/perl
#SBATCH --mem=40gb	 		# Job memory request

system("/gpfs/home/lmak/programs/jre1.8.0_181/bin/java -jar /gpfs/home/lmak/PHX_NonPerfect_Data3/programs/PHX8.2.jar /gpfs/home/lmak/PHX_NonPerfect_Data3/input/PHX.properties $ARGV[0] 20");
