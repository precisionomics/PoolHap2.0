#!/usr/bin/env python
#SBATCH --mem=10gb	 		# Job memory request

import sys
import datetime
import subprocess
import fileinput

master_dir = '/gpfs/home/lmak/PHX_Perfect_Data/'
input_dir = master_dir + 'input/'
program_dir = master_dir + 'programs/'
gs_dir = master_dir + 'gold_standard/'
inter_dir = master_dir + 'intermediate/'
out_dir = master_dir + 'output/'
java_cmd = '/gpfs/home/lmak/programs/jre1.8.0_181/bin/java'
s = str(sys.argv[1])
p = str(sys.argv[2])

# Step 1) Generate perfect paired-end read VEF files from coalescence-simulated haplotypes. 
# @input: input_dir/FS.properties, input_dir/ref_seq
# @output: inter_dir/*.ms.txt, gs_dir/gs.inter_freq_vars.txt, gs_dir/gs.intra_freq.txt (haps and vars), inter_dir/*.bwa.read1.fastq.gz, inter_dir/*.vef
	# (inter_dir is work_dir in *.properties)
print('Start of simulation ' + s + ': ' + str(datetime.datetime.now()))
subprocess.call([java_cmd, '-jar', program_dir + 'FS2.jar', input_dir + 'FS2.properties', s])
print('Step 1 of testing ' + s + ' finished at ' + str(datetime.datetime.now()))

# Step 2) Solve full-length haplotypes and their global frequencies from VEF information.
# @input: input_dir/PHX.properties, gs_dir/gs.intra_freq.txt (vars), inter_dir/*.vef
# @output: inter_dir/p.in.list, inter_dir/*.in, inter_dir/dc_plan_file.txt, inter_dir/RH.inter_freq_vars.txt (later, intra_freq.txt as well).
print('Start of reconstruction ' + s + ': ' + str(datetime.datetime.now()))
subprocess.call([java_cmd, '-jar', program_dir + 'PHX8.1.jar', input_dir + 'PHX.properties', s, p])
print('Step 2 of testing ' + s + ' finished at ' + str(datetime.datetime.now()))

# Step 3) Compare the accuracy of the variant composition and frequencies with the original data. 
# @input: inter_dir/RH.inter_freq_vars.txt, gs_dir/gs.inter_freq_vars.txt, gs_dir/gs.intra_freq.txt (vars, maybe?)
# @output: inter_dir/full_results.txt, out_dir/sim_results.txt
subprocess.call([java_cmd, '-jar', program_dir + 'O2R2.jar', gs_dir, out_dir, s, p])
print('Step 3 of testing ' + s + ' finished at ' + str(datetime.datetime.now()))