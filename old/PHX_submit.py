#!/usr/bin/env python
#SBATCH -J PHX_Submission	# The job name

import sys
import subprocess

for s in range(0, int(sys.argv[1])):
	subprocess.call(['sbatch', '-J', 'RunPD_' + str(s), '-o', 'RunPD_' + str(s), 'PHX_run.py', str(s), str(sys.argv[2])])