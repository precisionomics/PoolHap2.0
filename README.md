### How to run PHX8:
Suggestion: Open the source files in Eclipse and refactor line 17 as needed. Currently, the observed alternate allele frequency data is taken from ***gold-standard*** files, instead of empirically observed frequencies from the variant-calling process. The observed frequencies file is currently drawn from the directory called gp.gs_dir. Update this directory as necessary, either directly in the code here or in the parameters file (or both). However, if this is acceptable (i.e.: you're testing with simulated 'perfect' data as well), then this version is fine and you can skip straight using PHX8.jar. 


_Setup and Running Steps:_
1. Make your own version of the PHX.properties file to tell PHX where the input VEF and observed allele frequencies are. 
2. The command-line/configuration arguments are as follows: args[0] = absolute/path/to/PHX.properties, args[1] =  prefix_string, args[2] = number of pools.
- For example, if prefix_string = "test", PHX will look for work_dir/test_p0.vef, gs_dir/test_vars.intra_freq.txt, etc.
3. Run PHX8.jar with Java 8 (will not work with Java 7). 20 pools (33 unique haplotypes) x 100 segregating sites took about 5 hours.


### How to simulate 'perfect' VEFs to test PHX performance with:

**Program:** FullSimulator2. Can be found in JAR form as well as in src/MiscFunctions as a .java file. 

**Input:** Path to FS2.properties, the string prefix. 

**Output:** inter_freq_vars (global haplotype frequency and variant composition), haps.intra_freq (in-pool haplotype frequencies), and vars.intra_freq (in-pool alternate allele frequencies) text files. Pool VEF files.  

**Prerequisites:** Installed copies of ms and DWGSIM. Java 1.8+. 


_Setup and Running Steps:_
1. Install ms and DWGSIM and add the full directory path to the FS2.properties text file. 
2. Create a working directory and place FS2.properties in it. Add the path of the working directory to FS2.properties. 
3. Put a reference sequence in the working directory and add its name to FS2.properties. In this repo, it's HIV_HXB2.fa.
4. Set the rest of the toggle-able parameters in FS2.properties. 
5. Run the FullSimulator2 JAR executable with the path to FS2.properties and a prefix you want to apply to every gold-standard and VEF file. All VEF and gold-standard files will be in the working directory. 
- Generating 93 haplotypes x 100 variant positions x 100 pools VEFs took about 17 minutes. The bulk of the time was spent simulating the paired-end reads using DWGSIM and/or gunzip-ing the outputs.
- If you want to make multiple sets of perfect data at once, then use the PHX_submit.py and PHX_run.py scripts with slurm. PHX_run.py is capable of generating 'perfect' data, running PHX8.jar on the perfect data, and evaluating the accuracy of the output. Like an sbatch array job, PHX_submit.py submits PHX_run.py many times. Currently, the last two steps are commented out as PHX8.jar is finishing tests.


### How to measure PHX reconstruction accuracy:
I have written a program called Orig2Recons2 (in the MiscFunctions project subfolder) that compares original haplotypes to the closest reconstructed haplotype. It is also available in JAR form as O2R2.jar. For each (simulated) dataset, O2R2 will print...
1. One full-detail multi-pool results file containing the following columns: the pool ID, the original haplotype ID, the index of the closest reconstructed haplotype, the number of variants different between the closest RH and the OH, the difference in frequency between the closest RH and OH, and the number of haplotypes globally that are close enough to be 'quasispecies' (here, the cutoff is 5% difference of the total variant length). 
2. One line of multi-pool-aggregated summary statistics about PHX performance across the dataset. It contains the following columns: The fraction of accurately reconstructed OH, the average reconstruction error per OH, the average difference between the frequencies of the closest RH and the OH, and the proportion of the population that was accurately constructed (sum(quasispecies frequencies)).
- If there is more than one dataset being evaluated at once, each evaluation will generate the one line of aggregated summary statistics into the same file for easy comparison across the multiple datasets.


_Setup and Running Steps:_
1. Run the Orig2Recons2 JAR executable with the directory of the gold-standard files, the directory of the output files, the desired prefix, and the number of pools. 
 

=== DEPRECATED ===

# PoolHapX and Auxiliary Functions 

========================================================================================

## How to use PHX (and BAMFormatter) as a standalone program:

### Purpose: 
This is useful if you just want to run PHX on your own data. 

### Info about current version:  
The current version of PHX is called PoolHap2.0.2.jar. The formatting program is called BAMFormatter.jar.

### Planned or suggested updates:  
- [ ] It is in an extremely rudimentary form right now, as in it can only reconstruct haplotypes in a limited subset of the whole genome accurately (and possibly not even that). So higher-level algorithm design ideas to expand its capabilities from single-gene to whole-genome reconstruction would be welcome!
- [ ] BAMFormatter is currently designed to read standard BAM and VCF files. If you know yours has some peculiarities (ex. the BAM CIGAR string has customized characters), then variations of the program are very much needed and welcome. 

### How to use it: 
Read through steps 7 (BAMFormatter) and 8 (PHX) in Pipeline_RandSim_v0.2.pl, and refer to the parameter_v0.2 file for quick descriptions (and default values) to run PHX with. 

========================================================================================

## How to use PHX with simulated sequences:

### Purpose: 
This is useful if you want to simulated 'smaller', more easily handled versions of your real sequences. In other words, sequences with the properties that you want (ex. no indels, high mutation rate, specified hotspots, etc.) that are under your control. 'Simulated' as I have used it here means random number of patients, random number of haplotypes, mutations with the assumption of 'neutral' evolution on each of the haplotypes.

### Info about current version:  
The current version of the fully-simulated pipeline is called Pipeline_RandSim_v0.2.pl, and the parameters file, parameters_v0.2. The automator script is called Pipeline_Submitter.pl, and the automator-friendly parameter file is called parameters_v0.2_multi. 

### Planned or suggested updates:  
- [ ] Rewriting in Python
- [ ] Changing out DWGSIM for ART (more popular NGS read simulation program)
- [ ] Writing more modules for it (ex. adding your own analysis programs)
- [ ] Cleaning up the order of the parameters file

### How to use it: 
1) Go through the script file to make sure that the system Perl, Java, Picard, Samtools BWA, and GATK (see exceptions below) are where you expect them to be. Alter the script and parameters file as necessary. The following programs may need to be downloaded and added/installed to your own programs directory.
ms: Download from the [UChicago Box page for ms](https://uchicago.app.box.com/s/l3e5uf13tikfjm7e1il1eujitlsjdx13) and follow the super-simple compilation instructions. I used rand2. The ms executable needs to be called by absolute pathname i.e.: does not install in /usr/local/bin. 
DWGSIM: Download from the [Github page for DWGSIM](https://github.com/nh13/DWGSIM/). Need to download samtools-0.1.19.tar.bz2, specifically the version from the [Sourceforge page](https://sourceforge.net/projects/samtools/?source=typ_redirect) separately and put all files in a directory called samtools the in DWGSIM-master directory because the pre-packaged samtools in DWGSIM does not unzip properly. The dwgsim executable needs to be called by absolute pathname i.e.: does not install in /usr/local/bin. 
GATK: For BaseRecalibrator and AnalyzeCovariates, I use the the system defaults (version <4). However, gVCF mode and GenotypeGVCFs seems to be rather buggy in version <4, so I've switched to calling HaplotypeCaller, CombineGVCFs, and GenotypeGVCFs in version 4 (using files made from version <4 in >4 doesn't work either). I have downloaded GATK version 4.0.2.1 from the [GATK Downloads page](https://software.broadinstitute.org/gatk/download/) and since it's written in Java, the executable can either be called as a standalone JAR file or, in the latest version, without a Java call at all (`./gatk`).The structure of the GATK call changes radically from version <4 to >4. The GATK JAR file/executable needs to be called by absolute pathname i.e.: does not install in /usr/local/bin. 
2) Go through your parameters file to make sure all of the parameters you want to pass to the above programs and PHX are as you want them (ex. mutation rate for DWGSIM mutation generator). 
3) Add the FullSimulator, BAMFormatter, and PHX JAR files are in your own program directory. chmod as necessary to make them runnable.
	NOTE: HaploWriter and PatientSim are deprecated, as their function have been incorporated into FullSimulator. Please use FullSimulator instead.
4) Add the script, the parameters file, and the reference sequence to your home directory. 
5) Get rid of all DOS newline characters just to be safe by running `perl -pi -e 's/\r\n/\n/g' $filename` on all three files. 
6) Almost ready to run the script! To run the MiniSAT program that is called in PHX (see below for more details), you need to add some libraries to the environment's PATH variable using the following two commands:
	`LD_LIBRARY_PATH=/home/lauren.mak/programs/minisat-master/lib`
	`export LD_LIBRARY_PATH`
7) Submit the script as in the following command to the queue/slurm scheduler. The following example is for queue only, so the formatting for slurm may be different. The queue output and error files will be separate and appear in your home directory (NOT the directory you submitted in), and all script/PHX-generated output and result files will be in the directory specified in the parameter file.
	`echo "/home/lauren.mak/PoolHap_Testing/Pipeline_RandSim_helix.pl /home/lauren.mak/PoolHap_Testing/parameters_helix" | qsub -N PoolHapTest -l walltime=24:00:00 -l mem=32GB -v LD_LIBRARY_PATH -k oe -m abe -M lauren.mak@ucalgary.ca`
	NOTE: Intermediary results (ex. the initial haplotypes guessed by graph-colouring) will appear in the qsub output file (\*.o\*). I am open to changing this if there is interest!
8) If you want to simulate conditions, use the format of Pipeline_Submitter.pl and alter the variable parameters inside to whatever you want. I have tried to simulate variable mutation rates and numbers of haplotypes to test out the optimal datasets for PHX. Get rid of DOS newline characters with `perl -pi -e 's/\r\n/\n/g' $filename`, and chmod. Run the Perl script with...
	`perl Pipeline_Submitter.pl`
...and it will write and auto-submit scripts to the queue based on the parameters you have specified. Note that the library path commands are included in the submitter script. There's another little script called Results_Zipper.pl that'll help gather up the results directories if you submitted a few parameters.

### Runtime: 
Submitted with the described parameters (3 kb region with <30 variants, 10 patients x about 10 haplotypes), the **data pre-processing time took 40 minutes** and **running PHX took 2 minutes**.

========================================================================================

## How to evaluate PHX's performance:

### Purpose: 
This is useful if you want to see how closely PHX's reconstructions match the original input files. It provides i) an extended report comparing the original (simulated) haplotypes and the PHX-reconstructed haplotypes in a text file and ii) a few summary statistics to STDOUT that wraps up how accurate the reconstruction was. ii) especially is useful for comparing the accuracy of PHX on different datasets, which is primarily how I'm using it.

### Info about current version:  
The program is called Orig2Recons. 

### Planned or suggested updates:  
- [ ] Orig2Recons is currently customized to take mutation records from DWGSIM or ms + FullSimulator.jar, so you'll have to adjust it to take your own mutation record files (ex. VCF).

### How to use it: 
Coming soon! Currently running analyses so O2R is actively being worked on and altered. 
