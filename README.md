## Users’ Manual of PoolHapX (Version 1.0)
### Preamble
The PoolHapX program reconstructs haplotypes within-host from pooled-sequencing data by integrating population genetic models (statistical linkage disequilibrium) with genomics reads (physical linkage). It approximate the resolution of single-cell sequencing using only pooled sequencing data, enabling within-host evolution analyses.

The workflow of PoolHapX is briefed as follows: (A) PoolHapX first determines locations
of physical linkage uncertainty using sequencing reads, and then divides the full
genome into smaller regions. (B) Regional haplotypes are solved for and joined together
using a statistical model for a parsimonious global distribution of haplotypes. C) The
within-pool frequency of each haplotype is estimated by regularized regression to solve
for each within-pool haplotype distribution.

Installation and a simple example are described below. Users can get the final within-
host (or within-pool) frequencies of each haplotype by running the functions step by step:
“script”, “format”, “gc”, “aem”, “l0l1”. ScriptForPHX.jar generates all commands required
by PoolHapX in a script so users can run the script easily. More description and
benchmarking of PoolHapX can be found in our publication:

URL of the bioariv.

### Installation
PoolHapX is a batteries-included JAR executable. All needed external jar packages are included in the downloadable, PoolHapX.jar. However, as we used an R package L0Learn, the users have to install R and L0Learn (https://cran.r-project.org/web/packages/L0Learn/index.html). The versions of R and R package L0Learn that we have used on our platform are: version 1.2.0 for L0Learn and version 3.6.1 for R. Other versions are not tested, although they may work. Users are also
expected to have java (version: 1.8) on their platform.Longranger (version: 2.2.2) should be installed if processing 10x linked-reads.

Several other tools are prerequisites for running. PoolHapX. Users can download and install them from the websites:
* bwa (if using paired-end reads): https://github.com/lh3/bwa
* samtools 0.1.19+: http://www.htslib.org/download/
* GATK 4.1: https://software.broadinstitute.org/gatk/download/index
* longranger (if using 10x linked-reads): https://support.10xgenomics.com/genome-exome/software/pipelines/latest/installation

### Functions
* script: the “script” function generates a script which contains all commands. Users can run the project_name.cmd to get the final results calculated by PoolHapX from the initial FASTQ files.

* format: the “format” function generates the vef file (vef file contains variant sites linking information extracted from BAM files) and calculates the variant frequency at different sites for different pools from the VCF/SAM files.

* gc: the “gc” function generates the graph coloring result from the vef file.

* aem: the “aem” function first divides the variation sites into several regions from the graph coloring result; after that the aem function infers haplotypes in local regions for different levels using hierarchical AEM algorithm.

* l0l1: the “l0l1” function finalizes the identity and frequency of the global haplotypes using L0L1 regulated regression.

### Quick start with included example data

Example data is provided. After decompressing the downloadable, users can see the reference folder, fastq files folder, fastq_file.txt, and config.properties under the “Example” folder. After updating absolute paths of executable (such as bwa etc) and parent folder in the config.properties file, users can run PoolHapX by a simple commands:

Usage:

`java -jar PoolHapX.jar script config.properties`

Then go to the folder of /PATH/TO/Working_dir/cmd and run:

`./<project_name>.cmd`

Users will then generate the final haplotype results for each pool at the “output” folder
under their working directory.

### Full Manual
### Data Preparation
#### Check fastq_name file format
Please put all FASTQ files under the same directory. All paired-end fastq files have to be named as: sample_id.read1.fastq and sample_id.read2.fastq. Write the name of all your fastq files into a single file following the format below:

`———————————————————————————————————————————————————`<br>
*#header information such as project name<br>*
*sample1.read1.fastq　　　sample1..read2.fastq<br>*
*sample2.read1.fastq　　　sample2.read2.fastq<br>*
*sample3.read1.fastq　　　sample3.read2.fastq<br>*
*sample4.read1.fastq　　　sample4.read2.fastq<br>*
*sample5.read1.fastq　　　sample5.read2.fastq<br>*
*sample6.read1.fastq　　　sample6.read2.fastq<br>*
*sample7.read1.fastq　　　sample7.read2.fastq<br>*
*sample8.read1.fastq　　　sample8.read2.fastq<br>*
*sample9.read1.fastq　　　sample9.read2.fastq<br>*
*sample10.read1.fastq　　 sample10.read2.fastq<br>*
`———————————————————————————————————————————————————`<br>
Each row is an observation(sample), and each name is separated by tab.

For 10x linked reads, the folder name for each pool should be listed in the file:<br>
`———————————————————————————————————————————————————`<br>
*#header information such as project name<br>*
*sample1<br>*
*sample2<br>*
*sample3<br>*
*sample4<br>*
*sample5<br>*
*sample6<br>*
*sample7<br>*
*sample8<br>*
*sample9<br>*
*sample10<br>*
`———————————————————————————————————————————————————`

#### Check config file format (configure to your setting)
`———————————————————————————————————————————————————`<br>
*#config_file
*Main_Dir = /PATH/TO/PoolHapX_work_dir<br>*
*Project_Name = Test<br>*
*Java = /PATH/TO/java<br>*
*samtools = /PATH/TO/samtools<br>*
*gatk = /PATH/TO/gatk<br>*
*PHX_JAR = /PATH/TO/PoolHapX.jar<br>*
*#If Sequencing_Technology is 10x_linked-reads, users may leave the parameter "bwa" blank.<br>*
*bwa = /PATH/TO/bwa<br>*
*#If Sequencing_Technology is paired-end_reads, users may leave the parameter "longranger" blank.<br>*
*longranger = /PATH/TO/longranger<br>*
*#/PATH/TO/longranger_ref_folde is generated by the "mkref" function of longranger. If Sequencing_Technology is paired-end_reads, users may leave the parameter "Longranger_Ref_Folder" blank.<br>*
*Longranger_Ref_Folder = /PATH/TO/longranger_ref_folder<br>*
*Fastq_Path = /PATH/TO/fastq_file_dir<br>*
*Fastq_File = /PATH/TO/fastq_file.txt<br>*
*Ref_Path = /PATH/TO/HIV_HXB2.fa<br>*
*Rscript_Path = /PATH/TO/Rscript<br>*
*Sequencing_Technology= 10x_linked-reads or paired-end_reads<br>*
`———————————————————————————————————————————————————`<br>
config file included in the package is configured to sample set

#### Generate the script

Command:

`> java -jar PoolHapX.jar script <config file>`

The script for running PoolHapX will be generated under “/PATH/TO/Working_dir/cmd/”. The name of the file will be “<project_name>.cmd”

#### Run PoolHapX

`> ./<project_name>.cmd`

Output format:

The output for each step, i.e., graph coloring, aem and l0l1 regression, is generated under “/PATH/TO/Working_dir/<project_name>/intermediate/”. The final output is generated under “/PATH/TO/Working_dir/<project_name>/output/”. Under the folder “output”, each pool (host) has a folder. Within its folder, one can find the haplotype frequencies in this pool the file named as “final_freq_haps.txt”. The format of the frequencies is exampled below:
`———————————————————————————————————————————————————————————————`<br>
*Hap_ID　 　　　h0　　　　h1　　　　h2　　　　h3　　　　h4　　　　...<br>*
*Freq 　 　 　　0.028 　　  0.022 　   　0.031　 　  0.061　　　    0.124　  　 ...<br>* 
*0;211;211;0:1　 0　　　  　　1　　 　　1　　 　　1　　 　　0　  　 　　 ...<br>* 
*0;231;231;0:1　 1　　　  　　0　　 　　0　　 　　0　　 　　0　  　 　　 ...<br>* 
*0;310;310;0:1　 0　　　  　　1　　 　　1　　 　　1　　 　　0　  　 　　 ...<br>* 
*... ...* 
`———————————————————————————————————————————————————————————————`<br>
The first row lists the haplotype IDs. The second row lists the frequencies of each haplotype. The first column denotes the ID of the genetic variants in the format of chromosome-ID; start-position; end-position; alleles. In the event of some viruses that have not chromosome number, PoolHapX will use 0 to denote the chromosome ID. For SNPs, the start-position and end-position are the same. In the rest of the file, each column represents the composition of the haplotypes, i.e., the alleles at each location.


#### PoolHapX properties file
We provide default parameters for users (“/PATH/TO/Working_dir/<project_name>/input/PHX.properties”). Under most
circumstances, the default parameters work well when compared with other existing tools. However, in the event that users may want to make change to the parameters themselves, the properties file is located under input directory, named as
“PHX.properties”. We list the meaning for all the parameters below (explanation of each parameters as well as their ranges). Users can modify these parameters according to their needs.

**#PoolHapX Parameters**<br>
##########<br>
##The name of the project, will be the prefix of names of cross-pool files.<br>
**Proj_Name**=project_name<br>

##File locations: input directory; output files directory; intermediate files directory; gold standard files directory.<br>
**Input_Dir**=/PATH/TO/Input_Dir<br>
**Intermediate_Dir**=/PATH/TO/Intermediate_Dir<br>
**Output_Dir**=/PATH/TO/Output_Dir<br>
##If users do not have the gold standard files, please just leave the parameter blank<br>
**Gold_Dir**=/home/chencao/Desktop/PoolHap/freq_0_pool_25_dep_100/gold_standard<br>

##########<br>
###Graph-Colouring: link all reads to generate candidate global haplotypes based on physical linkage.
##Maximum number of positions in a window. [Default 20: Range: 1 - 100]<br>
**Num_Pos_Window**=20<br>
##Maximum number of gaps in a window. [Default 2: Range: 1 – 20]<br>
**Num_Gap_Window**=2<br>

##########<br>
###Divide-and-Conquer: divide the genome into multiple regions based on linkage uncertainty.<br>
##Proportion of raw GC-haplotypes that contain the gap in the pool. [Default: 0.6, Range: 0 - 1]<br>
**In-pool_Gap_Support_Min**=1<br>
##Proportion of raw GC-haplotypes that contain the gap across all pools. [Default: 0.1, Range: 0 - 1]<br>
**All-pool_Gap_Support_Min**=1<br>
##Minimum number of SNPs in a Level 1 region. [Default: 10, Range: 8 - 12]<br>
**Level_1_Region_Size_Min**=10<br>
##Maximum number of SNPs in a Level 1 region. [Default: 14, Range: 10 - 14]<br>
**Level_1_Region_Size_Max**=12<br>
##Minimum number of SNPs in a Level 1 tiling region. [Default: 10, Range: 8 - 12]<br>
**Level_1T_Region_Size_Min**=10<br>
##Maximum number of SNPs in a Level 1 tiling region. [Default: 14, Range: 10 - 14]<br>
**Level_1T_Region_Size_Max**=12<br>
##Estimated number of individuals in a pool. [Default: 1000000, Range: 1000 - 1000000]<br>
**Est_Ind_PerPool**=1000000<br>
##Number of maximum mismatch positions in constructing Level 2. [Default: 1: Range: 0 - 2]<br>
**Level_1_2_Region_Mismatch_Tolerance**=1<br>
##Number of maximum mismatch positions in constructing Level 3. [Default 2: Range: 1 - 3]<br>
**Level_2_3_Region_Mismatch_Tolerance**=2<br>
##Number of maximum mismatch positions in constructing Level 4. [Default 2: Range: 3 - 7]<br>
**Level_3_4_Region_Mismatch_Tolerance**=5<br>
##Number of AEM levels. [Default 4: Range: 1 - 4]<br>
**AEM_Maximum_Level**=4<br>
##Number of maximum mismatch positions in BFS. [Default 6: Range: 4 - 8]<br>
**BFS_Mismatch_Tolerance**=6<br>

##########<br>
###Approximate Expectation-Maximization: generate regional haplotype sets and their frequencies.<br>
##Maximum number of iterations regardless of convergence. [Default: 200, Range: 50 - 1000]<br>
**AEM_Iterations_Max**=200<br>
##The epsilon that controls the stop criteria of AEM (i.e. convergence). [Default: 0.00001, Range: 0 - 0.000001]<br>
**AEM_Convergence_Cutoff**=0.00001<br>
##For each iteration of AEM, some very rare haplotypes with frequencies below this parameter will be set to a frequency of zero. [Default: 0.00001, Range: 0.0 - 0.000001]<br>
**AEM_Zero_Cutoff**=0.00001<br>
##Initial value for regional cross-pool frequency cutoff immediately after AEM. [Default: 0.01, Range: 0.0 - 0.05]<br>
**AEM_Regional_Cross_Pool_Freq_Cutoff**=0.01<br>
##Maximum number of regional haplotypes in a region for AEM. [Default: 50, Range: 1-200]<br>
**AEM_Regional_HapSetSize_Max**=50<br>
##Minimum number of regional haplotypes in a region for AEM. [Default: 5, Range: 1-20]<br>
**AEM_Regional_HapSetSize_Min**=3<br>
##If both denominator and numerator are very close to zero, the Importance Factor (IF) value. [Default: 5.0, Range: 1.0-10.0]<br>
**IF_0_0=0.1<br>
##if the denominator is close to zero but the numerator is not, the IF value. [Default: 50.0, Range: 10.0-1000.0]<br>
**IF_Denominator_0 = 10.0**<br>








-----------------------------------------------------------------------------------


This is specifically for PoolHapX developers who are working directly with the code to expand its applicability on different types of data. For all applications to real data, see https://github.com/theLongLab/PoolHapX. The `TenSQR_Testing/` directory contains all of the programs needed to convert the TenSQR output format into the PoolHapX standard output format, and the `external_jars` directory contains all executables (mostly for the LASSO regression part) needed to compile PoolHapX.jar.

## Getting Started

### Installation
Create a PoolHapX directory and download the source code from repository. 

### Prerequisites
Java 1.8+: https://www.java.com/en/download/

If there is no data-specific variant-caller needed:

* bwa: https://github.com/lh3/bwa

* samtools 0.1.19+: http://www.htslib.org/download/

* GATK 4+: https://software.broadinstitute.org/gatk/download/index

* ms: https://uchicago.app.box.com/s/l3e5uf13tikfjm7e1il1eujitlsjdx13

* DWGSIM: https://github.com/nh13/DWGSIM (**Note**: do NOT use the Bioconda version, that one is currently bugged)

* A job scheduler capable of running job arrays and job dependencies (ex. Slurm). Makes it easy to run simulations in parallel. 

## How to Use

**PoolHapX input:** Allele-annotated reads and observed allele frequencies. For examples, see 0_0_p0.vef and 0_0_vars.intra_freq.txt in sample_io_files/. 

**PoolHapX output:** Within-pool haplotype distributions (allelic compositions, within-pool and across-pool frequencies). For examples, see 0_0.inter_freq_vars.txt and 0_0.intra_freq.txt in sample_io_files/. 

The basic command for running PoolHapX. if the sample prefix is "prefix" and there are 10 pools...

`java -jar PoolHapX.jar PHX.properties prefix 10`

Because PoolHapX requires specific input in the variant-encoded format (VEF), raw sequencing data must be processed one of the subsequently described pipelines, or a variation thereof. 

#### What is a VEF file?

The input data for PoolHapX is a modified version of the BAM file, termed the variant encoded (VEF) file. Each line of the VEF contains information about allele at each segregating site in sequence region spanned by the read. For example, if a read called ReadName spanned positions 100 to 200 relative to the reference sequence, and carried the reference allele at position 50 and the alternate allele at position 150, the VEF file would contain the line...

`ReadName 50=0;150=1; // 100 200`

One VEF file is generated per pool. In terms of biological applications, a single sample could be the NGS reads of a patient’s HIV population, or a single microbiome sample.

### What is your input data? 

#### Option 1:
You have real data, in the form of raw sequencing reads (FASTQ) or aligned sequencing reads (BAM) and segregating sites called from data-specific variant-caller. Please see https://github.com/theLongLab/PoolHapX for all instructions.

#### Option 2:
You want to simulate 'perfect' sequencing data and evaluate PoolHapX's performance on a completely idealized dataset (perfect annotation of alleles on reads, perfectly observed alternate alternate allele frequencies). 

#### Option 3:
You want to simulate 'non-perfect' sequencing data and evaluate PoolHapX's performance on a partially idealized dataset (variant-calling-based annotation of alleles on reads, variant-caller-observed alternate alternate allele frequencies). 

### Running Options 2 or 3:
0. Find the option-specific executables in the similarly-named directories.
1. Make the subdirectories `input/`, `programs/`, `intermediate/`, `output/` and `gold_standard/` in the main working directory. These will correspond to the organizational structures referred to in all of the scripts.  
2. PHX.properties is a JSON-like file setting parameters for PoolHapX.jar. Update the absolute paths of the directories. For an explanation of each parameter, see below. 
3. Similarly, FS\*.properties sets the simulation parameters for the simulation wrapper program FullSimulator (FS2.jar for 'perfect' or FSFQ.jar for 'non-perfect'). Update the absolute paths of the directories and programs. Place this in `input/`. 
4. Place your reference sequence in `input/`. Index it with BWA, Samtools, and Picard.
5. Place FS2/3.jar, O2R2/3.jar, and PoolHapX.jar in the `programs/` directory. If you're going with **Option 3**, also place BF3.jar in the `programs/` directory. 
6. If you only want to make one set of simulation data, use \*PD_Step2-4.\*. If you want to make multiple sets of simulation data, use \*PD_Step1.\* to submit \*PD_Step2-4.\* in parallel. 
	* All scripts will all require some degree of customization since I have tailored them to work on my HIV data. 
7. The result files (allele composition and within/across-pool frequencies) will show up in `output/`. The gold-standard files, which are generated by FS2/FQ.jar, will be in `gold_standard/`. 
8. When all of your simulations and PoolHapX runs have finished, compare the the result files in `output/` to the answers in `gold_standard/`. The script **O2R3_submit.pl** automates this comparison across many simulated datasets. The command inside can easily be applied to single sets of data. 
9. Use the output of O2R2/3.jar to analyze reconstruction accuracy. It will output...
	* One full-detail multi-pool results file containing the following columns: the pool ID, the original haplotype ID, the index of the closest reconstructed haplotype, the number of variants different between the closest RH and the OH, the difference in frequency between the closest RH and OH, and the number of haplotypes globally that are close enough to be 'quasispecies'. 
	* One line of multi-pool-aggregated summary statistics about PoolHapX performance across the dataset. It contains the following columns: The fraction of accurately reconstructed OH, the average reconstruction error per OH, the average difference between the frequencies of the closest RH and the OH, and the proportion of the population that was accurately constructed (sum(quasispecies frequencies)).
		* If there is more than one dataset being evaluated at once, each evaluation will generate the one line of aggregated summary statistics into the same file for easy comparison across the multiple datasets.
