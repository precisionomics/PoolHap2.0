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

*#header information such as project name<br>*
*sample1.read1.fastq\tsample1..read2.fastq<br>*
*sample2.read1.fastq\tsample2.read2.fastq<br>*
*sample3.read1.fastq     sample3.read2.fastq<br>*
sample4.read1.fastq     sample4.read2.fastq<br>
sample5.read1.fastq     sample5.read2.fastq<br>
sample6.read1.fastq     sample6.read2.fastq<br>
sample7.read1.fastq     sample7.read2.fastq<br>
sample8.read1.fastq     sample8.read2.fastq<br>
sample9.read1.fastq     sample9.read2.fastq<br>
sample10.read1.fastq    sample10.read2.fastq<br>










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
