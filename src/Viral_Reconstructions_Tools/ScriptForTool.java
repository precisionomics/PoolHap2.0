package Viral_Reconstructions_Tools;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;

import org.apache.hadoop.hive.metastore.api.ThriftHiveMetastore.Processor.partition_name_has_valid_characters;
import org.datanucleus.store.rdbms.identifier.IdentifierFactory;

import spire.optional.intervalGeometricPartialOrder;

public class ScriptForTool {
	
	int[] pools = new  int [] {50, 200};
//	double[] mut_rates = new  double [] { 1e-6, 1e-7, 1e-8, 1e-9};
	int[] depths = new  int  [] { 500, 1000, 2000};
//	int[] num_haps = new  int  [] { 15, 30};
	
	
	String slim_script; 
	
	int genome_len = 9719;
	String genome_path ="/export/home/jhe/project/"
			+ "Viral_reconstruction/SLiM/SLiM/Reference/HIV_HXB2.fa"; 
	
	
	public ScriptForTool(String prefix_folder, String phx_folder_prefix, String tool_name) throws IOException {
		
		new File(prefix_folder + "/cmd/").mkdir();
		for (int i =0; i< this.pools.length; i++) {
			for (int j =0; j< this.depths.length; j++) {
				for (int k=1;k<8;k++) {
				
				
				String project_name = "pool_"+ Integer.toString(this.pools[i])
				+ "_dep_"+ Integer.toString(this.depths[j])+"_"+k;
				
				new File(prefix_folder + "/"+project_name ).mkdir();
//				# O2R.properties
//				##########
//				# General:
//				Output_Dir = /export/home/jhe/project/Viral_reconstruction/PredictHaplo/testing/negative_fitness/50_pool/non_migration/50_loci/500_cov/0_3_p0/
//				Gold-Standard_Dir = /export/home/jhe/project/Viral_reconstruction/SLiM/SLiM/negative_fitness/50_pool/non_migration/50_loci/500_cov/gold_standard/
//				Final_Output_Dir = /export/home/jhe/project/Viral_reconstruction/PredictHaplo/testing/negative_fitness/50_pool/non_migration/50_loci/500_cov/0_3_p0/
//				Out_File_Name = 0_3_p0.out
//				Fasta_File_Name = 0_3_p0
//				Proj_Name = 0_3
//				Tool_Name = PredictHaplo
//				Ref_Seq_Path =  /export/home/jhe/project/Viral_reconstruction/SLiM/SLiM/Reference/HIV_HXB2.fa
//				Cutoff_Lowest_Freq = 0.0000001
				for(int p=0;p<this.pools[i];p++) {
					String pool_name = project_name +"_p"+p;
					new File(prefix_folder + "/"+ project_name +"/"+pool_name).mkdir();
					String pool_folder = prefix_folder + "/"+ project_name +"/"+pool_name+"/";
					BufferedWriter bw1 = new BufferedWriter(new FileWriter(pool_folder + pool_name+"_O2R.properties"));
					
					if(tool_name.equals("CliqueSNV")) {
						bw1.write("Output_Dir = "+ prefix_folder + "/"+ project_name +"/snv_output/"+"\n");
					}else {
						bw1.write("Output_Dir = "+ pool_folder+"\n");
					}

					bw1.write("Gold-Standard_Dir ="+ phx_folder_prefix +"/"+ project_name+ "/gold_standard/\n");
					bw1.write("Final_Output_Dir ="+ pool_folder+"\n");
					if(tool_name.equals("CliqueSNV")) {
						bw1.write("Out_File_Name ="+ pool_name+".rg.fasta");
					}else if(tool_name.equals("TenSQR")) {
						bw1.write("Out_File_Name ="+ pool_name+"_ViralSeq.txt");
					}else if(tool_name.equals("PredictHaplo")) {
						bw1.write("Out_File_Name ="+ pool_name+".out");
					}else if(tool_name.equals("QuasiRecomb")) {
						bw1.write("Out_File_Name ="+ "quasispecies.fasta");
					}else {
						System.out.println("This tools are not in the list. Please choose "
								+ "from CliqueSNV,QuasiRecomb,RegressHaplo,PredictHaplo "
								+ "and TenSQR.");
					}
					bw1.write("Fasta_File_Name = "+ pool_name+"\n");
					bw1.write("Proj_Name ="+ project_name+"\n");
					bw1.write("Tool_Name ="+ tool_name+"\n");
					bw1.write("Ref_Seq_Path ="+ this.genome_path+"\n");
					bw1.write("Cutoff_Lowest_Freq ="+ prefix_folder + "/pool_"+ Integer.toString(this.pools[i])
					+ "_dep_"+ Integer.toString(this.depths[j])+"_"+k+ "/gold_standard\n");
					bw1.close();
				}
				
				BufferedWriter bw = new BufferedWriter(new FileWriter(prefix_folder + "/cmd/"+project_name+".cmd"));
//				#!/bin/sh
//				#SBATCH --job-name=Qpne50_500c
//				#SBATCH --workdir=/export/home/jhe/project/Viral_reconstruction/QuasiRecomb/output/SLiM/negative_fitness/50_pool/non_migration/50_loci/500_cov
//				#SBATCH --error=Quasi.error
//				#SBATCH --output=Quasi.out
//				#SBATCH --array=0-49
//				#SBATCH --time=1-0
//				#SBATCH --mem=30gb
//				#SBATCH --ntasks=1
//				#SBATCH --cpus-per-task=48
//				#SBATCH --time=99-00:00:00
//				#SBATCH --nodes=1
//				#SBATCH --exclude=node[029-033]
				bw.write("#!/bin/bash\n");
				bw.write("#SBATCH --job-name="+ tool_name + project_name+"\n");
				bw.write("#SBATCH --workdir="+ prefix_folder + "/"+ project_name +"\n");
				bw.write("#SBATCH --error=p_"+Integer.toString(this.pools[i])
				+ "_d_"+ Integer.toString(this.depths[j])+"_"+k+".error\n" );
				bw.write("#SBATCH --output=p_"+Integer.toString(this.pools[i])
				+ "_d_"+ Integer.toString(this.depths[j])+"_"+k+".out\n" );
				bw.write("##SBATCH --mem=30gb\n");
				bw.write("##SBATCH --ntasks=1\n");
				bw.write("##SBATCH --cpus-per-task=8\n");
				bw.write("##SBATCH --time=99-00:00:00\n");
				bw.write("##SBATCH --nodes=1\n");
				
				bw.write("java=/export/home/jhe/download/java_jdk_8u201/jdk1.8.0_201/bin/java\n");
				bw.write("java7=/export/home/jhe/download/jdk1.7.0_80/bin/java\n");
				bw.write("QuasiRecomb=/export/home/jhe/project/Viral_reconstruction/QuasiRecomb/QuasiRecomb.jar\n");
				bw.write("python=/export/home/jhe/download/anaconda/anaconda3/envs/tensqr_env/bin/python\n");
				bw.write("ExtractMatrix=/export/home/jhe/download/TenSQR/TenSQR-master/ExtractMatrix\n");
				bw.write("bwa=/export/home/jhe/download/bwa-0.7.17/bwa\n");
				bw.write("ref_dir=/export/home/jhe/project/Viral_reconstruction/SLiM/SLiM/Reference\n");
				bw.write("ref="+ this.genome_path+"\n");
				bw.write("PredictHaplo_Paired=/export/home/jhe/project/Viral_reconstruction/PredictHaplo/PredictHaplo-Paired-0.4/PredictHaplo-Paired\n");
				
				bw.write("start=$SECONDS");
				for (int p=0;p < this.pools[i];p ++) {
					String pool_name = project_name +"_p"+p;
					String pool_folder = prefix_folder + "/"+ project_name +"/"+pool_name+"/";
					bw.write("inbam="+ phx_folder_prefix + "/"+ project_name+"/input/bam/"
						+ pool_name +".rg.bam\n");
					bw.write("fastq_fastq="+ phx_folder_prefix + "/"+ project_name+"/input/fastq/"
							+ pool_name);
					if(tool_name.equals("QuasiRecomb")) {
						bw.write("$java7 -XX:+UseParallelGC -XX:NewRatio=9 -Xms10G "
								+ "-Xmx20G -jar $QuasiRecomb -i "
								+ "$inbam -conservative");
					}else if(tool_name.equals("TenSQR")){
						BufferedWriter bw_config = new BufferedWriter(new FileWriter(pool_folder+pool_name+".config"));
						bw_config.write("filename of reference sequence (FASTA) :"+ this.genome_path+"\n");
						bw_config.write("filname of the aligned reads (sam format) :"+ pool_name+".sam\n");
						bw_config.write("SNV_thres : 0.01\n");
						bw_config.write("reconstruction_start : 1\n");
						bw_config.write("reconstruction_stop: 9719\n");
						bw_config.write("min_mapping_qual : 40\n");
						bw_config.write("min_read_length : 100\n");
						bw_config.write("max_insert_length : 400\n");
						bw_config.write("characteristic zone name :"+pool_name+"\n");
						bw_config.write("seq_err (assumed sequencing error rate(%)) : 0.1\n");
						bw_config.write("MEC improvement threshold : 0.0312\n");
						bw_config.write("initial population size :100\n");
						
//						bwa mem $ref ${fastq_prefix}.bwa.read1.fastq ${fastq_prefix}.bwa.read2.fastq > ${curr_prefix}.sam
//						$ExtractMatrix ${curr_prefix}.config
//						$python /export/home/jhe/download/TenSQR/TenSQR-master/TenSQR.py ${curr_prefix}.config

					}else if(tool_name.equals("CliqueSNV")) {
//						$java -Xmx20G -jar /export/home/jhe/project/Viral_reconstruction/CliqueSNV/clique-snv.jar -m snv-illumina -tf 0.000000001 -in ${bam_dir}${bam_prefix}.rg.bam -log

					}else if(tool_name.equals("PredictHaplo")) {
						BufferedWriter bw_config = new BufferedWriter(new FileWriter(pool_folder+pool_name+".config"));
						
//						% configuration file for the HIVhaplotyper
//						% prefix
//						0_1_p0_
//						% filename of reference sequence (FASTA)
//						/export/home/jhe/project/Viral_reconstruction/SLiM/SLiM/Reference/HIV_HXB2.fa
//						% do_visualize (1 = true, 0 = false)
//						0
//						% filname of the aligned reads (sam format)
//						0_1_p0.sam
//						% have_true_haplotypes  (1 = true, 0 = false)
//						0
//						% filname of the true haplotypes (MSA in FASTA format) (fill in any dummy filename if there is no "true" haplotypes)
//						dump
//						% configuration file for the HIVhaplotyper
//						% prefix
//						0_1_p0_
//						% filename of reference sequence (FASTA)
//						/export/home/jhe/project/Viral_reconstruction/SLiM/SLiM/Reference/HIV_HXB2.fa
//						% do_visualize (1 = true, 0 = false)
//						0
//						% filname of the aligned reads (sam format)
//						0_1_p0.sam
//						% have_true_haplotypes  (1 = true, 0 = false)
//						0
//						% filname of the true haplotypes (MSA in FASTA format) (fill in any dummy filename if there is no "true" haplotypes)
//						dump
//						% do_local_analysis  (1 = true, 0 = false) (must be 1 in the first run)
//						1
//						% max_reads_in_window;
//						5000
//						% entropy_threshold
//						4e-4
//						%reconstruction_start
//						0
//						%reconstruction_stop
//						9719
//						%min_mapping_qual
//						0
//						%min_readlength
//						50
//						%max_gap_fraction (relative to alignment length)
//						0.01
//						%min_align_score_fraction (relative to read length)
//						0.10
//						%alpha_MN_local (prior parameter for multinomial tables over the nucleotides)
//						20
//						%min_overlap_factor (reads must have an overlap with the local reconstruction window of at least this factor times the window size)
//						0.5
//						%local_window_size_factor (size of  local reconstruction window relative to the median of the read lengths)
//						0.7
//						% max number of clusters (in the truncated Dirichlet process)
//						200
//						% MCMC iterations
//						501
//						% include deletions (0 = no, 1 = yes)
//						0
						
						
//						sed 's/rep1/'"${curr_prefix}"'/g' /export/home/jhe/project/Viral_reconstruction/PredictHaplo/testing/SLiM/config_5V > ${curr_prefix}.config
//						bwa mem $ref ${fastq_prefix}.bwa.read1.fastq ${fastq_prefix}.bwa.read2.fastq > ${curr_prefix}.sam
//						$PredictHaplo_Paired ${curr_prefix}.config > ${curr_prefix}.out
					}
				}//end of p_for_loop
				bw.write("end=$SECONDS\n"); 
				bw.write("echo \"duration: $((end-start)) seconds.\"");
				bw.write("\n");
// For QuasiRecomb
//				echo '$1 = ' $1  ## for example: 10
//
//
//				bam_dir="/export/home/jhe/project/Viral_reconstruction/SLiM/SLiM/negative_fitness/50_pool/non_migration/50_loci/500_cov/input/"$1"/bam/"
//				curr_pool=$((${SLURM_ARRAY_TASK_ID}%50))
//				bam_prefix=""$1"_p${curr_pool}"
//				mkdir ${bam_prefix}
//				cd ${bam_prefix}
//
//				start=$SECONDS
//				java7=/export/home/jhe/download/jdk1.7.0_80/bin/java
//				$java7 -XX:+UseParallelGC -XX:NewRatio=9 -Xms10G -Xmx20G -jar /export/home/jhe/project/Viral_reconstruction/QuasiRecomb/QuasiRecomb.jar -i ${bam_dir}${bam_prefix}.rg.bam -conservative
//				end=$SECONDS
//				echo "duration: $((end-start)) seconds.
				
				bw.write("$java -jar $poolsim "+ 	prefix_folder + "/cmd/pool_"+ 
				Integer.toString(this.pools[i])+ "_dep_"+ 
						Integer.toString(this.depths[j])+"_"+k +".PoolSimulator.properties\n");
				
				
				for (int p=0;p < this.pools[i];p ++) {
					
					bw.write("prefix="+ prefix_folder + "/pool_"+ Integer.toString(this.pools[i])+ 
							"_dep_"+ Integer.toString(this.depths[j])+"_"+k+"/input/fastq/"
							+ "pool_"+ Integer.toString(this.pools[i])+ 
							"_dep_"+ Integer.toString(this.depths[j])+"_"+k+
							"_p"+ Integer.toString(p) +"\n");
					
					
					bw.write("prefix_bam="+ prefix_folder + "/pool_"+ Integer.toString(this.pools[i])+ 
								"_dep_"+ Integer.toString(this.depths[j])+"_"+k+"/input/bam/"
								+ "pool_"+ Integer.toString(this.pools[i])+ 
								"_dep_"+ Integer.toString(this.depths[j])+"_"+k+
								"_p"+ Integer.toString(p) +"\n");

// Step 2: For each pool, align the simulated reads to a reference sequence.	
					
					bw.write("gunzip   $prefix\\.bwa.read1.fastq\n");
					bw.write("gunzip   $prefix\\.bwa.read2.fastq\n");
					bw.write("$bwa mem $ref $prefix\\.bwa.read1.fastq $prefix\\.bwa.read2.fastq "
							+ "| samtools view -Shub - > $prefix_bam\\.bam\n");
					bw.write("samtools sort -o  " +	"$prefix_bam\\.srt.bam  $prefix_bam\\.bam\n");
					
					
// Step 3: For each pool, call variants using GATK HaplotypeCaller in gVCF mode.	
					

								
					bw.write("$java -jar $gatk  AddOrReplaceReadGroups -I  $prefix_bam\\.srt.bam -O $inbam"
							+ " -R $ref -ID " +   "pool_"+ Integer.toString(this.pools[i])+ 
							"_dep_"+ Integer.toString(this.depths[j])+"_"+k+ "_p"+ Integer.toString(p)
							+" -LB NPD -PL Illumina -PU NPD -SM pool_"+ Integer.toString(this.pools[i])+ 
							"_dep_"+ Integer.toString(this.depths[j])+"_"+k+ "_p"+ Integer.toString(p)+ "\n");
					
					bw.write("samtools index $inbam\n");
										
					bw.write("prefix_vcf="+ prefix_folder + "/pool_"+ Integer.toString(this.pools[i])+ 
							"_dep_"+ Integer.toString(this.depths[j])+"_"+k+"/input/vcf/"
							+ "pool_"+ Integer.toString(this.pools[i])+ 
							"_dep_"+ Integer.toString(this.depths[j])+"_"+k+
							"_p"+ Integer.toString(p) +"\n");
					
					
					bw.write("outgvcf=$prefix_vcf\\.raw.g.vcf\n");
					bw.write("$java -jar $gatk  HaplotypeCaller -R  $ref -I $inbam"
							+ " -ERC  GVCF -ploidy 8 --heterozygosity 0.01  --max-alternate-alleles 1 -O "
							+ "$outgvcf\n");
				}
//	Step 4: Join all pool-specific gVCFs into a joint gVCF file and convert to VCF.
				
				bw.write("prefix="+ prefix_folder + "/pool_"+ Integer.toString(this.pools[i])+ 
						"_dep_"+ Integer.toString(this.depths[j])+"_"+k+"/input/vcf/"
						+ "pool_"+ Integer.toString(this.pools[i])+ 
						"_dep_"+ Integer.toString(this.depths[j])+"_"+k +"\n");
				
				String tmp= "$java -jar $gatk CombineGVCFs -R $ref ";
				for (int p=0;p < this.pools[i];p ++) {
					tmp=tmp+" -V " + "$prefix\\"+
							"_p"+ Integer.toString(p)+".raw.g.vcf ";
				}
				
				
				tmp=tmp+ " -O   $prefix\\.g.vcf\n"  ;
				bw.write(tmp );
				
				bw.write("$java -jar $gatk   GenotypeGVCFs -R $ref -V  $prefix\\.g.vcf" + 
						" -ploidy 8 -O $prefix\\.raw.vcf\n");
				
				bw.write("prefix_vcf="+ prefix_folder + "/pool_"+ Integer.toString(this.pools[i])+ 
						"_dep_"+ Integer.toString(this.depths[j])+"_"+k+"/input/"
						+ "pool_"+ Integer.toString(this.pools[i])+ 
						"_dep_"+ Integer.toString(this.depths[j])+"_"+k +"\n");
				
				bw.write( "$java -jar $gatk  SelectVariants -R $ref -V  $prefix\\.raw.vcf" + 
						" -O $prefix_vcf\\.vcf\n");
				
//	Step 5. Convert each pool-specific BAM file to SAM (i.e.: text), then VEF files.	
				for (int p=0;p < this.pools[i];p ++) {
					
					bw.write("prefix_sam="+prefix_folder + "/pool_"+ Integer.toString(this.pools[i])+ 
							"_dep_"+ Integer.toString(this.depths[j])+"_"+k+"/input/sam/"
							+ "pool_"+ Integer.toString(this.pools[i])+ 
							"_dep_"+ Integer.toString(this.depths[j])+"_"+k+
							"_p"+ Integer.toString(p) +"\n");
					
					bw.write("prefix_bam="+prefix_folder + "/pool_"+ Integer.toString(this.pools[i])+ 
							"_dep_"+ Integer.toString(this.depths[j])+"_"+k+"/input/bam/"
							+ "pool_"+ Integer.toString(this.pools[i])+ 
							"_dep_"+ Integer.toString(this.depths[j])+"_"+k+
							"_p"+ Integer.toString(p) +"\n");
					
					bw.write("samtools view -ho $prefix_sam\\.sam $prefix_bam\\.srt.bam\n");
					
				}
				
				new File(prefix_folder + "/pool_"+ Integer.toString(this.pools[i])+ 
							"_dep_"+ Integer.toString(this.depths[j])+"_"+k+"/input/").mkdir();         
				
				BufferedWriter bw_properties = new BufferedWriter(new FileWriter(
					prefix_folder + "/pool_"+ Integer.toString(this.pools[i])+ 
					"_dep_"+ Integer.toString(this.depths[j])+"_"+k+"/input/PHX.properties"));
				
				bw_properties.write("Proj_Name = "+ "pool_"+ Integer.toString(this.pools[i])+ 
						"_dep_"+ Integer.toString(this.depths[j])+"_"+k+"\n" );

				bw_properties.write("Input_Dir = "+prefix_folder + "/pool_"+ Integer.toString(this.pools[i])+ 
						"_dep_"+ Integer.toString(this.depths[j])+"_"+k +"/input/\n" );
				
				bw_properties.write("Intermediate_Dir = "+prefix_folder + "/pool_"+ Integer.toString(this.pools[i])+ 
						"_dep_"+ Integer.toString(this.depths[j])+"_"+k +"/intermediate/\n" );
				
				bw_properties.write("Output_Dir = "+prefix_folder + "/pool_"+ Integer.toString(this.pools[i])+ 
						"_dep_"+ Integer.toString(this.depths[j])+"_"+k +"/output/\n" );
				
				bw_properties.write("Gold_Dir = "+prefix_folder + "/pool_"+ Integer.toString(this.pools[i])+ 
						"_dep_"+ Integer.toString(this.depths[j])+"_"+k +"/gold_standard/\n" );
				
				bw_properties.write("Num_Pos_Window = 1000\n");
				bw_properties.write("Num_Gap_Window = 2\n");
				bw_properties.write("Num_Pos_Job = 1000\n");
																
				bw_properties.write("In-pool_Gap_Support_Min = 1\n");
				bw_properties.write("All-pool_Gap_Support_Min = 1\n");
				
				bw_properties.write("Level_1_Region_Size_Min = 10\n");
				bw_properties.write("Level_1_Region_Size_Max = 12\n");
				bw_properties.write("Level_2_Region_Size_Min = 10\n");
				bw_properties.write("Level_2_Region_Size_Max = 12\n");
				
				bw_properties.write("Est_Ind_PerPool = 1000000\n");
				bw_properties.write("Level_3_4_Region_Mismatch_Tolerance = 1\n");
				bw_properties.write("Level_5_6_Region_Mismatch_Tolerance = 2\n");
				bw_properties.write("Level_7_8_Region_Mismatch_Tolerance = 5\n");
				bw_properties.write("AEM_Maximum_Level = 7\n");
				bw_properties.write("BFS_Mismatch_Tolerance = 8\n");
				bw_properties.write("AEM_Iterations_Max = 200\n");
				bw_properties.write("AEM_Convergence_Cutoff = 0.0\n");
				bw_properties.write("AEM_Zero_Cutoff = 0.0\n");
				bw_properties.write("AEM_Regional_Cross_Pool_Freq_Cutoff = 0.0\n");
				bw_properties.write("AEM_Regional_HapSetSize_Max = 50\n");
				bw_properties.write("AEM_Regional_HapSetSize_Min = 3\n");
				bw_properties.write("Virtual_Coverage_Link_GC = 1500\n");
				
				bw_properties.write("Hc_Similarity_Cutoff =0.95\n");
				bw_properties.write("MCC_Freq_Cutoff = 0.01\n");
				bw_properties.write("Rscript_path = /home/jingni.he1/anaconda3/envs/Regress_Haplo/bin/Rscript\n");
				bw_properties.write("Regression_Distance_Max_Weight = 2.5\n");
				bw_properties.write("Regression_Coverage_Weight = 1.0\n");
				bw_properties.write("Regression_Mismatch_Tolerance= 7\n");
				
				bw_properties.write("Regression_One_Vector_Weight = 5.0 \n");
				bw_properties.write("Regression_Hap_VC_Weight = 2.0 \n");
				bw_properties.write("Regression_Hap_11_Weight = 1.0 \n");
				
				bw_properties.write("Regression_Regional_HapSetSize_Max = 35\n");
				bw_properties.write("Regression_Regional_HapSetSize_Min = 40\n");
				bw_properties.write("Regression_Gamma_Min = 0.01\n");
				bw_properties.write("Regression_n_Gamma = 10\n");
				bw_properties.write("Regression_Gamma_Max = 0.2\n");
				bw_properties.write("Regression_Maximum_Regions = 3\n");
				bw_properties.write("Maximum_Selected_HapSet = 15\n");
				bw_properties.write("Sequencing_Technology = paired-end reads\n");
				bw_properties.write("Number_Threads = 5\n");
															
				bw_properties.close();
				
				bw.write("properties="+prefix_folder + "/pool_"+ Integer.toString(this.pools[i])+ 
				"_dep_"+ Integer.toString(this.depths[j])+"_"+k+"/input/PHX.properties\n");
				bw.write("/home/jingni.he1/project/Viral_reconstruction/SLiM/programs/PoolHapX.jar\n");
				bw.write("start=$SECONDS\n");
				bw.write("$java -jar $poolhapx format $properties\n");
				bw.write("$java -jar $rewrite_vars "+ prefix_folder + "/pool_"
						+ Integer.toString(this.pools[i]) + "_dep_"
						+ Integer.toString(this.depths[j])+"_"+k+ "/gold_standard/"
						+ "pool_"+ Integer.toString(this.pools[i])+ 
						"_dep_"+ Integer.toString(this.depths[j])+"_"+k
						+"_haps.inter_freq_vars.txt "+prefix_folder + "/pool_"
						+ Integer.toString(this.pools[i]) + "_dep_"
						+ Integer.toString(this.depths[j])+"_"+k+ "/intermediate/"
						+ "pool_"+ Integer.toString(this.pools[i])+ 
						"_dep_"+ Integer.toString(this.depths[j])+"_"+k
						+"_vars.intra_freq.txt\n");
				bw.write("$java -jar $poolhapx gc $properties\n");
				bw.write("$java -jar $poolhapx aem $properties\n");
				bw.write("$java -jar $poolhapx evaluate $properties\n");
				bw.write("end=$SECONDS\n" + 
						"echo \"duration: $((end-start)) seconds.\"");
				bw.write("\n");
				
	        	bw.close();
	        	
				}
			}
		}
	}

	public static void main(String[] args) throws IOException, InterruptedException {
		System.out.println("PoolHapX Comparison Simulation... ...");
		String prefix_folder= args[0];//"/export/qlong/chencao/Work/poolhapx/slim/sim/";
		String phx_folder_prefix = args[1];
		String tool_name = args[2];
		ScriptForTool cs = new ScriptForTool(prefix_folder,phx_folder_prefix,tool_name);
		System.out.println("Done, Enjoy!");
		
	}	
}