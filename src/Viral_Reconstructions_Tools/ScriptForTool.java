package Viral_Reconstructions_Tools;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;

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
						bw1.write("Out_File_Name ="+ pool_name+".rg.fasta\n");
					}else if(tool_name.equals("TenSQR")) {
						bw1.write("Out_File_Name ="+ pool_name+"_ViralSeq.txt\n");
					}else if(tool_name.equals("PredictHaplo")) {
						bw1.write("Out_File_Name ="+ pool_name+".out\n");
					}else if(tool_name.equals("QuasiRecomb")) {
						bw1.write("Out_File_Name ="+ "quasispecies.fasta\n");
					}
					
					bw1.write("Fasta_File_Name = "+ pool_name+"\n");
					bw1.write("Proj_Name ="+ project_name+"\n");
					bw1.write("Tool_Name ="+ tool_name+"\n");
					bw1.write("Ref_Seq_Path ="+ this.genome_path+"\n");
					bw1.write("Cutoff_Lowest_Freq ="+ 0.000001 +"\n");
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
				bw.write("#SBATCH --mem=30gb\n");
				bw.write("#SBATCH --ntasks=1\n");
				bw.write("#SBATCH --cpus-per-task=8\n");
				bw.write("#SBATCH --time=99-00:00:00\n");
				bw.write("#SBATCH --nodes=1\n");
				
				bw.write("java=/export/home/jhe/download/java_jdk_8u201/jdk1.8.0_201/bin/java\n");
				bw.write("java7=/export/home/jhe/download/jdk1.7.0_80/bin/java\n");
				bw.write("python=/export/home/jhe/download/anaconda/anaconda3/envs/tensqr_env/bin/python\n");
				bw.write("ExtractMatrix=/export/home/jhe/download/TenSQR/TenSQR-master/ExtractMatrix\n");
				bw.write("clustalo=/export/home/jhe/.local/bin/clustalo\n");
				bw.write("bwa=/export/home/jhe/download/bwa-0.7.17/bwa\n");
				bw.write("ref_dir=/export/home/jhe/project/Viral_reconstruction/SLiM/SLiM/Reference\n");
				bw.write("ref="+ this.genome_path+"\n");
				bw.write("PredictHaplo_Paired=/export/home/jhe/project/Viral_reconstruction/PredictHaplo/PredictHaplo-Paired-0.4/PredictHaplo-Paired\n");
				bw.write("TransferToFasta = /export/home/jhe/project/Viral_reconstruction/Other_Tools/programs/TransferToFasta.jar\n");
				bw.write("TransferOutput = /export/home/jhe/project/Viral_reconstruction/Other_Tools/programs/Transfer_Output.jar\n");
				bw.write("poolhapx=/export/qlong/PoolHapX/PoolHapX.jar\n");
				
				bw.write("start1=$SECONDS\n");
				for (int p=0;p < this.pools[i];p ++) {
					String pool_name = project_name +"_p"+p;
					String pool_folder = prefix_folder + "/"+ project_name +"/"+pool_name+"/";
					bw.write("prefix="+pool_folder+pool_name+"\n");
					bw.write("inbam="+ phx_folder_prefix + "/"+ project_name+"/input/bam/"
						+ pool_name +".rg.bam\n");
					bw.write("fastq_prefix="+ phx_folder_prefix + "/"+ project_name+"/input/fastq/"
							+ pool_name+"\n");

					if(tool_name.equals("QuasiRecomb")) {
						bw.write("cd "+ pool_folder +"\n");
						bw.write("$java7 -XX:+UseParallelGC -XX:NewRatio=9 -Xms10G "
								+ "-Xmx20G -jar /export/home/jhe/project/Viral_reconstruction/QuasiRecomb/QuasiRecomb.jar -i "
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
						bw_config.close();
						bw.write("cd "+ pool_folder +"\n");
						bw.write("$bwa mem $ref $fastq_prefix\\.bwa.read1.fastq $fastq_prefix\\.bwa.read2.fastq > $prefix\\.sam"+"\n");
						bw.write("$ExtractMatrix  $prefix\\.config\n");
						bw.write("$python /export/home/jhe/download/TenSQR/TenSQR-master/TenSQR.py $prefix\\.config\n");
						
					}else if(tool_name.equals("CliqueSNV")) {
						bw.write("$java -Xmx20G -jar /export/home/jhe/project/Viral_reconstruction/"
								+ "CliqueSNV/clique-snv.jar -m snv-illumina -tf "
								+ "0.000000001 -in $inbam -log\n");
					}else if(tool_name.equals("PredictHaplo")) {
						BufferedWriter bw_config = new BufferedWriter(new FileWriter(pool_folder+pool_name+".config"));
						bw_config.write("% configuration file for the HIVhaplotyper\n");
						bw_config.write("% prefix\n");
						bw_config.write(pool_name+"_\n");
						bw_config.write("% filename of reference sequence (FASTA)\n");
						bw_config.write(this.genome_path+"\n");
						bw_config.write("% do_visualize (1 = true, 0 = false)\n");
						bw_config.write("0\n");
						bw_config.write("% filname of the aligned reads (sam format)\n");
						bw_config.write(pool_name+".sam\n");
						bw_config.write("% have_true_haplotypes  (1 = true, 0 = false)\n");
						bw_config.write("0\n");
						bw_config.write("% filname of the true haplotypes (MSA in FASTA format) fill in any dummy filename if there is no true haplotypes)\n");
						bw_config.write("dump\n");
						bw_config.write("% do_local_analysis  (1 = true, 0 = false) (must be 1 in the first run)\n");
						bw_config.write("1\n");
						bw_config.write("% max_reads_in_window\n");
						bw_config.write("5000\n");
						bw_config.write("% entropy_threshold\n");
						bw_config.write("4e-4\n");
						bw_config.write("%reconstruction_start\n");
						bw_config.write("0\n");
						bw_config.write("%reconstruction_stop\n");
						bw_config.write("9719\n");
						bw_config.write("%min_mapping_qual\n");
						bw_config.write("0\n");
						bw_config.write("%min_readlength\n");
						bw_config.write("50\n");
						bw_config.write("%max_gap_fraction (relative to alignment length)\n");
						bw_config.write("0.01\n");
						bw_config.write("%min_align_score_fraction (relative to read length)\n");
						bw_config.write("0.10\n");
						bw_config.write("%alpha_MN_local (prior parameter for multinomial tables over the nucleotides)\n");
						bw_config.write("20\n");
						bw_config.write("%min_overlap_factor (reads must have an overlap with the local reconstruction window of at least this factor times the window size)\n");
						bw_config.write("0.5\n");
						bw_config.write("%local_window_size_factor (size of  local reconstruction window relative to the median of the read lengths)\n");
						bw_config.write("0.7\n");
						bw_config.write("% max number of clusters (in the truncated Dirichlet process)\n");
						bw_config.write("200\n");
						bw_config.write("% MCMC iterations\n");
						bw_config.write("501\n");
						bw_config.write("% include deletions (0 = no, 1 = yes)\n");
						bw_config.write("0\n");
						bw_config.close();
						
						bw.write("cd "+ pool_folder +"\n");
						bw.write("$bwa mem $ref $fastq_prefix\\.bwa.read1.fastq $fastq_prefix\\.bwa.read2.fastq > $prefix\\.sam"+"\n");
						bw.write("$PredictHaplo_Paired $prefix\\.config > $prefix\\.out"+"\n");
					}
				}//end of p_for_loop
				bw.write("end1=$SECONDS\n"); 
				bw.write("\n");
	        	
				bw.write("start2=$SECONDS\n");
				for(int p=0;p<this.pools[i];p++) {
					String pool_name = project_name +"_p"+p;
					String pool_folder = prefix_folder + "/"+ project_name +"/"+pool_name+"/";
						if(tool_name.equals("QuasiRecomb")) {
							BufferedReader br = new BufferedReader(new FileReader(pool_folder 
								+ "quasispecies.fasta"));
							String currline = br.readLine(); 
							int line_count=1;
							while(currline!=null) {
								currline = br.readLine(); 
								line_count++;
							}
							int num_file = line_count/100 + 1;
							br.close();
							bw.write("$java -jar $TransferToFasta "+ pool_folder+pool_name+"_O2R.properties"+"\n");
							for (int f=0;f<num_file;f++) {
								bw.write("$clustalo -i "+ pool_folder+pool_name+"_"+f+".fasta -o "+ 
										pool_folder+pool_name+"_"+f+".fa\n");
							}
							bw.write("$java -Xmx20g -jar $TransferOutput "+ pool_folder+pool_name+"_O2R.properties"+"\n");
						}else {
							bw.write("$java -jar $TransferToFasta "+ pool_folder+pool_name+"_O2R.properties"+"\n");
							bw.write("cp "+ pool_folder+pool_name+".fasta "+ pool_folder+pool_name+".fa"+"\n");
							bw.write("$java -Xmx20g -jar $TransferOutput "+ pool_folder+pool_name+"_O2R.properties"+"\n");
						}
					
					}//end of pool_for_loop
				
					BufferedWriter bw_properties = new BufferedWriter(new FileWriter(
						prefix_folder+"/"+project_name+"/PHX.properties"));
					
					bw_properties.write("Proj_Name = "+ project_name+"\n" );
					bw_properties.write("Input_Dir = "+phx_folder_prefix + "/"+ project_name +"/input/\n" );					
					bw_properties.write("Intermediate_Dir = "+phx_folder_prefix + "/"+ project_name +"/intermediate/\n" );					
					bw_properties.write("Output_Dir = "+prefix_folder + "/"+ project_name +"/\n" );			
					bw_properties.write("Gold_Dir = "+phx_folder_prefix + "/"+ project_name +"/gold_standard/\n" );
					
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
					bw_properties.write("Rscript_path = /export/home/jhe/download/anaconda/anaconda3/envs/Regress_Haplo/bin/Rscript\n");
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
					
					bw.write("properties="+prefix_folder + "/"+project_name+"/PHX.properties\n");
					bw.write("$java -jar $poolhapx evaluate $properties\n");
					
					bw.write("end2=$SECONDS\n"); 
					bw.write("echo \"duration_runing: $((end1-start1)) seconds.\"");
					bw.write("echo \"duration_calculating: $((end2-start2)) seconds.\"");
					bw.write("\n");
				}
			}
		}
	}

	public static void main(String[] args) throws IOException, InterruptedException {
		System.out.println("Viral Reconstruction Tool Script... ...");
		String prefix_folder= args[0];//
		String phx_folder_prefix = args[1];
		String tool_name = args[2];
		ScriptForTool cs = new ScriptForTool(prefix_folder,phx_folder_prefix,tool_name);
		System.out.println("Done, Enjoy!");
	}	
}