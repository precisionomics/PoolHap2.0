package Viral_Reconstructions_Tools;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;

public class ScriptForTool_arc_new_split {
	

	// project_name = "1_1,1_2..1_9;2_1,2_2...2_9...5_1,5_2...5_9"
	//int[] project_idx1= new int[] {5};
	//int[] project_idx2= new int[] {15};
	int num_pool = 25;
	int num_coverage = 5000;

	String slim_script; 
	
	int genome_len = 9719;
	String genome_path ="/home/jingni.he1/project/"
			+ "Viral_reconstruction/SLiM/Reference/HIV_HXB2.fa";
	
	
	public ScriptForTool_arc_new_split(String prefix_folder, String phx_folder_prefix, String tool_name) throws IOException {
		
		new File(prefix_folder + "/cmd/").mkdir();
		//for (int i =0; i< this.project_idx1.length; i++) {
		//	for (int j =0; j< this.project_idx2.length; j++) {
				
				
				//String project_name = this.project_idx1[i]+"_"+ this.project_idx2[j];
		String project_name = "5_15";
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
		for(int p=0;p<this.num_pool;p++) {
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
			bw1.write("Cutoff_Lowest_Freq ="+ 0.01 +"\n");
			bw1.close();
		}
		for (int p=0;p < this.num_pool;p ++) {
			String pool_name = project_name +"_p"+p;
			String pool_folder = prefix_folder + "/"+ project_name +"/"+pool_name+"/";
			BufferedWriter bw = new BufferedWriter(new FileWriter(prefix_folder + "/cmd/"+project_name+"_p"+p+".cmd"));
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
				bw.write("#SBATCH --job-name="+ tool_name + project_name+"_p"+p+"\n");
				bw.write("#SBATCH --workdir="+ pool_folder+"\n");
				bw.write("#SBATCH --error="+project_name+"_p"+p+".error\n" );
				bw.write("#SBATCH --output="+project_name+"_p"+p+".out\n" );
				bw.write("#SBATCH --mem=30gb\n");
				bw.write("#SBATCH --ntasks=1\n");
				bw.write("#SBATCH --cpus-per-task=6\n");
				bw.write("#SBATCH --time=7-00:00:00\n");
				bw.write("#SBATCH --nodes=1\n");
				bw.write("#SBATCH --partition=theia,parallel,cpu2019,cpu2013,lattice\n");
				
				bw.write("java=/home/jingni.he1/download/java_jdk_8u201/jdk1.8.0_201/bin/java\n");
				bw.write("java7=/home/jingni.he1/download/jdk1.7.0_80/bin/java\n");
				bw.write("python=/home/jingni.he1/anaconda3/envs/tensqr_env/bin/python\n");
				bw.write("ExtractMatrix=/home/jingni.he1/project/Viral_reconstruction/SLiM/programs/ExtractMatrix\n");
				bw.write("bwa=/home/jingni.he1/download/bwa-0.7.17/bwa\n");
				bw.write("ref_dir=/home/jingni.he1/project/Viral_reconstruction/SLiM/Reference\n");
				bw.write("ref="+ this.genome_path+"\n");
				bw.write("PredictHaplo_Paired=/home/jingni.he1/project/Viral_reconstruction/SLiM/PredictHaplo-Paired-0.4/PredictHaplo-Paired\n");
				
				bw.write("start=$SECONDS\n");
				//for (int p=0;p < this.num_pool;p ++) {
					
					bw.write("prefix="+pool_folder+pool_name+"\n");
					bw.write("inbam="+ phx_folder_prefix + "/"+ project_name+"/input/bam/"
						+ pool_name +".rg.bam\n");
					bw.write("fastq_prefix="+ phx_folder_prefix + "/"+ project_name+"/input/fastq/"
							+ pool_name+"\n");

					if(tool_name.equals("QuasiRecomb")) {
						bw.write("cd "+ pool_folder +"\n");
						bw.write("$java7 -XX:+UseParallelGC -XX:NewRatio=9 -Xms10G "
								+ "-Xmx20G -jar /home/jingni.he1/project/Viral_reconstruction/SLiM/programs/QuasiRecomb.jar -i "
								+ "$inbam -conservative"+"\n");
					}else if(tool_name.equals("TenSQR")){
						BufferedWriter bw_config = new BufferedWriter(new FileWriter(pool_folder+pool_name+".config"));
						bw_config.write("filename of reference sequence (FASTA) :"+ this.genome_path+"\n");
						bw_config.write("filname of the aligned reads (sam format) : "+ pool_name+".sam\n");
						bw_config.write("SNV_thres : 0.01\n");
						bw_config.write("reconstruction_start : 1\n");
						bw_config.write("reconstruction_stop: 9719\n");
						bw_config.write("min_mapping_qual : 40\n");
						bw_config.write("min_read_length : 100\n");
						bw_config.write("max_insert_length : 400\n");
						bw_config.write("characteristic zone name : "+pool_name+"\n");
						bw_config.write("seq_err (assumed sequencing error rate(%)) : 0.1\n");
						bw_config.write("MEC improvement threshold : 0.0312\n");
						bw_config.write("initial population size :100\n");
						bw_config.close();
						bw.write("cd "+ pool_folder +"\n");
						bw.write("$bwa mem $ref $fastq_prefix\\.bwa.read1.fastq $fastq_prefix\\.bwa.read2.fastq > $prefix\\.sam"+"\n");
						bw.write("$ExtractMatrix  $prefix\\.config\n");
						bw.write("$python /home/jingni.he1/project/Viral_reconstruction/SLiM/programs/TenSQR.py $prefix\\.config\n");
						
					}else if(tool_name.equals("CliqueSNV")) {
						bw.write("$java -Xmx30G -jar /home/jingni.he1/project/Viral_reconstruction/SLiM/programs/clique-snv.jar -m snv-illumina "
								//+ "-tf 0.01 "
								+"-in $inbam -log\n"
								);
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
						bw_config.write("10000\n");
						bw_config.write("% entropy_threshold\n");
						bw_config.write("4e-4\n");
						bw_config.write("%reconstruction_start\n");
						bw_config.write("0\n");
						bw_config.write("%reconstruction_stop\n");
						bw_config.write("9719\n");
						bw_config.write("%min_mapping_qual\n");
						bw_config.write("30\n");
						bw_config.write("%min_readlength\n");
						bw_config.write("100\n");
						bw_config.write("%max_gap_fraction (relative to alignment length)\n");
						bw_config.write("0.05\n");
						bw_config.write("%min_align_score_fraction (relative to read length)\n");
						bw_config.write("0.35\n");
						bw_config.write("%alpha_MN_local (prior parameter for multinomial tables over the nucleotides)\n");
						bw_config.write("25\n");
						bw_config.write("%min_overlap_factor (reads must have an overlap with the local reconstruction window of at least this factor times the window size)\n");
						bw_config.write("0.85\n");
						bw_config.write("%local_window_size_factor (size of  local reconstruction window relative to the median of the read lengths)\n");
						bw_config.write("0.7\n");
						bw_config.write("% max number of clusters (in the truncated Dirichlet process)\n");
						bw_config.write("25\n");
						bw_config.write("% MCMC iterations\n");
						bw_config.write("501\n");
						bw_config.write("% include deletions (0 = no, 1 = yes)\n");
						bw_config.write("0\n");
						bw_config.close();
						
						bw.write("cd "+ pool_folder +"\n");
						bw.write("$bwa mem $ref $fastq_prefix\\.bwa.read1.fastq $fastq_prefix\\.bwa.read2.fastq > $prefix\\.sam"+"\n");
						bw.write("$PredictHaplo_Paired $prefix\\.config > $prefix\\.out"+"\n");
					}
					bw.write("end=$SECONDS\n"); 
					bw.write("\n");
					bw.write("echo \"duration_runing: $((end-start)) seconds.\""+"\n");
					bw.close();
					
				}//end of p_for_loop			
		
	}

	public static void main(String[] args) throws IOException, InterruptedException {
		System.out.println("Viral Reconstruction Tool Script... ...");
		String prefix_folder= args[0];//
		String phx_folder_prefix = args[1];
		String tool_name = args[2];
		ScriptForTool_arc_new_split cs = new ScriptForTool_arc_new_split(prefix_folder,phx_folder_prefix,tool_name);
		System.out.println("Done, Enjoy!");
	}	
}