package Viral_Reconstructions_Tools;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;

public class ScriptCompareOut_arc_new {
	
	// project_name = "1_1,1_2..1_9;2_1,2_2...2_9...5_1,5_2...5_9"
	int[] project_idx1= new int[] {1,2,3,4,5};
	int[] project_idx2= new int[] {1,2,3,4,5,6,7,8,9};
	int num_pool = 25;
	int num_coverage = 5000;

	String slim_script; 
	
	int genome_len = 9719;
	String genome_path ="/home/jingni.he1/project/"
			+ "Viral_reconstruction/SLiM/Reference/HIV_HXB2.fa";
	 
	
	
	public ScriptCompareOut_arc_new(String prefix_folder, String phx_folder_prefix, String tool_name) throws IOException {
		new File(prefix_folder + "/cmd_out/").mkdir();
		for (int i =0; i< this.project_idx1.length; i++) {
			for (int j =0; j< this.project_idx2.length; j++) {
				
				
				String project_name = this.project_idx1[i]+"_"+ this.project_idx2[j];
				
				BufferedWriter bw = new BufferedWriter(new FileWriter(prefix_folder + "/cmd_out/"+project_name+".cmd"));
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
				bw.write("#SBATCH --job-name=CO"+ tool_name + project_name+"\n");
				bw.write("#SBATCH --workdir="+ prefix_folder + "/"+ project_name +"\n");
				bw.write("#SBATCH --error="+ project_name+"finalout.error\n" );
				bw.write("#SBATCH --output="+project_name+"finalout.out\n" );
				bw.write("##SBATCH --mem=20gb\n");
				bw.write("#SBATCH --ntasks=1\n");
				bw.write("#SBATCH --cpus-per-task=4\n");
				bw.write("#SBATCH --time=7-00:00:00\n");
				bw.write("#SBATCH --nodes=1\n");
				bw.write("#SBATCH --partition=theia,parallel,cpu2019,cpu2013\n");
				
				bw.write("java=/home/jingni.he1/download/java_jdk_8u201/jdk1.8.0_201/bin/java\n");
				bw.write("clustalo=/home/jingni.he1/.local/bin/clustalo\n");
				bw.write("TransferToFasta =/home/jingni.he1/project/Viral_reconstruction/SLiM/programs/TransferToFasta.jar\n");
				bw.write("TransferOutput = /home/jingni.he1/project/Viral_reconstruction/SLiM/programs/Transfer_Output.jar\n");
				bw.write("poolhapx=/home/jingni.he1/project/Viral_reconstruction/SLiM/programs/PoolHapX.jar\n");
				
				bw.write("start=$SECONDS\n");
				for(int p=0;p<this.num_pool;p++) {
					String pool_name = project_name +"_p"+p;
					String pool_folder = prefix_folder + "/"+ project_name +"/"+pool_name+"/";
						if(tool_name.equals("QuasiRecomb")) {
//							BufferedReader br = new BufferedReader(new FileReader(pool_folder 
//								+ "quasispecies.fasta"));
//							String currline = br.readLine(); 
//							int line_count=1;
//							while(currline!=null||line_count<=100) {
//								currline = br.readLine(); 
//								line_count++;
//							}
//							int num_file = line_count/100 + 1;
//							br.close();
							bw.write("$java -jar /home/jingni.he1/project/Viral_reconstruction/SLiM/programs/TransferToFasta.jar "+ pool_folder+pool_name+"_O2R.properties"+"\n");
							bw.write("/home/jingni.he1/.local/bin/clustalo -i "+ pool_folder+pool_name+".fasta -o "+ 
									pool_folder+pool_name+".fa --force"+"\n");
//							for (int f=0;f<num_file;f++) {
//								bw.write("$clustalo -i "+ pool_folder+pool_name+"_"+f+".fasta -o "+ 
//										pool_folder+pool_name+"_"+f+".fa\n");
//							}
							bw.write("$java -Xmx20g -jar /home/jingni.he1/project/Viral_reconstruction/SLiM/programs/Transfer_Output.jar "+ pool_folder+pool_name+"_O2R.properties"+"\n");
						}else {
							bw.write("$java -jar /home/jingni.he1/project/Viral_reconstruction/SLiM/programs/TransferToFasta.jar "+ pool_folder+pool_name+"_O2R.properties"+"\n");
							bw.write("cp "+ pool_folder+pool_name+".fasta "+ pool_folder+pool_name+".fa"+"\n");
							bw.write("$java -Xmx20g -jar /home/jingni.he1/project/Viral_reconstruction/SLiM/programs/Transfer_Output.jar "+ pool_folder+pool_name+"_O2R.properties"+"\n");
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
					bw_properties.write("Species = virus\n");
					bw_properties.close();
					
					bw.write("properties="+prefix_folder + "/"+project_name+"/PHX.properties\n");
					bw.write("$java -jar $poolhapx evaluate $properties\n");
					
					bw.write("end=$SECONDS\n"); 
					bw.write("echo \"duration_calculating: $((end-start)) seconds.\"");
					bw.write("\n");
					bw.close();
					
				}
			}
		
	}

	public static void main(String[] args) throws IOException, InterruptedException {
		System.out.println("Generating output for tools... ...");
		String prefix_folder= args[0];//
		String phx_folder_prefix = args[1];
		String tool_name = args[2];
		ScriptCompareOut_arc_new cs = new ScriptCompareOut_arc_new(prefix_folder,phx_folder_prefix,tool_name);
		System.out.println("Done, Enjoy!");
	}	
}
