package MiscFunctions;

import java.io.*;
import java.util.*;

import PoolHap.Entrance;
import breeze.macros.expand.args;
import shapeless.newtype;
import spire.math.UInt;
import spire.optional.intervalGeometricPartialOrder;



public class ComparisonSim_island_arc_new_phx {

	// project_name = "1_1,1_2..1_9;2_1,2_2...2_9...5_1,5_2...5_9"
	int[] project_idx1= new int[] {1,2,3,4,5};
	int[] project_idx2= new int[] {1,2,3,4,5,6,7,8,9};
	int num_pool = 25;
	int num_coverage = 5000;
	String slim_script; 
	
	int genome_len = 9719;
	String genome_path ="/home/jingni.he1/project/"
			+ "Viral_reconstruction/SLiM/Reference/HIV_HXB2.fa"; 
	
	
	public ComparisonSim_island_arc_new_phx(String prefix_folder, String slimout_folder) throws IOException {
		
		new File(prefix_folder + "/cmd/").mkdir();
		for (int i =0; i< this.project_idx1.length; i++) {
			for (int j =0; j< this.project_idx2.length; j++) {
				
				String project_name = this.project_idx1[i]+"_"+ this.project_idx2[j];
				BufferedWriter bw = new BufferedWriter(new FileWriter(prefix_folder + "/cmd/"+ project_name +"_phx.cmd"));
				
				bw.write("#!/bin/bash\n");
				bw.write("#SBATCH --job-name="+ project_name +"\n");
				bw.write("#SBATCH --workdir="+ prefix_folder + "/"+ project_name +"\n");
				bw.write("#SBATCH --error="+project_name+".error\n" );
				bw.write("#SBATCH --output="+project_name+".out\n" );
				bw.write("#SBATCH --mem=40gb\n");
				bw.write("#SBATCH --ntasks=1\n");
				bw.write("#SBATCH --cpus-per-task=3\n");
				bw.write("#SBATCH --time=7-00:00:00\n");
				bw.write("#SBATCH --nodes=1\n");
				bw.write("#SBATCH --partition=theia\n");
				
				
				bw.write("source /home/jingni.he1/anaconda3/bin/activate R\n");
				bw.write("java=/home/jingni.he1/download/java_jdk_8u201/jdk1.8.0_201/bin/java\n");
				bw.write("slim=/home/jingni.he1/download/SliM/build/slim\n");
				bw.write("poolhapx=/home/jingni.he1/project/Viral_reconstruction/SLiM/programs/PoolHapX.jar\n");
				bw.write("poolsim=/home/jingni.he1/project/Viral_reconstruction/SLiM/programs/PoolSimulator_SLiM.jar\n");
				bw.write("rewrite_vars=/home/jingni.he1/project/Viral_reconstruction/SLiM/programs/Rewrite_VarsFile.jar\n");
				BufferedWriter bw_properties = new BufferedWriter(new FileWriter(
						prefix_folder + "/"+ project_name+"/input/PHX.properties"));
					
					bw_properties.write("Proj_Name = "+ project_name+"\n" );

					bw_properties.write("Input_Dir = "+prefix_folder + "/"+ project_name +"/input/\n" );
					
					bw_properties.write("Intermediate_Dir = "+prefix_folder + "/"+ project_name +"/intermediate/\n" );
					
					bw_properties.write("Output_Dir = "+prefix_folder + "/"+ project_name +"/output/\n" );
					
					bw_properties.write("Gold_Dir = "+prefix_folder + "/"+ project_name +"/gold_standard/\n" );
					
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
					bw_properties.write("Rscript_path = /home/jingni.he1/anaconda3/envs/R/bin/Rscript\n");
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
					
					
				bw.write("properties="+prefix_folder + "/"+ project_name+"/input/PHX.properties\n");
				bw.write("poolhapx=/home/jingni.he1/project/Viral_reconstruction/SLiM/programs/PoolHapX.jar\n");				
				bw.write("$java -jar $poolhapx l0l1 $properties\n");
				bw.write("$java -jar $poolhapx evaluate $properties\n");
				bw.write("\n");
				
				
	        	bw.close();	        	
					
				}
			}
	}
	public static  void main(String[] args) throws IOException, InterruptedException {
		System.out.println("PoolHapX Comparison Simulation... ...");
		String prefix_folder= args[0];//"/export/qlong/chencao/Work/poolhapx/slim/sim/";
		String slimout_folder = args[1];
		ComparisonSim_island_arc_new_phx cs = new ComparisonSim_island_arc_new_phx(prefix_folder,slimout_folder);
		System.out.println("Done, Enjoy!");
		
	}		
}
				
				

				