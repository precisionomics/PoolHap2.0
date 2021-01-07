package MiscFunctions;

import java.io.*;
import java.util.*;

import PoolHap.Entrance;
import breeze.macros.expand.args;
import shapeless.newtype;
import spire.math.UInt;
import spire.optional.intervalGeometricPartialOrder;



public class ComparisonSim_island_arc_new {

	// project_name = "1_1,1_2..1_9;2_1,2_2...2_9...5_1,5_2...5_9"
	int[] project_idx1= new int[] {4};
	//int[] project_idx2= new int[] {1,2,3,4,5,6,7,8,9};
	int num_pool = 25;
	int num_coverage = 5000;
	String slim_script; 
	
	int genome_len = 9719;
	String genome_path ="/home/jingni.he1/project/"
			+ "Viral_reconstruction/SLiM/Reference/HIV_HXB2.fa"; 
	
	
	public ComparisonSim_island_arc_new(String prefix_folder, String slimout_folder) throws IOException {
		
		new File(prefix_folder + "/cmd/").mkdir();
		for (int i =0; i< this.project_idx1.length; i++) {
			for (int j =1; j< 44; j++) {
				
				String project_name = this.project_idx1[i]+"_"+ j;
						
				new File(prefix_folder +"/"+ project_name).mkdir();
				
				BufferedWriter bw1 = new BufferedWriter(new FileWriter(prefix_folder +
						"/cmd/"+project_name +".PoolSimulator.properties"));
				
				double hap_freq_cutoff = (double)0*0.005;
				
				bw1.write("Input_Dir = "+ prefix_folder +"/" + project_name+ "/input\n");
				bw1.write("Intermediate_Dir ="+ prefix_folder+"/" + project_name+ "/intermediate\n");
				bw1.write("Gold-Standard_Dir ="+ prefix_folder+"/" + project_name+ "/gold_standard\n");
				bw1.write("Slim_Output_Path ="+ prefix_folder+"/" + project_name+ "/gold_standard\n");
//				Slim_Model = panmictic_haploid
				bw1.write("Slim_Model = island_haploid\n");
				bw1.write("Proj_Name =" + project_name +"\n");
				
				bw1.write("Is_Single_Population = false\n");
				bw1.write("Is_Ms_Output = false\n");
				bw1.write("DWGSIM = /home/jingni.he1/download/DWGSIM-master/dwgsim\n");
				bw1.write("Num_Pools = "+this.num_pool+ "\n");
				bw1.write("Ref_Seq_Len = "+ Integer.toString(this.genome_len)+"\n");
				bw1.write("Reference_Seq= " + this.genome_path+"\n");
				bw1.write("Is_Perfect = false\n");
				bw1.write("Error_Rate_Per_Base = 0.001 \n");
				bw1.write("Hap_Freq_Cutoff = "+ Double.toString(hap_freq_cutoff)+  "\n");
				bw1.write("Coverage = " +this.num_coverage + "\n");
				bw1.write("Read_Len = 150\n");
				bw1.write("Outer_Dist = 400\n");
				bw1.write("Weak_Length = 550\n");
				bw1.write("Num_Haps_Pool = 30\n");
				bw1.close();

				BufferedWriter bw = new BufferedWriter(new FileWriter(prefix_folder + "/cmd/"+ project_name +".cmd"));
				
				bw.write("#!/bin/bash\n");
				bw.write("#SBATCH --job-name="+ project_name +"\n");
				bw.write("#SBATCH --workdir="+ prefix_folder + "/"+ project_name +"\n");
				bw.write("#SBATCH --error="+project_name+".error\n" );
				bw.write("#SBATCH --output="+project_name+".out\n" );
				bw.write("#SBATCH --mem=25gb\n");
				bw.write("#SBATCH --ntasks=1\n");
				bw.write("#SBATCH --cpus-per-task=2\n");
				bw.write("#SBATCH --time=7-00:00:00\n");
				bw.write("#SBATCH --nodes=1\n");
				bw.write("#SBATCH --partition=theia\n");
				
				
//				java=/home/jingni.he1/download/java_jdk_8u201/jdk1.8.0_201/bin/java
//				slim=/home/jingni.he1/download/SliM/build/slim
				bw.write("source /home/jingni.he1/anaconda3/bin/activate R\n");
				bw.write("java=/home/jingni.he1/download/java_jdk_8u201/jdk1.8.0_201/bin/java\n");
				bw.write("slim=/home/jingni.he1/download/SliM/build/slim\n");
				bw.write("poolhapx=/home/jingni.he1/project/Viral_reconstruction/SLiM/programs/PoolHapX.jar\n");
				bw.write("poolsim=/home/jingni.he1/project/Viral_reconstruction/SLiM/programs/PoolSimulator_SLiM.jar\n");
				bw.write("rewrite_vars=/home/jingni.he1/project/Viral_reconstruction/SLiM/programs/Rewrite_VarsFile.jar\n");
				bw.write("bwa=/home/jingni.he1/download/bwa-0.7.17/bwa\n");
				bw.write("samtools=/home/jingni.he1/download/samtools-1.9/samtools\n");
				bw.write("ref="+ this.genome_path+"\n");
				bw.write("gatk=/home/jingni.he1/download/gatk-4.0.0.0/"
						+ "gatk-package-4.0.0.0-local.jar\n");
			
				
				bw.write("mkdir "+ prefix_folder + "/"+ project_name+ "/input\n");
				bw.write("mkdir "+ prefix_folder + "/"+ project_name+ "/gold_standard\n");
				bw.write("mkdir "+ prefix_folder + "/"+ project_name+ "/output\n");
				bw.write("mkdir "+ prefix_folder + "/"+ project_name+ "/intermediate\n");
				
				bw.write("cp  "+ slimout_folder + "/"+project_name+ "_new.slim "
						+ prefix_folder + "/"+ project_name+ "/gold_standard/"+ 
						project_name + "_island_haploid.out\n");
				
// Step 1: Generate coalescence-simulated haplotypes, distribute to each of the pools, and simulate reads for each pool.
				
				bw.write("$java -jar $poolsim "+ 	prefix_folder + "/cmd/"+ 
						project_name +".PoolSimulator.properties\n");
			
				
				for (int p=0;p < this.num_pool;p ++) {
					
					bw.write("prefix="+ prefix_folder + "/"+ project_name+"/input/fastq/"+project_name+
							"_p"+ Integer.toString(p) +"\n");
					
					
					bw.write("prefix_bam="+ prefix_folder + "/"+ project_name+"/input/bam/"+project_name+
								"_p"+ Integer.toString(p) +"\n");
					

// Step 2: For each pool, align the simulated reads to a reference sequence.	
					
					bw.write("gunzip   $prefix\\.bwa.read1.fastq\n");
					bw.write("gunzip   $prefix\\.bwa.read2.fastq\n");
					bw.write("$bwa mem $ref $prefix\\.bwa.read1.fastq $prefix\\.bwa.read2.fastq "
							+ "| $samtools view -Shub - > $prefix_bam\\.bam\n");
					bw.write("$samtools sort -o  " +	"$prefix_bam\\.srt.bam  $prefix_bam\\.bam\n");
					
					
// Step 3: For each pool, call variants using GATK HaplotypeCaller in gVCF mode.	
					
					bw.write("inbam="+ prefix_folder + "/"+ project_name +"/input/bam/"+project_name+
									"_p"+ Integer.toString(p) +".rg.bam\n");
								
					bw.write("$java -jar $gatk  AddOrReplaceReadGroups -I  $prefix_bam\\.srt.bam -O $inbam"
							+ " -R $ref -ID " +  project_name + "_p"+ Integer.toString(p)
							+" -LB NPD -PL Illumina -PU NPD -SM "+ project_name + "_p"+ Integer.toString(p)+ "\n");
					
					bw.write("$samtools index $inbam\n");
										
					bw.write("prefix_vcf="+ prefix_folder + "/"+ project_name +"/input/vcf/"+project_name+
							"_p"+ Integer.toString(p) +"\n");
					
					
					bw.write("outgvcf=$prefix_vcf\\.raw.g.vcf\n");
					bw.write("$java -jar $gatk  HaplotypeCaller -R  $ref -I $inbam"
							+ " -ERC  GVCF -ploidy 8 --heterozygosity 0.01  --max-alternate-alleles 1 -O "
							+ "$outgvcf\n");
				}
//	Step 4: Join all pool-specific gVCFs into a joint gVCF file and convert to VCF.
				
				bw.write("prefix="+ prefix_folder + "/"+ project_name +"/input/vcf/"+project_name+"\n");
				
				String tmp= "$java -jar $gatk CombineGVCFs -R $ref ";
				for (int p=0;p < this.num_pool;p ++) {
					tmp=tmp+" -V " + "$prefix\\"+
							"_p"+ Integer.toString(p)+".raw.g.vcf ";
				}
				
				
				tmp=tmp+ " -O   $prefix\\.g.vcf\n"  ;
				bw.write(tmp );
				
				bw.write("$java -jar $gatk   GenotypeGVCFs -R $ref -V  $prefix\\.g.vcf" + 
						" -ploidy 8 -O $prefix\\.raw.vcf\n");
				
				bw.write("prefix_vcf="+ prefix_folder + "/"+ project_name +"/input/"
                        +project_name+"\n");

				
				bw.write( "$java -jar $gatk  SelectVariants -R $ref -V  $prefix\\.raw.vcf" + 
						" -O $prefix_vcf\\.vcf\n");
				
//	Step 5. Convert each pool-specific BAM file to SAM (i.e.: text), then VEF files.	
				for (int p=0;p < this.num_pool;p ++) {
					
					bw.write("prefix_sam="+prefix_folder + "/"+ project_name+"/input/sam/"
							+ project_name+
							"_p"+ Integer.toString(p) +"\n");
					
					bw.write("prefix_bam="+prefix_folder + "/"+ project_name+"/input/bam/"
							+ project_name +
							"_p"+ Integer.toString(p) +"\n");
					
					bw.write("$samtools view -ho $prefix_sam\\.sam $prefix_bam\\.srt.bam\n");
					
				}
				
				new File(prefix_folder + "/"+ project_name +"/input/").mkdir();         
				
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
				bw.write("start=$SECONDS\n");
				bw.write("$java -jar $poolhapx format $properties\n");
				bw.write("$java -jar $rewrite_vars "+ prefix_folder + "/"
						+ project_name+ "/gold_standard/"
						+ project_name
						+"_haps.inter_freq_vars.txt "+prefix_folder + "/"
						+ project_name+ "/intermediate/"
						+ project_name
						+"_vars.intra_freq.txt\n");
				bw.write("/usr/bin/python2.7  /home/jingni.he1/project/Viral_reconstruction/SLiM/programs/make_gc.py "
						+ prefix_folder + "/"+ project_name+ "/intermediate "+ project_name +"\n");
				//bw.write("$java -jar $poolhapx gc $properties\n");
				bw.write("$java -jar $poolhapx aem $properties\n");
				bw.write("$java -jar $poolhapx evaluate $properties\n");
				bw.write("end=$SECONDS\n" + 
						"echo \"duration: $((end-start)) seconds.\"");
				bw.write("\n");
	        	bw.close();	        	
					
				}
			}
		
	}

	public static  void main(String[] args) throws IOException, InterruptedException {
		System.out.println("PoolHapX Comparison Simulation... ...");
		String prefix_folder= args[0];//"/export/qlong/chencao/Work/poolhapx/slim/sim/";
		String slimout_folder = args[1];
		ComparisonSim_island_arc_new cs = new ComparisonSim_island_arc_new(prefix_folder,slimout_folder);
		System.out.println("Done, Enjoy!");
		
	}	
	
	
}

