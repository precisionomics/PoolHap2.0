package Viral_Reconstructions_Tools;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.io.InputStream;
import java.util.ArrayList;
import java.util.Properties;

public class ScriptForPHX {
	String main_dir;
	String project_name;
	String java;
	String phx_jar;
	String bwa;
	String samtools;
	String gatk;
	String fastq_path;
	String fastq_file;
	String ref_path;
	String Rscript_path;
	String species;
	ArrayList<String> sample_name_list= new ArrayList<String>();
	
	public ScriptForPHX(String parameter_file) throws IOException {
		InputStream is = new FileInputStream(parameter_file);
	    Properties prop = new Properties();
	    prop.load(is);
	    this.main_dir=prop.getProperty("Main_Dir");
	    this.project_name = prop.getProperty("Project_Name");
	    this.java=prop.getProperty("Java");
	    this.phx_jar=prop.getProperty("PHX_JAR");
	    this.bwa=prop.getProperty("bwa");
	    this.samtools=prop.getProperty("samtools");
	    this.gatk=prop.getProperty("gatk");
	    this.fastq_path=prop.getProperty("Fastq_Path");
	    this.ref_path=prop.getProperty("Ref_Path");
	    this.Rscript_path=prop.getProperty("Rscript_Path");
	    this.fastq_file=prop.getProperty("Fastq_Name");	    
	    is.close();
	    System.out.println(this.fastq_file);
	    System.out.println(this.fastq_file);
	    BufferedReader br = new BufferedReader(new FileReader(this.fastq_file)); 
	    String currline=br.readLine();
	    while(currline!=null) {
	    	String[] tmpcurrpos = currline.split("\t");
	    	String curr_sample_name=tmpcurrpos[0].split(".read1.fastq")[0];
	    	this.sample_name_list.add(curr_sample_name);
	    	currline=br.readLine();
	    }
		br.close();
	}
	
	public void ComparisonSim() throws IOException {
		new File(this.main_dir + "/cmd/").mkdir();
		String project_name=this.project_name;
		new File(this.main_dir +"/"+ project_name).mkdir();
		new File(this.main_dir+ "/"+ project_name +"/input/").mkdir();    
		BufferedWriter bw = new BufferedWriter(new FileWriter(this.main_dir + "/cmd/"+ project_name +".cmd"));
			bw.write("java="+this.java+"\n");
			bw.write("poolhapx="+this.phx_jar+"\n");
			bw.write("bwa="+this.bwa+"\n");
			bw.write("samtools="+this.samtools+"\n");
			bw.write("gatk="+this.gatk+"\n");				
			bw.write("ref="+ this.ref_path+"\n");
			bw.write("mkdir "+ this.main_dir + "/"+ project_name+ "/input/bam\n");
			bw.write("mkdir "+ this.main_dir + "/"+ project_name+ "/input/vcf\n");
			bw.write("mkdir "+ this.main_dir + "/"+ project_name+ "/input/sam\n");
			bw.write("mkdir "+ this.main_dir + "/"+ project_name+ "/output\n");
			bw.write("mkdir "+ this.main_dir + "/"+ project_name+ "/intermediate\n");
				
			for (int i =0; i< this.sample_name_list.size(); i++) {
				String sample_name=sample_name_list.get(i);
				bw.write("prefix_fastq="+this.fastq_path+"/"+sample_name+"\n");
				bw.write("prefix_bam="+ this.main_dir + "/"+ project_name+"/input/bam/"+sample_name +"\n");
				bw.write("inbam="+ this.main_dir + "/"+ project_name +"/input/bam/"+sample_name +".rg.bam\n");
				bw.write("prefix_vcf="+this.main_dir + "/"+ project_name +"/input/vcf/"+sample_name+"\n");
				bw.write("outgvcf=$prefix_vcf\\"+".raw.g.vcf\n");
// Step 1: For each pool, align the simulated reads to a reference sequence.	
				bw.write("$bwa mem $ref $prefix_fastq\\.read1.fastq $prefix_fastq\\.read2.fastq "
							+ "| $samtools view -Shub - > $prefix_bam\\.bam\n");
				bw.write("$samtools sort -o  " +	"$prefix_bam\\.srt.bam  $prefix_bam\\.bam\n");	
// Step 2: For each pool, call variants using GATK HaplotypeCaller in gVCF mode.												
				bw.write("$gatk  AddOrReplaceReadGroups -I  $prefix_bam\\.srt.bam -O $inbam"
							+ " -R $ref -ID " +  sample_name
							+" -LB NPD -PL Illumina -PU NPD -SM "+ sample_name+ "\n");
				bw.write("$samtools index $inbam\n");						
				bw.write("$gatk  HaplotypeCaller -R  $ref -I $inbam"
							+ " -ERC  GVCF -ploidy 8 --heterozygosity 0.01  "
							+ "--max-alternate-alleles 1 -O "
							+ "$outgvcf\n");
			}//end_of_for
//	Step 3: Join all pool-specific gVCFs into a joint gVCF file and convert to VCF.
			String tmp= "$gatk CombineGVCFs -R $ref ";
			for (int i=0;i < this.sample_name_list.size();i ++) {
				String sample_name=sample_name_list.get(i);
				tmp=tmp+" -V " + this.main_dir + "/"+ project_name +"/input/vcf/"+sample_name+".raw.g.vcf";
			}
			bw.write("prefix_project_vcf="+ this.main_dir + "/"+ project_name +"/input/vcf/"+project_name+"\n");
			tmp=tmp+ " -O $prefix_project_vcf\\.g.vcf\n"  ;
			bw.write(tmp);
			bw.write("$gatk GenotypeGVCFs -R $ref -V  $prefix_project_vcf\\.g.vcf" + 
						" -ploidy 8 -O $prefix_project_vcf\\.raw.vcf\n");
			bw.write("prefix_project="+ this.main_dir + "/"+ project_name +"/input/"+project_name+"\n");
			bw.write( "$gatk  SelectVariants -R $ref -V  $prefix_project_vcf\\.raw.vcf" + 
						" -O $prefix_project\\.vcf\n");				
//	Step 4. Convert each pool-specific BAM file to SAM (i.e.: text), then VEF files.	
			for (int i=0;i < this.sample_name_list.size();i ++) {
				String sample_name=sample_name_list.get(i);
				bw.write("prefix_sam="+this.main_dir + "/"+ project_name+"/input/sam/"+sample_name+"\n");
				bw.write("prefix_bam="+this.main_dir+ "/"+ project_name+"/input/bam/"+sample_name+"\n");
				bw.write("samtools view -ho $prefix_sam\\.sam $prefix_bam\\.srt.bam\n");
					
			}
//	Step 5. Write PoolHapX Properties file and run PoolHapX.	
			BufferedWriter bw_properties = new BufferedWriter(new FileWriter(
					this.main_dir + "/"+ project_name+"/input/PHX.properties"));
				bw_properties.write("Proj_Name = "+ project_name+"\n" );
				bw_properties.write("Input_Dir = "+this.main_dir + "/"+ project_name +"/input/\n" );				
				bw_properties.write("Intermediate_Dir = "+this.main_dir + "/"+ project_name +"/intermediate/\n" );				
				bw_properties.write("Output_Dir = "+this.main_dir + "/"+ project_name +"/output/\n" );								
				bw_properties.write("Num_Pos_Window = 20\n");
				bw_properties.write("Num_Gap_Window = 2\n");					
				bw_properties.write("In-pool_Gap_Support_Min = 1\n");
				bw_properties.write("All-pool_Gap_Support_Min = 1\n");
				bw_properties.write("Level_1_Region_Size_Min = 10\n");
				bw_properties.write("Level_1_Region_Size_Max = 12\n");
				bw_properties.write("Level_1T_Region_Size_Min = 10\n");
				bw_properties.write("Level_1T_Region_Size_Max = 12\n");
				bw_properties.write("Est_Ind_PerPool = 1000000\n");
				bw_properties.write("Level_1_2_Region_Mismatch_Tolerance = 1\n");
				bw_properties.write("Level_2_3_Region_Mismatch_Tolerance = 2\n");
				bw_properties.write("Level_3_4_Region_Mismatch_Tolerance = 5\n");
				bw_properties.write("AEM_Maximum_Level = 4\n");
				bw_properties.write("BFS_Mismatch_Tolerance = 6\n");
				bw_properties.write("AEM_Iterations_Max = 200\n");
				bw_properties.write("AEM_Convergence_Cutoff = 0.00001\n");
				bw_properties.write("AEM_Zero_Cutoff = 0.00001\n");
				bw_properties.write("AEM_Regional_Cross_Pool_Freq_Cutoff = 0.01\n");
				bw_properties.write("AEM_Regional_HapSetSize_Max = 50\n");
				bw_properties.write("AEM_Regional_HapSetSize_Min = 3\n");
				bw_properties.write("IF_0_0 = 0.1\n");
				bw_properties.write("IF_Denominator_0= 10.0\n");
				bw_properties.write("Rscript_path ="+this.Rscript_path+"\n");
				bw_properties.write("Regression_Distance_Max_Weight = 2.0\n");
				bw_properties.write("Regression_Coverage_Weight = 2.0\n");
				bw_properties.write("Regression_One_Vector_Weight = 5.0 \n");
				bw_properties.write("Regression_Hap_MAF_Weight = 2.0 \n");
				bw_properties.write("Regression_Hap_LD_Weight = 1.0 \n");	
				bw_properties.write("Regression_Mismatch_Tolerance = 7 \n");
				bw_properties.write("Maximum_Selected_HapSetSize = 25\n");
				bw_properties.write("Regression_Gamma_Min = 0.0001\n");
				bw_properties.write("Regression_Gamma_Max = 0.1\n");
				bw_properties.write("Regression_n_Gamma = 10\n");
				bw_properties.write("Regression_Maximum_Regions = 3\n");
				bw_properties.write("Sequencing_Technology = paired-end reads\n");
				bw_properties.write("Number_Threads = 3\n");															
				bw_properties.close();

				bw.write("properties="+this.main_dir + "/"+ project_name+"/input/PHX.properties\n");
				bw.write("poolhapx="+this.phx_jar+"\n");
				bw.write("$java -jar $poolhapx format $properties\n");
				bw.write("$java -jar $poolhapx gc $properties\n");
				bw.write("$java -jar $poolhapx aem $properties\n");
				bw.write("$java -jar $poolhapx l0l1 $properties\n");
	        	bw.close();	    
    	
	}

	public static void main(String[] args)throws IOException {
		String parameter= args[0];
		ScriptForPHX ps=new ScriptForPHX(parameter);
		ps.ComparisonSim();

	}

}
