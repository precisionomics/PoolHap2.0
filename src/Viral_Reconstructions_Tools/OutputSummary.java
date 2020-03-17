package Viral_Reconstructions_Tools;

import java.io.BufferedReader;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;


public class OutputSummary {
	
	public static void write_output(String project_name,  int num_index,int num_loci, String input_dir, 
			String output_dir,String tool_name) throws IOException {
		 String input_file_dir;
		 
		 if(tool_name.equals("PoolHapX")) {
			 input_file_dir= input_dir +"/"+ project_name+"_"+num_index+"/output";
		 }else {
			 input_file_dir= input_dir +"/"+ project_name+"_"+num_index;
		 }
				 
		 BufferedReader br_mcc = new BufferedReader(new FileReader(input_file_dir+"/MCC.result"));
		 BufferedReader br_jsd = new BufferedReader(new FileReader(input_file_dir+"/JSD.result"));
		 PrintWriter pw = new PrintWriter(new FileWriter(output_dir+"/"
		            +project_name + ".output_summary.txt", true)); 
		 //pw.append("Tool"+"\t"+"Loci"+"\t"+"MCC"+"\t"+"JSD\n");
		 String currLine_mcc = br_mcc.readLine(); 
		 String currLine_jsd = br_jsd.readLine(); 
		 while (currLine_mcc!=null || currLine_jsd!=null) {
			 String[] tmpVarPos_mcc = currLine_mcc.split("\t");
			 String[] tmpVarPos_jsd = currLine_jsd.split("\t");
			 if(!tmpVarPos_mcc[tmpVarPos_mcc.length-1].equals("0.0")||!tmpVarPos_jsd[tmpVarPos_jsd.length-1].equals("1.0")) {
				 pw.append(tool_name+"\t"+num_loci+"\t"+
						 tmpVarPos_mcc[tmpVarPos_mcc.length-1]+"\t"+
						 tmpVarPos_jsd[tmpVarPos_jsd.length-1]+"\n");
			 }
			 currLine_mcc = br_mcc.readLine(); 
			 currLine_jsd = br_jsd.readLine(); 
		 }
		 br_mcc.close();
		 br_jsd.close();
		 pw.close();
		
	}
	public static void write_output_loop(int num_index_50,int num_index_100,int num_index_200,String output_dir,String model_name) throws IOException {
		int[] pools = new  int [] {25,50};
		int[] depths = new  int  [] { 100, 250,1000};
		int[] freq_cutoff = new int [] {0,2};
		int[] loci=new int[] {50,100,200};
		String[] tools = new String[] {"PoolHapX", "CliqueSNV","PredictHaplo","TenSQR",};
		for (int i =0; i< pools.length; i++) {
			for (int j =0; j< depths.length; j++) {
				for(int h=0;h<freq_cutoff.length;h++) {
					String project_name = "freq_"+freq_cutoff[h]+"_pool_"+ Integer.toString(pools[i])
					+ "_dep_"+ Integer.toString(depths[j]);
					 PrintWriter pw = new PrintWriter(new FileWriter(output_dir+"/"
					            +project_name + ".output_summary.txt", true)); 
					 pw.append("Tool"+"\t"+"Loci"+"\t"+"MCC"+"\t"+"JSD\n");
					 pw.close();
					 for(int l=0;l<loci.length;l++) {
						 for(int t=0;t<tools.length;t++) {
							 String input_dir;
							 int num_index=0;
							 if(loci[l]==50) {
								 num_index=num_index_50;
							 }else if(loci[l]==100){
								 num_index=num_index_100;
							 }else if(loci[l]==200) {
								 num_index=num_index_200;
							 }
							 
							 if(model_name.equals("non_evolution")) {
								 if(tools[t].equals("PoolHapX")) {
									 input_dir = "/export/home/jhe/project/Viral_reconstruction/SLiM/SLiM"
											 +"/migration/sim/non_evolution/"+loci[l]+"_loci"; 
								 }else {
									 input_dir = "/export/home/jhe/project/Viral_reconstruction/Other_Tools/"
										 +tools[t]+"/migration/non_evolution/"+loci[l]+"_loci"; 
								 }
								 write_output(project_name, num_index, loci[l],input_dir,output_dir,tools[t]);
							 }else if(model_name.equals("negative_fitness")) {
								 if(tools[t].equals("PoolHapX")) {
									 input_dir = "/export/home/jhe/project/Viral_reconstruction/SLiM/SLiM"
											 +"/migration/sim/negative_fitness/"+loci[l]+"_loci"; 
								 }else {
									 input_dir = "/export/home/jhe/project/Viral_reconstruction/Other_Tools/"
										 +tools[t]+"/migration/negative_fitness/"+loci[l]+"_loci"; 
								 }
								 write_output(project_name, num_index, loci[l],input_dir,output_dir,tools[t]);
							 }else if(model_name.equals("positive_fitness")) {
								 if(tools[t].equals("PoolHapX")) {
									 input_dir = "/export/home/jhe/project/Viral_reconstruction/SLiM/SLiM"
											 +"/migration/sim/positive_fitness/"+loci[l]+"_loci"; 
								 }else {
									 input_dir = "/export/home/jhe/project/Viral_reconstruction/Other_Tools/"
										 +tools[t]+"/migration/positive_fitness/"+loci[l]+"_loci"; 
								 }
								 write_output(project_name, num_index, loci[l],input_dir,output_dir,tools[t]);
							 }else if(model_name.equals("sweep")) {
								 if(tools[t].equals("PoolHapX")) {
									 input_dir = "/export/home/jhe/project/Viral_reconstruction/SLiM/SLiM"
											 +"/sweep/sim/"+loci[l]+"_loci"; 
								 }else {
									 input_dir = "/export/home/jhe/project/Viral_reconstruction/Other_Tools/"
										 +tools[t]+"/sweep/"+loci[l]+"_loci"; 
								 }
								 write_output(project_name, num_index, loci[l],input_dir,output_dir,tools[t]);
							 }else if(model_name.equals("sweep2")) {
								 if(tools[t].equals("PoolHapX")) {
									 input_dir = "/export/home/jhe/project/Viral_reconstruction/SLiM/SLiM"
											 +"/sweep2/sim/"+loci[l]+"_loci"; 
								 }else {
									 input_dir = "/export/home/jhe/project/Viral_reconstruction/Other_Tools/"
										 +tools[t]+"/sweep2/"+loci[l]+"_loci"; 
								 }
								 write_output(project_name, num_index, loci[l],input_dir,output_dir,tools[t]);
							 }
							 
						 }
					 }
				}
			}
		}
	}
	

	public static void main(String[] args) throws IOException {
		//String input_dir = args[0];
		String output_dir =args[0];
		String model_name =args[1];
		//int num_loci = Integer.parseInt(args[3]);
		int num_index_50 = Integer.parseInt(args[2]);
		int num_index_100 = Integer.parseInt(args[3]);
		int num_index_200 = Integer.parseInt(args[4]);
		write_output_loop(num_index_50,num_index_100,num_index_200, output_dir,model_name);
	}

}
