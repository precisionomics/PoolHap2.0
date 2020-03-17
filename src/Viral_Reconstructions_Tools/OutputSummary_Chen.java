package Viral_Reconstructions_Tools;

import java.io.BufferedReader;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;


public class OutputSummary_Chen {
	
	public static void write_output(String main_dir,String project_name)throws IOException {
		int[] loci = new  int [] {50,100,200};
		PrintWriter pw = new PrintWriter(new FileWriter(main_dir+"/"
	            +project_name + ".output_summary.txt", true)); 
		for (int l =0; l<loci.length ; l++) {
			BufferedReader br = new BufferedReader(new FileReader(main_dir+"/snp_"+
														loci[l]+"/"+project_name+".txt"));
			String currline = br.readLine();
			while(currline!=null) {
				pw.append(currline+"\n");
				currline = br.readLine();
			}
			br.close();
		}
		pw.close();
	}

	public static void main(String[] args) throws IOException{
		String main_dir = args[0];
		int[] pools = new  int [] {25,50};
		int[] depths = new  int  [] { 50,100};
		int[] freq_cutoff = new int [] {0,2};
		for (int i =0; i< pools.length; i++) {
			for (int j =0; j< depths.length; j++) {
				for(int h=0;h<freq_cutoff.length;h++) {
					String project_name = "freq_"+freq_cutoff[h]+"_pool_"+ Integer.toString(pools[i])
					+ "_dep_"+ Integer.toString(depths[j]);
					 PrintWriter pw = new PrintWriter(new FileWriter(main_dir+"/"
					            +project_name + ".output_summary.txt", true)); 
					 pw.append("Tool"+"\t"+"Loci"+"\t"+"MCC"+"\t"+"JSD\n");
					 pw.close();
					 write_output(main_dir,project_name);
				}
			}	
		}
	}
}
