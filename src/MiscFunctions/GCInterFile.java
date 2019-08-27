package MiscFunctions;

import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;

import spire.optional.intervalGeometricPartialOrder;

public class GCInterFile {
	
	public static void generate_inter_file(String gcf_file,String vars_intra_file) 
			throws IOException, InterruptedException  {
		
		ArrayList<String> loci_string_list=new ArrayList<String>();
		//hap_full_string_list ignoring "?", and take into account all possible haploypes
		ArrayList<String> hap_full_string_list=new ArrayList<String>();	
		//hap_partial_string_list only takes those complete haplotypes into account(those without ?)
		ArrayList<String> hap_partial_string_list=new ArrayList<String>();	
		ArrayList<String> final_hap_string_list=new ArrayList<String>();
		
		BufferedReader br_gcf = new BufferedReader(new FileReader(gcf_file));
		BufferedReader br_var = new BufferedReader(new FileReader(vars_intra_file));
		String curr_var_line = br_var.readLine();
		curr_var_line = br_var.readLine(); 
		while(curr_var_line!=null) {
			String[] curr_pos_arr = curr_var_line.split("\t");
			loci_string_list.add(curr_pos_arr[0]);
			curr_var_line = br_var.readLine();
		}
		
		String curr_gcf_line = br_gcf.readLine();
		while(curr_gcf_line!=null) {
			String[] hap_array = curr_gcf_line.split(" ");
			String[] every_position = curr_gcf_line.split("");
			ArrayList<String> tmp_hap_list=new ArrayList<String>();	
			String new_full_hap = "";
			for(int i=0;i<every_position.length;i++) {
				tmp_hap_list.add(every_position[i]);
				if(every_position[i].equals("0")||every_position[i].equals("1")) {
					new_full_hap=new_full_hap+every_position[i];
				}
			}
			hap_full_string_list.add(new_full_hap); 
			if(!tmp_hap_list.contains("?")) {
				for(int i=0;i<tmp_hap_list.size();i++) {
					
				}
			}
			
			curr_gcf_line = br_gcf.readLine();
		}
		
		
		
	}
	
	

	public static void main(String[] args)
			throws IOException, InterruptedException {
		String project_name= args[0];//"0";
		String main_dir=args[1]; 
		int num_pools = Integer.parseInt(args[2]);
		String gs_dir = main_dir + "/gold_standard/";
    	String output_dir = main_dir + "/output/";
    	String inter_dir = main_dir + "/intermediate/";
    	String gcf_dir= main_dir + "/intermediate/gcf/";
    	String vars_intra_file = inter_dir+project_name+"_vars.intra_freq.txt";
    	for(int p=0;p<num_pools;p++) {
    		String gcf_file=gcf_dir+project_name+"_p"+p+".gcf";
    		generate_inter_file(gcf_file,vars_intra_file);
    	}
    	


	}

}
