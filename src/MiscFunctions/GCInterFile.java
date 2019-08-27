package MiscFunctions;

import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;

import org.apache.hadoop.yarn.webapp.hamlet.Hamlet.P;

import spire.optional.intervalGeometricPartialOrder;


public class GCInterFile {
	
	public static void generate_inter_file(String gcf_file,String vars_intra_file) 
			throws IOException, InterruptedException  {
		
		ArrayList<String> loci_string_list=new ArrayList<String>();
		//hap_full_string_list ignoring "?", and take into account all possible haploypes
		ArrayList<String> hap_full_string_list=new ArrayList<String>();	
		ArrayList<ArrayList<String>> full_hap_seq_listlist=new ArrayList<ArrayList<String>>();
		//patial_hap_seq_listlist only takes those complete haplotypes into account(those without ?)
		ArrayList<ArrayList<String>> patial_hap_seq_listlist=new ArrayList<ArrayList<String>>();
		ArrayList<String> final_hap_full_string_list=new ArrayList<String>();
		
		BufferedReader br_gcf = new BufferedReader(new FileReader(gcf_file));
		BufferedReader br_var = new BufferedReader(new FileReader(vars_intra_file));
		String curr_var_line = br_var.readLine();
		curr_var_line = br_var.readLine(); 
		while(curr_var_line!=null) {
			String[] curr_pos_arr = curr_var_line.split("\t");
			loci_string_list.add(curr_pos_arr[0]);
			ArrayList<String> new_var_list=new ArrayList<String>();
			patial_hap_seq_listlist.add(new_var_list); // num_var * num_hap
			full_hap_seq_listlist.add(new_var_list); // num_var * num_hap
			curr_var_line = br_var.readLine();
		}
		br_var.close();
		
		
		String curr_gcf_line = br_gcf.readLine();
		while(curr_gcf_line!=null) {
			String[] hap_array = curr_gcf_line.split("\t");
			String[] every_position = hap_array[0].split("");
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
				int var_pos=0;
				for(int i=0;i<tmp_hap_list.size();i++) {
					if(every_position[i].equals("0")||every_position[i].equals("1")) {
						patial_hap_seq_listlist.get(var_pos).add(every_position[i]);
						var_pos++;
					}
				}
			}
			curr_gcf_line = br_gcf.readLine();
		}
		System.out.println(patial_hap_seq_listlist);
		br_gcf.close();
		// get rid of the identical haplotype in the hap_full_string_list
		for(int h=0;h<hap_full_string_list.size();h++) {
			if(!final_hap_full_string_list.contains(hap_full_string_list.get(h))) {
				final_hap_full_string_list.add(hap_full_string_list.get(h));
			}
		}
		// hap_full_string_list transfer into full_hap_seq_listlist
		for(int h=0;h<final_hap_full_string_list.size();h++) {
			
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
    	
    	String gcf_file=gcf_dir+"0_p0.gcf";
    	generate_inter_file(gcf_file,vars_intra_file);
    	
//    	for(int p=0;p<num_pools;p++) {
//    		String gcf_file=gcf_dir+project_name+"_p"+p+".gcf";
//    		generate_inter_file(gcf_file,vars_intra_file);
//    	}
    	


	}

}
