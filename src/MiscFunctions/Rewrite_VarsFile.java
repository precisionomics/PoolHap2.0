package MiscFunctions;

import java.io.BufferedReader;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.HashMap;

public class Rewrite_VarsFile {
	
	public static void compare_loci(String gs_inter_file, String ori_vars_file
			) throws IOException, InterruptedException{
    	 
    	 BufferedReader br_ori_vars = new BufferedReader(new FileReader(
    			 ori_vars_file));
    	 BufferedReader br_gs_inter = new BufferedReader(new FileReader(
    			 gs_inter_file));
		 HashMap<String, String[]> ori_variant_file = new HashMap<>();
		 int num_pools;
		 String header_line;
    	 //Read original_vars_file, and generate ori_variant_position_list
		 String curr_ori_vars = br_ori_vars.readLine(); // read Pool_ID line
		 header_line = curr_ori_vars;
		 num_pools=curr_ori_vars.split("\t").length-1;
		 curr_ori_vars = br_ori_vars.readLine();
		 while(curr_ori_vars !=null){
			 String[] line_split = curr_ori_vars.split("\t");
			 String[] freq_per_pool = new String[num_pools];
			 for(int i=1;i<line_split.length;i++) {
				freq_per_pool[i-1]=line_split[i];
		     }
			 ori_variant_file.put(line_split[0],freq_per_pool );
			 curr_ori_vars = br_ori_vars.readLine();
		}
		PrintWriter pw_vars = new PrintWriter(new FileWriter(ori_vars_file, false)); 
		pw_vars.append(header_line);
		pw_vars.append("\n");
		String curr_gs_inter = br_gs_inter.readLine();// read hap_ID
		curr_gs_inter = br_gs_inter.readLine();// read freq_line
		curr_gs_inter = br_gs_inter.readLine();
		while(curr_gs_inter!=null) {
			String var_postion = curr_gs_inter.split("\t")[0];
			if(ori_variant_file.containsKey(var_postion)) {
				pw_vars.append(var_postion);
				for(int p=0;p<num_pools;p++) {
					pw_vars.append("\t"+ori_variant_file.get(var_postion)[p]);
				}
				pw_vars.append("\n");
			}else {
				pw_vars.append(var_postion);
				for(int p=0;p<num_pools;p++) {
					pw_vars.append("\t"+"0.0");
				}
				pw_vars.append("\n");
			}
			curr_gs_inter = br_gs_inter.readLine();
		}
		br_gs_inter.close();
		br_ori_vars.close();
		pw_vars.close();
		
	} 
	
	public static void main(String[] args)throws IOException, InterruptedException{
		String gs_inter_file = args[0];
		String ori_vars_file = args[1];
		compare_loci(gs_inter_file, ori_vars_file);
	}

}
