package PoolHap;

import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.HashMap;

import PoolHap.GraphColoring;
import PoolHap.HapConfig;
import PoolHap.DivideConquer_Testing;

import java.time.format.DateTimeFormatter;  
import java.time.LocalDateTime;    

public class MainTest {

	public static void main(String[] args) throws IOException {
		String work_dir = "/home/lmak/Documents/v0.7_test/";
		String input_dir = "/home/lmak/Dropbox/University of Calgary/PoolHap_Testing/FullSimulator2_Testing/Pool VEFs/";  
		// String input_gc_folder=working_folder0+"input/gc/";
		// String freq_real=working_folder+"/input/freq/simvars.intra_freq.txt";
		int num_pools = 20; 
		String observed_freq = work_dir + "test_vars.intra_freq.txt";
		String parameter_file = work_dir + "config_v0_6.properties";
		String dc_plan_outfile=work_dir + "dc_plan_file.txt";
		String hap_segment_outfile=work_dir + "hap_segment.txt";
		
		HashMap<Integer, Integer> pos_dict = new HashMap<Integer, Integer>();
		Integer pos_index = 0; 
	    BufferedReader br = new BufferedReader(new FileReader(observed_freq));
	    String currLine = br.readLine(); // Skip header. 
	    currLine = br.readLine();
	    while (currLine != null) {
	    	pos_dict.put(Integer.parseInt(currLine.split(";")[1]), pos_index); 
	    	pos_index++; 
	    	currLine = br.readLine();
	    }
	    br.close();
	    
	    DateTimeFormatter dtf = DateTimeFormatter.ofPattern("yyyy/MM/dd HH:mm:ss");  
	    
	    /* 
	     * System.out.println("PHX Initiation: " + dtf.format(LocalDateTime.now()));  
		PrintWriter in_list = new PrintWriter(work_dir + "p.in.list"); 
		for (int p = 0; p < num_pools; p++) {
			GraphColoring.gc(input_dir + "p" + p + ".vef", pos_dict, work_dir + "p" + p + ".in");
			System.out.println("Graph colouring for pool " + p + " is finished.");
			in_list.println(work_dir + "p" + p + ".in"); 
		}
		in_list.close();
	    System.out.println("\nGC Finished: " + dtf.format(LocalDateTime.now()) + "\n"); 
	    */  
		
		try{
			DivideConquer_Testing dc_tester = new DivideConquer_Testing(observed_freq, work_dir + "p.in.list", parameter_file, dc_plan_outfile);
		    System.out.println("DC Finished: " + dtf.format(LocalDateTime.now()) + "\n");  
			HapConfig[] level_I_config = dc_tester.analyze_regions(dc_tester.regions_level_I, parameter_file);
		    System.out.println("Solving Finished: " + dtf.format(LocalDateTime.now()) + "\n");  
			// HapConfig[] level_II_config = dc_tester.analyze_regions(dc_tester.regions_level_II, parameter_file);
			// HapConfig[] level_II_config = dc_tester.analyze_regions(dc_tester.regions_level_II, parameter_file);
			// dc_tester.combine_level_I_and_II(level_I_config, level_II_config, hap_segment_outfile);
		}catch(Exception e){e.printStackTrace();}		
	}

}
