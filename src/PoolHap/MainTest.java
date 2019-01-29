package PoolHap;

import PoolHap.HapConfig;

public class MainTest {

	public static void main(String[] args) {
		String working_folder0="/home/lmak/Documents/gc/";
		// String input_gc_folder=working_folder0+"input/gc/";
		// String freq_real=working_folder+"/input/freq/simvars.intra_freq.txt";
		String input_freq_gatk=working_folder0+"p.all.ct.txt";
		String parameter_file=working_folder0+"config_v0_6.properties";
		String dc_plan_outfile=working_folder0+"dc_plan_file.txt";
		String hap_segment_outfile=working_folder0+"hap_segment.txt";
		
		try{
			DivideConquer dc_tester = new DivideConquer(input_freq_gatk, working_folder0+"file_paths_linux.txt", parameter_file, dc_plan_outfile);
			HapConfig[] level_I_config = dc_tester.analyze_regions(dc_tester.regions_level_I, parameter_file);
			// HapConfig[] level_II_config = dc_tester.analyze_regions(dc_tester.regions_level_II, parameter_file);
			// HapConfig[] level_II_config = dc_tester.analyze_regions(dc_tester.regions_level_II, parameter_file);
			// dc_tester.combine_level_I_and_II(level_I_config, level_II_config, hap_segment_outfile);
		}catch(Exception e){e.printStackTrace();}		
	}

}
