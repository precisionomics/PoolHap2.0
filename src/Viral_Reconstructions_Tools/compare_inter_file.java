package Viral_Reconstructions_Tools;

import java.io.IOException;
import java.util.ArrayList;


import Viral_Reconstructions_Tools.HapConfig_inter_file;

public class compare_inter_file {
	   
	
	public void global_hap_evaluator(String orig_inter_file, 
			String recon_inter_file, double quasi_cutoff, String dir_prefix,
			int num_pools) throws IOException,InterruptedException  {
		
		HapConfig_inter_file orig_haps = new HapConfig_inter_file(orig_inter_file);
		HapConfig_inter_file recon_haps = new HapConfig_inter_file(recon_inter_file);
		int max_pos_diff = (int)Math.floor(orig_haps.num_loci*quasi_cutoff);
		double diff_ct = 0.0;
		double diff_abs = 0;
		int num_accurate = 0;
		if (orig_haps.num_loci != recon_haps.num_loci) {
			System.out.println("Error: orig_haps.num_loci!=recon_haps.num_loci: \n"
	                + "Returned without comparison!");
		}
		ArrayList<String> ori_hap_id = new ArrayList<String>();
	       for (int ho = 0; ho < orig_haps.num_global_hap; ho++) {
	    	   ori_hap_id.add(orig_haps.hap_IDs[ho]);
	       }
//		int num_inpool_ori = orig_haps.num_global_hap;
		
		int[] min_diff_hap = new int[orig_haps.num_global_hap];
	    // Index of the closest reconstructed haplotype.
	    String[] min_diff_ID = new String[orig_haps.num_global_hap];
	    // Storing the IDs of matched haps in reconstructed HapConfig
	    int[] min_diff_pos = new int[orig_haps.num_global_hap];
	    // Number of differing loci to closest reconstruction.
	    double[] min_diff_freq = new double[orig_haps.num_global_hap];
	    // Frequency difference (absolutes) between the number of haplotypes.

        int[] max_diff_haps = new int[orig_haps.num_global_hap];
        int[][] pairwise_diff_num = new int[orig_haps.num_global_hap][recon_haps.num_global_hap];		
		
		
		
	}
	
	public static void main(String[] args){
		
	}
	
	
	
}
