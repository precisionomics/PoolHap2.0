package Viral_Reconstructions_Tools;

import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.HashSet;
import Viral_Reconstructions_Tools.HapConfig_inter_file;

public class compare_inter_file {
	   
	
	public static double[] global_hap_evaluator(String orig_inter_file, 
			String recon_inter_file, double quasi_cutoff, String dir_prefix,
			String project_name) throws IOException, InterruptedException {
		
		HapConfig_inter_file orig_haps = new HapConfig_inter_file(orig_inter_file);
		HapConfig_inter_file recon_haps = new HapConfig_inter_file(recon_inter_file);
		int max_pos_diff = (int)Math.floor(orig_haps.num_loci*quasi_cutoff);
		double diff_ct = 0.0;
		double diff_abs = 0;
		int num_accurate = 0;
		HashSet<Integer> quasi_indices = new HashSet<Integer>();
		
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
		
        for (int ori_h_index = 0; ori_h_index <orig_haps.num_global_hap; ori_h_index++) {
            int h_ori = orig_haps.hapID2index.get(orig_haps.hap_IDs[ori_h_index]);
            int min_hap_index = 0;
            int min_diff = orig_haps.num_loci;
            for (int h_recon = 0; h_recon < recon_haps.num_global_hap; h_recon++) {
                int diff = 0;
                for (int p = 0; p < orig_haps.num_loci; p++) {
                    if (Double.compare(orig_haps.global_haps[h_ori][p],
                        recon_haps.global_haps[h_recon][p]) != 0)
                        diff++;
                }
                pairwise_diff_num[ori_h_index][h_recon] = diff;
                if (diff < min_diff) {
                    min_hap_index = h_recon;
                    min_diff = diff;
                }
            } // end_of_for(h_recon)
            if (min_diff <= max_pos_diff)
                num_accurate++; 
            min_diff_hap[ori_h_index] = min_hap_index;
            min_diff_ID[ori_h_index] = recon_haps.hap_IDs[min_hap_index];
            min_diff_pos[ori_h_index] = min_diff;
            diff_ct += min_diff; // Cumulative min_diff
            double hid_quasi_freq = 0.0; // combine "similar" haps for total hap-freq
            // TODO: [Quan] the same reconstructed hap may contribute to multiple hid_quasi_freq
            for (int h_recon = 0; h_recon < recon_haps.num_global_hap; h_recon++) {
                if (pairwise_diff_num[ori_h_index][h_recon] <= min_diff)
                	hid_quasi_freq += recon_haps.global_haps_freq[h_recon];
            }
            min_diff_freq[ori_h_index] =
                orig_haps.global_haps_freq[h_ori] - hid_quasi_freq;
            diff_abs += Math.abs(min_diff_freq[ori_h_index]);
        }// end_of_for(ori_h_index)
        for (int h_id = 0; h_id < orig_haps.num_global_hap; h_id++) {
            for (int hr = 0; hr < recon_haps.num_global_hap; hr++) {
                if (pairwise_diff_num[h_id][hr] <= max_pos_diff) {
                    max_diff_haps[h_id]++;
                    quasi_indices.add(hr); 
                }
            }
        }
        double freq_tot_wt = 0.0;
        for (int hr : quasi_indices)
            freq_tot_wt += recon_haps.global_haps_freq[hr];

        PrintWriter pw = new PrintWriter(
            new FileWriter(dir_prefix + "_" + quasi_cutoff + "_global_haplotypes.result.txt", true));
        pw.append("Project_Name"+"\t"+"Ori_Hap_ID"+"\t"+"Closest_Recon_Hap_ID"+"\t"
            +"Min_diff_Pos"+"\t"+"Min_freq_diff"+"\t"+"Num_of_Recon_Meet_Cutoff\n");
        for (int h_index = 0; h_index < orig_haps.num_global_hap; h_index++) {
            pw.append(project_name + "\t" + orig_haps.hap_IDs[h_index] + "\t"
                + min_diff_ID[h_index] + "\t" + min_diff_pos[h_index]
                + "\t" + min_diff_freq[h_index] + "\t" + max_diff_haps[h_index] + "\n");
            // also_min_diff[h_id] + "\t" + min_diff_freq_prop[h_id] + "\t" +
        }
        pw.close();
        return new double[] {
            (double) num_accurate / orig_haps.num_global_hap,
            diff_ct /orig_haps.num_global_hap,
            diff_abs / orig_haps.num_global_hap,
            freq_tot_wt,
            orig_haps.num_global_hap,
            recon_haps.num_global_hap
        };
        
	}
	
	
	
	public static void main(String[] args) throws IOException, InterruptedException{
    	String project_name="0_0";//args[0];
    	double quasi_cutoff=0.01;//Double.parseDouble(args[1]); // "0.01"
    	String gs_dir="D:\\PhD-Studying\\Informatics\\Project\\HIV project\\PoolHapX_testing\\gold_standard\\";//args[2];//
    	String output_dir="D:\\PhD-Studying\\Informatics\\Project\\HIV project\\PoolHapX_testing\\output\\";//args[3];// 
        String orig_inter_file=gs_dir+project_name+"_haps.inter_freq_vars.txt";
        String recon_inter_file=output_dir+project_name+".inter_freq_vars.txt";
        String output_files_prefix=output_dir+project_name;
		int num_of_pools = 20;//Integer.parseInt(args[1]);
		double[] multi_pool_record = new double[6];
		multi_pool_record = global_hap_evaluator(orig_inter_file, 
				recon_inter_file, quasi_cutoff, output_files_prefix,
				project_name);
		
		PrintWriter pw2 = new PrintWriter(
	            new FileWriter(output_files_prefix + "_" + quasi_cutoff + "_global_haps_average_results.txt", true));
	        pw2.append("## parameters: cut-off = " + quasi_cutoff + "\n");
	        pw2.append("#Project_Name\t"
	                + "Mean_OH_Recovery\t"
	                + "Mean_Mean_Min_Seq_Diff\t" // mean of means
	                + "Mean_Mean_Min_Freq_Diff\t"
	                + "Mean_Valid_QS_Freq_Sum\t"
	                + "Total_OH\t"
	                + "Total_RH\n");
	            pw2.append(project_name + "\t" + num_of_pools + "\t");
	            for (int i = 0; i < 6; i++) {
	                pw2.append(multi_pool_record[i] + "\t");
	            }
	            pw2.append("\n");
	        pw2.close();
	}

}
