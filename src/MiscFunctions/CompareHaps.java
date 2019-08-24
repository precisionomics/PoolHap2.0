package MiscFunctions;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;

import PoolHap.HapConfig;
import PoolHap.Parameters.GenParameters;
import spire.optional.intervalGeometricPartialOrder;

public class CompareHaps {

    /**
     * @author quanlong 2019-07
     * 
     * Compare two HapConfig objects, typically the original (true) haplotypes and reconstructed
     * ones.
     */

    public static HashSet<Integer> multipool_quasispecies;

    /**
     * 
     * @param orig_haps
     * @param recon_haps
     * @param quasi_cutoff
     * @param pool_index
     * @param dir_prefix
     * @return
     * @throws IOException
     * 
     * For each simulation, print 1) one multi-pool aggregated results file, 2) one line summarizing
     * the multi-pool aggregated results to be put together across all simulations. So, we will also
     * need the simulation. 1) For each haplotype, the pool, the OHID, the closest RHID, the variant
     * difference, the frequency difference, number of quasispecies. 2) Need to track the RHID of
     * quasispecies. Calculate the summed in-pool frequencies of all qualified quasispecies RHID.
     *
     * The number of segregating sites and their locations have to be the same between two
     * HapConfigs
     */
    public static void compare_loci(String ori_inter_file, 
    		String recon_inter_file) throws
    IOException, InterruptedException{
    	 BufferedReader br_ori_inter = new BufferedReader(new FileReader(
    			 ori_inter_file));
    	 BufferedReader br_recon_inter = new BufferedReader(new FileReader(
    			 recon_inter_file));
    	 ArrayList<String> ori_variant_position_list = new ArrayList<>();
    	 HashMap<String, String> var2hap = new HashMap<String, String>();
    	 String curr_ori_inter = br_ori_inter.readLine(); // read header line
    	 curr_ori_inter = br_ori_inter.readLine(); // read freq line
    	 curr_ori_inter = br_ori_inter.readLine(); // read third line
    	 while(curr_ori_inter != null) {
    		 String[] var_position_line = curr_ori_inter.split("\t");
    	 	 String[] var_position = var_position_line[0].split(";");
    	 	 ori_variant_position_list.add(var_position[1]);
    	 	 curr_ori_inter = br_ori_inter.readLine();
    	 }
    	 br_ori_inter.close();
    	 String curr_recon_inter = br_recon_inter.readLine();
    	 String[] var_position_line = curr_recon_inter.split("\t");
    	 //bw_new_recon_inter.write(curr_recon_inter);
    	 var2hap.put("header", curr_recon_inter);
    	 curr_recon_inter = br_recon_inter.readLine();
    	 int num_hap = var_position_line.length -1;
    	 //bw_new_recon_inter.write(curr_recon_inter);
    	 var2hap.put("freq", curr_recon_inter);
    	 curr_recon_inter = br_recon_inter.readLine(); //read the third line
    	 while(curr_recon_inter != null) {
    		 var_position_line = curr_recon_inter.split("\t");	 
    	 	 String[] var_position = var_position_line[0].split(";");
    	 	 var2hap.put(var_position[1], curr_recon_inter);
    	 	 curr_recon_inter = br_recon_inter.readLine();
    	 }
    	 br_recon_inter.close();
         PrintWriter pw = new PrintWriter(new FileWriter(recon_inter_file, false)); 
    	 pw.append(var2hap.get("header")+"\n");
    	 pw.append(var2hap.get("freq")+"\n");
    	 for(int p=0;p<ori_variant_position_list.size();p++) {
    		 if(var2hap.containsKey(ori_variant_position_list.get(p))) {
    			 pw.append(var2hap.get(ori_variant_position_list.get(p)));
    			 pw.append("\n");
    		 }else{
        		 pw.append("0;"+ori_variant_position_list.get(p)+";"+
        		    	 ori_variant_position_list.get(p)+";0:1");
        		 for(int h =0; h < num_hap; h++ ) {
        			 pw.append("\t"+"0");
        		 }
        		 pw.append("\n");
    		 }
    	 }
    	 pw.close();
    }
    
    public static double[] single_pool_evaluator(
        HapConfig orig_haps,
        HapConfig recon_haps,
        // int compare_loci, // total number of loci in the comparison
        // int[] o2r_indices, // index mapping between original and reconstructed
        double quasi_cutoff, // a ratio of acceptable differences.
        String pool_ID,
        String dir_prefix) throws IOException {
        if (orig_haps.num_loci != recon_haps.num_loci) {
            System.out.println("Error: orig_haps.num_loci!=recon_haps.num_loci: \n"
                + "Returned without comparison!");
        }
        int ori_p_index = -1;
        for (int id_index = 0; id_index < orig_haps.pool_IDs.length; id_index++) {
            if (orig_haps.pool_IDs[id_index].equals(pool_ID)) {
                ori_p_index = id_index;
                break;
            }
        }
        int recon_p_index = -1;
        for (int id_index = 0; id_index < recon_haps.pool_IDs.length; id_index++) {
            if (recon_haps.pool_IDs[id_index].equals(pool_ID)) {
                recon_p_index = id_index;
                break;
            }
        }
        if (recon_p_index == -1 || ori_p_index == -1) {
            System.out
                .println("Error: " + pool_ID + " can't be found. Returned without comparison!");
        }
        int max_pos_diff = (int) Math.floor(orig_haps.num_loci * quasi_cutoff);
        // Number of differing positions that still count as quasispecies.
        double diff_ct = 0.0;
        double diff_abs = 0;
        int num_accurate = 0;
        HashSet<Integer> quasi_indices = new HashSet<Integer>();

        // Figure out which original haplotypes are in the pool.
        ArrayList<String> none_0_ori_hap_id = new ArrayList<String>();
        for (int ho = 0; ho < orig_haps.num_global_hap; ho++)
            if (orig_haps.in_pool_haps_freq[ho][ori_p_index] > 0) {
                none_0_ori_hap_id.add(orig_haps.hap_IDs[ho]);
            }
        int num_inpool_ori = none_0_ori_hap_id.size();

        // Compare all original haplotypes to all reconstructed haplotypes to determine
        // i) how many OH were accurately reconstructed,
        // ii) the average error between the OH and the closest RH,
        // iii) the average frequency difference between the OH and the closest RH, and
        // iv) the proportion of frequencies in the pool that matched OHs (predicted proportion).

        // NOTE: The higher the i) fraction, the smaller the ii) number, the lower the iii) number,
        // and the higher the iv) number is, the better the reconstruction.

        int[] min_diff_hap = new int[num_inpool_ori];
        // Index of the closest reconstructed haplotype.
        String[] min_diff_ID = new String[num_inpool_ori];
        // Storing the IDs of matched haps in reconstructed HapConfig
        int[] min_diff_pos = new int[num_inpool_ori];
        // Number of differing loci to closest reconstruction.
        double[] min_diff_freq = new double[num_inpool_ori];
        // Frequency difference (absolutes) between the number of haplotypes.

        int[] max_diff_haps = new int[num_inpool_ori];
        int[][] pairwise_diff_num = new int[num_inpool_ori][recon_haps.num_global_hap];
        // Pairwise number of differing loci to all reconstruction.
        for (int ori_h_index = 0; ori_h_index < num_inpool_ori; ori_h_index++) {
            int h_ori = orig_haps.hapID2index.get(none_0_ori_hap_id.get(ori_h_index));
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
            }
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
                    hid_quasi_freq += recon_haps.in_pool_haps_freq[h_recon][recon_p_index];
            }
            min_diff_freq[ori_h_index] =
                orig_haps.in_pool_haps_freq[h_ori][ori_p_index] - hid_quasi_freq;
            diff_abs += Math.abs(min_diff_freq[ori_h_index]);
        }
        for (int h_id = 0; h_id < num_inpool_ori; h_id++) {
            for (int hr = 0; hr < recon_haps.num_global_hap; hr++) {
                if (pairwise_diff_num[h_id][hr] <= max_pos_diff) {
                    max_diff_haps[h_id]++;
                    quasi_indices.add(hr);
                }
            }
        }
        double freq_tot_wt = 0.0;
        for (int hr : quasi_indices)
            freq_tot_wt += recon_haps.in_pool_haps_freq[hr][recon_p_index];
        multipool_quasispecies.addAll(quasi_indices);

        PrintWriter pw = new PrintWriter(
            new FileWriter(dir_prefix + "_" + quasi_cutoff + "_single_pools.result.txt", true));
        pw.append("Pool_ID"+"\t"+"Ori_Hap_ID"+"\t"+"Closest_Recon_Hap_ID"+"\t"
                +"Min_diff_Pos"+"\t"+"Min_freq_diff"+"\t"+"Num_of_Recon_Meet_Cutoff"+"\n");
        for (int h_index = 0; h_index < num_inpool_ori; h_index++) {
            pw.append(pool_ID + "\t" + none_0_ori_hap_id.get(h_index) + "\t"
                + min_diff_ID[h_index] + "\t" + min_diff_pos[h_index]
                + "\t" + min_diff_freq[h_index] + "\t" + max_diff_haps[h_index] + "\n");
            // also_min_diff[h_id] + "\t" + min_diff_freq_prop[h_id] + "\t" +
        }
        // pw.append("\n");
        // pw.append("Pool\tAverage_Min_Diff\tAverage_Diff_Freq\tAverage_Diff_Freq_Prop\n");
        // double avg_diff_ct = diff_ct / num_inpool;
        // double avg_diff_freq = diff_abs / num_inpool;
        // double avg_diff_prop = diff_prop / num_inpool;
        // pw.append(pool + "\t" + avg_diff_ct + "\t" + avg_diff_freq + "\t" + avg_diff_prop +
        // "\n\n");
        pw.close();
        return new double[] {
            (double) num_accurate / num_inpool_ori,
            diff_ct / num_inpool_ori,
            diff_abs / num_inpool_ori,
            freq_tot_wt
        };
    }

    public static void multi_pool_summary(
        String project_name,
        double quasi_cutoff,
        String ori_inter_file,
        String ori_intra_file,
        String recon_inter_file,
        String recon_intra_file,
        String output_files_prefix) throws IOException{

        CompareHaps.multipool_quasispecies = new HashSet<Integer>();
        HapConfig orig_haps = new HapConfig(ori_inter_file, ori_intra_file);
        int num_ori_pools = orig_haps.num_pools;
        HapConfig recon_haps = new HapConfig(recon_inter_file, recon_intra_file);
        double[][] multi_pool_record = new double[num_ori_pools][4];
        for (int p = 0; p < num_ori_pools; p++) {
            multi_pool_record[p] = single_pool_evaluator(orig_haps,
                recon_haps,
                quasi_cutoff,
                orig_haps.pool_IDs[p],
                output_files_prefix);
        }
        PrintWriter pw1 = new PrintWriter(
            new FileWriter(output_files_prefix + "_" + quasi_cutoff + "_extended_results.txt", true));
        pw1.append("## parameters: cut-off = " + quasi_cutoff + "\n");
//        pw1.append("## orig_hap_files: " + ori_inter_file + "\t" + ori_intra_file + "\n");
//        pw1.append("## recon_hap_files: " + recon_inter_file + "\t" + recon_intra_file + "\n");
        pw1.append(
//            "# Project_name\tPool_ID\t" + "Prop_of_OH_recovered\t" + "Ave_dist_btw OH_closest_RH\t"
//                + "Ave_freq_diff_btw OH_closest_RH\t" + "Sum_in-pool_freq_valid_quasispecies\n");
        		"#Project_Name\t"
                + "Pool_ID\t"
                + "OH_Recovery\t"
                + "Mean_Min_Seq_Diff\t"
                + "Mean_Min_Freq_Diff\t"
                + "Valid_QS_Freq_Sum\n");
        double[] multi_pool_results = new double[4];
        for (int p = 0; p < num_ori_pools; p++) {
            pw1.append(project_name + "\t" + orig_haps.pool_IDs[p] + "\t");
            for (int i = 0; i < 4; i++) {
                pw1.append(multi_pool_record[p][i] + "\t");
                multi_pool_results[i] += multi_pool_record[p][i];
            }
            pw1.append("\n");
        }
        pw1.close();

        PrintWriter pw2 = new PrintWriter(
            new FileWriter(output_files_prefix + "_" + quasi_cutoff + "_aggregated_results.txt", true));
        pw2.append("## parameters: cut-off = " + quasi_cutoff + "\n");
//        pw1.append("## orig_hap_files: " + ori_inter_file + "\t" + ori_intra_file + "\n");
//        pw1.append("## recon_hap_files: " + recon_inter_file + "\t" + recon_intra_file + "\n");
//        pw1.append(
//            "# Project_name\tPool_ID\t" + "Prop_of_OH_recovered\t" + "Ave_dist_btw OH_closest_RH\t"
//                + "Ave_freq_diff_btw OH_closest_RH\t" + "Sum_in-pool_freq_valid_quasispecies\n");
        
        pw2.append("#Project_Name\t"
                + "Mean_OH_Recovery\t"
                + "Mean_Mean_Min_Seq_Diff\t" // mean of means
                + "Mean_Mean_Min_Freq_Diff\t"
                + "Mean_Valid_QS_Freq_Sum\t"
                + "Total_OH\t"
                + "Total_RH\t"
                + "Total_Valid_QS_Proportion\n");
        pw2.append(project_name+"\t");
        for (int i = 0; i < 4; i++) {
            pw2.append(multi_pool_results[i] / orig_haps.num_pools + "\t");
        }
        pw2.append(orig_haps.num_global_hap + "\t" + recon_haps.num_global_hap + "\t"
            + Double.toString((double) multipool_quasispecies.size() / recon_haps.num_global_hap)
            + "\t");

        pw2.append("\n");
        pw2.close();

    }
    
    public static void main(String[] args) throws IOException, InterruptedException {
    	String project_name= "0";//args[0];  //for example 0_1
    	double quasi_cutoff= 0.1;//Double.parseDouble(args[1]); // cutoff
    	String gs_dir= "D:\\PhD-Studying\\Informatics\\Project\\HIV project\\PoolHapX_testing\\PHX_Perfect_Data40_2\\gold_standard\\"; //args[2]; //gold_standard dir
    	String output_dir= "D:\\PhD-Studying\\Informatics\\Project\\HIV project\\PoolHapX_testing\\PHX_Perfect_Data40_2\\output\\" ;//args[3]; //out_put dir
        String ori_inter_file=gs_dir+project_name+"_haps.inter_freq_vars.txt";
        String ori_intra_file=gs_dir+project_name+"_haps.intra_freq.txt";
        String recon_inter_file=output_dir+project_name+".inter_freq_haps.txt";
        String recon_intra_file=output_dir+project_name+".intra_freq_haps.txt";
        String output_files_prefix=output_dir+project_name;
        compare_loci(ori_inter_file,recon_inter_file);
        recon_inter_file=output_dir+project_name+".inter_freq_haps.txt";
    	multi_pool_summary(
    	        project_name,
    	        quasi_cutoff,
    	        ori_inter_file,
    	        ori_intra_file,
    	        recon_inter_file,
    	        recon_intra_file,
    	        output_files_prefix);
	}
}
