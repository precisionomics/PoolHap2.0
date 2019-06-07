package MiscFunctions;

import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.HashSet;

import MiscFunctions.HapConfigTenSQR;

public class Orig2Recons2 {

    public static HashSet<Integer> multipool_quasispecies;

    public static void main(String[] args) throws IOException {
        String gs_dir = args[0] + "/";
        String out_dir = args[1] + "/";
        String prefix = args[2];    // This is also the simulation number.
        int num_pools = Integer.parseInt(args[3]);
        double quasi_cutoff = Double.parseDouble(args[4]);

        multipool_quasispecies = new HashSet<Integer>();
        HapConfigTenSQR orig_haps = new HapConfigTenSQR(gs_dir + prefix + "_haps.inter_freq_vars.txt", gs_dir + prefix + "_haps.intra_freq.txt");
        HapConfigTenSQR recon_haps = new HapConfigTenSQR(out_dir + prefix + ".inter_freq_vars.txt", out_dir + prefix + ".intra_freq.txt");
        double[][] multi_pool_record = new double[20][4];
        for (int p = 0; p < num_pools; p++) multi_pool_record[p] = pool_evaluator(orig_haps, recon_haps, quasi_cutoff, p, out_dir + prefix);

        PrintWriter pw1 = new PrintWriter(new FileWriter(out_dir + "PD_multipool_extended_results.txt", true));
        double[] multi_pool_results = new double[4];
        for (int p = 0; p < num_pools; p++) {
            pw1.append(prefix + "\t" + p + "\t");
            for (int i = 0; i < 4; i++) {
                pw1.append(multi_pool_record[p][i] + "\t");
                multi_pool_results[i] += multi_pool_record[p][i];
            }
            pw1.append("\n");
        }
        pw1.close();

        PrintWriter pw2 = new PrintWriter(new FileWriter(out_dir + "PD_multipool_aggregated_results.txt", true));
        pw2.append(prefix + "\t");
        for (int i = 0; i < 4; i++) pw2.append(multi_pool_results[i] / num_pools + "\t");
        pw2.append(orig_haps.num_global_hap + "\t" + recon_haps.num_global_hap + "\t" + Double.toString((double) multipool_quasispecies.size() / recon_haps.num_global_hap) + "\n");
        pw2.close();
    }

    // For each simulation, print 1) one multi-pool aggregated results file, 2) one line summarizing the multi-pool aggregated results to be put together across all simulations. So, we will also need the simulation.
    // 1) For each haplotype, the pool, the OHID, the closest RHID, the variant difference, the frequency difference, number of quasispecies.
    // 2) Need to track the RHID of quasispecies. Calculate the summed in-pool frequencies of all qualified quasispecies RHID.
    public static double[] pool_evaluator(HapConfigTenSQR orig_haps, HapConfigTenSQR recon_haps, double quasi_cutoff, int pool, String dir_prefix) throws IOException {
        int max_pos_diff = (int) Math.floor(orig_haps.num_loci * quasi_cutoff); // Number of differing positions that still count as quasispecies.
        double diff_ct = 0.0;
        double diff_abs = 0;
        // double diff_prop = 0;
        int num_accurate = 0;
        HashSet<Integer> quasi_indices = new HashSet<Integer>();

        // Figure out which original haplotypes are in the pool.
        ArrayList<Integer> inpool_hap_ids = new ArrayList<Integer>();
        for (int ho = 0; ho < orig_haps.num_global_hap; ho++)
            if (orig_haps.in_pool_haps_freq[ho][pool] != 0) inpool_hap_ids.add(ho);
        int num_inpool = inpool_hap_ids.size();

        // Compare all original haplotypes to all reconstructed haplotypes to determine i) how many OH were accurately reconstructed, ii) the average error between the OH and the closest RH,
        // iii) the average frequency difference between the OH and the closest RH, and iv) the proportion of frequencies in the pool that matched OHs (predicted proportion).
        // The higher the i) fraction, the smaller the ii) number, the lower the iii) number, and the higher the iv) number is, the better the reconstruction.
        int[] min_diff_hap = new int[num_inpool];	// Index of the closest reconstructed haplotype.
        String[] min_diff_ID = new String[num_inpool];
        int[] min_diff_pos = new int[num_inpool]; 	// Number of differing alleles to closest reconstruction.
        double[] min_diff_freq = new double[num_inpool];	// Frequency difference (absolutes) between the number of haplotypes.
        // double[] min_diff_freq_prop = new double[num_inpool];	// Frequency difference (relative to OH frequency) between the number of haplotypes.
        // int[] also_min_diff = new int[num_inpool];
        int[] max_diff_haps = new int[num_inpool];
        int[][] all_diff_num = new int[num_inpool][recon_haps.num_global_hap];  // Number of differing alleles to all reconstruction.
        // int[][] orig_compare = new int[num_inpool][num_inpool];   // Pairwise differences of original haplotypes, for reference.
        for (int h_id = 0; h_id < num_inpool; h_id++) {
            int ho = inpool_hap_ids.get(h_id);
            int min_hap = 0;
            int min_diff = orig_haps.num_loci;
            for (int hr = 0; hr < recon_haps.num_global_hap; hr++) {
                int diff = 0;
                // int unk = 0;
                for (int p = 0; p < orig_haps.num_loci; p++) {
                    if (orig_haps.global_haps[ho][p] != recon_haps.global_haps[hr][p]) {
                        diff++;
                    }
                }
                all_diff_num[h_id][hr] = diff;
                if (diff < min_diff) {
                    min_hap = hr;
                    min_diff = diff;
                }
            }
            if (min_diff <= max_pos_diff) num_accurate++;
            min_diff_hap[h_id] = min_hap;
            min_diff_ID[h_id] = recon_haps.hap_IDs[min_hap];
            min_diff_pos[h_id] = min_diff;
            diff_ct += min_diff;
            double hid_quasi_freq = 0.0;
            for (int hr = 0; hr < recon_haps.num_global_hap; hr++) if (all_diff_num[h_id][hr] <= min_diff) hid_quasi_freq += recon_haps.in_pool_haps_freq[hr][pool];
            min_diff_freq[h_id] = orig_haps.in_pool_haps_freq[ho][pool] - hid_quasi_freq;
            diff_abs += Math.abs(min_diff_freq[h_id]);
            // min_diff_freq_prop[h_id] = (orig_haps.in_pool_haps_freq[ho][pool] - recon_haps.hap_freq[min_hap]) / orig_haps.in_pool_haps_freq[ho][pool];
            // diff_prop += min_diff_freq_prop[h_id];
            //  min_diff_freq_prop[h_id] = (orig_haps.in_pool_haps_freq[ho][pool] - recon_haps.hap_freq[min_hap]) / orig_haps.in_pool_haps_freq[ho][pool];
        }
        for (int h_id = 0; h_id < num_inpool; h_id++) {
            for (int hr = 0; hr < recon_haps.num_global_hap; hr++) {
                // if (all_diff_num[h_id][hr] == min_diff_pos[h_id]) also_min_diff[h_id]++;
                if (all_diff_num[h_id][hr] <= max_pos_diff) {
                    max_diff_haps[h_id]++;
                    quasi_indices.add(hr);
                }
            }
        }
        double freq_tot_wt = 0.0;
        for (int hr : quasi_indices) freq_tot_wt += recon_haps.in_pool_haps_freq[hr][pool];
        multipool_quasispecies.addAll(quasi_indices);

        PrintWriter pw = new PrintWriter(new FileWriter(dir_prefix + "_results.txt", true));
        // pw.append("Original_ID\tClosest_Recon_ID\tNum_Allele_Diff\tFreq_Diff_Abs\tFreq_Diff_Prop\tNum_Also_MinDiff\tNum_Quasi\n"); // \tNum_Unk_Pos
        for (int h_id = 0; h_id < num_inpool; h_id++)
            pw.append(pool + "\t" + inpool_hap_ids.get(h_id) + "\t" + recon_haps.hapID2index.get(min_diff_ID[h_id]) + "\t" + min_diff_pos[h_id] + "\t" + min_diff_freq[h_id] + "\t" + max_diff_haps[h_id] + "\n");
            // also_min_diff[h_id] + "\t" + min_diff_freq_prop[h_id] + "\t" +
        // pw.append("\n");
        // pw.append("Pool\tAverage_Min_Diff\tAverage_Diff_Freq\tAverage_Diff_Freq_Prop\n");
        // double avg_diff_ct = diff_ct / num_inpool;
        // double avg_diff_freq = diff_abs / num_inpool;
        // double avg_diff_prop = diff_prop / num_inpool;
        // pw.append(pool + "\t" + avg_diff_ct + "\t" + avg_diff_freq + "\t" + avg_diff_prop + "\n\n");
        pw.close();
        return new double[] {(double) num_accurate / num_inpool, diff_ct / num_inpool, diff_abs / num_inpool, freq_tot_wt};
    }
}