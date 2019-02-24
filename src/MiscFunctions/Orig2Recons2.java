package MiscFunctions;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileOutputStream;
import java.io.FileReader;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.ArrayList;

public class Orig2Recons2 {

	public static void main(String[] args) throws IOException {
		String work_dir = "/home/lmak/Documents/v0.7_test/GC_AEM_Testing/"; // args[0]; //
		SimpleConfig raw_orig_haps = new SimpleConfig(work_dir + "test_haps.inter_freq_vars.txt"); 
		// Test on GS = new SimpleConfig(work_dir + "test_haps.inter_freq_vars.txt"); 
		// GC = new SimpleConfig(work_dir, num_pools, work_dir + "test_vars.intra_freq.txt");
		// rjMCMC = new SimpleConfig(work_dir + "region_0_haps_allpool.txt"); 
		
		ArrayList<Integer> lev2reg = new ArrayList<Integer>();
		BufferedReader dc_read = new BufferedReader(new FileReader(work_dir + "dc_plan_file.txt")); 
		String line = dc_read.readLine(); 
		while (line != null) {
			if (line.contains("Level"))	{
				lev2reg.add(dc_read.readLine().split("\t").length); 
			}
			line = dc_read.readLine(); 
		} dc_read.close();
		// int num_pools = 20; // Integer.parseInt(args[2]);
		System.out.println("Read in the divide-and-conquer plan."); 
		
		ArrayList<int[]> region_all_length = new ArrayList<int[]>();
		ArrayList<int[]> region_var_length = new ArrayList<int[]>();
		ArrayList<double[]> avg_region_diff_pos = new ArrayList<double[]>();
		ArrayList<double[]> avg_region_freq_abs = new ArrayList<double[]>();
		ArrayList<double[]> avg_region_freq_prop = new ArrayList<double[]>();
		SimpleConfig recon_haps = null;
		for (int l = 1; l <= lev2reg.size(); l++) {
			region_all_length.add(new int[lev2reg.get(l - 1)]);
			region_var_length.add(new int[lev2reg.get(l - 1)]);
			avg_region_diff_pos.add(new double[lev2reg.get(l - 1)]); 
			avg_region_freq_abs.add(new double[lev2reg.get(l - 1)]); 
			avg_region_freq_prop.add(new double[lev2reg.get(l - 1)]); 
			for (int r = 0; r < lev2reg.get(l - 1); r++) {
				int avg_region_diff_ct = 0;
				double avg_region_diff_abs = 0;
				double avg_region_diff_prop = 0;
				recon_haps = new SimpleConfig(work_dir + "level_" + l + "_region_" + r + "_RH.inter_freq_vars.txt"); 
				region_all_length.get(l - 1)[r] = recon_haps.locus_positions[recon_haps.num_loci - 1] - recon_haps.locus_positions[0] + 1;
				region_var_length.get(l - 1)[r] = recon_haps.num_loci;
				SimpleConfig orig_haps = raw_orig_haps; 
				if (orig_haps.num_loci != recon_haps.num_loci) { // But if not...
					int r_start = 0; // This is the index of the orig_haps variant position to start comparing at.
					int i = recon_haps.locus_positions[0]; 
					while (orig_haps.locus_positions[r_start] != i) r_start++; 
					int r_end = r_start; 
					int j = recon_haps.locus_positions[recon_haps.num_loci - 1]; 
					while (orig_haps.locus_positions[r_end] != j) r_end++; 
					orig_haps = raw_orig_haps.subsetter(r_start, r_end); 
				}
				int max_pos_diff = (int) (recon_haps.num_loci * 0.1); // Integer.parseInt(args[3]); Number of differing positions that still count as quasispecies.
				int[] min_diff_hap = new int[orig_haps.num_hap];	// Index of the closest reconstructed haplotype. 
				int[] min_diff_pos = new int[orig_haps.num_hap]; 	// Number of differing alleles to closest reconstruction. 
				double[] min_diff_freq = new double[orig_haps.num_hap];	// Frequency difference (absolutes) between the number of haplotypes. 
				double[] min_diff_freq_prop = new double[orig_haps.num_hap];	// Frequency difference (relative to OH frequency) between the number of haplotypes. 
				int[] also_min_diff = new int[orig_haps.num_hap]; 
				int[] max_diff_haps = new int[orig_haps.num_hap];
				int[][] all_diff_num = new int[orig_haps.num_hap][recon_haps.num_hap]; 	// Number of differing alleles to all reconstruction.
				int[][] orig_compare = new int[orig_haps.num_hap][orig_haps.num_hap]; 	// Pairwise differences of original haplotypes, for reference.
				for (int ho = 0; ho < orig_haps.num_hap; ho++) {
					int min_hap = 0;
					int min_diff = orig_haps.num_loci;
					for (int hr = 0; hr < recon_haps.num_hap; hr++) {
						int diff = 0;
						for (int p = 0; p < orig_haps.num_loci; p++)
							if (orig_haps.hap_var_comp[ho][p] != recon_haps.hap_var_comp[hr][p]) diff++;
						all_diff_num[ho][hr] = diff; 
						if (diff < min_diff) {
							min_hap = hr;
							min_diff = diff;
						} 
					}
					min_diff_hap[ho] = min_hap; 
					min_diff_pos[ho] = min_diff;
					avg_region_diff_ct += min_diff; 
					min_diff_freq[ho] = orig_haps.hap_freq[ho] - recon_haps.hap_freq[min_hap];
					min_diff_freq_prop[ho] = (orig_haps.hap_freq[ho] - recon_haps.hap_freq[min_hap]) / orig_haps.hap_freq[ho];
					avg_region_diff_abs += min_diff_freq[ho]; 
					avg_region_diff_prop += min_diff_freq_prop[ho]; 
					for (int hp = 0; hp < orig_haps.num_hap; hp++)
						for (int p = 0; p < orig_haps.num_loci; p++)
							if (orig_haps.hap_var_comp[ho][p] != orig_haps.hap_var_comp[hp][p]) orig_compare[ho][hp]++;
				}
				avg_region_diff_pos.get(l - 1)[r] = avg_region_diff_ct / (double) (orig_haps.num_hap * orig_haps.num_loci); 
				// System.out.println(avg_region_diff_ct + "\t" + orig_haps.num_hap);
				avg_region_freq_abs.get(l - 1)[r] = avg_region_diff_abs / (double) orig_haps.num_hap; 
				// System.out.printf("%.3f\t", avg_region_diff_abs);
				avg_region_freq_prop.get(l - 1)[r] = avg_region_diff_prop / (double) orig_haps.num_hap; 
				
				for (int ho = 0; ho < orig_haps.num_hap; ho++) {
					for (int hr = 0; hr < recon_haps.num_hap; hr++) {
						if (hr != min_diff_hap[ho]) {
							if (all_diff_num[ho][hr] == min_diff_pos[ho]) also_min_diff[ho]++; 
							else if (all_diff_num[ho][hr] <= min_diff_pos[ho] + max_pos_diff) max_diff_haps[ho]++;
						}
					}
				}
				
				// What do I want to report?
				// Orig hap num vs. closest recon hap num + num diff + freq diff + num also min diff + num quasi
				PrintWriter pw = new PrintWriter(new FileOutputStream(new File(work_dir + "level_region_all_summary_analysis.txt"), true));
				pw.append("Level\t" + l + "\tRegion\t" + r + "\tNum_Recon_Hap\t" + recon_haps.num_hap + "\n");
				pw.append("Original_ID\tClosest_Recon_ID\tNum_Allele_Diff\tFreq_Diff_Abs\tFreq_Diff_Prop\tNum_Also_MinDiff\tNum_Quasi\n");
				for (int ho = 0; ho < orig_haps.num_hap; ho++) {
					pw.append(ho + "\t" + min_diff_hap[ho] + "\t" + min_diff_pos[ho] + "\t" + min_diff_freq[ho] + "\t" + min_diff_freq_prop[ho] 
						+ "\t" + also_min_diff[ho] + "\t" + max_diff_haps[ho] + "\n");
				}
				/*
				 * pw.append("\nOriginal Haplotypes Pairwise Comparison\nOriginal_ID"); for (int
				 * ho = 0; ho < orig_haps.num_hap; ho++) pw.append("\t" + ho); for (int ho = 0;
				 * ho < orig_haps.num_hap; ho++) { pw.append("\n" + ho); for (int hp = 0; hp <
				 * orig_haps.num_hap; hp++) pw.append("\t" + orig_compare[ho][hp]); }
				 * pw.append("\n\nOriginal Haplotypes Global Frequency\nOriginal_ID"); for (int
				 * ho = 0; ho < orig_haps.num_hap; ho++) pw.append("\t" + ho);
				 * pw.append("\nFrequency"); for (int ho = 0; ho < orig_haps.num_hap; ho++)
				 * pw.append("\t" + orig_haps.hap_freq[ho]);
				 * 
				 * pw.append("\n\nPosition\t"); for (int p = 0; p < orig_haps.num_loci; p++)
				 * pw.append(orig_haps.locus_positions[p] + "\t"); for (int i = 0; i < 3; i++) {
				 * pw.append("\n" + i + "\t"); for (int p = 0; p < orig_haps.num_loci; p++)
				 * pw.append(orig_haps.hap_var_comp[i][p] + "\t"); pw.append("\n" +
				 * min_diff_hap[i] + "\t"); for (int p = 0; p < orig_haps.num_loci; p++)
				 * pw.append(recon_haps.hap_var_comp[min_diff_hap[i]][p] + "\t");
				 * pw.append("\n\n\n"); }
				 */
				pw.append("\n\n");
				pw.close();
				orig_haps.write_global_file_string(work_dir + "level_" + l + "_region_" + r + "_OH.inter_freq_vars.txt", false);
				System.out.println("Level " + l + " region " + r + " has been analyzed.");
			}
		}
		PrintWriter pw = new PrintWriter(new FileOutputStream(new File(work_dir + "level_region_all_summary_analysis.txt"), true));
		for (int r = 0; r < lev2reg.get(0); r++) pw.append("\t" + r); 
		for (int l = 1; l <= lev2reg.size(); l++) {
			pw.append("\n" + l); 		
			for (int r = 0; r < lev2reg.get(l - 1); r++) pw.append("\t" + region_all_length.get(l - 1)[r]); 
		}
		pw.append("\n\n");
		for (int r = 0; r < lev2reg.get(0); r++) pw.append("\t" + r); 
		for (int l = 1; l <= lev2reg.size(); l++) {
			pw.append("\n" + l); 		
			for (int r = 0; r < lev2reg.get(l - 1); r++) pw.append("\t" + region_var_length.get(l - 1)[r]); 
		}
		pw.append("\n\n");
		for (int r = 0; r < lev2reg.get(0); r++) pw.append("\t" + r); 
		for (int l = 1; l <= lev2reg.size(); l++) {
			pw.append("\n" + l); 		
			for (int r = 0; r < lev2reg.get(l - 1); r++) pw.format("\t%.4f", avg_region_diff_pos.get(l - 1)[r]); 
		}
		pw.append("\n\n");
		for (int r = 0; r < lev2reg.get(0); r++) pw.append("\t" + r); 
		for (int l = 1; l <= lev2reg.size(); l++) {
			pw.append("\n" + l); 		
			for (int r = 0; r < lev2reg.get(l - 1); r++) pw.format("\t%.4f", avg_region_freq_abs.get(l - 1)[r]); 
		}
		pw.append("\n\n");
		for (int r = 0; r < lev2reg.get(0); r++) pw.append("\t" + r); 
		for (int l = 1; l <= lev2reg.size(); l++) {
			pw.append("\n" + l); 		
			for (int r = 0; r < lev2reg.get(l - 1); r++) pw.format("\t%.4f", avg_region_freq_prop.get(l - 1)[r]); 
		}
		pw.close();
	}
} 