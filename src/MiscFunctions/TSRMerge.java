package MiscFunctions;

import MiscFunctions.HapConfig;

public class TSRMerge {

	public static void main(String[] args) {
		String inter_dir = args[0];
		String prefix = args[1];
		String out_dir = args[2];
		int pools = Integer.parseInt(args[3]);
		System.out.println("This is combination " + prefix.split("_")[0] + " simulation " + prefix.split("_")[0] + ".");
		
		HapConfig[] final_local_haps = new HapConfig[pools]; 
		for (int p = 0; p < pools; p++) {
			final_local_haps[p] = new HapConfig(inter_dir + prefix + "_p" + p + "_final_outcome.txt", null);
			System.out.println(final_local_haps[p].num_global_hap + " haplotypes from pool " + p + " have been loaded.");
		}
		HapConfig final_reconstruction = new HapConfig(final_local_haps); 
		System.out.println("Haplotype reconstruction finished. There are " + final_reconstruction.num_global_hap + " global haplotypes.");
		final_reconstruction.write2files(out_dir + prefix + ".inter_freq_vars.txt", out_dir + prefix + ".intra_freq.txt", out_dir + prefix + ".prop_recon.txt", "string", false);
	}
}
