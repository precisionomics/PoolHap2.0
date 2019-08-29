package MiscFunctions;

import PoolHap.HapConfig;

public class test {

	public static void main(String[] args) {
		String project_name= args[0];  
    	double quasi_cutoff= Double.parseDouble(args[1]);
    	String gs_dir= args[2]+"/"; 
    	String output_dir= args[3]+"/"; 
    	String recon_inter_file=output_dir+project_name+".inter_freq_haps.txt";
        String recon_intra_file=output_dir+project_name+".intra_freq_haps.txt";
        String ori_inter_file=gs_dir+project_name+"_haps.inter_freq_vars.txt";
        String ori_intra_file=gs_dir+project_name+"_haps.intra_freq.txt";
		HapConfig recon_haps = new HapConfig(recon_inter_file, recon_intra_file);
//		HapConfig ori_haps = new HapConfig(ori_inter_file, ori_intra_file);

	}

}
