package PoolHap;

import java.io.FileInputStream;
import java.io.IOException;
import java.io.InputStream;
import java.util.Properties;

public class Parameters {

	public static class DivideParameters extends Parameters {
		public double gap_inpool_cutoff;  // a ratio
		public double gap_all_pool_cutoff;// a ratio
		public int min_level_I_region_size;
		public int max_level_I_region_size;
		public int min_level_II_region_size;
		public int max_level_II_region_size;
		public int max_num_rounds_forming_initial_haps; // because of trimming, we may not be able to form all full haplotypes. Cutoff to halt that.
		public int max_num_haps;	// the max number of haps that will be allowed in the AEM algorithm. 
									//	(When there are too many haps from GC, this parameter will be used to filter low-support haps out.)
		public int est_ind_pool; 
		public int level_I_and_II_alignment_cutoff; 

		DivideParameters(String propFilePath) throws IOException { 
			InputStream is = null; 
			try {
				Properties prop = new Properties();
				is = new FileInputStream(propFilePath);
				prop.load(is);
				this.gap_inpool_cutoff = Double.parseDouble(prop.getProperty("In-pool_Gap_Support_Min"));
				this.gap_all_pool_cutoff = Double.parseDouble(prop.getProperty("All-pool_Gap_Support_Min"));
				this.min_level_I_region_size = Integer.parseInt(prop.getProperty("Level_1_Region_Size_Min"));
				this.max_level_I_region_size = Integer.parseInt(prop.getProperty("Level_1_Region_Size_Max"));
				this.min_level_II_region_size = Integer.parseInt(prop.getProperty("Level_2_Region_Size_Min"));
				this.max_level_II_region_size = Integer.parseInt(prop.getProperty("Level_2_Region_Size_Max"));
				this.max_num_rounds_forming_initial_haps = Integer.parseInt(prop.getProperty("Rounds_Haps4Regions_Max"));
				this.max_num_haps = Integer.parseInt(prop.getProperty("Haps_Per_Region_Max"));
				this.est_ind_pool = Integer.parseInt(prop.getProperty("Est_Ind_PerPool"));
				this.level_I_and_II_alignment_cutoff = Integer.parseInt(prop.getProperty("Region_Alignment_Cutoff"));
			} catch (Exception e) {
				System.out.println("Exception: " + e);
			} finally {
				is.close();
			}
		} 
	}

	public static class AemParameters extends Parameters {
		public int max_iteration; 
		public int est_ind_pool; 
		public double epsilon; 
		public double rare_cutoff; 
		
		AemParameters(String propFilePath) throws IOException { 
			InputStream is = null; 
			try {
				Properties prop = new Properties();
				is = new FileInputStream(propFilePath);
				prop.load(is);
				this.max_iteration = Integer.parseInt(prop.getProperty("Iterations_AEM_Max"));
				this.est_ind_pool = Integer.parseInt(prop.getProperty("Est_Ind_PerPool"));
				this.epsilon = Double.parseDouble(prop.getProperty("Difference_Cutoff"));
				this.rare_cutoff = Double.parseDouble(prop.getProperty("Hap_All-Pool_Freq_Min"));
			} catch (Exception e) {
				System.out.println("Exception: " + e);
			} finally {
				is.close();
			}
		}
	}

	public static class McmcParameters extends Parameters {
		
		public int num_round;
		public double freqs_sum; // How much frequency to assign the randomly generated haplotypes to fill the matrix.
		public int burn_in;
		public int max_iteration; 
		public double alpha; 
		public double beta_a;
		public double beta_c;
		public double gamma;
		public double c_old;
		public double c_new;
		public double p_add;
		public int coalescing_mismatch;
		public double rare_cutoff;

		McmcParameters(String propFilePath) throws IOException { 
			InputStream is = null; 
			try {
				Properties prop = new Properties();
				is = new FileInputStream(propFilePath);
				prop.load(is);
				this.num_round = Integer.parseInt(prop.getProperty("Rounds_rjMCMC_Max"));
				this.freqs_sum = Double.parseDouble(prop.getProperty("Freq2Assign_RandHaps"));
				this.burn_in = Integer.parseInt(prop.getProperty("Burnin_rjMCMC"));
				this.max_iteration = Integer.parseInt(prop.getProperty("Iterations_rjMCMC"));
				this.alpha = Double.parseDouble(prop.getProperty("Shape_Param_Alpha"));
				this.beta_a = Double.parseDouble(prop.getProperty("Shape_Param_Beta_A"));
				this.beta_c = Double.parseDouble(prop.getProperty("Shape_Param_Beta_C"));
				this.gamma = Double.parseDouble(prop.getProperty("Shape_Param_Gamma"));
				this.c_old = Double.parseDouble(prop.getProperty("Shape_Param_C_Old"));
				this.c_new = Double.parseDouble(prop.getProperty("Shape_Param_C_New"));
				this.p_add = Double.parseDouble(prop.getProperty("Probability_AddHap"));
				this.coalescing_mismatch = Integer.parseInt(prop.getProperty("Coalesce_Mismatch_Max"));
				this.rare_cutoff = Double.parseDouble(prop.getProperty("Hap_In-Pool_Freq_Min"));
			} catch (Exception e) {
				System.out.println("Exception: " + e);
			} finally {
				is.close();
			}
		}
	}

}
