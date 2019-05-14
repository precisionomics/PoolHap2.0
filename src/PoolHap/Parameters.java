package PoolHap;

import java.io.FileInputStream;
import java.io.IOException;
import java.io.InputStream;
import java.util.Properties;

public class Parameters {
	public static class GenParameters extends Parameters {
		public String inter_dir; 
		public String gs_dir; 
		public String out_dir;
		public int fragments; 
		public double final_cutoff; 
		public double lambda;
		public double[] lasso_weights; 
		
		GenParameters(String propFilePath) throws IOException { 
			InputStream is = null; 
			try {
				Properties prop = new Properties();
				is = new FileInputStream(propFilePath);
				prop.load(is);
				this.inter_dir = prop.getProperty("Intermediate_Dir");
				this.gs_dir = prop.getProperty("Gold-Standard_Dir");
				this.out_dir = prop.getProperty("Output_Dir");
				this.fragments = Integer.parseInt(prop.getProperty("Fragments"));
				this.final_cutoff = Double.parseDouble(prop.getProperty("FullLength_Local_Freq_Min"));
				this.lambda = Double.parseDouble(prop.getProperty("LASSO_Lambda_Penalty"));
				this.lasso_weights = new double[] {Double.parseDouble(prop.getProperty("One_Vector_Weight")), 
						Double.parseDouble(prop.getProperty("Hap_VC_Weight")), Double.parseDouble(prop.getProperty("Hap_11_Weight"))};
			} catch (Exception e) {
				System.out.println("Exception: " + e);
			} finally {
				is.close();
			}
		}
	}

	public static class DivideParameters extends Parameters {
		public double gap_inpool_cutoff;  // a ratio
		public double gap_all_pool_cutoff;// a ratio
		public int min_level_I_region_size;
		public int max_level_I_region_size;
		public int min_level_II_region_size;
		public int max_level_II_region_size;
		public int est_ind_pool; 
		public double final_cutoff; 
		public double lambda;
		public double[] lasso_weights; 

		// public int max_num_rounds_forming_initial_haps; // because of trimming, we may not be able to form all full haplotypes. Cutoff to halt that.
		// public int max_num_haps;	// the max number of haps that will be allowed in the AEM algorithm. 
									//	(When there are too many haps from GC, this parameter will be used to filter low-support haps out.)
		// public int est_ind_pool; 
		// public int level_I_and_II_alignment_cutoff; 

		public DivideParameters(String propFilePath) throws IOException { 
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
				this.est_ind_pool = Integer.parseInt(prop.getProperty("Est_Ind_PerPool"));
				this.final_cutoff = Double.parseDouble(prop.getProperty("Regional_Global_Freq_Min"));
				this.lambda = Double.parseDouble(prop.getProperty("DC_Lambda_Penalty"));
				this.lasso_weights = new double[] {Double.parseDouble(prop.getProperty("One_Vector_Weight")), 
						Double.parseDouble(prop.getProperty("Hap_VC_Weight")), Double.parseDouble(prop.getProperty("Hap_11_Weight"))};
				// this.max_num_rounds_forming_initial_haps = Integer.parseInt(prop.getProperty("Rounds_Haps4Regions_Max"));
				// this.max_num_haps = Integer.parseInt(prop.getProperty("Haps_Per_Region_Max"));
				// this.est_ind_pool = Integer.parseInt(prop.getProperty("Est_Ind_PerPool"));
				// this.level_I_and_II_alignment_cutoff = Integer.parseInt(prop.getProperty("Region_Alignment_Cutoff"));
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
		public double final_cutoff; 
		
		AemParameters(String propFilePath) throws IOException { 
			InputStream is = null; 
			try {
				Properties prop = new Properties();
				is = new FileInputStream(propFilePath);
				prop.load(is);
				this.max_iteration = Integer.parseInt(prop.getProperty("Iterations_AEM_Max"));
				this.est_ind_pool = Integer.parseInt(prop.getProperty("Est_Ind_PerPool"));
				this.epsilon = Double.parseDouble(prop.getProperty("Difference_Cutoff"));
				this.rare_cutoff = Double.parseDouble(prop.getProperty("Running_Freq_Min"));
				this.final_cutoff = Double.parseDouble(prop.getProperty("Regional_Global_Freq_Min"));
			} catch (Exception e) {
				System.out.println("Exception: " + e);
			} finally {
				is.close();
			}
		}
	}

}
