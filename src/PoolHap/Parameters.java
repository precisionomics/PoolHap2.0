package PoolHap;

import java.io.FileInputStream;
import java.io.IOException;
import java.io.InputStream;
import java.util.HashSet;
import java.util.Properties;

// chen test change


/*
 *  A typical property file looks like:
 *
 *  # FullSimulator Parameters
 *  ##########
 *  # General: Commands and file locations
 *  Input_Dir = /gpfs/home/lmak/PHX_Perfect_Data/input/
 *  Intermediate_Dir = /gpfs/home/lmak/PHX_Perfect_Data/intermediate/
 *  Gold-Standard_Dir = /gpfs/home/lmak/PHX_Perfect_Data/gold_standard/
 *  ms = /gpfs/home/lmak/programs/msdir/ms
 *  DWGSIM = /gpfs/home/lmak/programs/DWGSIM-master/dwgsim
 *  ##########
 *  # ms: Generates populations of genotypes under a variety of neutral models to investige their
 *  # stati
 *  Haps_Per_Pool = 10
 *  Num_Pools = 20
 *  Est_Ind_Per_Pool = 1000000
 *  Mutaton_Rate_Per_Base = 0.00001
 *  Segregating_Sites = 80
 *  Ref_Seq_Len = 9718
 *  ##########
 *  # dwgsim: Simulating a variety of next- and third-generation sequencing reads from input
 *  # genetic se
 *  Reference_Seq = HIV_HXB2.fa
 *  Error_Rate_Per_Base = 0
 *  Coverage = 100
 *  Read_Len = 100
 *  Outer_Dist = 400
 *  ##########
 */


public class Parameters {
    // TODO: [ReconEP]:: add any new parameters as needed.

    public static class GenParameters extends Parameters {
        // TODO: [ReconEP]:: refactor variable names to be more sensible (e.g. camelCase for Java).

        // Added by Quan Long 2019-06-27.
        // gc = Graph Coloring; aem = Divide and Conquer and then AEM;
        public static String[] supported_functions_array = {"format", "gc", "aem", "lasso"};
        public String function; // which module to run: format, gc, aem or lasso.

        /*
         *  General parameter set.
         */

        public String project_name;
        // input files: SAM files and VCF files. Added by Quan Long 2019-07-01
        public String input_dir;
        public String inter_dir; // intermediate directory, including the following files:

        // TODO: LEFTOVER ML 20190702
        // public String gs_dir; // gold standard directory

        public String out_dir; // output directory
        public int fragments; // TODO: [Question]:: what is this?
        public double final_cutoff;
        public double lambda;

        // TODO: [Question]:: what are one_vector, Hap_VC, and Hap_11?
        public double[] lasso_weights;

        public double min_r2;
        public double lasso_penalty_step;


        /**
         *  General parameters object constructor.
         *
         *  @param propFilePath (required) general parameters properties file path string.
         *  @throws IOException On input error.
         */
        // TODO: [Question]:: why isn't this public?
        public GenParameters(String propFilePath) throws IOException {
            // TODO: [Question]:: why does this have to be initialized to a value? There are other
            // instances of the variable name being set but initialized.
            // Would it affect things if we change this to:
            // InputStream is;
            InputStream is = null; // initialize input stream to null
            HashSet<String> supported_functions = new HashSet<String>();
            for (int k = 0; k < supported_functions_array.length; k++) {
                supported_functions.add(supported_functions_array[k]);
            }

            try {
                /*
                 *  Load parameters from properties file.
                 */
                Properties prop = new Properties(); // properties object from properties file
                is = new FileInputStream(propFilePath);
                prop.load(is);


                /*
                 *  Extract parameters to general parameter object variables from properties object.
                 */
                this.function = prop.getProperty("Function");
                if (!supported_functions.contains(this.function)) {
                    System.out.println("Function "+this.function+" is not supported. A typo?");
                    System.exit(0);
                }

                this.project_name= prop.getProperty("Proj_Name");
                this.input_dir= prop.getProperty("Input_Dir");
                this.inter_dir = prop.getProperty("Intermediate_Dir");

                // TODO: LEFTOVER Removed by Quan Long, 2019-07-01
                // this.gs_dir = prop.getProperty("Gold-Standard_Dir");

                this.out_dir = prop.getProperty("Output_Dir");
                this.fragments = Integer.parseInt(prop.getProperty("Fragments"));
                this.final_cutoff = Double.parseDouble(
                    prop.getProperty("FullLength_Local_Freq_Min"));

                this.lambda = Double.parseDouble(prop.getProperty("Lambda_Penalty"));
                this.lasso_weights = new double[] {
                    Double.parseDouble(prop.getProperty("One_Vector_Weight")),
                    Double.parseDouble(prop.getProperty("Hap_VC_Weight")),
                    Double.parseDouble(prop.getProperty("Hap_11_Weight"))};

                this.min_r2 = Double.parseDouble(prop.getProperty("Minimum_R2_Fit"));
                this.lasso_penalty_step = Double.parseDouble(prop.getProperty("Penalty_Step_Size"));

            } catch (Exception e) {
                e.printStackTrace();
            } finally {
                is.close(); // close input stream
            }

        }

    }


    public static class DivideParameters extends Parameters {
        // TODO: [ReconEP] refactor variable names to be more sensible (e.g. camelCase for Java).

        /*
         *  Divide and conquer parameter set.
         */
        public String project_name;
        public double gap_inpool_cutoff; // a ratio
        public double gap_all_pool_cutoff; // a ratio
        public double gap_support_step; // a ratio
        public int min_level_I_region_size;
        public int max_level_I_region_size;
        public int min_level_I_last_size;
        public int min_level_II_region_size;
        public int max_level_II_region_size;
        public int est_ind_pool;
        public double final_cutoff;
        public double lambda;
        public double[] lasso_weights;
        public double min_r2;
        public double lasso_penalty_step;
        public int hapset_size_max;
        public int hapset_size_min;
        public double hapset_size_rand;

        // TODO: [LEFTOVER]
        // // Because of trimming, we may not be able to form all full haplotypes.
        // // Cutoff to halt that.
        // public int max_num_rounds_forming_initial_haps;
        //
        // // The max number of haps that will be allowed in the AEM algorithm.
        // // (When there are too many haps from GC, this parameter will be used to filter
        // // low-support haps out.)
        // public int max_num_haps;
        //
        // public int est_ind_pool;
        // public int level_I_and_II_alignment_cutoff;

        /**
         *  Divide and conquer parameters object constructor.
         *
         *  @param propFilePath (required) DC parameters properties file path string.
         *  @throws IOException on input error.
         */
        public DivideParameters(String propFilePath) throws IOException {
            // TODO: [Question]:: same as above.
            InputStream is = null; // initialize input stream to null

            try {
                /*
                 *  Load parameters from properties file.
                 *
                 *  TODO: [ReconEP]:: extract the following into a static method (maybe in utils) as
                 *  this is the second time appearing.
                 */
                Properties prop = new Properties(); // properties object from properties file
                is = new FileInputStream(propFilePath);
                prop.load(is);


                /*
                 *  Extract parameters to divide and conquer parameter object variables from
                 *  properties object.
                 */
                this.project_name= prop.getProperty("Proj_Name");
                this.gap_inpool_cutoff = Double.parseDouble(
                    prop.getProperty("In-pool_Gap_Support_Min"));

                this.gap_all_pool_cutoff =  Double.parseDouble(
                    prop.getProperty("All-pool_Gap_Support_Min"));

                this.gap_support_step = Double.parseDouble(
                    prop.getProperty("Gap_Support_Step_Size"));

                this.min_level_I_region_size = Integer.parseInt(
                    prop.getProperty("Level_1_Region_Size_Min"));

                this.max_level_I_region_size = Integer.parseInt(
                    prop.getProperty("Level_1_Region_Size_Max"));

                this.min_level_I_last_size = Integer.parseInt(
                    prop.getProperty("Level_1_Last_Region_Min"));

                this.min_level_II_region_size = Integer.parseInt(
                    prop.getProperty("Level_2_Region_Size_Min"));

                this.max_level_II_region_size = Integer.parseInt(
                    prop.getProperty("Level_2_Region_Size_Max"));

                this.est_ind_pool = Integer.parseInt(prop.getProperty("Est_Ind_PerPool"));
                this.final_cutoff = Double.parseDouble(
                    prop.getProperty("Regional_Global_Freq_Min"));

                this.lambda = Double.parseDouble(prop.getProperty("Lambda_Penalty"));
                this.lasso_weights = new double[] {
                    Double.parseDouble(prop.getProperty("One_Vector_Weight")),
                    Double.parseDouble(prop.getProperty("Hap_VC_Weight")),
                    Double.parseDouble(prop.getProperty("Hap_11_Weight"))};

                this.min_r2 = Double.parseDouble(prop.getProperty("Minimum_R2_Fit"));
                this.lasso_penalty_step = Double.parseDouble(prop.getProperty("Penalty_Step_Size"));
                this.hapset_size_max = Integer.parseInt(
                    prop.getProperty("Regional_HapSetSize_Max"));

                this.hapset_size_min = Integer.parseInt(
                    prop.getProperty("Regional_HapSetSize_Min"));

                this.hapset_size_rand = Double.parseDouble(prop.getProperty("DC_HapSetSize_Rand"));

            } catch (Exception e) {
                e.printStackTrace();
            } finally {
                is.close();
            }

        }

    }


    public static class AemParameters extends Parameters {
        // TODO: [ReconEP] refactor variable names to be more sensible (e.g. camelCase for Java).

        /*
         *  Approximate expectation-maximization parameter set.
         */
        public int max_iteration;
        public int est_ind_pool;
        public double epsilon;
        public double rare_cutoff;
        public double final_cutoff;
        public int hapset_size_max;
        public int adhoc_freq_cutoff;
        public int hapset_size_min;
        public double hapset_size_rand;


        /**
         *  Approximate expectation-maximization parameters object constructor.
         *
         *  @param propFilePath (required) AEM parameters properties file path string.
         *  @throws IOException on input error.
         */
        AemParameters(String propFilePath) throws IOException {
             // TODO: [Question]:: same as above.
            InputStream is = null; // initialize input stream to null

            try {
                /*
                 *  Load parameters from properties file.
                 *
                 *  TODO: [ReconEP]:: extract the following into a static method (maybe in utils) as
                 *  this is the third time appearing.
                 */
                Properties prop = new Properties();
                is = new FileInputStream(propFilePath);
                prop.load(is);

                /*
                 *  Extract parameters to AEM parameter object variables from properties object.
                 */
                this.max_iteration = Integer.parseInt(prop.getProperty("Iterations_AEM_Max"));
                this.est_ind_pool = Integer.parseInt(prop.getProperty("Est_Ind_PerPool"));
                this.epsilon = Double.parseDouble(prop.getProperty("Difference_Cutoff"));
                this.rare_cutoff = Double.parseDouble(prop.getProperty("Running_Freq_Min"));
                this.final_cutoff = Double.parseDouble(
                    prop.getProperty("Regional_Global_Freq_Min"));

                this.hapset_size_max = Integer.parseInt(
                    prop.getProperty("Regional_HapSetSize_Max"));

                this.adhoc_freq_cutoff = Integer.parseInt(prop.getProperty("Adhoc_Freq_Cutoff"));
                this.hapset_size_min = Integer.parseInt(
                    prop.getProperty("Regional_HapSetSize_Min"));

                this.hapset_size_rand = Integer.parseInt(prop.getProperty("AEM_HapSetSize_Rand"));

            } catch (Exception e) {
                // TODO: [Question]:: same as above.
                // TODO: [ReconEP]:: the try-catch statement colour is different from above on my
                // machine, double check to make sure nothing is wrong on a different machine.
                System.out.println("Exception: " + e);

            } finally {
                is.close();
            }

        }

    }

}
