package PoolHap;

import java.io.PrintWriter;
import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.InputStreamReader;

import PoolHap.Parameters;
import java.time.format.DateTimeFormatter;
import java.util.Arrays;
import java.util.HashMap;
import java.util.HashSet;

import MiscFunctions.BAMFormatterGATK;

import java.time.LocalDateTime;

public class Entrance {
    // TODO: [ReconEP]:: refactor variable names to be more sensible (e.g. camelCase for Java).
    // TODO: [Question]:: confusion over the HapLASSO portion.
	
    /**
     * Reconstructed by @author quanlong 2019-07. (1) split into three modules (2) The users need to
     * put input files into directories named as sam, vef, and vcf, instead of specified prefix plus
     * "_p" and number. Variables [name2index] and [names_array] will take over the numbering and
     * indexing. The possible suffixes are specified in the variable [suffixes]
     */

    public static String[] suffixes = {
        "sam", "vcf", "vef", "gcf"
    };
    
    public static int num_pools = -1;
    public static HashMap<String, Integer> name2index = new HashMap<String, Integer>();
    public static String[] names_array;
    
    /**
     * @author Quan Long 2019-07
     * 
     * Get all file paths from a folder. All the files should have the same suffix and the number of
     * files must be the same as the number of pools. All files should have the same names
     * 
     * If @param initialize_names ==true, the three static fields will be initialized. Otherwise,
     * the established values in the three fields will be compared with the outcome for a
     * double-check.
     */

    public static String[] get_filepaths(
        String name_file,
        String folder,
        String suffix,
        boolean write_name_file) {
        String[] filepaths = null;
        try {
            File[] files = new File(folder).listFiles();
            if (write_name_file) { // write the name_file and initiate structures.
                filepaths = new String[files.length];
                for (int k = 0; k < files.length; k++) {
                    filepaths[k] = files[k].toString();
                }
                Arrays.sort(filepaths); // sort by names.
                Entrance.num_pools = files.length;
                Entrance.names_array = new String[Entrance.num_pools];
                for (int index = 0; index < Entrance.num_pools; index++) {
                    String separat_dirs[] = filepaths[index].split("/");
                    String[] separated = separat_dirs[separat_dirs.length - 1].split("\\.");
                    if (!separated[separated.length - 1].equals(suffix)) {
                        System.out.println(
                            "ERROR: The suffix of " + filepaths[index] + " is not " + suffix);
                        System.exit(0);
                    } else {
                        String name = separated[0]; // get the string before the suffix, i.e., name
                        if (separated.length >= 3) { // except for [0] and the suffix, there are
                                                     // elements
                            for (int i = 1; i < separated.length - 1; i++) {
                                name = name + "." + separated[i];
                            }
                        }
                        Entrance.name2index.put(name, index);
                        Entrance.names_array[index] = name;
                    }
                }
                // write the names to the document for future functions to load
                BufferedWriter bw = new BufferedWriter(new FileWriter(name_file));
                for (int p = 0; p < Entrance.num_pools; p++) {
                    bw.write(names_array[p] + "\n");
                }
                bw.close();
            } else { // load names from the "name_file" and check whether they are correct.
                BufferedReader br = new BufferedReader(new FileReader(name_file));
                String line = br.readLine();
                int p_index = 0;
                while (line != null) {
                    Entrance.name2index.put(line, p_index++);
                    line = br.readLine();
                }
                br.close();
                Entrance.num_pools = p_index;
                Entrance.names_array = new String[Entrance.num_pools];
                for (String name : Entrance.name2index.keySet()) {
                    Entrance.names_array[Entrance.name2index.get(name)] = name;
                }
                // check the correctness in the target folder.
                if (Entrance.num_pools != files.length) { // first, number must be correct.
                    System.out.println("ERROR: No. of files in " + folder + " is " + files.length
                        + ", NOT consistent to the No. of lines in the " + name_file + ": "
                        + Entrance.num_pools);
                    System.exit(0);
                }
                filepaths = new String[files.length];
                // observed_files are the files in the folder, may or may not align the names
                // and order specified by name_file. So need to check and return ordered paths
                String[] observed_files = new String[files.length];
                for (int k = 0; k < files.length; k++) {
                    observed_files[k] = files[k].toString();
                }
                for (int index = 0; index < Entrance.num_pools; index++) {
                    String[] separated = observed_files[index]
                        .split("/")[observed_files[index].split("/").length - 1].split("\\.");
                    if (!separated[separated.length - 1].equals(suffix)) {
                        System.out.println(
                            "ERROR: The suffix of " + observed_files[index] + " is not " + suffix);
                        System.exit(0);
                    } else {
                        String name = separated[0]; // get the string before the suffix, i.e., name
                        if (separated.length >= 3) { // except for [0] and the suffix, there are
                                                     // elements
                            for (int i = 1; i < separated.length - 1; i++) {
                                name = name + "." + separated[i];
                            }
                        }
                        if (!Entrance.name2index.containsKey(name)) {
                            // System.out.println(name);
                            System.out.println(
                                "ERROR: The name of " + observed_files[index] + " is not in "
                                    + "existing name2index table specified by previsous file sets");
                            System.exit(0);
                        } else {
                            int prespecified_index = name2index.get(name);
                            filepaths[prespecified_index] = observed_files[index];
                        }
                    }
                }
            }
        } catch (Exception e) {
            e.printStackTrace();
        }
        return filepaths;
    }

    
    public static void main(String[] args) throws Exception {
        /*
         * Input arguments.
         */
    	String[] supported_functions_array = {"non_perfect", "format", "gc", "aem", "lasso", 
    			"split" ,"clustering", "evaluate", "analysis", "l0l1"};
    	HashSet<String> supported_functions = new HashSet<String>();
        for (int k = 0; k < supported_functions_array.length; k++) {
             supported_functions.add(supported_functions_array[k]);
        }
    	String function  = args[0]; // parameter file path
    	
    	if (!supported_functions.contains(function)) {
            System.out.println("Function "+ function+" is not supported. A typo?");
            System.exit(0);
        }
    	
        String parameter_file = args[1]; // parameter file path
        
        // int num_pools = Integer.parseInt(args[2]); // number of pools present
        // removed by Quan Long 2019-07-03. This will be identified by counting how many files in a
        // folder
        
        /*
         * Program initialization.
         */
        // New general parameter set object from parameter file.
        // TODO: [Review]:: GenParameters(args[0]) changed to GenParameters(parameter_file)
        
        Parameters gp = new Parameters(parameter_file);

        // Gold standard variant positions file path (intra pool).
        String gs_var_pos = gp.inter_dir + gp.project_name + "_vars.intra_freq.txt"; // TODO: change
                                                                             // name
        String name_file = gp.inter_dir + gp.project_name + "_sample_names.txt";
        
        // Print start time.
        DateTimeFormatter dtf = DateTimeFormatter.ofPattern("yyyy/MM/dd HH:mm:ss");
        System.out
            .println("PoolHapX " + gp.function + " Initiated: " + dtf.format(LocalDateTime.now()));

        /*
         * Local haplotype configuration.
         */
        Evaluate  eva = new Evaluate ();
        // By design, HapConfig contains haps of all pools.
        if (function.equals("format")) {
        	
        	
        	 
        	 
            // SAM files in input_dir/sam initiates pool-IDs and orders, and write to the
            // "name_file".
            // The function below generates vef files and inpool-var-freq file.
            String[] sam_files =
                Entrance.get_filepaths(name_file, gp.input_dir + "/sam/", "sam", true);
            // genrate VEF files (based on SAM files) and var-freq file (based on joint VCF file)
//            System.out.println(gp.sequencing_technology);
            BAMFormatterGATK.generate_vef_varfreq(gp.input_dir,
                gp.inter_dir,
                gp.project_name,
                sam_files,
                Entrance.name2index,
                gp.sequencing_technology);
            
            File file = new File(gp.gold_dir+"/"+gp.project_name+"_vars.intra_freq.txt");
			
			if (file.exists()) {
				BAMFormatterGATK.rewrite_inter_vars(gp.gold_dir+"/"+gp.project_name+"_vars.intra_freq.txt"
	             		, gs_var_pos);
			}
            
            
        } else if (function.equals("gc")) {
            // PrintWriter in_list = new PrintWriter( // writer object for graph colored pools
            // new FileWriter(new File(gp.inter_dir + prefix + "_p.in.list")));
            // the folder "vef/" should exist under gp.inter_dir.
            String[] vef_files =
                Entrance.get_filepaths(name_file, gp.inter_dir + "vef", "vef", false);
            // make a new dir "in/" to store GC-outcome files.
            new File(gp.inter_dir + "/gcf/").mkdir();
            // Apply graph coloring for each pool. Note, initial "global" haplotypes are roughly
            // linked by GC with gaps between VARIANTS (not regions).
            for (int p = 0; p < num_pools; p++) {
                // TODO: [Question]:: this "pool_in" variable not used anywhere.
                // Previously added to initial_local_haps but now commented out. ML 20190702
            	try {
	                GraphColoring pool_in = new GraphColoring(vef_files[p],
	                    gs_var_pos,
	                    gp.inter_dir + "gcf/" + Entrance.names_array[p] + ".gcf" ,
	                    gp.num_pos_window, gp.num_gap_window);
            	}  catch (Exception e) {
            		GraphColoring pool_in = new GraphColoring(vef_files[p],
    	                    gs_var_pos, gp.inter_dir + "gcf/" + Entrance.names_array[p] + ".gcf"  );
            	}
                
                System.out.println("Graph colouring for pool " + p + ":" + Entrance.names_array[p]
                    + " is finished.");
            }
            System.out
                .println("\nGraph-Coloring Finished: " + dtf.format(LocalDateTime.now()) + "\n");
        } else if (function.equals("aem_previous")) { // Note: aem includes DC and AEM
            // Apply divide and conquer across all pools.
        	

//        	System.exit(0);
            String dc_out_file = gp.inter_dir + gp.project_name + "_dc_plan.txt"; // dc output
                                                                                  // filepath string
            String[] vef_files = 
                Entrance.get_filepaths(name_file, gp.inter_dir + "vef", "vef", false);
            String[] gcf_files = 
                Entrance.get_filepaths(name_file, gp.inter_dir + "gcf", "gcf", false);
            DivideConquer dc_maker = new DivideConquer(gs_var_pos,
                gcf_files,
                // gp.inter_dir + prefix + "_p.in.list",
                parameter_file,
                dc_out_file);
            
            
            new File(gp.inter_dir + "/aem/").mkdir();            
         
            HapConfig[] level_I_config = dc_maker.regional_AEM(
                Entrance.names_array,
                vef_files,
                dc_maker.regions_level_I,
                parameter_file,
                gp.inter_dir + "/aem/" + gp.project_name,
                gp.inter_dir + "/aem_fail_regional_lasso/",
                1);
            
            System.out.println("Level 1 AEM Finished: " + dtf.format(LocalDateTime.now()) + "\n");

            // Resolve level 2 regions with the same strategy as above.
            HapConfig[] level_II_config = dc_maker.regional_AEM(
                Entrance.names_array,
                vef_files,
                dc_maker.regions_level_II,
                parameter_file,
                gp.inter_dir + "/aem/" + gp.project_name,
                gp.inter_dir + "/aem_fail_regional_lasso/",
                2);

            System.out.println("Level 2 AEM Finished: " + dtf.format(LocalDateTime.now()) + "\n");
            
            eva.AemEvaluate(gp.project_name, dc_out_file, 
            		"/home/chencao/Desktop/sim001/gold_standard/sim001_haps.txt", 
        			"/home/chencao/Desktop/sim001/intermediate/aem/");
            
//            System.exit(0);

            // Link regions by applying graph coloring across the level 1 and level 2 regional
            // haplotype configurations.
            GraphColoring region_linker =
                new GraphColoring(level_I_config, level_II_config, gs_var_pos, gp.virtual_cov_link_gc);

            // Write final global haplotype configurations (inter pool) to output.
            HapConfig final_global_haps; // final global haplotype configuration object
            final_global_haps = region_linker.hapOut(Entrance.names_array); // HapConfig from
                                                                            // GC-linked regions
            final_global_haps.recode_HapIDs_to_base16();
            final_global_haps.write_global_file_string( // write to output
                gp.out_dir + gp.project_name + "_gc.inter_freq_haps.txt");

            System.out.println("\nPost-AEM GC Finished: " + dtf.format(LocalDateTime.now()) + "\n");
            eva.GcAemEvaluate("/home/chencao/Desktop/sim001/gold_standard/sim001_haps.txt", 
        			"/home/chencao/Desktop/sim001/output/sim001_gc.inter_freq_haps.txt");
        	
        } else if (function.equals("lasso")) {
            /*
             * Final intra pool haplotype configuration via LASSO selection. TODO: [Question]:: are
             * these global intra pool haplotypes?
             */
        	
        	
        	
	            new File(gp.inter_dir + "/lasso/").mkdir();
	            Entrance.get_filepaths(name_file, gp.inter_dir + "/gcf/", "gcf", false);
	            
	            HapConfig final_global_haps =
	                new HapConfig(gp.out_dir + gp.project_name + "_hc.inter_freq_haps.txt", null);
	            
	
	            SiteInPoolFreqAnno siteInPoolFreqAnno = new SiteInPoolFreqAnno(gs_var_pos);
	            
	            final_global_haps.inpool_site_freqs = siteInPoolFreqAnno.inpool_freqs;
	            final_global_haps.num_pools = Entrance.num_pools;
	            final_global_haps.pool_IDs = siteInPoolFreqAnno.pool_IDs;
	            final_global_haps.in_pool_haps_freq =
	                new double[final_global_haps.num_global_hap][final_global_haps.num_pools];
	
	            HapConfig[] final_inpool_haps = new HapConfig[num_pools]; // final intra pool HapConfigs
	            
	
        	
            for (int pool_index = 0; pool_index < num_pools; pool_index++) {
                HapLASSO inpool_lasso = new HapLASSO(pool_index,
                    Entrance.names_array[pool_index],
                    gp.lasso_global_lambda,
                    final_global_haps,
                    gp.lasso_full_hap_freq_cutoff,
                    gp.lasso_global_memory,
                    gp.inter_dir + "/lasso/" + Entrance.names_array[pool_index],
                	gp.lasso_coverage_weight,
                	gp.lasso_distance_max_weight	);

                inpool_lasso.estimate_frequencies_lasso(gp.inter_dir + "/vef/"
                    + Entrance.names_array[pool_index] + ".vef", null, gp.lasso_weights);

                // Learn optimal LASSO lambda value.
                //double new_penalty = gp.global_lambda;
                //while (inpool_lasso.r2 == gp.min_r2) {something deleted}
                
//                String[] single_pool_IDs=new String[] {names_array[pool_index]};
                //     add final intra pool
//                final_inpool_haps[pool_index] = inpool_lasso.hapOut(single_pool_IDs); 
                // HapConfig to array
//                System.out.println(final_inpool_haps[pool_index].num_global_hap
//                    + " haplotypes were generated for pool " + pool_index + ":"
//                    + Entrance.names_array[pool_index] + ".");
            }

            System.out.println("\nGlobal LASSO Finished: " + dtf.format(LocalDateTime.now())+"\n");

            // Writing final outputs to file.
//            HapConfig final_reconstruction = new HapConfig(final_inpool_haps);
//            System.out.println("Haplotype reconstruction finished. There are "
//                + final_reconstruction.num_global_hap + " global haplotypes.");
//
//            final_reconstruction.write2files(gp.out_dir + gp.project_name + ".inter_freq_haps.txt",
//                gp.out_dir + gp.project_name + ".intra_freq_haps.txt",
//                "string");
            
        	
        	
        	
        }else if (function.equals("clustering")) {
            /*
             * Chen: cluster the potential haplotypes using  hierarchical clustering
             * 
             */

	           	HapConfig clustering_global_haps =
	                new HapConfig(gp.out_dir + gp.project_name + "_gc.inter_freq_haps.txt", null);
//	           	HierarchicalClustering hc= new HierarchicalClustering(clustering_global_haps.hap_IDs, 
//	           			clustering_global_haps.global_haps_string, clustering_global_haps.hapID2index, 
//	           			gp.out_dir + gp.project_name + "_hc.haps.txt", gp.hc_similarity_cutoff);
//	           	hc.gc_solver(gs_var_pos);
//	           	Entrance.get_filepaths(name_file, gp.inter_dir + "/vef/", "vef", false);
//	           	
//	           	HapConfig hc_global_haps; // final global haplotype configuration object
//	           	hc_global_haps = hc.hapOut(Entrance.names_array); // HapConfig from Hierarchical Clustering                                                           
//	           	hc_global_haps.recode_HapIDs_to_base16();
//	           	hc_global_haps.write_global_file_string( // write to output
//                    gp.out_dir + gp.project_name + "_hc.inter_freq_haps.txt");


	           	
	           	
	           	
        }else if (function.equals("aem"))   {
       
        	
            String dc_out_file = gp.inter_dir + gp.project_name + "_dc_plan.txt"; // dc output
                        
            // filepath string       	
            String[] vef_files = 
                    Entrance.get_filepaths(name_file, gp.inter_dir + "vef", "vef", false);
            String[] gcf_files = 
                    Entrance.get_filepaths(name_file, gp.inter_dir + "gcf", "gcf", false);  
            
            
            DivideConquer dc_maker = new DivideConquer(gs_var_pos,
                    gcf_files,
                    // gp.inter_dir + prefix + "_p.in.list",
                    parameter_file,
                    dc_out_file);
                        
            
            new File(gp.inter_dir + "/aem/").mkdir();   
            HashMap<Integer, String> index_var_prefix_dict = dc_maker.gs_map(gs_var_pos);
            
//            if (gp.max_level_I_region_size >= index_var_prefix_dict.size()) {
//            	gp.max_level_I_region_size = index_var_prefix_dict.size()-1;
//            	gp.min_level_I_region_size = index_var_prefix_dict.size()-1;
//            	gp.max_level_II_region_size = index_var_prefix_dict.size()-1;
//            	gp.min_level_II_region_size = index_var_prefix_dict.size()-1;
//            }
            
            String reg_dc_out_file = gp.inter_dir + 
            		gp.project_name + "_regression__dc_plan.txt"; // regression dc output
            
            HashMap<Integer, Integer> index_var_pos_dict = dc_maker.gs_map_pos(gs_var_pos);
            
            DivideConquer reg_dc_maker = new DivideConquer(dc_maker.regions_level_final, 
            		dc_maker.regions_level_final_link, reg_dc_out_file, gp.regression_maximum_regions, 
        			index_var_pos_dict);  
            
            System.out.println("DC Finished Time:\t" + dtf.format(LocalDateTime.now()) + "\n");
            
            
            for (int p = 0; p < num_pools; p++) {
                // TODO: [Question]:: this "pool_in" variable not used anywhere.
                // Previously added to initial_local_haps but now commented out. ML 20190702
                GraphColoring pool_in = new GraphColoring(vef_files[p],
                    gs_var_pos,
                    gp.inter_dir + "gcf/" + Entrance.names_array[p] + ".gcf",
                    dc_maker.regions_level_I );
                dc_maker.loci_link_freq[p]= pool_in.loci_link_freq;
                System.out.println("Linking information extraction for pool " + p + ":" + Entrance.names_array[p]
                    + " is finished.");
            }
            
            HapConfig[] level_I_config = null;
            HapConfig[] level_II_config = null;
            HapConfig[] level_III_config = null;
            HapConfig[] level_IV_config = null;
            HapConfig[] level_V_config = null;
            HapConfig[] level_VI_config = null;
            HapConfig[] level_VII_config = null;
            HapConfig[] level_VIII_config = null;
            
            if (dc_maker.final_level >=1 ) {
	            level_I_config = dc_maker.regional_AEM(
	                    Entrance.names_array,
	                    vef_files,
	                    dc_maker.regions_level_I,
	                    parameter_file,
	                    gp.inter_dir + "/aem/" + gp.project_name,
	                    gp.inter_dir + "/aem_fail_regional_lasso/",
	                    1, Entrance.name2index);
	            dc_maker.level_I_config= level_I_config;
	            	System.out.println("Level 1 AEM Finished: " + dtf.format(LocalDateTime.now()) + "\n");
	            	
	            level_II_config = dc_maker.regional_AEM(
	                    Entrance.names_array,
	                    vef_files,
	                    dc_maker.regions_level_II,
	                    parameter_file,
	                    gp.inter_dir + "/aem/" + gp.project_name,
	                    gp.inter_dir + "/aem_fail_regional_lasso/",
	                    2,  Entrance.name2index);
	            dc_maker.level_II_config= level_II_config;
	            System.out.println("Level 2 AEM Finished: " + dtf.format(LocalDateTime.now()) + "\n");
            
	          eva.AemEvaluate(gp.project_name, dc_out_file, 
	          gp.gold_dir+"/"+gp.project_name+"_haps.inter_freq_vars.txt", 
	       		gp.inter_dir+"/aem/");
           }
            
           if (dc_maker.final_level >=3 ) {
	            level_III_config = dc_maker.level_III_regional_AEM( 
	            		Entrance.names_array,
	                    vef_files,
	                    dc_maker.regions_level_III,
	                    parameter_file,
	                    gp.inter_dir + "/aem/" + gp.project_name,
	                    gp.inter_dir + "/aem_fail_regional_lasso/",
	                    3,  Entrance.name2index,
	                    gp.level_III_IV_region_mismatch_tolerance);
	            dc_maker.level_III_config= level_III_config;
	            System.out.println("Level 3 AEM Finished: " + dtf.format(LocalDateTime.now()) + "\n");
	          
	            level_IV_config = dc_maker.level_IV_regional_AEM( 
	            		Entrance.names_array,
	                    vef_files,
	                    dc_maker.regions_level_IV,
	                    parameter_file,
	                    gp.inter_dir + "/aem/" + gp.project_name,
	                    gp.inter_dir + "/aem_fail_regional_lasso/",
	                    4,  Entrance.name2index,
	                    gp.level_III_IV_region_mismatch_tolerance);
	            dc_maker.level_IV_config= level_IV_config;
	            System.out.println("Level 4 AEM Finished: " + dtf.format(LocalDateTime.now()) + "\n");
            
            
	            eva.Aem_III_IV_Evaluate(gp.project_name, dc_out_file, 
                    gp.gold_dir+"/"+gp.project_name+"_haps.inter_freq_vars.txt", 
                 		gp.inter_dir+"/aem/");
            }
            if (dc_maker.final_level >=5 ) {
            	level_V_config = dc_maker.level_V_regional_AEM( 
            		Entrance.names_array,
                    vef_files,
                    dc_maker.regions_level_V,
                    parameter_file,
                    gp.inter_dir + "/aem/" + gp.project_name,
                    gp.inter_dir + "/aem_fail_regional_lasso/",
                    5,  Entrance.name2index,
                    gp.level_V_VI_region_mismatch_tolerance);
            	dc_maker.level_V_config= level_V_config;
            	System.out.println("Level 5 AEM Finished: " + dtf.format(LocalDateTime.now()) + "\n");
            
            level_VI_config = dc_maker.level_VI_regional_AEM( 
            		Entrance.names_array,
                    vef_files,
                    dc_maker.regions_level_VI,
                    parameter_file,
                    gp.inter_dir + "/aem/" + gp.project_name,
                    gp.inter_dir + "/aem_fail_regional_lasso/",
                    6,  Entrance.name2index,
                    gp.level_V_VI_region_mismatch_tolerance);
	            dc_maker.level_VI_config= level_VI_config;
	            System.out.println("Level 6 AEM Finished: " + dtf.format(LocalDateTime.now()) + "\n");
            
	            eva.Aem_V_VI_Evaluate(gp.project_name, dc_out_file, 
	                    gp.gold_dir+"/"+gp.project_name+"_haps.inter_freq_vars.txt", 
	                 		gp.inter_dir+"/aem/");
            }
            if (dc_maker.final_level >=7 ) {
	            level_VII_config = dc_maker.level_VII_regional_AEM( 
	            		Entrance.names_array,
	                    vef_files,
	                    dc_maker.regions_level_VII,
	                    parameter_file,
	                    gp.inter_dir + "/aem/" + gp.project_name,
	                    gp.inter_dir + "/aem_fail_regional_lasso/",
	                    7,  Entrance.name2index,
	                    gp.level_VII_VIII_region_mismatch_tolerance);
	            dc_maker.level_VII_config= level_VII_config;
	            System.out.println("Level 7 AEM Finished: " + dtf.format(LocalDateTime.now()) + "\n");
	            
	            level_VIII_config = dc_maker.level_VIII_regional_AEM( 
	            		Entrance.names_array,
	                    vef_files,
	                    dc_maker.regions_level_VIII,
	                    parameter_file,
	                    gp.inter_dir + "/aem/" + gp.project_name,
	                    gp.inter_dir + "/aem_fail_regional_lasso/",
	                    8,  Entrance.name2index,
	                    gp.level_VII_VIII_region_mismatch_tolerance);
	            dc_maker.level_VIII_config= level_VIII_config;
	            System.out.println("Level 8 AEM Finished: " + dtf.format(LocalDateTime.now()) + "\n");
	            
	            eva.Aem_VII_VIII_Evaluate(gp.project_name, dc_out_file, 
	                    gp.gold_dir+"/"+gp.project_name+"_haps.inter_freq_vars.txt", 
	                 		gp.inter_dir+"/aem/");
            }
            
            dc_maker.level_final_config= null;
            dc_maker.level_final_link_config =null;
            
            
            if (dc_maker.final_level == 7  ) {
            	dc_maker.level_final_config= level_VII_config;
            	dc_maker.level_final_link_config= level_VIII_config;
            } else if (dc_maker.final_level == 5  ) {
            	dc_maker.level_final_config= level_V_config;
            	dc_maker.level_final_link_config= level_VI_config;
            } else if (dc_maker.final_level == 3  ) {
            	dc_maker.level_final_config= level_III_config;
            	dc_maker.level_final_link_config= level_IV_config;
            } else if (dc_maker.final_level == 1  ) {
            	dc_maker.level_final_config= level_I_config;
            	dc_maker.level_final_link_config= level_II_config;
            }
            
            System.out.println("AEM Finished Time:\t" + dtf.format(LocalDateTime.now()) + "\n");
            System.out
                 .println("\nAEM Finished and Start regional L0L1 regression: " + dtf.format(LocalDateTime.now()) + "\n");
            
            String reg_prefix=  gp.inter_dir + "/regression/";
//------------------------------------------------------------
            try {  
                String shpath="rm -rf  "+ reg_prefix;  
                Process ps = Runtime.getRuntime().exec(shpath, null);  
                ps.waitFor();  
                BufferedReader br = new BufferedReader(new InputStreamReader(ps.getInputStream()));  
                StringBuffer sb = new StringBuffer();  
                String line;  
                while ((line = br.readLine()) != null) {  
                    sb.append(line).append("\n");  
                }  
                String result = sb.toString();  
                br.close();
            }   catch (Exception e) {  
                e.printStackTrace();  
            }  
           
            
            new File(reg_prefix).mkdir();
            new File(reg_prefix+"/Rscript").mkdir();
            
            for (int pool_index = 0; pool_index < num_pools; pool_index++) {
            	new File(gp.inter_dir + "/regression/"+ Entrance.names_array[pool_index]).mkdir();
	            for (int r=0; r< reg_dc_maker.regession_level[0].length;r++) {
	            	
	            	GraphColoring region_linker =
	                        new GraphColoring( reg_dc_maker.regession_level[0][r][0],  
	                        		reg_dc_maker.regession_level[0][r][1],
	                        		dc_maker.level_final_config, dc_maker.level_final_link_config, 
	                        		gs_var_pos, dc_maker.regions_level_final,
	                        		dc_maker.regions_level_final_link, gp.bfs_mismatch_tolerance, 
	                        		r+1 ,index_var_prefix_dict,  reg_prefix+ Entrance.names_array[pool_index],
	                        		gp.regression_maximum_regions);
	            }
            }
            
            System.out
            .println("Have generated initialization regional haptypes from AEM level "+ dc_maker.final_level+
            		" regional haplotypes. " + dtf.format(LocalDateTime.now()) + "\n");
            
            for (int regression_level=0; regression_level< reg_dc_maker.regession_level.length; 
            		regression_level++) {
            	for (int r=0; r< reg_dc_maker.regession_level[regression_level].length;r++) {
            		for (int pool_index = 0; pool_index < num_pools; pool_index++) {
		            	System.out
		                .println("Generating level "+Integer.toString(regression_level+1) + 
		                		" region "+ Integer.toString(r+1)+ "  potential haplotypes for pool "
		                		+Entrance.names_array[pool_index] + ". "
		                		 + dtf.format(LocalDateTime.now()) );
		            	if ((reg_dc_maker.regession_level[regression_level][r][0]!=-1 ) && 
		            			(reg_dc_maker.regession_level[regression_level][r][1]!=-1 )) {
		            		if (regression_level!=0) {
				            	GraphColoring regression_region_linker =
				            	    new GraphColoring( reg_dc_maker.regession_level[regression_level][r][0],
				            	    	reg_dc_maker.regession_level[regression_level][r][1],
					                    reg_dc_maker.regession_level[regression_level-1],  
					                    dc_maker.level_final_link_config, 
					                    gs_var_pos,
					                    dc_maker.regions_level_final_link, 
					                    gp.regression_link_mismatch_tolerance, 
					                    r+1 ,regression_level+1, index_var_prefix_dict,  
					                    reg_prefix+ Entrance.names_array[pool_index], 
					                    gp.regression_maximum_regions, 
					                    gp.regression_maximum_selected_haplotypes);
		            		}
		            		
		            		HapConfig final_global_haps =
			    	                new HapConfig(gp.inter_dir + "/regression/" + Entrance.names_array[pool_index] +
			    	                		"/regression_level_"+   Integer.toString(regression_level+1)+
			    	                		"_region_"+ Integer.toString(r+1)+".potential.haps", null);
			            	SiteInPoolFreqAnno siteInPoolFreqAnno = new SiteInPoolFreqAnno(gs_var_pos, 
			            			reg_dc_maker.regession_level[regression_level][r][0],
			            			reg_dc_maker.regession_level[regression_level][r][1], pool_index);
			            	final_global_haps.inpool_site_freqs = siteInPoolFreqAnno.inpool_freqs;
			            	final_global_haps.num_pools = Entrance.num_pools;
			            	final_global_haps.pool_IDs = siteInPoolFreqAnno.pool_IDs;
			            	final_global_haps.in_pool_haps_freq =
			    	                new double[final_global_haps.num_global_hap][final_global_haps.num_pools];
			            	
			            	HapLASSO inpool_lasso = new HapLASSO(pool_index,
			                        Entrance.names_array[pool_index],
			                        gp.lasso_global_lambda,
			                        final_global_haps,
			                        gp.lasso_full_hap_freq_cutoff,
			                        gp.lasso_global_memory,
			                        gp.inter_dir + "/regression/" + Entrance.names_array[pool_index]+
			                        "/regression_level_"+Integer.toString(regression_level+1)+
			                        "_region_"+ Integer.toString(r+1) ,
			                    	gp.lasso_coverage_weight,
			                    	gp.lasso_distance_max_weight);
			            	
			            	
			            	inpool_lasso.estimate_frequencies_lasso(gp.inter_dir + "/vef/"
			                        + Entrance.names_array[pool_index] + ".vef", null, gp.lasso_weights ,
			                        gp.sequencing_technology, 
			                        reg_dc_maker.regession_level_pos_come_from.get(regression_level).get(r),
			                        reg_dc_maker.regession_level_num_come_from[regression_level][r], 
			                        index_var_pos_dict, 
			                        regression_level);
		            	}
		            }
	            }
            	
	            System.out
                .println("\nRunning L0L1 regression for all regions in level "+Integer.toString(regression_level+1) +
                		". "+ dtf.format(LocalDateTime.now()) + "\n");
	            
	            HapLASSO regression_rscript  = new HapLASSO(Entrance.names_array, 
	            		reg_dc_maker.regession_level[regression_level],
	            		gp.inter_dir + "/regression/Rscript/" , 
	            		gp.inter_dir + "/regression/" , 
	            		regression_level+1,
	            		gp.regression_hapset_size_max,
	            		gp.regression_hapset_size_min,
	            		gp.regression_gamma_min,
	            		gp.regression_gamma_max,
	            		gp.ngammas, 
	            		gp.lasso_weights[1], 
	            		gp.species);
	            
	            for (int r=0; r< reg_dc_maker.regession_level[regression_level].length;r++) {
	            	
	            	if ((reg_dc_maker.regession_level[regression_level][r][0]!=-1 ) && 
	            			(reg_dc_maker.regession_level[regression_level][r][1]!=-1 )) {
	            		regression_rscript.run_L0L1learn( Entrance.names_array,
	            				gp.rscript_path, gp.inter_dir + "/regression/Rscript/regression_level_"+
	            				Integer.toString(regression_level+1)+ "_region_"+Integer.toString(r+1), 
	            				 gp.inter_dir + "/regression/",
	            				regression_level+1, 
	            				r+1, 
	            				gp.num_threads
	            				);
	            		System.out
		                .println("Level "+Integer.toString(regression_level+1) +" region "+ 
		                		Integer.toString(r+1)+" L0L1 regression finished.\t"+ 
		                		dtf.format(LocalDateTime.now()) + "\n");
	            	}
	            }
            }

            int regression_last_level= reg_dc_maker.regession_level.length-1;
            
            for (int pool_index = 0; pool_index < num_pools; pool_index++) {
            	new File(gp.out_dir + Entrance.names_array[pool_index]).mkdir();
            	HapConfig final_haps = new HapConfig (gp.inter_dir + "/regression/" + Entrance.names_array[pool_index] +
                		"/regression_level_"+   Integer.toString(regression_last_level+1)+
                		"_region_1.regression_out", 
                		gp.out_dir  + Entrance.names_array[pool_index]+ "/final_freq_haps.txt" ,
                		index_var_prefix_dict);
            }
            System.out.println("L0L1 Regression Finished Time:\t" + 
            		dtf.format(LocalDateTime.now()) + "\n");
            
            System.out.println("PoolHapX Successfully Finished, Enjoy!\n");
            
            
        }    else if (function.equals("evaluate")) {
        	
//        	eva.GcAemEvaluate("/home/chencao/Desktop/sim001/gold_standard/sim001_haps.txt", 
//        			"/home/chencao/Desktop/sim001/intermediate/regression/sim001_p5/"
//        			+ "regression_level_3_region_1.potential.haps");
        	
        	String[] vef_files =
                    Entrance.get_filepaths(name_file, gp.inter_dir + "vef", "vef", false);
//        	System.out.println(Entrance.num_pools); 
        	
        	BufferedWriter bw_mcc = new BufferedWriter(new FileWriter(gp.out_dir+"/MCC.result"));
            
            
        	double mcc_total =0;
        	for (int pool_index = 0; pool_index < Entrance.num_pools; pool_index++) {
//            	new File(gp.out_dir + Entrance.names_array[pool_index]).mkdir();
        		eva.MCCEvaluate(gp.gold_dir +"/"+ gp.project_name+"_haps.inter_freq_vars.txt",
        				gp.gold_dir +"/"+ gp.project_name+"_haps.intra_freq.txt",
        			gp.out_dir+"/"+ Entrance.names_array[pool_index]+"/final_freq_haps.txt",
        			gp.out_dir+"/"+Entrance.names_array[pool_index]+"/MCC.txt",
        			gp.mcc_freq_cutoff, 
        			Entrance.names_array[pool_index]);
        		bw_mcc.write("MCC for " +Entrance.names_array[pool_index] + 
        				" is:\t" + eva.mcc_value +"\n");
        		mcc_total+= eva.mcc_value;
            }
        	
        	System.out.println("Average MCC for all "+ Entrance.num_pools
        		+ " pools:\t"+ mcc_total/ (double) Entrance.num_pools );
        	bw_mcc.write("Average MCC for all "+ Entrance.num_pools
            		+ " pools:\t"+ mcc_total/ (double) Entrance.num_pools+"\n");
        	bw_mcc.close();
        	
//        	System.exit(0);
//        	   eva.GcAemEvaluate("/home/chencao/Desktop/sim001/gold_standard/sim001_haps.txt", 
//           			"/home/chencao/Desktop/sim001/output/sim001_hc.inter_freq_haps.txt");
        	// Chen:	MCC
//        		eva.MCCEvaluate(gp.gold_dir +"/"+ gp.project_name+"_haps.inter_freq_haps.txt",
//        				gp.gold_dir +"/"+ gp.project_name+"_haps.intra_freq.txt",
//        			gp.out_dir+"/"+gp.project_name+".inter_freq_haps.txt",
//        			gp.out_dir+"/"+gp.project_name+".intra_freq_haps.txt",
//        			gp.out_dir+"/"+gp.project_name+".MCC.txt",
//        			gp.mcc_freq_cutoff);
        	// Chen:	JSD
        	
        	BufferedWriter bw_jsd = new BufferedWriter(new FileWriter(gp.out_dir+"/JSD.result"));
        	
        	double jsd_total =0;
        	for (int pool_index = 0; pool_index < Entrance.num_pools; pool_index++) {
	        	eva.JSDEvaluate(gp.gold_dir +"/"+ gp.project_name+"_haps.inter_freq_vars.txt",
	        				gp.gold_dir +"/"+ gp.project_name+"_haps.intra_freq.txt",
	        			gp.out_dir+"/"+ Entrance.names_array[pool_index]+"/final_freq_haps.txt",
	        			Entrance.names_array[pool_index],
	        			gp.out_dir+"/"+Entrance.names_array[pool_index]+"/JSD.txt");
	        	jsd_total+= eva.jsd_value;
	        	bw_jsd.write("JSD for " +Entrance.names_array[pool_index] + 
        				" is:\t" + eva.jsd_value +"\n");
        	}
        	
        	System.out.println("Average JSD for all "+ Entrance.num_pools
        			+ " pools:\t"+ jsd_total/ (double) Entrance.num_pools );
        	
        	bw_jsd.write("Average JSD for all "+ Entrance.num_pools
        			+ " pools:\t"+ jsd_total/ (double) Entrance.num_pools+"\n" );;
        	
        	bw_jsd.close();
        	System.out
                .println("\nEvaluation Finished: " + dtf.format(LocalDateTime.now()) + "\n");
        	
        	
        }else if (function.equals("l0l1")) { 
        	
            String dc_out_file = gp.inter_dir + gp.project_name + "_dc_plan.txt"; // dc output
            
            // filepath string       	
            String[] vef_files = 
                    Entrance.get_filepaths(name_file, gp.inter_dir + "vef", "vef", false);
            String[] gcf_files = 
                    Entrance.get_filepaths(name_file, gp.inter_dir + "gcf", "gcf", false);     
            DivideConquer dc_maker = new DivideConquer(gs_var_pos,
                    gcf_files,
                    parameter_file,
                    dc_out_file);
            
            new File(gp.inter_dir + "/aem/").mkdir();   
            HashMap<Integer, String> index_var_prefix_dict = dc_maker.gs_map(gs_var_pos);
            
            String reg_dc_out_file = gp.inter_dir + 
            		gp.project_name + "_regression__dc_plan.txt"; // regression dc output
            
            HashMap<Integer, Integer> index_var_pos_dict = dc_maker.gs_map_pos(gs_var_pos);
            
            DivideConquer reg_dc_maker = new DivideConquer(dc_maker.regions_level_final, 
            		dc_maker.regions_level_final_link, reg_dc_out_file, gp.regression_maximum_regions, 
        			index_var_pos_dict);  
            
            System.out.println("DC Finished:\t" + dtf.format(LocalDateTime.now()) + "\n");
            
          
            dc_maker.level_final_config =  new HapConfig [dc_maker.num_regions_level_final ];
            dc_maker.level_final_link_config= new HapConfig [dc_maker.num_regions_level_final_link];
            for (int r=0 ; r< dc_maker.num_regions_level_final;r++) {
	            HapConfig final_aem_config =
		                new HapConfig(gp.inter_dir +"aem/"
		                		+gp.project_name
		                        + "_level_"
		                        + Integer.toString(dc_maker.final_level)
		                        + "_region_"
		                        + r
		                        + ".inter_freq_haps.txt", null);
	            dc_maker.level_final_config[r] = final_aem_config;
            }
            
            for (int r=0 ; r< dc_maker.num_regions_level_final_link;r++) {
            	HapConfig final_link_aem_config =
		                new HapConfig(gp.inter_dir +"aem/" 
		                		+ gp.project_name
		                        + "_level_"
		                        + Integer.toString(dc_maker.final_level+1)
		                        + "_region_"
		                        + r
		                        + ".inter_freq_haps.txt", null);
            	dc_maker.level_final_link_config[r] = final_link_aem_config;
            }
            
//            eva.Aem_V_VI_Evaluate(gp.project_name, dc_out_file, 
//                    gp.gold_dir+"/"+gp.project_name+"_haps.inter_freq_vars.txt", 
//                 		gp.inter_dir+"/aem/");
//            
            
            System.out
            .println("\nAEM Finished and Start regional L0L1 regression: " + dtf.format(LocalDateTime.now()) + "\n");
       
       String reg_prefix=  gp.inter_dir + "/regression/";
//------------------------------------------------------------
       try {  
           String shpath="rm -rf  "+ reg_prefix;  
           Process ps = Runtime.getRuntime().exec(shpath, null);  
           ps.waitFor();  
           BufferedReader br = new BufferedReader(new InputStreamReader(ps.getInputStream()));  
           StringBuffer sb = new StringBuffer();  
           String line;  
           while ((line = br.readLine()) != null) {  
               sb.append(line).append("\n");  
           }  
           String result = sb.toString();  
           br.close();
       }   catch (Exception e) {  
           e.printStackTrace();  
       }  
      
       
       new File(reg_prefix).mkdir();
       new File(reg_prefix+"/Rscript").mkdir();
       
       for (int pool_index = 0; pool_index < num_pools; pool_index++) {
       	new File(gp.inter_dir + "/regression/"+ Entrance.names_array[pool_index]).mkdir();
           for (int r=0; r< reg_dc_maker.regession_level[0].length;r++) {
           	
           	GraphColoring region_linker =
                       new GraphColoring( reg_dc_maker.regession_level[0][r][0],  
                       		reg_dc_maker.regession_level[0][r][1],
                       		dc_maker.level_final_config, dc_maker.level_final_link_config, 
                       		gs_var_pos, dc_maker.regions_level_final,
                       		dc_maker.regions_level_final_link, gp.bfs_mismatch_tolerance, 
                       		r+1 ,index_var_prefix_dict,  reg_prefix+ Entrance.names_array[pool_index],
                       		gp.regression_maximum_regions);
           	
           }
       }
       
       System.out
       .println("Have generated initialization regional haptypes from AEM level "+ dc_maker.final_level+
       		" regional haplotypes. " + dtf.format(LocalDateTime.now()) + "\n");
       
       for (int regression_level=0; regression_level< reg_dc_maker.regession_level.length; 
       		regression_level++) {
       	for (int r=0; r< reg_dc_maker.regession_level[regression_level].length;r++) {
       		for (int pool_index = 0; pool_index < num_pools; pool_index++) {
	            	System.out
	                .println("Generating level "+Integer.toString(regression_level+1) + 
	                		" region "+ Integer.toString(r+1)+ "  potential haplotypes for pool "
	                		+Entrance.names_array[pool_index] + ". "
	                		 + dtf.format(LocalDateTime.now()) );
	            	if ((reg_dc_maker.regession_level[regression_level][r][0]!=-1 ) && 
	            			(reg_dc_maker.regession_level[regression_level][r][1]!=-1 )) {
	            		if (regression_level!=0) {
			            	GraphColoring regression_region_linker =
			            	    new GraphColoring( reg_dc_maker.regession_level[regression_level][r][0],
			            	    	reg_dc_maker.regession_level[regression_level][r][1],
				                    reg_dc_maker.regession_level[regression_level-1],  
				                    dc_maker.level_final_link_config, 
				                    gs_var_pos,
				                    dc_maker.regions_level_final_link, 
				                    gp.regression_link_mismatch_tolerance, 
				                    r+1 ,regression_level+1, index_var_prefix_dict,  
				                    reg_prefix+ Entrance.names_array[pool_index], 
				                    gp.regression_maximum_regions, 
				                    gp.regression_maximum_selected_haplotypes);
	            		}
	            		
	            		HapConfig final_global_haps =
		    	                new HapConfig(gp.inter_dir + "/regression/" + Entrance.names_array[pool_index] +
		    	                		"/regression_level_"+   Integer.toString(regression_level+1)+
		    	                		"_region_"+ Integer.toString(r+1)+".potential.haps", null);
		            	SiteInPoolFreqAnno siteInPoolFreqAnno = new SiteInPoolFreqAnno(gs_var_pos, 
		            			reg_dc_maker.regession_level[regression_level][r][0],
		            			reg_dc_maker.regession_level[regression_level][r][1], pool_index);
		            	final_global_haps.inpool_site_freqs = siteInPoolFreqAnno.inpool_freqs;
		            	final_global_haps.num_pools = Entrance.num_pools;
		            	final_global_haps.pool_IDs = siteInPoolFreqAnno.pool_IDs;
		            	final_global_haps.in_pool_haps_freq =
		    	                new double[final_global_haps.num_global_hap][final_global_haps.num_pools];
		            	
		            	HapLASSO inpool_lasso = new HapLASSO(pool_index,
		                        Entrance.names_array[pool_index],
		                        gp.lasso_global_lambda,
		                        final_global_haps,
		                        gp.lasso_full_hap_freq_cutoff,
		                        gp.lasso_global_memory,
		                        gp.inter_dir + "/regression/" + Entrance.names_array[pool_index]+
		                        "/regression_level_"+Integer.toString(regression_level+1)+
		                        "_region_"+ Integer.toString(r+1) ,
		                    	gp.lasso_coverage_weight,
		                    	gp.lasso_distance_max_weight);
		            	
		            	
		            	inpool_lasso.estimate_frequencies_lasso(gp.inter_dir + "/vef/"
		                        + Entrance.names_array[pool_index] + ".vef", null, gp.lasso_weights ,
		                        gp.sequencing_technology, 
		                        reg_dc_maker.regession_level_pos_come_from.get(regression_level).get(r),
		                        reg_dc_maker.regession_level_num_come_from[regression_level][r], 
		                        index_var_pos_dict, 
		                        regression_level);
	            	}
	            }
           }
       	
           System.out
           .println("\nRunning L0L1 regression for all regions in level "+Integer.toString(regression_level+1) +
           		". "+ dtf.format(LocalDateTime.now()) + "\n");
           
           HapLASSO regression_rscript  = new HapLASSO(Entrance.names_array, 
           		reg_dc_maker.regession_level[regression_level],
           		gp.inter_dir + "/regression/Rscript/" , 
           		gp.inter_dir + "/regression/" , 
           		regression_level+1,
           		gp.regression_hapset_size_max,
           		gp.regression_hapset_size_min,
           		gp.regression_gamma_min,
           		gp.regression_gamma_max,
           		gp.ngammas, 
           		gp.lasso_weights[1], gp.species);
           
       for (int r=0; r< reg_dc_maker.regession_level[regression_level].length;r++) {
           	
           	if ((reg_dc_maker.regession_level[regression_level][r][0]!=-1 ) && 
           			(reg_dc_maker.regession_level[regression_level][r][1]!=-1 )) {
           		regression_rscript.run_L0L1learn( Entrance.names_array,
           				gp.rscript_path, gp.inter_dir + "/regression/Rscript/regression_level_"+
           				Integer.toString(regression_level+1)+ "_region_"+Integer.toString(r+1), 
           				 gp.inter_dir + "/regression/",
           				regression_level+1, 
           				r+1, 
           				gp.num_threads
           				);
           		System.out
	                .println("Level "+Integer.toString(regression_level+1) +" region "+ 
	                		Integer.toString(r+1)+" L0L1 regression finished.\t"+ 
	                		dtf.format(LocalDateTime.now()) + "\n");
           	}
           }
       }
       

       int regression_last_level= reg_dc_maker.regession_level.length-1;
       
       for (int pool_index = 0; pool_index < num_pools; pool_index++) {
       	new File(gp.out_dir + Entrance.names_array[pool_index]).mkdir();
       	HapConfig final_haps = new HapConfig (gp.inter_dir + "/regression/" + Entrance.names_array[pool_index] +
           		"/regression_level_"+   Integer.toString(regression_last_level+1)+
           		"_region_1.regression_out", 
           		gp.out_dir  + Entrance.names_array[pool_index]+ "/final_freq_haps.txt" ,
           		index_var_prefix_dict);
       }
       
       System.out.println("PoolHapX Successfully Finished, Enjoy!\n");
       
      }
    }
}
