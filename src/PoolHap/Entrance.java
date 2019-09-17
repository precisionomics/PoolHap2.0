package PoolHap;

import java.io.PrintWriter;
import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
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
    	String[] supported_functions_array = {"format", "gc", "aem", "lasso", 
    			"split" ,"clustering"};
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
            BAMFormatterGATK.generate_vef_varfreq(gp.input_dir,
                gp.inter_dir,
                gp.project_name,
                sam_files,
                Entrance.name2index);
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
                GraphColoring pool_in = new GraphColoring(vef_files[p],
                    gs_var_pos,
                    gp.inter_dir + "gcf/" + Entrance.names_array[p] + ".gcf" ,
                    gp.num_pos_window, gp.num_gap_window);

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
            
            System.out.println("DC Finished: " + dtf.format(LocalDateTime.now()) + "\n");
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
	            
	
	            // Apply LASSO for each pool from final global haplotypes to estimate frequencies of
	            // each
	            // in each pool.
	            // String prefix=null;
        	
            for (int pool_index = 0; pool_index < num_pools; pool_index++) {
                HapLASSO inpool_lasso = new HapLASSO(pool_index,
                    Entrance.names_array[pool_index],
                    gp.lasso_global_lambda,
                    final_global_haps,
                    gp.lasso_full_hap_freq_cutoff,
                    gp.lasso_global_memory,
                    gp.inter_dir + "/lasso/" + Entrance.names_array[pool_index],
                	gp.lasso_coverage_weight,
                	gp.lasso_distance_max_weight);

                inpool_lasso.estimate_frequencies_lasso(gp.inter_dir + "/vef/"
                    + Entrance.names_array[pool_index] + ".vef", null, gp.lasso_weights);

                // Learn optimal LASSO lambda value.
                //double new_penalty = gp.global_lambda;
                //while (inpool_lasso.r2 == gp.min_r2) {something deleted}
                String[] single_pool_IDs=new String[] {names_array[pool_index]};
                //     add final intra pool
                final_inpool_haps[pool_index] = inpool_lasso.hapOut(single_pool_IDs); 
                // HapConfig to array
                System.out.println(final_inpool_haps[pool_index].num_global_hap
                    + " haplotypes were generated for pool " + pool_index + ":"
                    + Entrance.names_array[pool_index] + ".");
            }

            System.out.println("\nGlobal LASSO Finished: " + dtf.format(LocalDateTime.now()) + "\n");

            // Writing final outputs to file.
            HapConfig final_reconstruction = new HapConfig(final_inpool_haps);
            System.out.println("Haplotype reconstruction finished. There are "
                + final_reconstruction.num_global_hap + " global haplotypes.");

            final_reconstruction.write2files(gp.out_dir + gp.project_name + ".inter_freq_haps.txt",
                gp.out_dir + gp.project_name + ".intra_freq_haps.txt",
                "string");
            
        	eva.MCCEvaluate(gp.inter_dir +"/../gold_standard/"+ gp.project_name+"_haps.inter_freq_haps.txt",
        			gp.out_dir+"/"+gp.project_name+".inter_freq_haps.txt", gp.mcc_freq_cutoff);
        	
        	
        }else if (function.equals("clustering")) {
            /*
             * Chen: cluster the potential haplotypes using  hierarchical clustering
             * 
             */
	           	HapConfig clustering_global_haps =
	                new HapConfig(gp.out_dir + gp.project_name + "_gc.inter_freq_haps.txt", null);
	           	HierarchicalClustering hc= new HierarchicalClustering(clustering_global_haps.hap_IDs, 
	           			clustering_global_haps.global_haps_string, clustering_global_haps.hapID2index, 
	           			gp.out_dir + gp.project_name + "_hc.haps.txt", gp.hc_similarity_cutoff);
	           	hc.gc_solver(gs_var_pos);
	           	Entrance.get_filepaths(name_file, gp.inter_dir + "/vef/", "vef", false);
	           	
	           	HapConfig hc_global_haps; // final global haplotype configuration object
	           	hc_global_haps = hc.hapOut(Entrance.names_array); // HapConfig from Hierarchical Clustering                                                           
	           	hc_global_haps.recode_HapIDs_to_base16();
	           	hc_global_haps.write_global_file_string( // write to output
                    gp.out_dir + gp.project_name + "_hc.inter_freq_haps.txt");

        }else if (function.equals("aem")) {
        	
//        	eva.GcAemEvaluate("/home/chencao/Desktop/sim001/gold_standard/sim001_haps.txt", 
//        			"/home/chencao/Desktop/sim001/output/sim001_gc.inter_freq_haps.txt");
//        	
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
            
            for (int p = 0; p < num_pools; p++) {
                // TODO: [Question]:: this "pool_in" variable not used anywhere.
                // Previously added to initial_local_haps but now commented out. ML 20190702
                GraphColoring pool_in = new GraphColoring(vef_files[p],
                    gs_var_pos,
                    gp.inter_dir + "gcf/" + Entrance.names_array[p] + ".gcf",
                    dc_maker.regions_level_I );
                dc_maker.loci_link_freq[p]= pool_in.loci_link_freq;
                System.out.println("Graph colouring for pool " + p + ":" + Entrance.names_array[p]
                    + " is finished.");
            }
            
//          new File(gp.inter_dir + "/gcaem/").mkdir();
            
            boolean has_run_bfs= false;
            
            HapConfig[] level_I_config = dc_maker.regional_AEM(
                    Entrance.names_array,
                    vef_files,
                    dc_maker.regions_level_I,
                    parameter_file,
                    gp.inter_dir + "/aem/" + gp.project_name,
                    gp.inter_dir + "/aem_fail_regional_lasso/",
                    1, Entrance.name2index);
            dc_maker.level_I_config= level_I_config;
            	System.out.println("Level 1 AEM Finished: " + dtf.format(LocalDateTime.now()) + "\n");
            	
            HapConfig[] level_II_config = dc_maker.regional_AEM(
                    Entrance.names_array,
                    vef_files,
                    dc_maker.regions_level_II,
                    parameter_file,
                    gp.inter_dir + "/aem/" + gp.project_name,
                    gp.inter_dir + "/aem_fail_regional_lasso/",
                    2,  Entrance.name2index);
            dc_maker.level_II_config= level_II_config;
            System.out.println("Level 2 AEM Finished: " + dtf.format(LocalDateTime.now()) + "\n");
            
            if ((dc_maker.regions_level_III.length==0) || (dc_maker.regions_level_IV.length==0)) {
            	has_run_bfs= true;
            	GraphColoring region_linker =
                        new GraphColoring(level_I_config, level_II_config, 
                        		gs_var_pos, dc_maker.regions_level_I,
                        		dc_maker.regions_level_II, gp.bfs_mismatch_tolerance);
            	HapConfig final_global_haps; // final global haplotype configuration object
                
                final_global_haps = region_linker.hapOut(Entrance.names_array); // HapConfig from
                                                                                // GC-linked regions
                final_global_haps.recode_HapIDs_to_base16();
                final_global_haps.write_global_file_string(
                		gp.out_dir + gp.project_name + "_gc.inter_freq_haps.txt");
            }
            if (!has_run_bfs) {
	            HapConfig[] level_III_config = dc_maker.level_III_regional_AEM( 
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
	            
	            HapConfig[] level_IV_config = dc_maker.level_IV_regional_AEM( 
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
	            if ((dc_maker.regions_level_V.length==0) || (dc_maker.regions_level_VI.length==0)) {
	            	has_run_bfs= true;
	            	GraphColoring region_linker =
	                        new GraphColoring(level_III_config, level_IV_config, 
	                        		gs_var_pos, dc_maker.regions_level_III,
	                        		dc_maker.regions_level_IV, gp.bfs_mismatch_tolerance);
	            	HapConfig final_global_haps; // final global haplotype configuration object
	                
	                final_global_haps = region_linker.hapOut(Entrance.names_array); // HapConfig from
	                                                                                // GC-linked regions
	                final_global_haps.recode_HapIDs_to_base16();
	                final_global_haps.write_global_file_string( // write to output
	                		
	                		
	                gp.out_dir + gp.project_name + "_gc.inter_freq_haps.txt");
	            }
	            
            }

                
//            eva.AemEvaluate(gp.project_name, dc_out_file, 
//               "/home/chencao/Desktop/sim001/gold_standard/sim001_haps.txt", 
//            		"/home/chencao/Desktop/sim001/intermediate/aem/");
//            eva.Level_III_AemEvaluate(gp.project_name, dc_out_file, 
//                    "/home/chencao/Desktop/sim001/gold_standard/sim001_haps.txt", 
//                 		"/home/chencao/Desktop/sim001/intermediate/aem/"); 
//            eva.Level_IV_AemEvaluate(gp.project_name, dc_out_file, 
//                    "/home/chencao/Desktop/sim001/gold_standard/sim001_haps.txt", 
//                 		"/home/chencao/Desktop/sim001/intermediate/aem/");
            if (!has_run_bfs) {
	            HapConfig[] level_V_config = dc_maker.level_V_regional_AEM( 
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
	           
	//            eva.Level_V_AemEvaluate(gp.project_name, dc_out_file, 
	//                    "/home/chencao/Desktop/sim001/gold_standard/sim001_haps.txt", 
	//                 		"/home/chencao/Desktop/sim001/intermediate/aem/");
	            
	            HapConfig[] level_VI_config = dc_maker.level_VI_regional_AEM( 
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
	            
	//            eva.Level_VI_AemEvaluate(gp.project_name, dc_out_file, 
	//                    "/home/chencao/Desktop/sim001/gold_standard/sim001_haps.txt", 
	//                 		"/home/chencao/Desktop/sim001/intermediate/aem/");
	            
	            GraphColoring region_linker =
	                    new GraphColoring(level_V_config, level_VI_config, 
	                    		gs_var_pos, dc_maker.regions_level_V,
	                    		dc_maker.regions_level_VI, gp.bfs_mismatch_tolerance);
	            HapConfig final_global_haps; // final global haplotype configuration object
                
                final_global_haps = region_linker.hapOut(Entrance.names_array); // HapConfig from
                                                                                // GC-linked regions
                final_global_haps.recode_HapIDs_to_base16();
                final_global_haps.write_global_file_string( 
                		gp.out_dir + gp.project_name + "_gc.inter_freq_haps.txt");
            }
            
            
                // Link regions by applying Breadth-First-Search across the level 1 and level 2 regional
                // haplotype configurations (Pruning strategy, Chen).
//                
//            GraphColoring region_linker =
//                    new GraphColoring(level_I_config, level_II_config, gs_var_pos, dc_maker.regions_level_I,
//                    		dc_maker.regions_level_II);
                
//                eva.StringEvaluate("/home/chencao/Desktop/sim001/gold_standard/sim001_haps.txt", 
//                		"/home/chencao/Desktop/ccc.txt");
                
//                eva.AemEvaluate("/home/chencao/Desktop/sim001/gold_standard/sim001_haps.txt", 
//            			"/home/chencao/Desktop/sim001/output/sim001.inter_freq_haps.txt");
//                System.exit(0);
                
              
//                GraphColoring region_linker =
//                        new GraphColoring(level_I_config, level_II_config, gs_var_pos, gp.virtual_cov_link_gc);

                // Write final global haplotype configurations (inter pool) to output.
                
                
//                eva.GcAemEvaluate("/home/chencao/Desktop/sim001/gold_standard/sim001_haps.txt", 
//            			"/home/chencao/Desktop/sim001/output/sim001_gc.inter_freq_haps.txt");
            
            System.out
                .println("\nAEM+BFS Finished: " + dtf.format(LocalDateTime.now()) + "\n");
            
        }    
    }
}
