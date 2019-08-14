package PoolHap;

import java.io.PrintWriter;
import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import PoolHap.Parameters.GenParameters;
import java.time.format.DateTimeFormatter;
import java.util.Arrays;
import java.util.HashMap;
import java.util.HashSet;

import MiscFunctions.BAMFormatterGATK;

import java.time.LocalDateTime;

public class Entrance_pf {
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
                Entrance_pf.num_pools = files.length;
                Entrance_pf.names_array = new String[Entrance_pf.num_pools];
                for (int index = 0; index < Entrance_pf.num_pools; index++) {
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
                        Entrance_pf.name2index.put(name, index);
                        Entrance_pf.names_array[index] = name;
                    }
                }
                // write the names to the document for future functions to load
                BufferedWriter bw = new BufferedWriter(new FileWriter(name_file));
                for (int p = 0; p < Entrance_pf.num_pools; p++) {
                    bw.write(names_array[p] + "\n");
                }
                bw.close();
            } else { // load names from the "name_file" and check whether they are correct.
                BufferedReader br = new BufferedReader(new FileReader(name_file));
                String line = br.readLine();
                int p_index = 0;
                while (line != null) {
                    Entrance_pf.name2index.put(line, p_index++);
                    line = br.readLine();
                }
                br.close();
                Entrance_pf.num_pools = p_index;
                Entrance_pf.names_array = new String[Entrance_pf.num_pools];
                for (String name : Entrance_pf.name2index.keySet()) {
                    Entrance_pf.names_array[Entrance_pf.name2index.get(name)] = name;
                }
                // check the correctness in the target folder.
                if (Entrance_pf.num_pools != files.length) { // first, number must be correct.
                    System.out.println("ERROR: No. of files in " + folder + " is " + files.length
                        + ", NOT consistent to the No. of lines in the " + name_file + ": "
                        + Entrance_pf.num_pools);
                    System.exit(0);
                }
                filepaths = new String[files.length];
                // observed_files are the files in the folder, may or may not align the names
                // and order specified by name_file. So need to check and return ordered paths
                String[] observed_files = new String[files.length];
                for (int k = 0; k < files.length; k++) {
                    observed_files[k] = files[k].toString();
                }
                for (int index = 0; index < Entrance_pf.num_pools; index++) {
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
                        if (!Entrance_pf.name2index.containsKey(name)) {
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
        String parameter_file = args[0]; // parameter file path
        // int num_pools = Integer.parseInt(args[2]); // number of pools present
        // removed by Quan Long 2019-07-03. This will be identified by counting how many files in a
        // folder

        /*
         * Program initialization.
         */
        // New general parameter set object from parameter file.
        // TODO: [Review]:: GenParameters(args[0]) changed to GenParameters(parameter_file)

        GenParameters gp = new GenParameters(parameter_file);

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
        // TODO: LEFTOVER ML 20190702
        // HapConfig[] initial_local_haps = new HapConfig[num_pools]; // HapConfig array; 1 per pool

        // By design, HapConfig contains haps of all pools.
        if (gp.function.equals("format")) {
            // SAM files in input_dir/sam initiates pool-IDs and orders, and write to the
            // "name_file".
            // The function below generates vef files and inpool-var-freq file.
            String[] sam_files =
                Entrance_pf.get_filepaths(name_file, gp.input_dir + "/sam/", "sam", true);
            BAMFormatterGATK.generate_vef_varfreq(gp.input_dir,
                gp.inter_dir,
                gp.project_name,
                sam_files,
                Entrance_pf.name2index);
        } else if (gp.function.equals("gc")) {
            // PrintWriter in_list = new PrintWriter( // writer object for graph colored pools
            // new FileWriter(new File(gp.inter_dir + prefix + "_p.in.list")));
            // the folder "vef/" should exist under gp.inter_dir.
            String[] vef_files =
                Entrance_pf.get_filepaths(name_file, gp.inter_dir + "vef", "vef", true);
            // make a new dir "in/" to store GC-outcome files.
            new File(gp.inter_dir + "/gcf/").mkdir();
            // Apply graph coloring for each pool. Note, initial "global" haplotypes are roughly
            // linked by GC with gaps between VARIANTS (not regions).
            for (int p = 0; p < num_pools; p++) {
                // TODO: [Question]:: this "pool_in" variable not used anywhere.
                // Previously added to initial_local_haps but now commented out. ML 20190702
                GraphColoring pool_in = new GraphColoring(vef_files[p],
                    gs_var_pos,
                    gp.inter_dir + "gcf/" + Entrance_pf.names_array[p] + ".gcf");

                System.out.println("Graph colouring for pool " + p + ":" + Entrance_pf.names_array[p]
                    + " is finished.");
                // in_list.println(gp.inter_dir + prefix + "_p" + p + ".in");
                // write GC pool output

                // TODO: LEFTOVER
                // initial_local_haps[p] = pool_in.hapOut(); // add pool HapConfig to array
            }

            // in_list.close();
            System.out
                .println("\nGraph-Coloring Finished: " + dtf.format(LocalDateTime.now()) + "\n");
        } else if (gp.function.equals("aem")) { // Note: aem includes DC and AEM
            // Apply divide and conquer across all pools.
            String dc_out_file = gp.inter_dir + gp.project_name + "_dc_plan.txt"; // dc output
                                                                                  // filepath string
            String[] vef_files =
                Entrance_pf.get_filepaths(name_file, gp.inter_dir + "vef", "vef", false);
            String[] gcf_files = 
                Entrance_pf.get_filepaths(name_file, gp.inter_dir + "gcf", "gcf", false);
            DivideConquer dc_maker = new DivideConquer(gs_var_pos,
                gcf_files,
                // gp.inter_dir + prefix + "_p.in.list",
                parameter_file,
                dc_out_file);
            
            System.out.println("DC Finished: " + dtf.format(LocalDateTime.now()) + "\n");
            new File(gp.inter_dir + "/aem/").mkdir();
            /*
             * Global haplotype configuration.
             */

            // TODO: LEFTOVER ML 20190702
            // // If there are enough gaps to run local reconstruction...
            // if (dc_maker.region_solve) {

            // Resolve level 1 regions with AEM; regional LASSO if AEM fails.
            
         
            HapConfig[] level_I_config = dc_maker.regional_AEM(
                Entrance_pf.names_array,
                vef_files,
                dc_maker.regions_level_I,
                parameter_file,
                gp.inter_dir + "/aem/" + gp.project_name,
                gp.inter_dir + "/aem_fail_regional_lasso/",
                1);

            System.out.println("Level 1 AEM Finished: " + dtf.format(LocalDateTime.now()) + "\n");

            // Resolve level 2 regions with the same strategy as above.
            HapConfig[] level_II_config = dc_maker.regional_AEM(
                Entrance_pf.names_array,
                vef_files,
                dc_maker.regions_level_II,
                parameter_file,
                gp.inter_dir + "/aem/" + gp.project_name,
                gp.inter_dir + "/aem_fail_regional_lasso/",
                2);

            System.out.println("Level 2 AEM Finished: " + dtf.format(LocalDateTime.now()) + "\n");

            // Link regions by applying graph coloring across the level 1 and level 2 regional
            // haplotype configurations.
            GraphColoring region_linker =
                new GraphColoring(level_I_config, level_II_config, gs_var_pos, gp.fragments);

            // Write final global haplotype configurations (inter pool) to output.
            HapConfig final_global_haps; // final global haplotype configuration object
            final_global_haps = region_linker.hapOut(Entrance_pf.names_array); // HapConfig from
                                                                            // GC-linked regions
            final_global_haps.recode_HapIDs_to_base16();
            final_global_haps.write_global_file_string( // write to output
                gp.out_dir + gp.project_name + "_gc.inter_freq_haps.txt");

            System.out.println("\nPost-AEM GC Finished: " + dtf.format(LocalDateTime.now()) + "\n");

            // TODO: LEFTOVER ML 20190702
            // // Else, not enough gaps to run regional construction; existing haplotype
            // // configurations deemed "global" enough already.
            // } else {
            // final_global_haps = new HapConfig(
            // gp.inter_dir + prefix + "_p",
            // gs_var_pos,
            // num_pools);

            // // final_global_haps = new HapConfig(initial_local_haps);
            // final_global_haps.write_global_file_string( // write to output
            // gp.out_dir + prefix + "_gc.inter_freq_vars.txt",
            // false);
            // }
        } else if (gp.function.equals("lasso")) {
            /*
             * Final intra pool haplotype configuration via LASSO selection. TODO: [Question]:: are
             * these global intra pool haplotypes?
             */
            new File(gp.inter_dir + "/lasso/").mkdir();
            Entrance_pf.get_filepaths(name_file, gp.inter_dir + "/gcf/", "gcf", false);
            HapConfig final_global_haps =
                new HapConfig(gp.out_dir + gp.project_name + "_gc.inter_freq_haps.txt", null);

            SiteInPoolFreqAnno siteInPoolFreqAnno = new SiteInPoolFreqAnno(gs_var_pos);
            final_global_haps.inpool_site_freqs = siteInPoolFreqAnno.inpool_freqs;
            final_global_haps.num_pools = Entrance_pf.num_pools;
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
                    Entrance_pf.names_array[pool_index],
                    gp.lambda,
                    final_global_haps,
                    gp.final_cutoff,
                    "500000000",
                    gp.inter_dir + "/lasso/" + Entrance_pf.names_array[pool_index]);

                inpool_lasso.estimate_frequencies_lasso(gp.inter_dir + "/vef/"
                    + Entrance_pf.names_array[pool_index] + ".vef", null, gp.lasso_weights);

                // Learn optimal LASSO lambda value.
                double new_penalty = gp.lambda;
                while (inpool_lasso.r2 == gp.min_r2) {
                    new_penalty -= gp.lasso_penalty_step;
                    System.out
                        .println("Local full-length LASSO has failed. Adjust the lambda penalty to "
                            + new_penalty + " and trying again.");

                    inpool_lasso.lambda = new_penalty;
                    inpool_lasso.estimate_frequencies_lasso(gp.inter_dir + "/vef/"
                        + Entrance_pf.names_array[pool_index] + ".vef", null, gp.lasso_weights);

                }
                String[] single_pool_IDs=new String[] {names_array[pool_index]};
                //     add final intra pool
                final_inpool_haps[pool_index] = inpool_lasso.hapOut(single_pool_IDs); 
                // HapConfig to array
                System.out.println(final_inpool_haps[pool_index].num_global_hap
                    + " haplotypes were generated for pool " + pool_index + ":"
                    + Entrance_pf.names_array[pool_index] + ".");

            }

            System.out.println("\nLASSO Finished: " + dtf.format(LocalDateTime.now()) + "\n");

            // Writing final outputs to file.
            HapConfig final_reconstruction = new HapConfig(final_inpool_haps);
            System.out.println("Haplotype reconstruction finished. There are "
                + final_reconstruction.num_global_hap + " global haplotypes.");

            final_reconstruction.write2files(gp.out_dir + gp.project_name + ".inter_freq_haps.txt",
                gp.out_dir + gp.project_name + ".intra_freq_haps.txt",
                "string");
        }
    }
}
