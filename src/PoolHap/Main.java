package PoolHap;

import java.io.PrintWriter;
import java.io.File;
import java.io.FileWriter;
import PoolHap.Parameters.GenParameters;
import java.time.format.DateTimeFormatter;
import java.util.HashMap;
import java.time.LocalDateTime;

public class Main {
    // TODO: [ReconEP]:: refactor variable names to be more sensible (e.g. camelCase for Java).
    // TODO: [Question]:: confusion over the HapLASSO portion.

    public HashMap<String, Integer> name2index = new HashMap<String, Integer>();
    public String[] names_array;

    public static void main(String[] args) throws Exception {
        /*
         *  Input arguments.
         */
        String parameter_file = args[0]; // parameter file path
        String prefix = args[1]; // output prefix
        int num_pools = Integer.parseInt(args[2]); // number of pools present


        /*
         *  Program initialization.
         */
        // New general parameter set object from parameter file.
        // TODO: [Review]:: GenParameters(args[0]) changed to GenParameters(parameter_file)
        GenParameters gp = new GenParameters(parameter_file);

        // Gold standard variant positions file path (intra pool).
        String gs_var_pos = gp.inter_dir + prefix + "_vars.intra_freq.txt";
        String dc_out_file = gp.inter_dir + prefix + "_dc_plan.txt"; // dc output filepath string

        // Print start time.
        DateTimeFormatter dtf = DateTimeFormatter.ofPattern("yyyy/MM/dd HH:mm:ss");
        System.out.println(
            "PoolHapX " + gp.function + " Initiated: " + dtf.format(LocalDateTime.now()));

        /*
         *  Local haplotype configuration.
         */
        // TODO: LEFTOVER ML 20190702
        // HapConfig[] initial_local_haps = new HapConfig[num_pools]; // HapConfig array; 1 per pool

        // By design, HapConfig contains haps of all pools.
        if (gp.function.equals("gc")) {
            PrintWriter in_list = new PrintWriter( // writer object for graph colored pools
                new FileWriter(new File(gp.inter_dir + prefix + "_p.in.list")));

            // Apply graph coloring for each pool. Note, initial "global" haplotypes are roughly
            // linked by GC with gaps between VARIANTS (not regions).
            for (int p = 0; p < num_pools; p++) {
                // TODO: [Question]:: this "pool_in" variable not used anywhere.
                // Previously added to initial_local_haps but now commented out. ML 20190702
                GraphColoring pool_in = new GraphColoring(
                    gp.inter_dir + prefix + "_p" + p + ".vef",
                    gs_var_pos,
                    gp.inter_dir + prefix + "_p" + p + ".in");

                System.out.println("Graph colouring for pool " + p + " is finished.");
                in_list.println(gp.inter_dir + prefix + "_p" + p + ".in"); // write GC pool output

                // TODO: LEFTOVER
                // initial_local_haps[p] = pool_in.hapOut(); // add pool HapConfig to array

            }

            in_list.close();
            System.out.println(
                "\nGraph-Coloring Finished: " + dtf.format(LocalDateTime.now()) + "\n");

        } else if (gp.function.equals("aem")) { // Note: aem includes DC and AEM
            // Apply divide and conquer across all pools.
            DivideConquer dc_maker = new DivideConquer(
                gs_var_pos,
                gp.inter_dir + prefix + "_p.in.list",
                parameter_file,
                dc_out_file);

            System.out.println("DC Finished: " + dtf.format(LocalDateTime.now()) + "\n");

            /*
             *  Global haplotype configuration.
             */

            // TODO: LEFTOVER ML 20190702
            // // If there are enough gaps to run local reconstruction...
            // if (dc_maker.region_solve) {

            // Resolve level 1 regions with AEM; regional LASSO if AEM fails.
            HapConfig[] level_I_config = dc_maker.analyze_regions(
                dc_maker.regions_level_I,
                parameter_file,
                gp.inter_dir + prefix,
                1);

            System.out.println("Level 1 AEM Finished: "
                + dtf.format(LocalDateTime.now())
                + "\n");

            // Resolve level 2 regions with the same strategy as above.
            HapConfig[] level_II_config = dc_maker.analyze_regions(
                dc_maker.regions_level_II,
                parameter_file,
                gp.inter_dir + prefix,
                2);

            System.out.println("Level 2 AEM Finished: "
                + dtf.format(LocalDateTime.now())
                + "\n");

            // Link regions by applying graph coloring across the level 1 and level 2 regional
            // haplotype configurations.
            GraphColoring region_linker = new GraphColoring(
                level_I_config,
                level_II_config,
                gs_var_pos,
                gp.fragments);

            // Write final global haplotype configurations (inter pool) to output.
            HapConfig final_global_haps; // final global haplotype configuration object
            final_global_haps = region_linker.hapOut(); // HapConfig from GC-linked regions
            final_global_haps.write_global_file_string( // write to output
                gp.out_dir + prefix + "_gc.inter_freq_vars.txt",
                false);

            System.out.println("\nPost-AEM GC Finished: " + dtf.format(LocalDateTime.now()) + "\n");

            // TODO: LEFTOVER ML 20190702
            // // Else, not enough gaps to run regional construction; existing haplotype
            // // configurations deemed "global" enough already.
            // } else {
            //     final_global_haps = new HapConfig(
            //         gp.inter_dir + prefix + "_p",
            //         gs_var_pos,
            //         num_pools);

            //     // final_global_haps = new HapConfig(initial_local_haps);
            //     final_global_haps.write_global_file_string( // write to output
            //             gp.out_dir + prefix + "_gc.inter_freq_vars.txt",
            //             false);
            // }

        } else if(gp.function.equals("lasso")) {
            /*
             *  Final intra pool haplotype configuration via LASSO selection.
             *  TODO: [Question]:: are these global intra pool haplotypes?
             */
            HapConfig final_global_haps = new HapConfig(
                gp.out_dir + prefix + "_gc.inter_freq_vars.txt",
                null);

            SiteInPoolFreqAnno siteInPoolFreqAnno = new SiteInPoolFreqAnno(gs_var_pos);
            final_global_haps.inpool_site_freqs = siteInPoolFreqAnno.inpool_freqs;
            final_global_haps.num_pools = num_pools;
            final_global_haps.pool_IDs = siteInPoolFreqAnno.pool_IDs;
            final_global_haps.in_pool_haps_freq =
                    new double[final_global_haps.num_global_hap][final_global_haps.num_pools];

            HapConfig[] final_local_haps = new HapConfig[num_pools]; // final intra pool HapConfigs

           // Apply LASSO for each pool from final global haplotypes to estimate frequencies of each
           // in each pool.
           for (int p = 0; p < num_pools; p++) {
               HapLASSO inpool_lasso = new HapLASSO(
                   p,
                   gp.lambda,
                   final_global_haps,
                   gp.final_cutoff,
                   "500000000",
                   gp.inter_dir + prefix + "_p" + p);

               inpool_lasso.estimate_frequencies_lasso(
                   gp.inter_dir + prefix + "_p" + p + ".vef",
                   null,
                   gp.lasso_weights);

               // Learn optimal LASSO lambda value.
               double new_penalty = gp.lambda;
               while (inpool_lasso.r2 == gp.min_r2) {
                   new_penalty -= gp.lasso_penalty_step;
                   System.out.println(
                       "Local full-length LASSO has failed. Adjust the lambda penalty to "
                       + new_penalty
                       + " and trying again.");

                   inpool_lasso.lambda = new_penalty;
                   inpool_lasso.estimate_frequencies_lasso(
                       gp.inter_dir + prefix + "_p" + p + ".vef",
                       null,
                       gp.lasso_weights);

               }

               final_local_haps[p] = inpool_lasso.hapOut(); // add final intra pool HapConfig to array
               System.out.println(final_local_haps[p].num_global_hap
                   + " haplotypes were generated for pool "
                   + p
                   + ".");

           }

           System.out.println("\nLASSO Finished: " + dtf.format(LocalDateTime.now()) + "\n");

           // Writing final outputs to file.
           HapConfig final_reconstruction = new HapConfig(final_local_haps);
           System.out.println("Haplotype reconstruction finished. There are "
               + final_reconstruction.num_global_hap
               + " global haplotypes.");

           final_reconstruction.write2files(
               gp.out_dir + prefix + ".inter_freq_vars.txt",
               gp.out_dir + prefix + ".intra_freq.txt",
               "string",
               false);
        }
    }
}
