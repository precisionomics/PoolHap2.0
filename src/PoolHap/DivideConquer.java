package PoolHap;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.HashSet;
import java.util.concurrent.ThreadLocalRandom;

import PoolHap.Parameters;
import PoolHap.HapConfig;
import PoolHap.SiteInPoolFreqAnno;
import PoolHap.LocusAnnotation;
import PoolHap.RegionEMSolver;

public class DivideConquer {
    /**
     *  @author  Quan Long. Oct 08, 2018
     *
     *  This class divides the full-length haplotype into multiple regions based on linkage support
     *  from reads measured by GC, solves for the regional haplotypes, and merges them. Recursively,
     *  the above procedure will be done multiple times.
     *
     *  There are two strategies on how to divide the genome:
     *  (1) Specify a fixed number of genetic markers (SNPs) for each region
     *  (2) Calculate the number of genetic markers in each region based on estimated LD structure
     *      The LD structure can be estimated by:
     *      (2.1) Another reference panel
     *      (2.2) Information from sequencing reads.
     *  Different layers of divide-and-Conquer may adopt different strategies.
     *
     *  All parameters are in the Object "parameters" to facilitate Machine Learning-based training.
     */

    public Parameters dp;

    public int num_pools;
    public int num_sites;
    public int num_regions_level_I;
    public int[][] regions_level_I;	// num_regsion_level_I x 2 (start and end)
    public int num_regions_level_II;
    public int[][] regions_level_II; // num_regsion_level_II x 2 (start and end)
    public String[][] global_haps_gc;
    public double[] global_gc_freq;
    public int num_haps_gc;

    // Note that we do NOT generate the fields for HapConfigs of each regions. These will be
    // generated on-the-fly during the calculation.

    // The in-pool frequencies and annotations for all sites.
    public SiteInPoolFreqAnno in_pool_sites_freq_anno;
    public String[][] gc_outcome;
    public boolean region_solve; // TODO: remove this when the debugging ends. Quan Long 2019-06-29


    /**
     *  Constructor that generates regions (i.e.: the dividing plan) based on GC outcomes.
     *
     *
     *  @param frequency_file
     *  @param gc_input_list
     *  @param parameter_file
     *  @param dc_out_file
     */
    public DivideConquer(
        String frequency_file,
        String[] gc_input_list,
        String parameter_file,
        String dc_out_file) {

        try {
            /*
             *  Load files.
             */
            this.dp = new Parameters(parameter_file);
            System.out.println("Finished loading the PoolHapX parameter file from " + parameter_file);
            //load_gc_outcome(parse_gc_input(gc_input_list)); 
            // above was removed and replaced by the line below by Quan Long 2019-07-07
            this.load_gc_outcome((gc_input_list));
            System.out.println("Finished loading graph-coloring files from " + dp.inter_dir+"gcf");
            System.out.println("Number of pools = "
                + this.num_pools
                + "\tNumber of segregating sites = "
                + this.num_sites);

            /*
             *  Division?
             */
            this.generate_dividing_plan_two_level(this.gc_outcome);

            // TODO: LEFTOVER ML 20190702
            // // If can be divided...
            // if (this.region_solve) {

            this.output_current_DC_plan(dc_out_file);
            System.out.println("The current dividing plan: ");
            System.out.println(
                "\nFinished generating divide-and-conquer plan. The plan has been written to "
                + dc_out_file);

            // TODO: LEFTOVER ML 20190702
            // // If not enough gaps...
            // }
            // else {
            //     System.out.println(
            //         "There are not enough locations of linkage uncertainty to run regional "
            //         + "reconstruction. Skipping straight to linkage-informed LASSO regression.");
            // }

            this.in_pool_sites_freq_anno = new SiteInPoolFreqAnno(frequency_file);

        } catch(Exception e) {
            e.printStackTrace();
        }
    }


    /**
     *  Read files recording the outcome of GC:
     *
     *  0?0-0-0?1-0-0-0-0-0-0-1-0-0-0-0-0?0-0-0?1?0?0-0-0-1-0-0?0-0\t25
     *  0?0-0-0?1-0-0-0-0-0-0-1-0-1-1-0-1?0-0-0?1?0?0-0-0-0-0-0?1-1\t10
     *  1?1-0-1?1-0-0-0-0-0-0-1-0-1-1-0-1?0-1-0?0?0?0-1-0-0-0-0?1-1\t8
     *
     *  It will return GC outcomes in the form of strings.
     *
     *  @param graph_coloring_outcome_files
     */
    public void load_gc_outcome(String[] graph_coloring_outcome_files) {
        this.num_pools = graph_coloring_outcome_files.length;
        this.gc_outcome = new String[num_pools][];
        try {
            for (int p_index = 0; p_index < num_pools; p_index++) {
                ArrayList<String> haps_list = new ArrayList<String>();
                BufferedReader br = new BufferedReader(
                    new FileReader(new File(graph_coloring_outcome_files[p_index])));

                String line = br.readLine();

                // TODO: [ReconEP]:: extract to static utils method.
                while (line.startsWith("#")) {
                    line = br.readLine();
                }

                // Changed from original code (+1 -> -1) because at the end of the draft haplotype
                // variant.
                this.num_sites = (line.split("\t")[0].length() + 1) / 2;

                // Composition, the output_ref_arr.get(x), which seems to be the count of that
                // draft.
                while (line != null) {
                    // Haplotype in the pool, is also outputted. Need to NOT read this into the
                    // variant composition.
                    haps_list.add(line);
                    line = br.readLine();
                }

                this.gc_outcome[p_index] = new String[haps_list.size()];
                for (int h = 0; h < haps_list.size(); h++) {
                    this.gc_outcome[p_index][h] = haps_list.get(h);
                }

                br.close();
            }

        } catch(Exception e) {
            e.printStackTrace();
        }
    }


    /**
     *
     *  @param gc_pathes_file
     *  @return
     */
    public String[] parse_gc_input(String gc_pathes_file) {
        ArrayList<String> pathes = new ArrayList<String>();
        try{
            BufferedReader br = new BufferedReader(new FileReader(gc_pathes_file));
            String line = br.readLine();

            // TODO: [ReconEP]:: extract to static utils method.
            while (line.startsWith("#")) { // skip headers, if any
                line = br.readLine();
            }

            while (line != null) {
                pathes.add(line);
                line = br.readLine();
            }

            this.num_pools = pathes.size();
            br.close();

        } catch(Exception e) {
            e.printStackTrace();
        }

        String[] gc_input_files = new String[this.num_pools];
        for (int p = 0; p < this.num_pools; p++) {
            gc_input_files[p] = pathes.get(p);
        }
        return gc_input_files;
    }


    /**
     *
     * @param graph_coloring_outcome
     */
    public void generate_dividing_plan_two_level(String[][] graph_coloring_outcome){
        // Step 1) Generate list of positions where there are linkage uncertainties i.e.: gaps.
        int[] gap_positions = identify_gaps(graph_coloring_outcome);

        // Step 2) Make the level 1 regions i.e.: windows of variants based on the generated gaps.
        // These are the boundaries of the regions/windows of variants.
        // TODO: [ReconEP]:: again, each step should be helpers?
        ArrayList<Integer> region_cuts = new ArrayList<Integer>();

        // This is the length of (number of variants in) the region/windows of variants so far.
        int curr_unmatched_len=0;
        for (int gap = 0; gap < gap_positions.length; gap++) { // For each gap...
            // The start position of the region currently being made.
            int previous_gap_position = (gap == 0)?0:gap_positions[gap - 1]; // No linting?

            // The number of variants in the region currently being made.
            int gap_len = gap_positions[gap] - previous_gap_position;
            // regardless whether the above is executed, add the gap_len to curr_unmatched_len
            curr_unmatched_len += gap_len;
            // If the newly updated curr_unmatched_len contains enough variants to fall within the
            // acceptable range, form the next region.
            if (inbetween(curr_unmatched_len, 
                dp.min_level_I_region_size,  dp.max_level_I_region_size)) {
                // Create the region.
                curr_unmatched_len = 0;
                region_cuts.add(gap_positions[gap]);
            } // ...otherwise, if the region currently being made is larger than the max allowable
            // size, split it up so that everything after the max size is added to the next region.
            else if (curr_unmatched_len > dp.max_level_I_region_size) {
                curr_unmatched_len = curr_unmatched_len - dp.max_level_I_region_size;
                int last_cut= (region_cuts.size()==0)?0:(region_cuts.get(region_cuts.size()-1));
                region_cuts.add(last_cut + dp.max_level_I_region_size);
                // the current_unmatched_len may still be larger than max_level_I_region_size
                // because of the gap_len just added is very large. 
                while(curr_unmatched_len > dp.max_level_I_region_size) {
                    curr_unmatched_len = curr_unmatched_len - dp.max_level_I_region_size;
                    last_cut= region_cuts.get(region_cuts.size()-1);
                    region_cuts.add(last_cut + dp.max_level_I_region_size);
                }// If the remained curr_unmatched_len contains enough variants to fall within the
                // acceptable range, form the next region.
                if (inbetween(curr_unmatched_len, 
                    dp.min_level_I_region_size,  dp.max_level_I_region_size)) {
                    // Create the region.
                    last_cut= region_cuts.get(region_cuts.size()-1);
                    region_cuts.add(last_cut+curr_unmatched_len);
                    curr_unmatched_len = 0;
                }
            } else { // curr_unmatched_len < parameters.min_level_I_region_size
                // do NOT form a region;
            }
        }

        this.num_regions_level_I = region_cuts.size() + 1;

        // The start and end positions of each region.
        this.regions_level_I = new int[num_regions_level_I][2];
        this.regions_level_I[0][0] = 0;

        // Note that because there are max region sizes, not every region will start and/or end at a
        // gap position.

        // Gap positions are merely guidelines to determine linkage regions.
        for (int r = 0; r < num_regions_level_I - 1; r++) {
            this.regions_level_I[r][1] = region_cuts.get(r) - 1;
            this.regions_level_I[r + 1][0] = region_cuts.get(r);
        }

        this.regions_level_I[num_regions_level_I - 1][1] = this.num_sites - 1;

        // Step 3) Make the level 2 regions i.e.: windows of variants that overlap with the end of
        // the ith level 1 window and the start of the (i + 1)th level 1 window.
        this.num_regions_level_II = this.num_regions_level_I - 1;
        this.regions_level_II = new int[num_regions_level_II][2];

        // !!! This is the midpoint of the 0th level 1 region.
        this.regions_level_II[0][0] = (this.regions_level_I[0][1] - this.regions_level_I[0][0]) / 2;
        for (int r = 0; r < num_regions_level_II - 1; r++) {
            int tmp_end = (this.regions_level_I[r+1][0] + this.regions_level_I[r + 1][1]) / 2;
            int tmp_reg_2_len = tmp_end - this.regions_level_II[r][0] + 1;
            if (inbetween(
                tmp_reg_2_len,
                dp.min_level_II_region_size,
                dp.max_level_II_region_size)) {

                this.regions_level_II[r][1] = tmp_end;

            } else if (tmp_reg_2_len > dp.max_level_II_region_size) {
                this.regions_level_II[r][1] = this.regions_level_II[r][0]
                    + dp.max_level_II_region_size;

            } else {
                this.regions_level_II[r][1] = tmp_end;
                while (tmp_reg_2_len < dp.min_level_II_region_size) {
                    this.regions_level_II[r][1]++;
                    tmp_reg_2_len = this.regions_level_II[r][1] - this.regions_level_II[r][0] + 1;

                }
            }

            this.regions_level_II[r + 1][0] = this.regions_level_II[r][1] + 1;
        }

        this.regions_level_II[num_regions_level_II - 1][1] = (this.num_sites
            + this.regions_level_II[num_regions_level_II - 1][0])
            / 2; // !!!
    }


    /**
     *
     *  @param value
     *  @param min
     *  @param max
     *  @return
     */
    public boolean inbetween(int value, int min, int max) {
        return (value >= min && value <= max);
    }


    /** TODO: remove this funtion when everything is stablized. 
     *  Identify gaps from the GC outcome.
     *  The GC outcome file is composed of inferred haplotypes in the format of:
     *  		0-1-0-1-0?1-0-1\t5
     *  where 0/1 stands for alleles, "-" stands for linked, and "?" stands for gap. The number
     *  after \t stands for the count
     *
     *  There are two criteria:
     *      (1) the gap in a single pool, largely due to the sequencing gap
     *      (2) the shared gap in all the pools, largely due to the long physical distance between
     *          two sites.
     *
     *  Parameters involved:
     *      double gap_all_pool_cutoff
     *      double gap_inpool_cutoff
     *
     *  @param graph_coloring_outcome
     *  @return
     */
    /**
    public int[] identify_gaps_Lauren_iterative_tobe_deleted_later(String[][] graph_coloring_outcome) {
        ArrayList<Integer> gap_indexes = new ArrayList<>();

        // Step 1) Count up the number of gaps within each and between all of the raw GC haplotypes
        // and all of the pools.
        // TODO: [ReconEP]:: separate helpers for each step?
        // The count at each potential gap in each pool. "this.num_sites-1" potential gaps.
        int[][] gap_counts_inpool = new int[this.num_pools][this.num_sites - 1];

        // The cumulative count at each potential gap position.
        int[] gap_counts_all = new int[this.num_sites - 1];

        // The number of types of raw GC haplotypes in each pool.
        double[] num_haps_inpool = new double[this.num_pools];
        double num_haps_all = 0.0;
        HashMap<String,Integer> hap_tracker = new HashMap<String,Integer>();
        int tot_ct = 0;
        for (int p = 0; p < this.num_pools; p++) { // for each pool...
            String[] haps = graph_coloring_outcome[p]; // the list of raw GC haplotypes
            num_haps_all += haps.length;
            num_haps_inpool[p] = haps.length;
            for (int h = 0; h < haps.length; h++){ // ...for each raw GC haplotype...
                String curr_vc = "";
                int hap_ct = Integer.parseInt(haps[h].split("\t")[1]);
                for(int k = 0; k < this.num_sites - 1; k++) { // ...for each potential gap...
                    curr_vc += haps[h].charAt(k * 2);

                    // If the linkage between the two variant positions is uncertain
                    // i.e.: gap is present...
                    if (haps[h].charAt(k * 2 + 1) == '?'){
                        gap_counts_inpool[p][k] += hap_ct; // increment the count within-pool and globally
                        gap_counts_all[k] += hap_ct;
                    }
                }

                curr_vc += haps[h].charAt((this.num_sites - 1) * 2);
                if (!hap_tracker.containsKey(curr_vc)) {
                    hap_tracker.put(curr_vc, hap_ct);

                } else {
                    int new_ct = hap_tracker.get(curr_vc) + hap_ct;
                    hap_tracker.put(curr_vc, new_ct);
                }

                tot_ct += hap_ct;
            }
        }

        // Step 2) Check if the potential gaps meet the within- OR between-pool frequency
        // thresholds.
        HashSet<Integer> gap_indexes_set = new HashSet<Integer>();
        double avg_region_size = this.num_sites;
        double across_cutoff = 1;
        double local_cutoff = 1;

        // While there aren't enough gaps to divide up the GC haplotypes...
        while (across_cutoff >= this.dp.gap_all_pool_cutoff
            || local_cutoff >= this.dp.gap_inpool_cutoff) {

            for (int k = 0; k <this.num_sites - 1; k++) {
                if((double) gap_counts_all[k] / num_haps_all >= across_cutoff
                    && !gap_indexes_set.contains(k)) {

                    gap_indexes.add(k);
                    gap_indexes_set.add(k);
                }

                for (int p = 0; p < this.num_pools; p++) {
                    if ((double) gap_counts_inpool[p][k] / num_haps_inpool[p] >= local_cutoff
                        && !gap_indexes_set.contains(k)){

                        gap_indexes.add(k);
                        gap_indexes_set.add(k);
                    }
                }
            }

            avg_region_size = (double) this.num_sites / (double) (gap_indexes.size() + 1);
            if (avg_region_size > this.dp.max_level_I_region_size) {
                break;
            }

            // Relax the cutoffs so that there many be more gaps.
            across_cutoff -= this.dp.gap_support_step;

            // The local cutoff is identical in case there are rare haplotypes with uncommon linkage
            // patterns.
            local_cutoff -= this.dp.gap_support_step;
        }

        System.out.println("The final across-pool cutoff is "
            + across_cutoff
            + " and the within-pool cutoff is "
            + local_cutoff
            + ". The average region size is "
            + avg_region_size
            + " segregating sites long.");

        int[] gap_indexes_array;

        // If there are enough gaps to run regional reconstruction....
        if (avg_region_size <= this.dp.max_level_I_region_size) {
            gap_indexes_array = new int[gap_indexes.size()];

            // Add the last site index as the final "gap" so that all the sites are within gaps.
            gap_indexes.add(this.num_sites - 1);

            // This is not useful for the moment, but keep the data integrity for potential future
            // use.
            gap_indexes_set.add(this.num_sites - 1);

            // Clean the outcome and return.
            for (int i = 0; i < gap_indexes.size(); i++) {
                gap_indexes_array[i] = gap_indexes.get(i);
            }
            Arrays.sort(gap_indexes_array); // sort the indices

            this.global_haps_gc = new String[hap_tracker.size()][this.num_sites];
            this.global_gc_freq = new double[hap_tracker.size()];
            this.num_haps_gc = hap_tracker.size();
            int hap_index = 0;
            for (String curr_vc : hap_tracker.keySet()) {
                String[] tmp = curr_vc.split("");
                for (int l = 0; l < this.num_sites; l++) {
                    this.global_haps_gc[hap_index][l] = tmp[l];
                }

                this.global_gc_freq[hap_index] = (double) hap_tracker.get(curr_vc)
                    / (double) tot_ct;

                hap_index++;
            }

            this.region_solve = true;

        } else { // ...otherwise, we can just reduce the candidates using LASSO
            gap_indexes_array = new int[1];
            gap_indexes_array[0] = -1;
            this.region_solve = false;
        }

        return gap_indexes_array;
    }
    */

    /*
     *  Identify gaps from the GC outcome.
     *  The GC outcome file is composed of inferred haplotypes in the format of:
     * 		0-1-0-1-0?1-0-1\t30
     *  where 0/1 stands for alleles, "-" stands for linked, and "?" stands for gap. The number
     *  after \t stands for the count of the hap
     *
     *  There are two criteria:
     * 	 (1) the gap in a single pool, largely due to the sequencing gap
     * 	 (2) the shared gap in all the pools, largely due to the long physical distance between two
     *       sites.
     *
     *  Parameters involved:
     * 	 double gap_all_pool_cutoff
     * 	 double gap_inpool_cutoff
     */
    public int[] identify_gaps(String[][] graph_coloring_outcome) {
        ArrayList<Integer> gap_indexes = new ArrayList<>();
        // Step 1) Count up the number of gaps within each and between all of the raw GC haplotypes
        // and all of the pools.
        // The count at each potential gap in each pool. "this.num_sites-1" potential gaps.
        int[][] gap_counts_inpool = new int[this.num_pools][this.num_sites-1];
        // The cumulative count at each potential gap position.
        int[] gap_counts_all=new int[this.num_sites-1];
        // The number of types of raw GC haplotypes in each pool.
        double[] num_haps_inpool=new double[this.num_pools];
  
        double num_haps_all = 0.0;
        HashMap<String,Integer> hap_tracker = new HashMap<String,Integer>();
        int tot_ct = 0;
        for (int p = 0; p < this.num_pools; p++) { // for each pool...
            String[] haps = graph_coloring_outcome[p]; // the list of raw GC haplotypes
            num_haps_all += haps.length;
            num_haps_inpool[p] = haps.length;
            for (int h = 0; h < haps.length; h++){ // ...for each raw GC haplotype...
                String curr_vc = ""; // the hap_string
                int hap_ct = Integer.parseInt(haps[h].split("\t")[1]);
                for(int k = 0; k < this.num_sites - 1; k++) { // ...for each potential gap...
                    curr_vc += haps[h].charAt(k * 2);
                    // If the linkage between the two variant positions is uncertain
                    // i.e.: gap is present...
                    if (haps[h].charAt(k * 2 + 1) == '?'){
                        gap_counts_inpool[p][k] += hap_ct; // increment the count within-pool and globally
                        gap_counts_all[k] += hap_ct;
                    }
                }
                curr_vc += haps[h].charAt((this.num_sites - 1) * 2);
                // put the haplotype "curr_vc" and its count to the table "hap_tracker"
                if (!hap_tracker.containsKey(curr_vc)) {
                    hap_tracker.put(curr_vc, hap_ct);
                } else {
                    int new_ct = hap_tracker.get(curr_vc) + hap_ct;
                    hap_tracker.put(curr_vc, new_ct);
                }
                tot_ct += hap_ct;
            }
        }
        //Step 2) Check if the potential gaps meet the within- OR between-pool frequency thresholds.
        HashSet<Integer> gap_indexes_set = new HashSet<Integer>();
        for (int k = 0; k < this.num_sites - 1; k++) {
            if ((double) gap_counts_all[k] / num_haps_all >= this.dp.gap_all_pool_cutoff
                && !gap_indexes_set.contains(k)) {

                gap_indexes.add(k);
                gap_indexes_set.add(k);
            }

            for (int p = 0; p < this.num_pools; p++) {
                if ((double) gap_counts_inpool[p][k] / num_haps_inpool[p] >= this.dp.gap_inpool_cutoff
                    && !gap_indexes_set.contains(k)) {

                    gap_indexes.add(k);
                    gap_indexes_set.add(k);
                }
            }
        }
        // Add the last site index as the final "gap" so that all the sites are within gaps.
        gap_indexes.add(this.num_sites - 1);
        // This is not useful for the moment, but keep the data integrity for potential future use.
        gap_indexes_set.add(this.num_sites - 1);

        // Clean the outcome and return.
        int[] gap_indexes_array=new int[gap_indexes.size()];
        for (int i = 0; i < gap_indexes.size(); i++) {
            gap_indexes_array[i] = gap_indexes.get(i);
        }
        Arrays.sort(gap_indexes_array); // sort the indices

        this.global_haps_gc = new String[hap_tracker.size()][this.num_sites];
        this.global_gc_freq = new double[hap_tracker.size()];
        this.num_haps_gc = hap_tracker.size();
        int hap_index = 0;
        for (String curr_vc : hap_tracker.keySet()) {
            String[] tmp = curr_vc.split("");
            for (int l = 0; l < this.num_sites; l++) {
                this.global_haps_gc[hap_index][l] = tmp[l];
            }
            this.global_gc_freq[hap_index] = ((double)hap_tracker.get(curr_vc))/ ((double)tot_ct);
            hap_index++;
        }
        return gap_indexes_array;
    }

    /**
     *  Output the current dividing plan to a file
     *
     *  @param dc_plan_outfile
     */
    public void output_current_DC_plan(String dc_plan_outfile) {
        try {
            BufferedWriter bw = new BufferedWriter(new FileWriter(dc_plan_outfile));

            // Output Level I:
            bw.write("Level_I\n");
            System.out.println("Level_I");
            for (int r = 0; r < this.num_regions_level_I; r++) {
                bw.write(this.regions_level_I[r][0] + ":" + this.regions_level_I[r][1] + "\t");
                System.out.print(this.regions_level_I[r][0]
                    + ":"
                    + this.regions_level_I[r][1]
                    + "\t");

            }

            // Output Level II:
            bw.write("\nLevel_II\n");
            System.out.println("\nLevel_II");
            for (int r = 0; r < this.num_regions_level_II; r++) {
                bw.write(this.regions_level_II[r][0] + ":" + this.regions_level_II[r][1] + "\t");
                System.out.print(this.regions_level_II[r][0]
                    + ":"
                    + this.regions_level_II[r][1]
                    + "\t");

            }

            bw.write("\n");
            bw.close();
            System.out.println();

        } catch(Exception e) {
            e.printStackTrace();
        }
    }


    /**
     *  Analyze a set of regions based on the dividing plan.
     *
     *  This function can be used for either level I or level II, or a subset of them
     *
     *  @param regions
     *  @param parameter_file
     *  @param dir_prefix
     *  @param level_index
     *  @return
     *  @throws Exception
     */
    public HapConfig[] regional_AEM(
    	String[] pool_IDs,
    	String[] vef_files,
        int[][] regions,
        String parameter_file,
        String dir_prefix,
        String aem_fail_lasso_path,
        int level_index) throws Exception {

        HapConfig[] region_haps = new HapConfig[regions.length];
        for (int r_index = 0; r_index < regions.length; r_index++) {
            System.out.print("AEM for level " + level_index + " region " + r_index + "... ");
            HapConfig hap_config = generate_hapconfig_2n(regions, r_index, pool_IDs);
            RegionEMSolver hap_solver = new RegionEMSolver(hap_config, parameter_file);
            if (hap_solver.failure) {
                System.out.println("AEM failed to converge. Initiating regional LASSO... ");

                // If AEM fails, at least some of the frequencies will be NaN. In that case, use
                // sub-optimal GC results. TODO [Quan] 
                region_haps[r_index] = generate_hapconfig_gc_regional_lasso(pool_IDs, vef_files,
                    regions[r_index], level_index, r_index, aem_fail_lasso_path);
            } else {
                region_haps[r_index] = hap_solver.final_Haps;
            }
            region_haps[r_index].recode_HapIDs_to_base16();
            region_haps[r_index].write_global_file_string(dir_prefix
                + "_level_"
                + level_index
                + "_region_"
                + r_index
                + ".inter_freq_haps.txt");

            System.out.print("Done. AEM ");
            if (hap_solver.failure) {
                System.out.print("failed to converge. ");
            } else {
                System.out.print("successfully converged. ");
            }

            System.out.println(region_haps[r_index].num_global_hap
                + " regional haplotypes have been generated.");
        }

        return region_haps;
    }


    /**
     *	generate the HapConfig with 2^n haplotypes (n is the number of sites).
     *
     *  @param the_region
     *  @param region_index
     *  @param pool_IDs
     *  @return
     *  @throws IOException
     */
    public HapConfig generate_hapconfig_2n(
        int[][] the_region, int region_index, String[] pool_IDs) throws IOException {

        int region_start = the_region[region_index][0];
        int region_end = the_region[region_index][1];
        int num_site_regional = region_end-region_start + 1;
        double[][] inpool_site_freqs = new double[num_site_regional][];
        LocusAnnotation[] locusInfo = new LocusAnnotation[num_site_regional];

        for (int l = 0; l < num_site_regional; l++) {
            locusInfo[l] = this.in_pool_sites_freq_anno.loci_annotations[l + region_start];
            inpool_site_freqs[l] = this.in_pool_sites_freq_anno.inpool_freqs[l + region_start];
        }
        int haps_2n = (int) Math.pow(2, num_site_regional);
        String[][] global_haps_string = new String[haps_2n][num_site_regional];
        String[] hap_IDs = new String[haps_2n];
        for (int h = 0; h < haps_2n; h++) {
            String curr_ID = "";
            String vc_str = Integer.toBinaryString(h);
            String[] vc_arr  = vc_str.split("");
            // the length of vc_str may not reach num_site_regional; so put zeros in.
            for (int locus = 0; locus < num_site_regional - vc_arr.length; locus++) {
                global_haps_string[h][locus] = "0";
                curr_ID += "0";
            }
            for (int locus = num_site_regional - vc_arr.length; locus<num_site_regional; locus++) {
                global_haps_string[h][locus] = vc_arr[locus - num_site_regional + vc_arr.length];
                curr_ID += vc_arr[locus - num_site_regional + vc_arr.length];
            }
            hap_IDs[h] = curr_ID;
        }

        double[] global_haps_freq = new double[haps_2n];
        Arrays.fill(global_haps_freq, 1.0 / (double) haps_2n);
        return new HapConfig(
            global_haps_string,
            global_haps_freq,
            null,
            inpool_site_freqs,
            locusInfo,
            this.num_pools,
            hap_IDs,
            pool_IDs,
            this.dp.est_ind_pool);
    }


    /**
     *
     *  @param region
     *  @param level
     *  @param r_index
     *  @param dir_prefix intermediate folder
     *  @return
     *  @throws FileNotFoundException
     */
    public HapConfig generate_hapconfig_gc_regional_lasso( 
        String[] pool_IDs,
        String[] vef_files,
        int[] region,
        int level,
        int r_index,
        String aem_fail_lasso_path)  // i 
            throws FileNotFoundException {

        new File(aem_fail_lasso_path).mkdir();
        int start = region[0]; 
        int end = region[1];
        int num_site_regional = end - start + 1;
        // collect all regional haps from raw GC output. 
        HashMap<String, Double> hap_tracker = new HashMap<String, Double>();
        for (int h = 0; h < this.num_haps_gc; h++) {
            String curr_vc = String.join("",Arrays.copyOfRange(this.global_haps_gc[h],start,end+1));
            if (!hap_tracker.containsKey(curr_vc)) {
                hap_tracker.put(curr_vc, this.global_gc_freq[h]); 
            } else {
                double new_freq = hap_tracker.get(curr_vc) + this.global_gc_freq[h];
                hap_tracker.put(curr_vc, new_freq);
            }
        }
        int num_hap_regional = hap_tracker.size();
        String[][] reg_haps_string = new String[num_hap_regional][num_site_regional];
        int hap_index = 0;
        String[] hap_IDs = new String[num_hap_regional];
        double[] global_haps_freq = new double[num_hap_regional];
        for (String curr_vc : hap_tracker.keySet()) {
            String[] tmp = curr_vc.split("");
            for (int s = 0; s < num_site_regional; s++) reg_haps_string[hap_index][s] = tmp[s];
            global_haps_freq[hap_index] = hap_tracker.get(curr_vc) / num_pools;
            hap_IDs[hap_index] = "h"+Integer.toString(hap_index);
            hap_index++;
        }

        double[][] inpool_site_freqs = new double[num_site_regional][];
        LocusAnnotation[] locusInfo = new LocusAnnotation[num_site_regional];
        for (int s = 0; s < num_site_regional; s++) {
            locusInfo[s] = this.in_pool_sites_freq_anno.loci_annotations[s + start];
            inpool_site_freqs[s] = this.in_pool_sites_freq_anno.inpool_freqs[s + start];
        }

        HapConfig gc_regional_haps = new HapConfig(
            reg_haps_string,
            global_haps_freq,
            null,
            inpool_site_freqs,
            locusInfo,
            this.num_pools,
            hap_IDs,
            pool_IDs,
            this.dp.est_ind_pool);

        HapLASSO regional_lasso = new HapLASSO( 
            -1,    // Specified by pool_index==-1, this is the *multi_pool* implementation of LASSO.
            null,  // pool_IDs is null as it is the multi_pool LASSO.
            dp.lasso_regional_lambda,
            gc_regional_haps,
            0,
            dp.lasso_regional_memory,
            aem_fail_lasso_path + "/"+dp.project_name+"_level_" + level + "_region_" + r_index + ".lasso_in");

        regional_lasso.estimate_frequencies_lasso(null, vef_files, dp.lasso_weights);

        HapConfig final_reg_haps = regional_lasso.hapOut(pool_IDs);

        // filter out low-frequent haps 
        boolean[] list_remove_haps = new boolean[final_reg_haps.num_global_hap];      
        int num_remove_hap = 0;
        for (int h = 0; h < final_reg_haps.num_global_hap; h++) {
            if (final_reg_haps.global_haps_freq[h] < dp.lasso_regional_cross_pool_cutoff) {
                list_remove_haps[h] = true;
                num_remove_hap++;
            }
        }
        double actual_cutoff = dp.lasso_regional_cross_pool_cutoff;
        
        // If too many of them are below the regional frequency minimum,
        // i.e., too few haps are remained 
        if (dp.lasso_hapset_size_min > final_reg_haps.num_global_hap - num_remove_hap) {
            if(final_reg_haps.num_global_hap >= dp.lasso_hapset_size_min) {
                list_remove_haps = Algebra.permute_sort_and_remove(
                    final_reg_haps.global_haps_freq.clone(), 
                    dp.lasso_hapset_size_min);
            }            
        }
        // Or, if too many haps are remained...
        else if (dp.lasso_hapset_size_max < final_reg_haps.num_global_hap - num_remove_hap) {
            list_remove_haps = Algebra.permute_sort_and_remove(
                final_reg_haps.global_haps_freq.clone(), 
                dp.lasso_hapset_size_max);
        }
        num_remove_hap = 0;
        actual_cutoff=1;
        for(int h=0; h < final_reg_haps.num_global_hap; h++) {
            if(list_remove_haps[h]) {
                num_remove_hap++;
            } else {
                //this.initial_Haps.global_haps_freq[h] = freq[h];
                if(actual_cutoff > final_reg_haps.global_haps_freq[h]) { 
                    actual_cutoff = final_reg_haps.global_haps_freq[h];
                }
            }
        }     
        final_reg_haps.remHaps(list_remove_haps, num_remove_hap);
        System.out.println("Of the "
            + final_reg_haps.num_global_hap
            + " regional haplotypes, "
            + num_remove_hap
            + " were removed. "
            + "The frequency cutoff was " + actual_cutoff + ".");        
//        if (Double.isNaN(actual_cutoff)) {
//            double[] tmp = new double[final_reg_haps.num_global_hap];
//            for (int h = 0; h < final_reg_haps.num_global_hap; h++) {
//                tmp[h] = 1.0 / ((double) final_reg_haps.num_global_hap);
//            }
//            final_reg_haps.global_haps_freq = tmp;
//        }
        return final_reg_haps;
    }
}
