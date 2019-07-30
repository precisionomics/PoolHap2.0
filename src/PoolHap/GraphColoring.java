package PoolHap;

import java.io.BufferedReader;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Iterator;
import java.util.Vector;
import java.util.Random;

public class GraphColoring {
    // TODO: [ReconEP]:: refactor variable names to be more sensible (e.g. camelCase for Java).
    // TODO: [Question]:: why are these Vectors instead of ArrayLists in the first place?

    // The index of the genotype when it was added to readinfo_arr_tmp.
    public Vector<Integer> readindex_arr_tmp;

    // The list of all possible genotypes, ArrayList<site=allele; ... site=allele;>
    public Vector<String> readinfo_arr_tmp;

    public HashMap<String, Integer> output_ref_arr;
    public HashMap<String, String> conf_ref_arr;
    public int num_loci;
    public int num_pools;
    public LocusAnnotation[] locusInfo; // # of loci, note that a locus can be a SNP or a region
    public double[][] inpool_site_freqs; // # of loci x # of pools; Added by Quan Dec. 2018.
    public int num_loci_window;
    public int max_num_gap;
    public String vef_file;

    // TODO: [LEFTOVER]
    // new GraphColoring(gp.inter_dir + prefix + "_p" + p + ".vef",
    //     gs_var_pos,
    //     gp.inter_dir + prefix + "_p" + p + ".in");

    /**
     *  Graph coloring object constructor with input files. For first round of graph coloring.
     *  TODO: [Javadoc]:: improve description of the constructor for GC round 1.
     *
     *  @param vef variant encoded file file path string.
     *  @param gs_var_pos (required) gold standard variant positions file path string.
     *  @param out_file output file path string.
     *  @throws IOException on input error.
     */
    public GraphColoring(String vef, String gs_var_pos, String out_file, int num_pos_window, 
    		int num_gap_window) throws IOException {
        /*
         *  Initialize read indices, read information, genotype dictionary, and reader variables.
         */
        // TODO: [Question]:: why are the variables named "arr" when they're Vector or HashMap
        // objects?
        this.readindex_arr_tmp = new Vector<Integer>();
        this.readinfo_arr_tmp = new Vector<String>();
        this.num_loci_window = num_pos_window;
        this.max_num_gap = num_gap_window;
        this.vef_file = vef;
        
        int count = 0;

        // HashMap<pos=allele;..., count>
        HashMap<String, Integer> geno_dict = new HashMap<String, Integer>();

        // I/O for VEF file.
        BufferedReader bufferedreader = new BufferedReader(new FileReader(vef));
        String line = "";

        // The maximum number of times a genotype can be counted in a single pool.
        // TODO: (old) [Question]::  Why does this exist?
        int max_num_geno = 32878;


        /*
         *  Read through VEF file.
         */
        while ((line = bufferedreader.readLine()) != null) {
            line = line.replace("\r", ""); // remove newline characters
            String[] line_arr = line.split("\t"); // Read_Name\tPos=Allele;...\t//\tStart\tEnd

            // TODO: [ReconEP]:: extract the if-else logic below into a separate method.
            // If the read contains a segregating site (i.e.: has a distinguishing genotype)...
            if (line_arr[1].contains("=")) {
                String tmp_geno = line_arr[1]; // the segregating site information column

                // If this combination of alleles hasn't been recorded yet, add to dictionary.
                if (!geno_dict.containsKey(tmp_geno)) {
                    geno_dict.put(tmp_geno, 1);

                    // The index of readinfo_arr_tmp that corresponds to this genotype.
                    this.readindex_arr_tmp.add(count);
                    count = count + 1; // TODO: (minor) [Question]:: change to += for consistency?
                    this.readinfo_arr_tmp.add(line_arr[1]);

                // Else, it has already been recorded, remove and re-add with count + 1.
                } else {
                    int tmp_num = geno_dict.get(tmp_geno); // current count of site combination

                    // TODO: [Question]:: do we need to remove and re-add? Wouldn't the following
                    // auto-update and eliminate most of the if-else outside of guarding against the
                    // maximum?
                    // map.put(key, map.getOrDefault(key, 0) + 1)
                    geno_dict.remove(tmp_geno);
                    geno_dict.put(tmp_geno, (tmp_num + 1));

                    // If site combination count is less than the stated max, re-add to dictionary.
                    if (geno_dict.get(tmp_geno) < max_num_geno) {
                        this.readindex_arr_tmp.add(count);
                        count = count + 1;
                        this.readinfo_arr_tmp.add(line_arr[1]);
                    }
                }
            }
        }
        bufferedreader.close();
        /*
         *  Solve and produce haplotype configurations.
         */
        this.gc_solver(gs_var_pos);
        this.fileOut(out_file);
    }


    /**
     *  2nd round of graph coloring for global haplotype through region linking.
     *  TODO: [Javadoc]:: improve description of this 2nd GC round constructor.
     *
     *  @param level_1 all level 1 haplotype configurations
     *  @param level_2 all level 2 haplotype configurations
     *  @param gs_var_pos (required) gold standard variant positions file path string
     *  @param virtual_cov_link_gc // See property file for the explanation of this parameter
     *  @throws IOException on input error.
     */
    public GraphColoring(HapConfig[] level_1,  HapConfig[] level_2,  String gs_var_pos,
        int virtual_cov_link_gc) throws IOException {
        // Initialize read indices and read information variables.
        this.readindex_arr_tmp = new Vector<Integer>();
        this.readinfo_arr_tmp = new Vector<String>();
        int count = 0;
        // Covert regional haplotypes to VEF format.
        // For each region in all level 1 regions...
        for (int r = 0; r < level_1.length; r++) {
            // For each haplotype in all global haplotypes in region...
            for (int h = 0; h < level_1[r].num_global_hap; h++) {
                String curr_vc = ""; // TODO: [Question]:: what is vc?
                // For each locus in all loci in region
                for (int l = 0; l < level_1[r].num_loci; l++) {
                    curr_vc += (level_1[r].locusInfo[l].start_loc // format read info
                        + "="
                        + level_1[r].global_haps_string[h][l]
                        + ";");
                }
                // TODO: [Question]:: what is ct? Count?
                int hap_ct = (int) (level_1[r].global_haps_freq[h] * virtual_cov_link_gc); // haplotype count(?)
                for (int c = 0; c < hap_ct; c++) { // add reads and corresponding indices to vectors
                    // The index of readinfo_arr_tmp that corresponds to this genotype.
                    this.readindex_arr_tmp.add(count);  // the number of line in the file
                    count++;
                    this.readinfo_arr_tmp.add(curr_vc);
                }
            }
        }
        // For each region in all level 2 regions...
        for (int r = 0; r < level_2.length; r++) {
            // For each haplotype in all global haplotypes in region...
            for (int h = 0; h < level_2[r].num_global_hap; h++) {
                String curr_vc = "";
                // For each locus in all loci in region...
                for (int l = 0; l < level_2[r].num_loci; l++) {
                    curr_vc += (level_2[r].locusInfo[l].start_loc // format read info
                    + "="
                    + level_2[r].global_haps_string[h][l]
                    + ";");
                }
                int hap_ct = (int) (level_2[r].global_haps_freq[h] * virtual_cov_link_gc); // haplotype count(?)
                for (int c = 0; c < hap_ct; c++) { // add reads and corresponding indices to vectors
                    // The index of readinfo_arr_tmp that corresponds to this genotype.
                    this.readindex_arr_tmp.add(count);
                    count++;
                    this.readinfo_arr_tmp.add(curr_vc);
                }
            }
        }
        System.out.println("There are " + count
            + " individual fragments (reads or regional haplotypes) in the dataset.");

        this.gc_solver(gs_var_pos);
    }


    /**
     *  TODO: [Question]:: what does this mean?
     *  Source_Path or source_path
     *
     *  Connects variants with graph colouring.
     *  TODO: [Javadoc]:: improve description of solver method.
     *  TODO: [ReconEP]:: this method is way too long, simply not maintainable or testable, refactor
     *  into multiple smaller helper methods.
     *
     *  @param gs_var_pos (required) gold standard variant positions file path string.
     *  @throws IOException on input error.
     */
    public void gc_solver(String gs_var_pos) throws IOException {
        /*
         *  Initialize variables.
         */
        // HashMap<Seg_Site, Index>
        HashMap<Integer, Integer> pos_dict = new HashMap<Integer, Integer>();
        Integer pos_index = 0;

        /*
         *  Read through gold standard variant positions file to load into dictionary and count
         *  total number of loci.
         *
         *  also this needs to be done in SiteInPoolFreqAnno as well, refactor to a static method?
         *  retain number of loci as a global variable?
         */
        // Read into file.
        BufferedReader br = new BufferedReader(new FileReader(gs_var_pos));
        String currLine = br.readLine(); // skip header
        currLine = br.readLine();

        // Read each line into a dictionary.
        while (currLine != null) {
            pos_dict.put(Integer.parseInt(currLine.split(";")[1]), pos_index);

            // TODO: [LEFTOVER]
            // System.out.println(currLine.split(";")[1] + "\t" +  pos_index);

            pos_index++; // move to next variant position index
            currLine = br.readLine(); // read next line
        }

        br.close();
        this.num_loci = pos_dict.size(); // number of loci in gold standard variants file


        /*
         *  Read through gold standard variant positions file again to load loci info into a matrix.
         */
        // Read into file.
        br = new BufferedReader(new FileReader(gs_var_pos));
        currLine = br.readLine(); // skip header
        currLine = br.readLine();

        // Initialize variables.
        int loci_index = 0;
        this.num_pools = currLine.split("\t").length - 1;
        this.locusInfo = new LocusAnnotation[this.num_loci];
        this.inpool_site_freqs = new double[this.num_loci][this.num_pools];

        // Read through file.
        while (currLine != null) {
            String[] tmp = currLine.split("\t");

            // The first column is the locus-info.
            this.locusInfo[loci_index] = new LocusAnnotation(tmp[0]);

            // For each pool, load locus frequency into matrix.
            for (int p = 0; p < this.num_pools; p++) {
                this.inpool_site_freqs[loci_index][p] = Double.parseDouble(tmp[p + 1]);
            }

            loci_index++; // move to next locus index
            currLine = br.readLine(); // read next line
        }

        br.close();

        // TODO: [LEFTOVER]
        // System.out.println("There are " + this.locusInfo.length + " positions.");

        /*
         *  Randomize
         *  TODO: [IMPORTANT]:: this entire next section makes no sense to me, it feels like there
         *  are a lot of inefficiencies and overhead that could be avoided by not using Vectors
         *  in the first place? We could literally replace the next 30 lines with like 4 lines of
         *  code if some variables were initialized as ArrayLists.
         */
        // ArrayList<pos=allele;> in random order
        Vector<Integer> readindex_arr = new Vector<Integer>();

        // ArrayList<original index of pos=allele;> TODO: (old) [Review] Confirm!
        Vector<String> readinfo_arr = new Vector<String>();

        // TODO: [Question]:: why not .toArray()?
        // Create new array from read index Vector.
        int[] index_arr_tmp = new int[this.readindex_arr_tmp.size()];
        for (int k = 0; k < index_arr_tmp.length; k++) {
            index_arr_tmp[k] = k;
        }

        // TODO: [Question]:: why not convert from array to arraylist?
        // TODO: [ReconEP]:: rename variables to something that intuitively makes sense...
        // e.g. ArrayList<Integer> list = ArrayList<Integer>(Arrays.asList(index_arr_tmp));
        // Create new ArrayList from read index array.
        ArrayList<Integer> list = new ArrayList<Integer>();
        for (int i = 0; i < index_arr_tmp.length; i++) {
            list.add(index_arr_tmp[i]);
        }

        // A list of the indices corresponding to each read, in random order.
        // Create new array the length of the read index Vector.
        int[] index_arr = new int[this.readindex_arr_tmp.size()];

        // index_arr_tmp, list, and index_arr are copies of this.readindex_arr_tmp so far.
        // Shuffle read index ArrayList.
        Collections.shuffle(list);

        // Now, index_arr_tmp, index_arr are copies of this.readindex_arr_tmp. List is now a
        // randomized version of this.readindex_arr_tmp.
        // TODO: [Question]:: why can't we clone the the shuffled list as index_arr? Is iterating
        // through it faster?
        // Honestly
        // e.g. int[] index_arr = list.toArray(new int[list.size()]);
        // Iterate through randomized read index ArrayList to create new read index array that is
        // randomized.
        Iterator<Integer> ite = list.iterator();
        int tmp_i = 0;
        while (ite.hasNext()) {

            // TODO: [LEFTOVER]
            // System.out.println(ite.next().toString()+", ");

            index_arr[tmp_i] = ite.next();
            tmp_i++;
        }

        // Now, index_arr_tmp is a copy of this.readindex_arr_tmp. List, index_arr are now
        // randomized versions of this.readindex_arr_tmp.

        // TODO: [LEFTOVER]
        // for (int i = 0; i < index_arr.length; i++) {
        //     System.out.println(index_arr[i]);
        // }
        // index_arr is a randomized list of indices corresponding to the genotypes from the VEF
        // file.

        // Iteratively add the randomized read index and non-randomized read information arrays
        // to the previously defined new Vectors.
        // TODO: [Question]:: I still don't understand why we need Vectors.
        for (int i = 0; i < this.readindex_arr_tmp.size(); i++) {
            // TODO: [Question]:: why not this.readinfo_arr_tmp?
            readinfo_arr.add(readinfo_arr_tmp.get(index_arr[i]));
            readindex_arr.add(this.readindex_arr_tmp.get(index_arr[i]));
        }
        // readinfo_arr contains the randomized order list of the genotypes from the VEF file.
        // readindex_arr contains the original indices corresponding to the genotypes from the VEF
        // file i.e.: the order they were read in.

        // TODO: [LEFTOVER]
        // for (int i = 0; i < readindex_arr.size(); i++) {
        // 	System.out.println(readinfo_arr.get(i));
        // }


        /*
         *  Set up matrices?
         *  TODO: [Question]:: clarify/understand this section better.
         */
        // Genotype index: ArrayList<ArrayList<Index_SS>>
        ArrayList<ArrayList<Integer>> read_pos_2D_arr= new ArrayList<ArrayList<Integer>>();

        // Genotype allele: ArrayList<ArrayList<Allele_SS>>
        ArrayList<ArrayList<String>> read_geno_2D_arr= new ArrayList<ArrayList<String>>();

        // TODO: [LEFTOVER]
        // System.out.println(readinfo_arr);

        // For each read in the read info array....
        for (int i = 0; i < readinfo_arr.size(); i++) {
            String tmp_str = readinfo_arr.get(i);
            String [] tmp_str_arr= tmp_str.split(";"); // split loci in the read
            ArrayList<String> i_geno_arr = new ArrayList<String>();
            ArrayList<Integer> i_pos_arr = new ArrayList<Integer>();

            // For each locus in the read...
            for (int k = 0; k < tmp_str_arr.length; k++) {
                String [] tmp2_arr = tmp_str_arr [k].split("=");
                int pos = Integer.parseInt(tmp2_arr[0]);
                i_pos_arr.add(pos_dict.get(pos)); // add position index to position row array
                i_geno_arr.add(tmp2_arr[1]); // add genotype to genotype row array
            }

            read_pos_2D_arr.add(i_pos_arr); // add row to read variant positions matrix
            read_geno_2D_arr.add(i_geno_arr); // add row to read genotype matrix

        }

        // TODO: [LEFTOVER]
        // System.out.println(read_pos_2D_arr);
        // System.out.println(read_geno_2D_arr);

        // ArrayList<Genotype, ArrayList<Index_Genotype>>
        ArrayList<ArrayList<Integer>> readadj_arr = new ArrayList<ArrayList<Integer>>();

        // For each possible genotype...
        for (int i = 0; i < read_pos_2D_arr.size(); i++) {

            // TODO: [LEFTOVER]
            // if (i % 500 == 0) {
            //     System.out.println(i + " fragments have been processed.");
            // }

            // TODO: [Question]:: what's this, adjacent?
            ArrayList<Integer> adj_arr = new ArrayList<Integer>();

            // ...comparing it to all other genotypes...
            for (int j = 0; j < read_pos_2D_arr.size(); j++) {
                if (i != j) {
                    boolean IsConnect = false;
                    for (int k = 0; k < read_pos_2D_arr.get(i).size(); k++) {
                        for (int l = 0; l < read_pos_2D_arr.get(j).size(); l++) {
                            // If pos(geno_i, index_k) == pos(geno_j, index_l) and the alleles are
                            // the same, the two reads can connect.
                            if ((read_pos_2D_arr.get(i).get(k) == read_pos_2D_arr.get(j).get(l))
                                && (!read_geno_2D_arr.get(i)
                                    .get(k)
                                    .equals(read_geno_2D_arr.get(j).get(l)))) {

                                IsConnect = true;
                            }
                        }
                    }

                    if (IsConnect){
                        // Add index of geno_j to the list of possible connects for geno_i.
                        adj_arr.add(readindex_arr.get(j));
                    }
                }
            }

            readadj_arr.add(adj_arr);

        }

        // TODO: [LEFTOVER]
        // System.out.println(
        //     "Finished identifying all possible conflicts between genotype fragments.");
        //
        // System.out.println(readadj_arr);

        int max_color = readindex_arr.size(); // only if there are all 2^loci genotypes present
        ArrayList<Integer> read_color_arr = new ArrayList<Integer>();
        ArrayList<HashSet<Integer>> nb_color_arr = new ArrayList<HashSet<Integer>>();

        // TODO: [LEFTOVER]
        // nb_color_arr.add(new HashSet());

        for (int i = 0; i < readadj_arr.size(); i++) {
            nb_color_arr.add(new HashSet<Integer>());
            read_color_arr.add(-1);
        }

        read_color_arr.set(0, 0); // set the first colour as the genotype at index 0

        for (int i = 0; i < readadj_arr.get(0).size(); i++) {
            int index = readadj_arr.get(0).get(i);

            // Add colour 0 as a possible colour to all genotypes that can connect with the genotype
            // index 0.
            nb_color_arr.get(index).add(0);
        }

        int real_max_color = 0;
        ArrayList<HashSet<String>> color_geno_set_arr = new ArrayList<HashSet<String>>();
        color_geno_set_arr.add(new HashSet<String>());

        // This is the full genotype of a colour.
        color_geno_set_arr.get(0).add(readinfo_arr.get(0));

        // TODO: [LEFTOVER]
        // System.out.println(color_geno_set_arr);

        // Make the conflict graph. Basically, assign colours to any genotype fragments that can't
        // go together.
        while (true) {
            int max_nb_color = -1;
            int index = -1;

            for (int i = 0; i < readadj_arr.size(); i++) { // for each genotype...
                // If there hasn't been a colour assigned to that genotype...
                if ((read_color_arr.get(i) == -1)
                    // ...and there are other genotypes that can go with it of the same colour...
                    && (nb_color_arr.get(i).size() > max_nb_color)) {

                    index = i;
                    max_nb_color = nb_color_arr.get(i).size();
                }
            }

            if (index == -1) { // If there are no colours left (?) end this step.
                break;
            }
            int color = -1;

            for (int i = 0; i < max_color; i++) { // for each genotype...
                // If the possible colours list for that genotype doesn't contain colour i...
                if (!nb_color_arr.get(index).contains(i)) {

                    // If colour i is smaller than the number of available colours.
                    if (i < color_geno_set_arr.size()) {
                        if (!color_geno_set_arr.get(i).contains(readinfo_arr.get(index))) { // and
                            color = i;
                            color_geno_set_arr.get(i).add(readinfo_arr.get(index));
                            break;
                        }
                    } else {
                        color_geno_set_arr.add(new HashSet<String>());
                        color_geno_set_arr.get(color_geno_set_arr.size() - 1)
                            .add(readinfo_arr.get(index));

                        color = i;
                        break;
                    }
                }
            }

            if (color > real_max_color) {
                real_max_color = color;
            }

            read_color_arr.set(index, color);
            for (int i = 0; i < readadj_arr.get(index).size(); i++) {
                nb_color_arr.get(readadj_arr.get(index).get(i)).add(color);
            }

            // TODO: [LEFTOVER]
            // System.out.println(max_nb_color);

        }

        // TODO: [LEFTOVER]
        // System.out.println(
        //     "Finished identifying potential full-genome genotypes i.e.: haplotypes.");
        //
        // System.out.println(read_color_arr);

        //
        String null_ref= "*";
        String conf_ref= "";
        for (int i = 1; i < num_loci; i++ ) {
            null_ref = null_ref + "*";
            conf_ref = conf_ref + "?";
        }

        //
        String[] ref_arr = new String[real_max_color+1];
        String[] conf_arr = new String[real_max_color+1];
        for (int i = 0; i <= real_max_color; i++) {
            ref_arr[i] = null_ref;
            conf_arr[i] = conf_ref;
        }

        // TODO: [LEFTOVER]
        // String ss = "1234567";
        // System.out.println(ss.substring(1, ss.length()));

        //
        for (int i = 0; i < read_color_arr.size(); i++) {
            int i_color = read_color_arr.get(i);
            for (int j = 0; j < read_pos_2D_arr.get(i).size(); j++) {

                // TODO: [LEFTOVER]
                   // System.out.println(read_pos_2D_arr.get(i).get(j).toString());

                int p = read_pos_2D_arr.get(i).get(j);
                p++;
                ref_arr[i_color] = ref_arr[i_color].substring(0, (p - 1))
                    + read_geno_2D_arr.get(i).get(j).toString()
                    + ref_arr[i_color].substring(p, ref_arr[i_color].length());
            }

            //
            if (read_pos_2D_arr.get(i).size() > 1) {
                for (int j = 0; j < (read_pos_2D_arr.get(i).size() - 1); j++) {
                    if (Math.abs(read_pos_2D_arr.get(i).get(j + 1) - read_pos_2D_arr.get(i).get(j))
                        == 1) {

                        int p = read_pos_2D_arr.get(i).get(j);
                        p++;
                        try {
                            conf_arr[i_color] = conf_arr[i_color].substring(0, (p - 1))
                                + "-"
                                + conf_arr[i_color].substring(p, conf_arr[i_color].length());

                        } catch (StringIndexOutOfBoundsException e) {
                            System.out.println(conf_arr[i_color].substring(0,(p - 1))
                                + "\t"
                                + p
                                + "\t"
                                + conf_arr[i_color].length());

                        }
                    }
                }
            }

        }

        // TODO: [LEFTOVER]
        // for (int i = 0; i < conf_arr.length; i++) {
        //     System.out.println(ref_arr[i]);
        // }

        this.output_ref_arr = new HashMap<String,Integer>();
        this.conf_ref_arr = new HashMap<String,String>();
        
        int num_window =  this.num_loci/ this.num_loci_window +1;
        if ((this.num_loci % this.num_loci_window )==0) {
        	num_window--;
        }
        ArrayList<ArrayList<String>> ref_reg_2D_arr= new ArrayList<ArrayList<String>>();
        ArrayList<ArrayList<String>> conf_reg_2D_arr= new ArrayList<ArrayList<String>>();
        int completeness_cutoff = this.max_num_gap;
        if (num_window >1 ) {
        	for (int i=0; i <(num_window); i++ ) {
        		ArrayList<String> tmp_ref_arr = new ArrayList<String>();
        		ArrayList<String> tmp_conf_arr = new ArrayList<String>();
        		int index_max = (i+1)* this.num_loci_window;
        		if (((i+1)* this.num_loci_window) > this.num_loci) {
        			index_max = this.num_loci;
        		}
        		
        		for (int j = 0; j <= real_max_color; j++) {
        			if (count(ref_arr[j].substring(this.num_loci_window*i ,index_max )  ) <= completeness_cutoff) {
        				tmp_ref_arr.add( ref_arr[j].substring(this.num_loci_window*i ,index_max )  );
        				tmp_conf_arr.add( conf_arr[j].substring(this.num_loci_window*i,index_max-1 )  );
        			}
        		}
        		ref_reg_2D_arr.add(tmp_ref_arr);
        		conf_reg_2D_arr.add(tmp_conf_arr);
        	}
        	
        	int max_size =0;
        	for (int i=0; i <ref_reg_2D_arr.size(); i++ ) {
        		if (ref_reg_2D_arr.get(i).size()>  max_size) {
        			max_size = ref_reg_2D_arr.get(i).size();
        		}
        	}
        	
        	for (int i=0; i <ref_reg_2D_arr.size(); i++ ) {
        		for (int j=ref_reg_2D_arr.get(i).size();j< max_size;j++) {
        			int max_index = ref_reg_2D_arr.get(i).size();
        			if (max_index>0 ) {
		        	    Random random = new Random();
		        	    int s = random.nextInt(max_index)%(max_index+1);
		        	    ref_reg_2D_arr.get(i).add(ref_reg_2D_arr.get(i).get(s));
		        	    conf_reg_2D_arr.get(i).add(conf_reg_2D_arr.get(i).get(s));
        			}else {
        				System.out.println("ERROR: Can not reconstruct haplotype for" + this.vef_file+ "using graph "
        						+ "coloring.\nPlease decrease the Num_Pos_Window or increase Num_Gap_Window!");
        	        	System.exit(0);
        			}
        			
        		}
        	}
        	
//        	for (int i =0; i< ref_reg_2D_arr.size();i++)
//        		System.out.println(ref_reg_2D_arr.get(i));
        	
        	for (int j=0; j< ref_reg_2D_arr.get(0).size(); j++) {
        		String tmp_ref = ref_reg_2D_arr.get(0).get(j); 
        		String tmp_conf =conf_reg_2D_arr.get(0).get(j); 
        		for (int i=1; i <ref_reg_2D_arr.size(); i++ ) {
        			tmp_ref=tmp_ref + ref_reg_2D_arr.get(i).get(j);
        			tmp_conf =tmp_conf+ "?" + conf_reg_2D_arr.get(i).get(j); 
        		}
        		if (this.output_ref_arr.containsKey(tmp_ref)) {
                    this.output_ref_arr.put(tmp_ref, this.output_ref_arr.get(tmp_ref)+1);
                    conf_ref_arr.put(tmp_ref, tmp_conf);

                //
                } else {
                    this.output_ref_arr.put(tmp_ref, 1);
                    this.conf_ref_arr.put(tmp_ref, tmp_conf);
                }
        	}
        	
        	
        }else {
        
	        for (int i = 0; i <= real_max_color; i++) {
	            if (count(ref_arr[i]) <= completeness_cutoff) {	// may implement this in the future
	
	                // TODO: [LEFTOVER]
	                // System.out.println(ref_arr[i]);
	
	                if (this.output_ref_arr.containsKey(ref_arr[i])) {
	                    this.output_ref_arr.put(ref_arr[i], this.output_ref_arr.get(ref_arr[i])+1);
	                    conf_ref_arr.put(ref_arr[i], conf_arr[i]);
	
	                //
	                } else {
	                    this.output_ref_arr.put(ref_arr[i], 1);
	                    this.conf_ref_arr.put(ref_arr[i], conf_arr[i]);
	                }
	
	            } else {
	
	                // TODO: [LEFTOVER]
	                // System.out.println(ref_arr[i] + "\t" + count(ref_arr[i]));
	
	            }
	        }
        }

        // TODO: [LEFTOVER]
        // if (this.conf_ref_arr.isEmpty()) {
        //
        // }

    }

    public void fileOut(String out_file) throws IOException {
        FileWriter mydata = new FileWriter(out_file,false);
        PrintWriter pw = new PrintWriter(mydata);

        // TODO: [LEFTOVER]
        // System.out.println(output_ref_arr );

        // Changed iterator from Map<K, V> -> K because just getting the key requires less memory.
        for (String entry : this.output_ref_arr.keySet()) {
            String b = this.conf_ref_arr.get(entry);
            String c = "";
            String tmp_str= "";
           
            for (int i = 0; i < b.length(); i++) {
            	tmp_str = entry.substring(i, i + 1);
            	if (tmp_str.equals("*")) {
            		tmp_str="0";
            	}
                c = c + tmp_str + b.substring(i, i + 1); // TODO: (minor) c+=?
            }
            tmp_str =  entry.substring(entry.length() - 1, entry.length()) ;
        	if (tmp_str.equals("*")) {
        		tmp_str="0";
        	}
            c = c+ tmp_str ; // TODO: (minor) c+=?
            pw.write(c + "\t" + this.output_ref_arr.get(entry).toString() + "\n");

            // TODO: [LEFTOVER]
            // System.out.println(c + "\t" + output_ref_arr.get(x).toString());

        }

        pw.flush();
        pw.close();

        // TODO: [LEFTOVER]
        // for (int i = 0; i< real_max_color + 1; i++) {
        //     System.out.println(ref_arr[i]);
        // }

        return;
    }

    /**
     *
     * @return
     */
    public HapConfig hapOut(String[] pool_IDs) {
        //
        int num_global_hap = this.output_ref_arr.size();
        String[][] global_haps_string = new String[num_global_hap][num_loci];
        int[] global_haps_ct = new int[num_global_hap];
        int tot_hap_ct = 0;
        int hap_index = 0;

        //
        for (String entry : this.output_ref_arr.keySet()) {
            String[] var_comp = entry.split("");
            for (int v = 0; v < num_loci; v++) {
                global_haps_string[hap_index][v] = var_comp[v];
            }

            global_haps_ct[hap_index] = this.output_ref_arr.get(entry);
            tot_hap_ct += this.output_ref_arr.get(entry);
            hap_index++;
        }

        //
        double[] global_haps_freq = new double[num_global_hap];
        System.out.println("There are " + num_global_hap + " haplotypes generated by GC.");
        for (int h = 0; h < num_global_hap; h++) {
            global_haps_freq[h] = (double) global_haps_ct[h] / (double) tot_hap_ct;
        }

        return new HapConfig(
            global_haps_string,
            global_haps_freq,
            null,
            this.inpool_site_freqs,
            this.locusInfo,
            this.num_pools,
            null,
            pool_IDs,
            0);

    }


    /**
     *
     * @param gc_hap
     * @return
     */
    int count(String gc_hap) {
        int unknown = 0;
        String[] allele_comp = gc_hap.split("");
        for (String a : allele_comp) {
            if (a.equals("*")) {
                unknown++;
            }
        }

        return unknown;
    }
}
