package PoolHap;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.FileReader;
import java.io.FileWriter;
import java.math.RoundingMode;
import java.text.DecimalFormat;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.concurrent.ThreadLocalRandom;

import org.apache.commons.math3.distribution.BetaDistribution;
import org.apache.commons.math3.exception.NumberIsTooSmallException;
import org.apache.commons.math3.linear.MatrixUtils;
import org.apache.commons.math3.linear.SingularValueDecomposition;

import PoolHap.LocusAnnotation;

public class HapConfig {
	
	/*
	 * @author  Quan Long. Oct 12, 2018
	 * 
	 * Recording the configuration of the haplotypes, both in-pool and global.
	 * 
	 * The objects can be a small region or a whole chromosome.  
	 * The object of this class will be input or output parameters of other relaying algorithms.
	 * 
	 * Corresponding to the data structure in this class, there are two files:
	 * (1) global_hap_file: 
	 * 		(1.1) Each row represents a locus; each column represents a haplotype. 
	 * 	   	(1.2) The first row is the header of the IDs of all haps (usually are indexes); 
	 * 			  the second row is the global frequency (or NaN if unknown).
	 * 		(1.3) The first column is the IDs of all loci.
	 * 				(1.3.1) chr_index;start_loc;end_loc;alleles (the alleles are separated by ":"). chr_index starts with zero. 
	 * 				(1.3.2) NOTE: if it is an indel, start_loc=end_loc; Here start_loc!=end_loc only if it is a region, instead of a primitive locus.
	 * 		(1.3) Use '\t' to separate columns.   
	 * (2) local_hap_file: local_haplotypes
	 * 		(2.1) Each row represents the in-pool frequencies of all haplotypes (many of them can be zero)
	 * 			the first row is the header of all hap IDs (in the global file)
	 * 			the first element is the pool ID 
	 * 			the rest elements are the frequencies of this haplotype in all pools: zero indicates the absence.
	 * 		(2.2) Each column represents a haplotype's frequencies in all pools.
	 * 			the first column is the IDs of all pools.
	 * 		(2.3) Use '\t' to separate columns. 
	 * 		(2.4) In general, it is recommended that the order of haplotypes in this file is the same as the order of 
	 * 			the haplotypes in the global hap file. But if it is not the case, the program can handle it.
	 * 		(Note that this matrix representation is convenient for analysis of relationship between pools, e.g., transmission or evolution.) 
	 * 
	 * Note that the information of how to encode alleles to numbers will be generated on-the-fly
	 * in the future, we may add a function to read them from a file.  
	 * 
	 */
	
	// Redundant variables for convenience
	public int num_global_hap;				// number of global haplotype in this region
	public int num_loci;					// number of loci in this region
	public int num_pools;					// number of pools under study	
	public HashMap<String, Integer> hapID2index;  // map the pool IDs to their indexes in this.hap_IDs. 
	// Main structures
	public String[] hap_IDs;				// #global_hap
	public String[][] global_haps_string; 	// #global_hap x #loci: haplotypes that ever show up in the total population. String coded 
	public double[][] global_haps; 			// #global_hap x #loci: haplotypes that ever show up in the total population. floating number coded. 
	public double[] global_haps_freq; 		// #global_hap
	public String[] pool_IDs;				// #pools 		
	public double[][] in_pool_haps_freq;	// #global_hap x #pools  
	public LocusAnnotation[] locusInfo; 	// #loci. Note that a locus can be a SNP or a region. 
	public double[][] inpool_site_freqs;	// #loci x #pools; Added by Quan Dec. 2018.
	public int[] region=null;				// 2: region_start, region_end (these are indexes, only assigned when a regional HapConfig is formed) 			
	// Summary statistics
	public double[] mu;						// The 'average' allele at each position i.e. average haplotype Average frequency of alternate allele.
	public double[][] sigma;				// The variance-covariance matrix for each position.
	public double logL; 					// The log-likellihood of this HapConfig explaining the observed variant data.
	public int est_ind_pool; 				// Estimated number of individuals per pool (for now, assumed same per pool but easily extended). Used to calculate logL.
	// Link to the solving method
	PoolSolver2_Testing solver;
	
	// Constructors
	
	/*
	 * Constructor, forming the object by the variables in the memory.
	 * This constructor is useful in the divide & conquer algorithm.
	 */
	public HapConfig(String[][] global_haps_string, double[] global_haps_freq, // *** int[][] in_poll_haps deleted, int num_pools added.
			double[][] in_pool_haps_freq, double[][] inpool_site_freqs, LocusAnnotation[] locusInput, 
			int num_pools, String[] hap_IDs, String[] pool_IDs, int est_ind_pool){
		this.num_global_hap=global_haps_freq.length;
		this.num_loci=locusInput.length;
		this.global_haps_freq=global_haps_freq.clone();
		Boolean setID = false;
		if(hap_IDs!=null){ // *** Changed from this.hap_IDs (global) to hap_IDs (parameter)  
			this.hap_IDs=hap_IDs.clone();
		}else{ // if no IDs assigned, use indexes.
			this.hap_IDs = new String[this.num_global_hap];
			setID = true;
		}
		this.global_haps_string=global_haps_string.clone();
		for(int k=0;k<this.num_global_hap;k++) {
			this.global_haps_string[k]=global_haps_string[k].clone();
			if (setID)
				this.hap_IDs[k] = Arrays.toString(global_haps_string[k]);
			// System.out.println(this.hap_IDs[k]);
		}
		if(in_pool_haps_freq!=null){ // *** Matched format as above. 
			this.in_pool_haps_freq=in_pool_haps_freq.clone();
			for(int j=0;j<this.num_global_hap;j++)
				for(int k=0;k<this.num_pools;k++)
					this.in_pool_haps_freq[j][k]=in_pool_haps_freq[j][k];
		}else{ // if no intra-pool frequencies provided, start at 0. 
			this.in_pool_haps_freq = new double[this.num_global_hap][this.num_pools]; 
		}
		if(inpool_site_freqs!=null){
			this.inpool_site_freqs=inpool_site_freqs.clone();
		}
		this.locusInfo=locusInput;  // Originally, clone(), which is not a deep-clone. But we assume that the locusInfo won't be changed in the algorithm.
		// System.out.print("\n" + locusInput.length + "\t" + this.num_loci + "\t" + this.locusInfo.length);
		/*
		/for (int v = 0; v < this.num_loci; v++) {
			System.out.println(v);
			System.out.println("HapConfig constructor input:\t" + locusInput[v].alleles_coding);
			this.locusInfo[v].alleles_coding = new HashMap<String, Float>();
			System.out.println("HapConfig constructor pre:\t" + this.locusInfo[v].alleles_coding.size());
			for (String k : locusInput[v].alleles_coding.keySet()) {
				System.out.println("setting:\t" + locusInput[v].alleles_coding.get(k));
				this.locusInfo[v].alleles_coding.put(k, locusInput[v].alleles_coding.get(k));
			}
			System.out.println("HapConfig constructor post:\t" + this.locusInfo[v].alleles_coding.size());
		} // Need to specifically copy the mappings over to the new locus storage object. Shallow copy alone doesn't work.
		*/
		this.construct_hapID2index_map();
		this.encoding_haps(); // initialize this.global_haps;
		this.num_pools=num_pools;
		if(pool_IDs!=null){ // *** Matched format as above. 
			this.pool_IDs=pool_IDs.clone();
		}else{ // if no IDs assigned, use indexes.
			this.pool_IDs = new String[this.num_pools];
			for(int k=0;k<this.num_pools;k++)
				this.pool_IDs[k]=k+"";
		}
		this.est_ind_pool = est_ind_pool; 
	}
	
	/*
	 * Constructor that reads information from files. 
	 * This constructor is for the connection between modules (such as AEM, rjMCMC, etc.)
	 */
	public HapConfig(String global_hap_input_file, String in_pool_hap_input_file){
		try{
			// parse global_hap_input_file
			BufferedReader br=new BufferedReader(new FileReader(global_hap_input_file));
			String[] header_ids=br.readLine().split("\t");
			String[] header_freqs=br.readLine().split("\t");
			this.num_global_hap=header_ids.length-1;
			this.hap_IDs=new String[this.num_global_hap];
			this.global_haps_freq=new double[this.num_global_hap];
			for(int h=0;h<num_global_hap;h++){
				this.hap_IDs[h]=header_ids[h+1];
				this.global_haps_freq[h]=Double.parseDouble(header_freqs[h+1]);
			}
			String line=br.readLine();
			while(line!=null){
				this.num_loci++;
				line=br.readLine();
			}br.close();
			this.global_haps_string=new String[this.num_global_hap][this.num_loci];
			this.locusInfo=new LocusAnnotation[this.num_loci];
			br=new BufferedReader(new FileReader(global_hap_input_file));
			line=br.readLine();line=br.readLine(); // skip two headers
			line=br.readLine();
			int loci_index=0;
			while(line!=null){
				String[] tmp=line.split("\t");
				this.locusInfo[loci_index]=new LocusAnnotation(tmp[0]); // the first column is the locus-info
				for(int h_index=0;h_index<this.num_global_hap;h_index++)
					this.global_haps_string[h_index][loci_index]=tmp[h_index+1];
				// System.out.println(this.locusInfo[loci_index].alleles_coding.get("1"));
				loci_index++;
				line=br.readLine();
			}br.close();
			this.construct_hapID2index_map();
			this.encoding_haps();  // initialize this.global_haps;
			// then read the local-hap file
			br=new BufferedReader(new FileReader(in_pool_hap_input_file));
			String[] header_hap_IDs=br.readLine().split("\t"); // the first row is the hap IDs.
			if (header_hap_IDs[0].equals("Hap_ID")) {	// A *_haps.intra_freq.txt (in-pool frequencies of haplotypes) has been provided.
				if(this.num_global_hap!=header_hap_IDs.length-1)
					System.out.println("WRONG: this.num_global_hap is not the same between two files!");
				line=br.readLine();
				while(line!=null){ // count how many pools
					this.num_pools++;
					line=br.readLine();
				}br.close();
				this.pool_IDs=new String[this.num_pools];
				this.in_pool_haps_freq=new double[this.num_global_hap][this.num_pools];
				br=new BufferedReader(new FileReader(in_pool_hap_input_file));
				line=br.readLine(); //read the file again, skip the header 
				line=br.readLine();
				int pool_index=0;
				while(line!=null){
					String[] tmp=line.split("\t");
					this.pool_IDs[pool_index]=tmp[0]; // assign pool IDs
					for(int h=1;h<header_hap_IDs.length;h++){
						int h_index=this.hapID2index.get(tmp[h]); // based on the hap-ID, find out the hap-index
						this.in_pool_haps_freq[h_index][pool_index]=Double.parseDouble(tmp[h]);
					}
					pool_index++;
					line=br.readLine();
				}
				br.close();
			} else { // A *_vars.intra_freq.txt (in-pool frequencies of alternate alleles) has been provided.
				this.num_pools = header_hap_IDs.length - 1; 
				this.pool_IDs=new String[this.num_pools];
				for (int p = 0; p < this.num_pools; p++) this.pool_IDs[p] = header_hap_IDs[p + 1];
				this.inpool_site_freqs = new double[this.num_loci][this.num_pools];
				int locus_index = 0; 
				line = br.readLine();
				while(line != null){
					String tmp[] = line.split("\t"); 
					for (int p = 0; p < this.num_pools; p++) this.inpool_site_freqs[locus_index][p] = Double.parseDouble(tmp[p + 1]);
					locus_index++;
					line = br.readLine();
				} br.close();
			}
		}catch(Exception e){e.printStackTrace();}
	}
	
	public void construct_hapID2index_map(){
		this.hapID2index = new HashMap<String, Integer>();
		for(int h=0;h<this.num_global_hap;h++)
			this.hapID2index.put(this.hap_IDs[h], h);
	}
	
	/*
	 * Maps global_haps_string to global_haps (that is coded by integers). 
	 * 
	 * It will be done by invoking the function encoding_alleles in the class LocusAnnotation. 
	 * The key algorithm of how to design the encoding is implemented in the   
	 *    
	 */
	public void encoding_haps(){
		this.global_haps=new double[this.num_global_hap][this.num_loci];
		// System.out.println("\n" + this.num_global_hap  +"\t" + this.num_loci);
		// System.out.println("\nEncoded:" + this.locusInfo[0].alleles_coding.size());
		for(int h=0;h<this.num_global_hap;h++){
			for(int l=0;l<this.num_loci;l++){
				// if (h == this.num_global_hap - 1) System.out.println(h + "\t" + l + "\t" + this.global_haps_string[h][l] + "\t" + this.locusInfo[l].alleles_coding.get(this.global_haps_string[h][l]));
				this.global_haps[h][l]=
						this.locusInfo[l].
						alleles_coding.
						get(this.global_haps_string[h][l]);
			}
		}
	}
	
	// Write (STDOUT and text file) functions 
	
	/*
	 * Prints intermediate haplotype (single-region) guesses to STDOUT. 
	 */
	public void write_global_stdout() {
		DecimalFormat df = new DecimalFormat("#.####");
		df.setRoundingMode(RoundingMode.CEILING);    
		System.out.println("Hap.\tVar. Comp.\tInter. Freq.");
		for(int h=0;h<this.num_global_hap;h++) {
			System.out.print(h + "\t");
			for(int l=0;l<this.num_loci;l++){
				System.out.print(this.global_haps_string[h][l]);  
			}
			System.out.println("\t" + df.format(this.global_haps_freq[h]));
		}
	}	
	
	/*
	 * Output the in global_haplotypes using string alleles. 
	 */
	public void write_global_file_string(String global_hap_output_file, boolean append){
		try{
			BufferedWriter bw= new BufferedWriter(new FileWriter(global_hap_output_file, append));
			bw.write("Hap_ID");
			for(int h=0;h<this.num_global_hap;h++)
				bw.write("\t"+this.hap_IDs[h]);
			bw.write("\nFreq");
			for(int h=0;h<this.num_global_hap;h++)
				bw.write("\t"+this.global_haps_freq[h]);
			bw.write("\n");
			for(int l=0;l<this.num_loci;l++){
				bw.write(this.locusInfo[l].output2string());
				for(int h=0;h<this.num_global_hap;h++)
					bw.write("\t"+this.global_haps_string[h][l]);
				bw.write("\n");
			}bw.close();
		}catch(Exception e){e.printStackTrace();}
	}
	
	/*
	 * Output the in global_haplotypes using coded alleles. 
	 * 
	 * During the analysis, e.g., divide-and-conquer, when read by the program, these coded alleles 
	 * will become strings and then the next round of coding will be performed. 
	 */
	public void write_global_file_code(String global_hap_output_file, boolean append){
		try{
			BufferedWriter bw= new BufferedWriter(new FileWriter(global_hap_output_file, append));
			bw.write("Hap_ID");
			for(int h=0;h<this.num_global_hap;h++)
				bw.write("\t"+this.hap_IDs[h]);
			bw.write("\nFreq");
			for(int h=0;h<this.num_global_hap;h++)
				bw.write("\t"+this.global_haps_freq[h]);
			bw.write("\n");
			for(int l=0;l<this.num_loci;l++){
				bw.write(this.locusInfo[l].output2string());
				for(int h=0;h<this.num_global_hap;h++)
					// the only difference between this method and "write_global_file_string" is the line below: 
					bw.write("\t"+this.global_haps[h][l]);  
				bw.write("\n");
			}bw.close();
		}catch(Exception e){e.printStackTrace();}
	}
	
	/*
	 * Output the in pool frequencies. 
	 */
	public void write_inpool(String inpool_hap_output_file, boolean append){
		try{
			BufferedWriter bw= new BufferedWriter(new FileWriter(inpool_hap_output_file, append));
			bw.write("Hap_ID");
			for(int h=0;h<this.num_global_hap;h++)
				bw.write("\t"+this.hap_IDs[h]);
			bw.write("\n");
			// System.out.println(this.in_pool_haps_freq.length + "\t" + this.in_pool_haps_freq[0].length);
			for(int p=0;p<this.num_pools;p++){
				bw.write(this.pool_IDs[p]);
				for(int h=0;h<this.num_global_hap;h++)	// TODO Report error! Formerly, this.num_pools.
					if (this.in_pool_haps_freq[h].length == 0) bw.write("\t0");
					else bw.write("\t"+this.in_pool_haps_freq[h][p]);
				bw.write("\n");
			}
			bw.close();
		}catch(Exception e){e.printStackTrace();}
	}
	
	public void write2files(String global_hap_output_file, String in_pool_hap_output_file, String type, boolean append){
		write_inpool(in_pool_hap_output_file, append);
		if(type.equals("string")){
			write_global_file_string(global_hap_output_file, append);
		}else if(type.equals("code")){
			write_global_file_code(global_hap_output_file, append);
		}else{
			System.out.println("Error: the type has to be string or code!");
		}
	}

	// rjMCMC
	
	/*
	 * Updates the composition of the main global variables when a set of new haplotypes are added.
	 * @param The variant composition of new haplotypes and how much frequency to allot them.  
	 */
	public void addHaps(ArrayList<double[]> list_add_haps, double freqs_sum, int old_hap_index) {
		// System.out.println(this.num_global_hap + "\t" + this.num_loci); 
		int num_new_haps = list_add_haps.size();
		int tot_new_haps = this.num_global_hap + num_new_haps; 
		String[][] tmp_global_haps_string = new String[tot_new_haps][this.num_loci];
		String[] tmp_hap_IDs = new String[tot_new_haps];
		for(int h=0;h<this.num_global_hap;h++) {
			tmp_global_haps_string[h]=this.global_haps_string[h].clone();
			tmp_hap_IDs[h]=this.hap_IDs[h];
		}
		for(int h=this.num_global_hap;h<tot_new_haps;h++) {
			double[] tmp_new_hap = list_add_haps.get(h - this.num_global_hap); 
			String tmp_id = ""; 
			for (int l = 0; l < this.num_loci; l++) {
				tmp_global_haps_string[h][l] = Integer.toString((int)(tmp_new_hap[l]));
				// System.out.print(tmp_global_haps_string[h][l] + " "); 
				tmp_id += (int)(tmp_new_hap[l]); 
			}
			// System.out.println();
			tmp_hap_IDs[h] = tmp_id;
		}
		this.global_haps_string = new String[tot_new_haps][this.num_loci];
		for(int h=0;h<tot_new_haps;h++) {
			for (int l = 0; l < this.num_loci; l++) this.global_haps_string[h][l] = tmp_global_haps_string[h][l];
		}
		this.hap_IDs = tmp_hap_IDs.clone(); 
		
		double[] tmp_global_haps_freq = new double[tot_new_haps];
		if (old_hap_index == -1) {
			for(int h=0;h<this.num_global_hap;h++)
				tmp_global_haps_freq[h]=this.global_haps_freq[h] * (1 - freqs_sum);
			for(int h=this.num_global_hap;h<tot_new_haps;h++)
				tmp_global_haps_freq[h] = (double) freqs_sum / num_new_haps;
		} else {
			for(int h=0;h<this.num_global_hap;h++) {
				if (h == old_hap_index)
					tmp_global_haps_freq[h]=this.global_haps_freq[h] - freqs_sum; 	
				else 
					tmp_global_haps_freq[h]=this.global_haps_freq[h];
			}
			tmp_global_haps_freq[this.num_global_hap] = freqs_sum;	// The last one is the global frequency of the new haplotype.
		}
		this.global_haps_freq = tmp_global_haps_freq.clone(); 
		this.num_global_hap=this.global_haps_freq.length;

		this.construct_hapID2index_map();
		this.encoding_haps(); // initialize this.global_haps;
	}

	
	/*
	 * Checks the rank of the initial haplotypes to make sure that every primitive locus is represented.
	 * @param How much frequency to allot new haplotypes if they're needed.  
	 */
	public void checkrank_and_fullfill(double freqs_sum){
		ArrayList<Integer> noSNP = new ArrayList<Integer>();
		if(this.num_global_hap<this.num_loci) {
			for(int l=0;l<num_loci;l++) {
				int varPresent = 0; 
				for(int h=0;h<num_global_hap;h++){
					varPresent += this.global_haps[h][l]; // Checks to see if any haplotype has an alternative allele in that variant position. 
				}
				if (varPresent == 0) noSNP.add(l); // If it doesn't then that rank will be unfilled.
			}
			int rank = this.num_loci - noSNP.size();
			this.fulfill(freqs_sum, rank, noSNP);
		}
		else{
			double[][] the_H_array=new double[this.num_global_hap][this.num_loci];
			for(int h=0;h<num_global_hap;h++){
				for(int l=0;l<num_loci;l++)the_H_array[h][l]=this.global_haps[h][l];
			}
			SingularValueDecomposition svd=new SingularValueDecomposition(MatrixUtils.createRealMatrix(the_H_array));
			double[] sv=svd.getSingularValues();
			int rank=0;
			for(int i=0;i<sv.length;i++){
				if(sv[i]!=0)rank++;
				if(sv[i]==0) noSNP.add((int) sv[i]); 
			}if(rank<num_loci)this.fulfill(freqs_sum,rank,noSNP);
		}
		this.update_sigma_mu_logL();
	}
	
	/*
	 * Generates semi-random haplotypes to make sure that every primitive locus is represented.
	 * @param How much frequency to allot new haplotypes, the number of, and which loci to cover.  
	 */
	public void fulfill(double freqs_sum, int rank, ArrayList<Integer> noSNP){
		 
		for(int k=0;k<this.num_global_hap;k++) this.global_haps_freq[k] *= (1-freqs_sum);
		int fillRank = this.num_loci - rank;
		ArrayList<double[]> list_add_haps = new ArrayList<double[]>();
		for (int p = 0; p < fillRank; p++) {
			double[] new_hap=new double[this.num_loci];
			new_hap[noSNP.get(p)]=1;	// Fill a 'rank'. 
			int addVars = ThreadLocalRandom.current().nextInt(0,(int) (num_loci / 2));	// To 'thin out' the new haplotypes being created. 
			for (int i = 0; i < addVars; i++) new_hap[ThreadLocalRandom.current().nextInt(0,num_loci)]=1;
			if(search_a_hap(new_hap, list_add_haps) == -1) // If this proposed haplotype doesn't exist yet.
				list_add_haps.add(new_hap); 
		} 
		addHaps(list_add_haps, freqs_sum, -1); 
	}

	/*
	 * Updates the composition of the main global variables when existing haplotypes are removed.
	 * @param Whether to remove the haplotype or not, and the total number to remove.  
	 */
	public void remHaps(boolean[] list_rem_haps, int num_rem_haps) {
		int tot_new_haps = this.num_global_hap - num_rem_haps; 
		// System.out.println(tot_new_haps + "\t" + this.num_global_hap + "\t" + num_rem_haps); 
		String[][] tmp_global_haps_string = new String[tot_new_haps][this.num_loci];
		String[] tmp_hap_IDs = new String[tot_new_haps];
		double tmp_tot_freq = 0; 
		int tmp_count = 0; 
		for(int h=0;h<this.num_global_hap;h++) {
			if (list_rem_haps[h]) {
				// System.out.print(h);
				continue;
			}
			tmp_global_haps_string[tmp_count]=this.global_haps_string[h].clone();
			tmp_hap_IDs[tmp_count]=this.hap_IDs[h];
			tmp_tot_freq += this.global_haps_freq[h];
			tmp_count++; 
		}
		this.global_haps_string = tmp_global_haps_string.clone(); 
		this.hap_IDs = tmp_hap_IDs.clone();
		double[] tmp_global_haps_freq = new double[tot_new_haps];
		int new_hap = 0;
		for(int h=0;h<this.num_global_hap;h++) {
			if (list_rem_haps[h]) continue;
			tmp_global_haps_freq[new_hap] = this.global_haps_freq[h] / tmp_tot_freq;
			new_hap++; 
		}
		this.global_haps_freq = tmp_global_haps_freq.clone(); 
		this.num_global_hap=this.global_haps_freq.length;
		this.construct_hapID2index_map();
		this.encoding_haps(); // initialize this.global_haps;
	}

	/*
	 * Looks through the list of proposed haplotypes to see if the proposed one is there already. 
	 * @param Variant composition of the proposed haplotype. The list of proposed haplotypes.
	 * @return The index of an identical-match haplotype. -1 if no match.  
	 */
	public int search_a_hap(double[] hap, ArrayList<double[]> list_add_haps){
		for(int h=0;h<list_add_haps.size(); h++){
			boolean match=true;
			double[] the_candidate=list_add_haps.get(h);
			for(int l=0;l<this.num_loci;l++){
				if(hap[l]!=the_candidate[l]){
					match=false; break;
				}
			}if (match) return h;
		}return -1;
	}
	
	/*
	 * Looks through existing haplotypes to see if the proposed one is there already. 
	 * @param Variant composition of the proposed haplotype. If the mutant has been created already, 'avoid' this index. avoid = -1 if the haplotype has not been added yet. 
	 * @return The index of an identical-match haplotype. -1 if no match.  
	 */
	public int search_a_hap(double[] hap, int avoid){
		for(int h=0;h<this.num_global_hap;h++){
			if (h == avoid) continue;
			boolean match=true;
			double[] the_candidate=this.global_haps[h];
			for(int l=0;l<this.num_loci;l++){
				if(hap[l]!=the_candidate[l]){
					match=false; break;
				}
			}if (match) return h;
		}return -1;
	}
	
	/*
	 * Calculates the average variant frequency, var.-covar. matrix for different primitive loci
	 * Calculates the log-likelihood of this HapConfig object explaining the data. 
	 * Updates the global variables mu, sigma, and logL respectively for this HapConfig object. 
	 */
	public void update_sigma_mu_logL(){
		this.mu=new double[this.num_loci];
	    // System.out.println("\nmu:"); 
		for(int l=0;l<this.num_loci;l++){
	    	for(int h=0;h<this.num_global_hap;h++){		    
	    		this.mu[l]=this.mu[l]+this.global_haps[h][l]*this.global_haps_freq[h];
	    	}
    		// System.out.print(this.mu[l] + "\t");
	    } 
	    double[][] eta=new double[this.num_loci][this.num_loci];
	    for(int l1=0;l1<this.num_loci;l1++){
	    	 for(int l2=0;l2<this.num_loci;l2++){
	    		 for(int h=0;h<this.num_global_hap;h++){
	 		    	eta[l1][l2]+=(this.global_haps[h][l1]*this.global_haps[h][l2]*this.global_haps_freq[h]);
	 		     }
			 }
	    }
	    // System.out.println("\n\nsigma:"); 
	    this.sigma=new double[this.num_loci][this.num_loci]; // TODO Update the sigma function to use LDx
	    for(int q1=0;q1<this.num_loci;q1++){
	    	 for(int q2=0;q2<this.num_loci;q2++){
	    		 this.sigma[q1][q2]=eta[q1][q2]-this.mu[q1]*this.mu[q2];
	    		 // System.out.print(this.sigma[q1][q2] + "\t");
	    	 }
	    	 // System.out.println();
	    }
	    this.logL=Algebra.logL_aems(this.sigma, this.mu, this.inpool_site_freqs);
	    // Got rid of Algebra.times(this.sigma, this.est_ind_pool), Algebra.times(this.mu, this.est_ind_pool) because everything is in frequencies.
	}

	/*
	 * Returns haplotypes that are above a certain frequency cutoff.
	 * @param the frequency cutoff.
	 */
	public HapConfig clone(double freq_cutoff){
		if (freq_cutoff == 0) return new HapConfig(this.global_haps_string, this.global_haps_freq, this.in_pool_haps_freq, this.inpool_site_freqs, this.locusInfo, this.num_pools, this.hap_IDs, this.pool_IDs, this.est_ind_pool);
		boolean[] list_rem_haps = new boolean[this.num_global_hap];
		int num_rem_haps = 0; 
		for (int h = 0; h < this.num_global_hap; h++) {
			if (this.global_haps_freq[h] < freq_cutoff) {
				list_rem_haps[h] = true;
				num_rem_haps++;
			}
		}
		// for (int h = 0; h < this.num_global_hap; h++)
			// System.out.println(this.hap_IDs[h]); 
		remHaps(list_rem_haps, num_rem_haps);
		// for (int h = 0; h < this.num_global_hap; h++)
			// System.out.println(this.hap_IDs[h]); 
		return new HapConfig(this.global_haps_string, this.global_haps_freq, this.in_pool_haps_freq, this.inpool_site_freqs, this.locusInfo, this.num_pools, this.hap_IDs, this.pool_IDs, this.est_ind_pool);
	}

	/*
	 * Sample two haplotypes completely randomly from the global set.
	 * @return The indices of two randomly selected global haplotypes. Donor is first.
	 */
	int[] sample_two_haps(){
		int i1=(int)(ThreadLocalRandom.current().nextDouble()*this.num_global_hap);
	    int i2=(int)(ThreadLocalRandom.current().nextDouble()*this.num_global_hap);
	    if(i1==i2){
	    	if(i2!=this.num_global_hap-1)i2++;
	    	else i1--;
	    }
	    return new int[] {i1, i2};
	}

	/*
	 * Sample two haplotypes according to their frequencies. 
	 * Higher frequency means more likely to be picked as the 'donor' haplotype, i1.
	 * The second haplotype is selected randomly from the rest of the haplotypes. 
	 * @return The indices of frequency-selected global haplotypes. Donor is first.
	 */
	int[] sample_two_haps_prob(){	// !!!
		double cumul_lim1 = ThreadLocalRandom.current().nextDouble();	// Real number between 0 and 1; 'floor' of maximum cumulative frequency.
		int i1 = 0, i2 = 0; 
		double cumul_freq = this.global_haps_freq[i1];
		while ((cumul_lim1 > cumul_freq) && (i1 < (this.num_global_hap - 1))) {
			i1++; 
			cumul_freq += this.global_haps_freq[i1];		
		}
		double cumul_lim2 = ThreadLocalRandom.current().nextDouble() * (1 - this.global_haps_freq[i1]);
		if (i1 == 0) i2 = 1;
		cumul_freq = this.global_haps_freq[i2]; 
		while ((cumul_lim2 > cumul_freq) && (i2 < (this.num_global_hap - 1))) {
			i2++; 
			if (i2 == i1) i2++; 
			if (i2 == this.num_global_hap) { 
				i2--;
				i1--; 
			}
			cumul_freq += this.global_haps_freq[i2]; 
		}
		return new int[] {i1, i2};
	}
	
	/*
	 * Returns either a successful or failed update to the existing global haplotype frequencies.
	 * @param First three are beta distribution (transition probability) parameters. iter determines whether or not to sample haplotypes randomly. The last two are for mu/sigma/logL updates.
	 * @return 1 for a successful update, -1 for a failed update. 
	 */
	int update_freqs(double beta_a, double beta_c, double alpha, int iter){	// !!!
		int[] indices = new int[2]; 
		if ((iter % 3) == 0) {	// Every third iteration, sample two haplotypes completely randomly according to index.
			indices = this.sample_two_haps();	
		} else {	// Otherwise, sample two haplotypes according to their frequencies (higher = more likely to be picked).
			indices = this.sample_two_haps_prob();
		}
		double ori_freq_1 = this.global_haps_freq[indices[0]];
		double ori_freq_2 = this.global_haps_freq[indices[1]]; // Backup for returning to the original of rejected. 
	    BetaDistribution beta_dist=new BetaDistribution(beta_a+ori_freq_1*beta_c, beta_a+ ori_freq_2*beta_c);	// TODO Why is the transition matrix a beta distribution? Why do the parameters change with the frequencies selected?   
	    double proportion = beta_dist.sample(); 		//gsl_ran_beta(rng,beta_a+p[i1]*beta_c,beta_a+p[i2]*beta_c);
	    this.global_haps_freq[indices[0]] = ori_freq_1 * proportion;	// Haplotype 1 donates part of its frequency to haplotype 2. 
	    this.global_haps_freq[indices[1]] = ori_freq_2 + ori_freq_1 * (1 - proportion);	// Note that this is different than hippo. Which for some reason uses (ori_freq_1 + ori_freq_2) * proportion. 
	    double new_freq_1 = this.global_haps_freq[indices[0]]; 
	    double new_freq_2 = this.global_haps_freq[indices[1]]; 
	    double tmp_logl = this.logL; 
	    this.update_sigma_mu_logL();	// v0.5 this.logL is <v0.4 logl2.
	    double logHR=this.logL - tmp_logl;	
	    BetaDistribution beta_dist2=new BetaDistribution(beta_a+new_freq_1*beta_c,beta_a+new_freq_2*beta_c);  
	    double prop_dens = 0;	// !!!
	    try {	// This was included because when the original frequency of haplotype 2 was too small, beta_dist2.logDensity returned NaN because it couldn't use a number so close to 0.
	    	double below_min = 0; 
	    	if (ori_freq_2 > Math.pow(10, -17)) below_min = ori_freq_2; 
	    	prop_dens = beta_dist2.logDensity(ori_freq_1/(ori_freq_1+below_min))-beta_dist.logDensity(new_freq_1/(new_freq_1+new_freq_2));
	    } catch (NumberIsTooSmallException e) {
	    }
	    logHR=logHR+prop_dens;    
	    logHR+=(alpha-1)*(Math.log(new_freq_1)+Math.log(new_freq_2)-Math.log(ori_freq_1)-Math.log(ori_freq_2)); // TODO Why is there an additional scaling parameter? 
	    if ((iter % 3) != 0) { // TODO Why is there an additional scaling parameter when the haplotypes are not chosen according to frequency?
	    	logHR+=Math.log(new_freq_1)+Math.log(new_freq_2)+Math.log(1-ori_freq_2)+Math.log(1-ori_freq_1)-(Math.log(1-new_freq_1)+Math.log(1-new_freq_2)+Math.log(ori_freq_2)+Math.log(ori_freq_1));
	    }
		if( (!(Double.isNaN(this.logL))) && Math.log(ThreadLocalRandom.current().nextDouble())<=logHR){//proposal accepted
			if(this.logL > solver.logL_best) solver.logL_best = this.logL;
			return 1;
		}
		else{//proposal rejected
		    this.global_haps_freq[indices[0]] = ori_freq_1;	
		    this.global_haps_freq[indices[1]] = ori_freq_2;	
		    this.update_sigma_mu_logL();	// Restore old values of mu, sigma, and logL. This might slow the MCMC down a lot. 
		    return 0;
		}		
	}
	
	/*
	 * Returns either a successful or failed update to the haplotype composition by mutating one locus of one haplotype.
	 * @param The last two are for mu/sigma/logL updates.
	 * @return 1 for a successful update, -1 for a failed update. 
	 */
	int mutate_a_hap(){
		int hap_index=this.sample_one_hap(); // TODO Does this happen anywhere else? If not, get rid of it!
		int loc_index = this.generate_mutant(hap_index, -1);
		double[] tmp_hap = new double[this.num_loci]; 
		for (int l = 0; l < this.num_loci; l++) tmp_hap[l] = this.global_haps[hap_index][l]; 
		if(search_a_hap(tmp_hap, loc_index)!=-1)return 0;
	    double tmp_logl = this.logL; 
	    this.update_sigma_mu_logL();	// v0.5 this.logL is <v0.4 logl2.

		double logHR = this.logL - tmp_logl; // This is the log-likehood only because the transition probabilities are the same.
		if( (!(Double.isNaN(this.logL))) && Math.log(ThreadLocalRandom.current().nextDouble())<=logHR){//proposal accepted
			if(this.logL > solver.logL_best) solver.logL_best = this.logL;
			return 1;
		}
		else{//proposal rejected
			this.generate_mutant(hap_index, loc_index);
		    this.update_sigma_mu_logL();	// Restore old values of mu, sigma, and logL. This might slow the MCMC down a lot. 
			return 0;
		}	
	}

	/*
	 * Sample a haplotype completely randomly from the global set.
	 * @return The index of one randomly selected global haplotype.
	 */
	int sample_one_hap(){
		return (int)(ThreadLocalRandom.current().nextDouble()*this.num_global_hap);
	}

	/*
	 * Generates a mutant at a single locus in a single haplotype from the existing global haplotype set. 
	 * The locus can be randomly selected (proposal) or specified (upon rejection).
	 * @return The index of the mutated locus.
	 */
	int generate_mutant(int hap_index, int revert){
		int loc_index = 0; 
		if (revert == -1) {	// Make a random mutant.
			loc_index = (int) (ThreadLocalRandom.current().nextDouble() * this.num_loci);			
		} else {	// If the mutation proposal fails, revert to the original haplotype.
			loc_index = revert;
		}
		double mut_allele = Math.abs(this.global_haps[hap_index][loc_index] - 1); 
		this.global_haps_string[hap_index][loc_index] = Long.toString(Math.round(mut_allele));
		this.global_haps[hap_index][loc_index] = mut_allele;	// Don't need to re-encode all haplotypes for a single locus change.
		return loc_index; 
	}

	double[] generate_mutant(int hap_index){
		int loc_index = (int) (ThreadLocalRandom.current().nextDouble() * this.num_loci);			
		double mut_allele = Math.abs(this.global_haps[hap_index][loc_index] - 1); 
		double[] tmp_hap = this.global_haps[hap_index].clone();
		tmp_hap[loc_index] = mut_allele;
		return tmp_hap;
	}

	/*
	 * Add one hap by mutating an existing hap, or remove one hap by coalescing it with a close one.
	 * @param A scaling parameter, the prior likelihood of adding a new haplotype, another scaling parameter, and the number of locus mismatches acceptable for coalescence.
	 * @return 0 if nothing done, 1 if mutated, -1 if mutant rejected, 2 if coalesced, -2 if coalescence rejected.  
	 */
	// TODO Figure out all of the derivations here.
	int add_mutant_or_coalesce(double alpha, double p_add, double gamma, int coalescing_mismatch){
		int hap_index = this.sample_one_hap();
		boolean proceed_mut;
		if(this.num_global_hap<=2)proceed_mut=true;
		else proceed_mut = ThreadLocalRandom.current().nextDouble() < p_add;
		if(proceed_mut){
			double proportion = solver.new_old_haps_beta.sample(); // Move some frequency from the old hap to the new one. 
			double ori_freq = this.global_haps_freq[hap_index];
			double[] tmp_hap = generate_mutant(hap_index); 
			if(search_a_hap(tmp_hap, -1)!=-1){return 0;}
			ArrayList<double[]> tmp_list = new ArrayList<double[]>();	// This might slow it down by a lot by creating a new ArrayList of size 1.
			tmp_list.add(tmp_hap); 
			this.addHaps(tmp_list, ori_freq * proportion, hap_index);

			double tmp_logl = this.logL; 
		    this.update_sigma_mu_logL();	// v0.5 this.logL is <v0.4 logl2.
			double logHR = this.logL - tmp_logl;
		    logHR-=gamma;	// TODO why is this yet another scaling constant? 

		    logHR+=(alpha-1)*(Math.log(this.global_haps_freq[hap_index])+Math.log(this.global_haps_freq[this.num_global_hap - 1])-Math.log(ori_freq)); // TODO Seems to be part of the prior for p, along w/ log-likelihood
		    logHR+=-Math.log(this.num_global_hap)+Math.log(1-p_add)+Math.log(coalescing_probability(hap_index, this.num_global_hap-1, coalescing_mismatch));//inverse transition probability
		    logHR+=Math.log((this.num_global_hap-1)*this.num_loci*1.0)-Math.log(solver.new_old_haps_beta.density(proportion)); //transition probability
		    logHR+=Math.log(ori_freq); //|Jacobian|    
		    
			if( (!(Double.isNaN(this.logL))) && Math.log(ThreadLocalRandom.current().nextDouble())<=logHR){//proposal accepted
				if(this.logL > solver.logL_best) solver.logL_best = this.logL;
				return 1;
			}
			else{//proposal rejected
				this.global_haps_freq[hap_index]=ori_freq;
				boolean[] list_rem_hap = new boolean[this.num_global_hap];
				list_rem_hap[this.num_global_hap - 1] = true;
				remHaps(list_rem_hap, 1);
				return -1;
			}	
		}else{ //coalesce index with the closest hap.
			double[] sampling_prob=new double[1];
			int coal_index=seek_coalescence(hap_index, coalescing_mismatch, sampling_prob);
			if(coal_index==-1)return 0; // no hap selected, return 0.
			double ori_freq_1=this.global_haps_freq[hap_index];
			double ori_freq_2=this.global_haps_freq[coal_index];
			this.global_haps_freq[hap_index]=0;
			this.global_haps_freq[coal_index]+=ori_freq_1;
			// calculating new logL:
			double tmp_logl = this.logL; 
		    this.update_sigma_mu_logL();	// v0.5 this.logL is <v0.4 logl2.
			double logHR = this.logL - tmp_logl;

		    logHR+=gamma;
		    logHR+=(alpha-1)*(Math.log(this.global_haps_freq[coal_index])-Math.log(ori_freq_1)-Math.log(ori_freq_2)); //prior for p
		    logHR-=-Math.log(this.num_global_hap)+Math.log(1-p_add)+Math.log(sampling_prob[0]);//transition probability
		    logHR+=-Math.log((this.num_global_hap-1)*this.num_loci*1.0)+
		    		Math.log(solver.new_old_haps_beta.density(ori_freq_2/this.global_haps_freq[coal_index])); //inverse transition probability
			logHR+=-Math.log(ori_freq_1+ori_freq_2); //|Jacobian|
			
			if( (!(Double.isNaN(this.logL))) && Math.log(ThreadLocalRandom.current().nextDouble())<=logHR){//proposal accepted
				if(this.logL > solver.logL_best) solver.logL_best = this.logL;
				boolean[] list_rem_hap = new boolean[this.num_global_hap];
				list_rem_hap[hap_index] = true;
				remHaps(list_rem_hap, 1);
				return 2;
			}
			else{//proposal rejected
				this.global_haps_freq[hap_index]=ori_freq_1;
				this.global_haps_freq[coal_index]=ori_freq_2;
				return -2;
			}	
		}					
	}

	/*
	 * The inverse transition probability for add_mutant. Calculates the likelihood that a haplotype coalesces toward the original, i1.
	 * Is calculated according to the total frequency of haplotypes that are less than the similarity limit away from the proposed haplotype.
	 * @param The indices of the two haplotypes, and the similarity limit.
	 * @return The likelihood of i2 coalescing towards i1 expressed as the frequency of i2 over total acceptable haplotypes.
	 */
	double coalescing_probability(int i1, int i2, int coalescing_mismatch){
		double sum=0;	// Frequencies of the haplotypes that are 1 position away from the original haplotype i1.
		for(int h=0;h<this.num_global_hap - 1;h++){	// The last one is the newest haplotype, so don't look at it. 
			int diff=0;
			for(int l=0;l<this.num_loci;l++){
				if(this.global_haps[i1][l]!=this.global_haps[h][l]) diff++;
				if(diff>coalescing_mismatch) break;
			}
			if(diff<=coalescing_mismatch)sum+=this.global_haps_freq[h];
			if(h==i2 && diff!=1)System.out.println("ERROR: in coalescing probability, d(h(i1),h(12))>1\n");
	    }
		return this.global_haps_freq[i2]/sum;
	}
	
	/*
	 * Job 1: The index picker for or_coalesce. From the haplotypes that are less than the similarity limit away from the selected haplotype, pick one for the selected haplotype to coalesce to. 
	 * Job 2: The transition probability for or_coalesce. Basically does what coalescing_probability does while scanning through the global set for potentials.
	 * @param The selected haplotype, the similarity limit, and a one-level-up container for the transition probability.
	 * @return The randomly picked haplotype for the selected haplotype to coalesce to.  
	 */
	int seek_coalescence(int hap_index, int coalescing_mismatch, double[] sampling_prob){
		int unacceptable=coalescing_mismatch+1;
		// double sum_freq_closest=0;
		int[] distance=new int[this.num_global_hap];
		int smallest_found=unacceptable; // initialize
		distance[hap_index]=unacceptable;  // itself will be excluded.
		double[] the_target=this.global_haps[hap_index].clone();
		for(int h=0;h<this.num_global_hap;h++){
			if(h==hap_index)continue;
			for(int l=0;l<this.num_loci;l++){				
				if(the_target[l]!=this.global_haps[h][l]){
					distance[h]++;
					if(distance[h]==unacceptable)break;
				}
			}if(distance[h]<smallest_found) smallest_found=distance[h];
		}if(smallest_found==unacceptable)return -1; // No haplotype found that is variant-wise close enough to coalesce with the chosen haplotype.
		ArrayList<Integer> the_closest_ones=new ArrayList<>();		
		for(int h=0;h<this.num_global_hap;h++){
			if(smallest_found==distance[h]){
				the_closest_ones.add(h);
				// sum_freq_closest+=nodes.get(h).freq;
			}
		}// randomly choose an index from the closest ones:
		int the_chosen=(int)(ThreadLocalRandom.current().nextDouble()*the_closest_ones.size());
		sampling_prob[0]=1.0/the_closest_ones.size();
		return the_closest_ones.get(the_chosen);
	}

}