package MiscFunctions;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.util.HashMap;

public class SimpleConfig {

	// Redundant variables for convenience
	public int num_hap;						// number of global haplotype in this region
	public int num_loci;					// number of loci in this region
	// Main structures
	public int[][] hap_var_comp; 			// #global_hap x #loci: haplotypes that ever show up in the total population. floating number coded. 
	public double[] hap_freq; 				// #global_hap
	public int[] locus_positions;			// #loci 

	public SimpleConfig(String hap_file) throws IOException {	// This constructor is for haplotypes in standard output files.
		BufferedReader br = new BufferedReader(new FileReader(hap_file));
		br.readLine();	// Skip IDs. 
		String[] header_freqs = br.readLine().split("\t");
		this.num_hap = header_freqs.length - 1;
		this.hap_freq = new double[this.num_hap];
		for(int h = 0; h < num_hap; h++)
			this.hap_freq[h] = Double.parseDouble(header_freqs[h + 1]);
		while(br.readLine() != null)
			this.num_loci++;
		br.close();
		this.hap_var_comp = new int[this.num_hap][this.num_loci];
		this.locus_positions = new int[this.num_loci];
		br = new BufferedReader(new FileReader(hap_file));
		String line = br.readLine(); line=br.readLine(); // skip two headers
		line = br.readLine();
		int loci_index = 0;
		while (line != null){
			String[] tmp = line.split("\t");
			this.locus_positions[loci_index] = Integer.parseInt(tmp[0].split(";")[1]); // the second column is the variant position
			for(int h = 0; h < this.num_hap; h++)
				this.hap_var_comp[h][loci_index] = Integer.parseInt(tmp[h + 1]);
			loci_index++; 
			line = br.readLine();
		} br.close();
	}

	public SimpleConfig(String work_dir, int num_pools, String var_file) throws IOException {	// This constructor is for haplotypes in GC output files.
		double hap_ct_all = 0; 
		HashMap<String, Integer> hap_tracker = new HashMap<String, Integer>();
		String[] tmp = null;
		for (int p = 0; p < num_pools; p++) {
			BufferedReader br = new BufferedReader(new FileReader(work_dir + "p" + p + ".in"));
			String line = br.readLine(); 
			while(line != null) {
				tmp = line.split("\\t|-|\\?"); 
				String hap = "";
				for (int l = 0; l < tmp.length - 1; l++) hap += tmp[l];
				int hap_ct = Integer.parseInt(tmp[tmp.length - 1]); 
				if (!hap_tracker.containsKey(hap)) hap_tracker.put(hap, hap_ct); 
				else {
					int prev_ct = hap_tracker.get(hap); 
					hap_tracker.put(hap, prev_ct + hap_ct);
				}
				hap_ct_all += hap_ct;
				line = br.readLine();
			} br.close();
		}
		this.num_loci = tmp.length - 1;
		this.num_hap = hap_tracker.size(); 
		this.hap_var_comp = new int[this.num_hap][this.num_loci];
		this.locus_positions = new int[this.num_loci];
		int[] hap_ct = new int[this.num_hap];
		int hap_index = 0;
		for (String hap : hap_tracker.keySet()) {
			String[] tmp2 = hap.split("");
			for (int l = 0; l < this.num_loci; l++)
				this.hap_var_comp[hap_index][l] = Integer.parseInt(tmp2[l]); 
			hap_ct[hap_index] = hap_tracker.get(hap);
			hap_index++;
		}		
		this.hap_freq = new double[this.num_hap];
		for(int h = 0; h < this.num_hap; h++)
			this.hap_freq[h] = (double) hap_ct[h] / hap_ct_all; 		
		this.locus_positions = new int[this.num_loci];
		BufferedReader br = new BufferedReader(new FileReader(var_file));
		int pos_index = 0; 
		String line = br.readLine(); // Skip header. 
		line = br.readLine();
		while (line != null) {
			this.locus_positions[pos_index] = Integer.parseInt(line.split(";")[1]); 
			pos_index++; 
			line = br.readLine();
		} br.close();
	}

	public SimpleConfig(int new_hap, int new_loci, int[][] new_var_comp, double[] new_freq, int[] new_positions) {
		this.num_hap = new_hap;
		this.num_loci = new_loci;
		this.hap_var_comp = new_var_comp;
		this.hap_freq = new_freq;
		this.locus_positions = new_positions;  
	}

	public SimpleConfig subsetter(int r_start, int r_end) {
		HashMap<String, Double> hap_tracker = new HashMap<String, Double>();
		int new_loci = r_end - r_start + 1;
		int[] new_positions = new int[new_loci]; 
		for (int h = 0; h < this.num_hap; h++) {
			String hap = ""; 
			for (int l = r_start; l <= r_end; l++) {
				hap += this.hap_var_comp[h][l]; 
				new_positions[l - r_start] = this.locus_positions[l]; 
			}
			if (!hap_tracker.containsKey(hap)) hap_tracker.put(hap, this.hap_freq[h]); 
			else {
				double prev_freq = hap_tracker.get(hap); 
				hap_tracker.put(hap, prev_freq + this.hap_freq[h]);
			}
		}
		int new_hap = hap_tracker.size(); 
		int[][] new_var_comp = new int[new_hap][new_loci];
		double[] new_freq = new double[new_hap]; 
		int hap_index = 0; 
		for (String hap : hap_tracker.keySet()) {
			String[] tmp2 = hap.split("");
			for (int l = 0; l < new_loci; l++) new_var_comp[hap_index][l] = Integer.parseInt(tmp2[l]); 
			new_freq[hap_index] = hap_tracker.get(hap);
			hap_index++;
		}
		return new SimpleConfig(new_hap, new_loci, new_var_comp, new_freq, new_positions);
	}
	
	public void write_global_file_string(String global_hap_output_file, boolean append){
		try{
			BufferedWriter bw= new BufferedWriter(new FileWriter(global_hap_output_file, append));
			bw.write("Hap_ID");
			for(int h=0;h<this.num_hap;h++)
				bw.write("\t"+h);
			bw.write("\nFreq");
			for(int h=0;h<this.num_hap;h++)
				bw.write("\t"+this.hap_freq[h]);
			bw.write("\n");
			for(int l=0;l<this.num_loci;l++){
				bw.write("0;" + locus_positions[l] + ";" + locus_positions[l] + ";0:1");
				for(int h=0;h<this.num_hap;h++)
					bw.write("\t"+this.hap_var_comp[h][l]);
				bw.write("\n");
			}bw.close();
		}catch(Exception e){e.printStackTrace();}
	}

}
