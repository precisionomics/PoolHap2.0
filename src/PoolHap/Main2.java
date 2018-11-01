package PoolHap;

import java.io.BufferedReader;
import java.io.FileReader;
import java.io.PrintWriter;
import java.sql.Date;
import java.text.SimpleDateFormat;
import java.util.ArrayList;

public class Main2 {

	public static void main(String[] args) throws Exception{

		String working_folder = args[0];
		String file_veflist = args[1];
		String file_data = args[2]; 
		String file_pos = args[19];	// Should be args[3]
		int num_pools = Integer.parseInt(args[18]);	// Should be args[4]
		int est_ind =Integer.parseInt(args[5]);
		String minisat_path = args[20];	// Should be args[6]

		int max_iterations=Integer.parseInt(args[14]), burn_in=Integer.parseInt(args[15]), coalescing_mismatch=Integer.parseInt(args[16]), 
				num_round=Integer.parseInt(args[17]);
		double alpha=Double.parseDouble(args[7]), gamma=Double.parseDouble(args[8]), beta_a=Double.parseDouble(args[9]), beta_c=Double.parseDouble(args[10]), 
			c_old=Double.parseDouble(args[11]), c_new=Double.parseDouble(args[12]), p_add=Double.parseDouble(args[13]);
		double rare_cutoff=Double.parseDouble(args[4]);
		int max_iteration_aem = Integer.parseInt(args[6]); 
		double epsilon=Double.parseDouble(args[3]); 
		
		System.out.println("************************************************************");
		System.out.println("\t\t\tPoolHapX");
		System.out.println("De novo Haplotype Reconstruction from Pooled Sequencing Data\n");
		System.out.println("************************************************************\n");

		Date startDate = new Date(System.currentTimeMillis());
		SimpleDateFormat sdf = new SimpleDateFormat("HH:mm:ss");
		String time = sdf.format(startDate);
		System.out.println("The current time is " + time + ".\n");

		System.out.println("Reading intra-pool allele frequencies from file_data...\n");
		BufferedReader br1 = new BufferedReader(new FileReader(file_data)); 
		String currLine = br1.readLine();	// Skip the pool ID line.
		ArrayList<double[]> tmp_intra_freq = new ArrayList<double[]>();
		currLine = br1.readLine();	// The first locus' alternate allele frequency.
		int currLocus = 0; 
		while(currLine != null) {	
			String[] tmp_locus = currLine.split("\t"); 
			tmp_intra_freq.add(new double[num_pools]); 
			for (int p = 0; p < num_pools; p++) tmp_intra_freq.get(currLocus)[p] = Double.valueOf(tmp_locus[p + 1]);
			currLine = br1.readLine();	
			currLocus++;
		} 
		br1.close();
		int numVars = currLocus;
		double[][] data = new double[num_pools][numVars]; 
		for (int v = 0; v < numVars; v++) 
			for (int p = 0; p < num_pools; p++) 
				data[p][v] = tmp_intra_freq.get(v)[p]; 

		System.out.println("Reading variant positions from file_pos...\n");
		BufferedReader br2 = new BufferedReader(new FileReader(file_pos));
		String[] tmpArr = br2.readLine().split("\t");
		int[] posArray = new int[tmpArr.length];
		for (int i = 0; i < tmpArr.length; i++) {
			posArray[i] = Integer.parseInt(tmpArr[i]);
		}
		br2.close();

		System.out.println("Guessing initial haplotypes from the VEF files using the graph-colouring algorithm...\n");
		PoolSolver pool=new PoolSolver(data, posArray, est_ind, file_veflist, minisat_path);
		System.out.println("There are " + pool.initial_Haps.num_global_hap + " initial haplotypes spanning " + pool.initial_Haps.num_loci + " variant positions.");
		pool.initial_Haps.write_global_stdout(); 
		Date gcDate = new Date(System.currentTimeMillis());
		long finSection = gcDate.getTime() - startDate.getTime(); 
		int seconds = (int) (finSection / 1000) % 60 ;
		int minutes = (int) ((finSection / (1000*60)) % 60);
		int hours = (int) (((finSection / (1000*60*60)) % 60));
		System.out.println("Graph colouring has taken " + hours + " hours, " + minutes + " minutes, and " + seconds + " seconds.");

		System.out.println("\nRefining initial haplotypes and their frequencies using the rjMCMC algorithm..."); 
		pool.population_freq_rjmcmc_multirun(num_round, max_iterations, burn_in, alpha, beta_a, beta_c, gamma, 
				c_old, c_new, p_add, coalescing_mismatch, rare_cutoff);
		System.out.println("There are " + pool.best_Haps.num_global_hap + " refined haplotypes spanning " + pool.best_Haps.num_loci + " variant positions. These haplotypes pass the frequency cutoff of " + rare_cutoff + ".");
		pool.best_Haps.write_global_stdout(); 
		Date rjDate = new Date(System.currentTimeMillis());
		finSection = rjDate.getTime() - gcDate.getTime(); 
		seconds = (int) (finSection / 1000) % 60 ;
		minutes = (int) ((finSection / (1000*60)) % 60);
		hours = (int) (((finSection / (1000*60*60)) % 60));
		System.out.println("rjMCMC has taken " + hours + " hours, " + minutes + " minutes, and " + seconds + " seconds.");

		System.out.println("\nCalculating the inter- and intra-pool frequencies of refined haplotypes using the AEM algorithm... "); 
		pool.population_freq_aem(epsilon, rare_cutoff, max_iteration_aem);
		Date emDate = new Date(System.currentTimeMillis());
		finSection = emDate.getTime() - rjDate.getTime(); 
		seconds = (int) (finSection / 1000) % 60 ;
		minutes = (int) ((finSection / (1000*60)) % 60);
		hours = (int) (((finSection / (1000*60*60)) % 60));
		System.out.println("AEM has taken " + hours + " hours, " + minutes + " minutes, and " + seconds + " seconds.");
		
		System.out.print("\nPrinting results. ");
		PrintWriter pw = new PrintWriter(working_folder + "/p.all.results");
		
		pw.append("************************************************************\n");
		pw.append("\t\t\tPoolHapX\n");
		pw.append("De novo Haplotype Reconstruction from Pooled Sequencing Data\n\n");
		pw.append("************************************************************\n\n");
		
		pw.append("Input Files and Parameters: \n");
		pw.append("VEF List: " + file_veflist + "\n");
		pw.append("Variant Count File: " + file_data + "\n");
		pw.append("Variant Position File: " + file_pos + "\n");
		
		pw.append("num pools = " + num_pools + "\t");
		pw.append("est num ind per pool = " + est_ind + "\t");
		pw.append("epsilon = " + epsilon + "\t");
		pw.append("rarity cutoff = " + rare_cutoff + "\n\n");

		pw.append("rjMCMC Algorithm Parameters: \n");		
		pw.append("num rounds rjmcmc = " + num_round + "\t");
		pw.append("num iter rjmcmc = " + max_iterations + "\t");
		pw.append("num burnin rjmcmc = " + burn_in + "\t");
		pw.append("alpha = " + alpha + "\t");
		pw.append("gamma = " + gamma + "\t");
		pw.append("beta a = " + beta_a + "\t");
		pw.append("beta c = " + beta_c + "\t");
		pw.append("c old = " + c_old + "\t");
		pw.append("c new = " + c_new + "\t");
		pw.append("p add = " + p_add + "\t");
		pw.append("coal mismatch = " + coalescing_mismatch + "\n\n");

		pw.append("AEM Algorithm Parameters: \n");		
		pw.append("num iter aem = " + max_iteration_aem + "\n\n");
				
		pw.append("+======BEST_INTERPOOL=======+\n");
		pool.best_Haps.write_global_file_string(working_folder + "/p.all.results", true);
		
		pw.append("\n");
		
		pw.append("+======BEST_INTRAPOOL=======+\n");
		pool.best_Haps.write_inpool(working_folder + "/p.all.results", true);
		pw.append("\n");
		
		pw.close();
		startDate = new Date(System.currentTimeMillis());
		time = sdf.format(startDate);
		System.out.print("Finished. The current time is " + time + ".");
		pw.close();
	}
	
}
