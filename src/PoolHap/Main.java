package PoolHap; 

import java.io.BufferedReader;
import java.io.FileReader;
import java.io.PrintWriter;
import java.sql.Date;
import java.text.SimpleDateFormat;
import java.util.ArrayList;

public class Main {

	public static void main(String[] args) throws Exception{

		// for (String s : args) {
			// System.out.println(s);
		// }
		String working_folder = args[0];
		String file_veflist = args[1];
		String file_data = args[2]; 
		String file_pos = args[19];
		int num_pools = Integer.parseInt(args[18]);

		double epsilon=Double.parseDouble(args[3]); 
		double rare_cutoff=Double.parseDouble(args[4]);
		int N=Integer.parseInt(args[5]);
		int max_iteration_aem = Integer.parseInt(args[6]); 
		double alpha=Double.parseDouble(args[7]), gamma=Double.parseDouble(args[8]), beta_a=Double.parseDouble(args[9]), beta_c=Double.parseDouble(args[10]), 
			c_old=Double.parseDouble(args[11]), c_new=Double.parseDouble(args[12]), p_add=Double.parseDouble(args[13]);
		int max_iterations=Integer.parseInt(args[14]), burnin=Integer.parseInt(args[15]), coalescing_mismatch=Integer.parseInt(args[16]), 
			round=Integer.parseInt(args[17]);
		String minisat_path = args[20];

		System.out.println("************************************************************");
		System.out.println("\t\t\tPoolHap2.0");
		System.out.println("De novo Haplotype Reconstruction from Pooled Sequencing Data\n");
		System.out.println("************************************************************\n");

		Date startDate = new Date(System.currentTimeMillis());
		SimpleDateFormat sdf = new SimpleDateFormat("HH:mm:ss");
		String time = sdf.format(startDate);
		System.out.println("The current time is " + time + ".");
		
		System.out.println("Reading allele counts from file_data...\n");
		BufferedReader br1 = new BufferedReader(new FileReader(file_data)); 
		String currLine = br1.readLine();
		String[] tmpVar = currLine.split("\t");
		ArrayList<Double> tmpPos = new ArrayList<Double>();
		for(String s : tmpVar) tmpPos.add(Double.parseDouble(s));
		int numVars = tmpPos.size(); 
		double[][] data = new double[num_pools][numVars]; 
		for (int i = 0; i < numVars; i++) data[0][i] = tmpPos.get(i) * N;
		int currPool = 0;
		currLine = br1.readLine();
		while(currLine != null) {
			currPool++; 
			String[] tmpVar2 = currLine.split("\t");
			// System.out.println(currPool + "\t" + tmpVar2[0]);
			for (int i = 0; i < numVars; i++) {
				data[currPool][i] = Double.parseDouble(tmpVar2[i]) * N;
			}
			currLine = br1.readLine();
		}
		br1.close();
		/*
		System.out.println("Normalized allele counts as pool x variant position...");
		for (int p = 0; p < num_pools; p++) {
			for (int i = 0; i < numVars; i++) System.out.print(data[p][i] + " ");
			System.out.println(); 
		}
		*/

		System.out.println("Reading variant positions from file_pos...");
		BufferedReader br2 = new BufferedReader(new FileReader(file_pos));
		String[] tmpArr = br2.readLine().split("\t");
		int[] posArray = new int[tmpArr.length];
		for (int i = 0; i < tmpArr.length; i++) {
			posArray[i] = Integer.parseInt(tmpArr[i]);
		}
		br2.close();
		System.out.println();
		
		System.out.println("Guessing initial haplotypes from the VEF files using the graph-colouring algorithm...\n"); 
		PoolAnalyzer pool2=new PoolAnalyzer(data, posArray, N, file_veflist, 1.0, minisat_path);
		System.out.println("There are " + pool2.curr_Haps.num_curr_H + " initial haplotypes spanning " + pool2.curr_Haps.num_snps + " variant positions.");
		pool2.curr_Haps.print_stdout();
		Date gcDate = new Date(System.currentTimeMillis());
		long finSection = gcDate.getTime() - startDate.getTime(); 
		int seconds = (int) (finSection / 1000) % 60 ;
		int minutes = (int) ((finSection / (1000*60)) % 60);
		int hours = (int) (((finSection / (1000*60*60)) % 60));
		System.out.println("Graph colouring has taken " + hours + " hours, " + minutes + " minutes, and " + seconds + " seconds.");
		
		System.out.println("\nRefining initial haplotypes and their frequencies using the rjMCMC algorithm..."); 
		pool2.population_freq_rjmcmc_multirun(alpha, beta_a, beta_c,  gamma, c_old, c_new, 
			p_add, max_iterations, burnin, coalescing_mismatch, round, working_folder, rare_cutoff);
		System.out.println("There are " + pool2.best_Haps.num_curr_H + " refined haplotypes spanning " + pool2.best_Haps.num_snps + " variant positions. These haplotypes pass the frequency cutoff of " + rare_cutoff + ".");
		pool2.best_Haps.print_stdout(rare_cutoff);
		Date rjDate = new Date(System.currentTimeMillis());
		finSection = rjDate.getTime() - gcDate.getTime(); 
		seconds = (int) (finSection / 1000) % 60 ;
		minutes = (int) ((finSection / (1000*60)) % 60);
		hours = (int) (((finSection / (1000*60*60)) % 60));
		System.out.println("rjMCMC has taken " + hours + " hours, " + minutes + " minutes, and " + seconds + " seconds.");

		System.out.println("\nCalculating the inter- and intra-pool frequencies of refined haplotypes using the AEM algorithm... "); 
		pool2.population_freq_aem(data, epsilon, rare_cutoff, max_iteration_aem);
		Date emDate = new Date(System.currentTimeMillis());
		finSection = emDate.getTime() - rjDate.getTime(); 
		seconds = (int) (finSection / 1000) % 60 ;
		minutes = (int) ((finSection / (1000*60)) % 60);
		hours = (int) (((finSection / (1000*60*60)) % 60));
		System.out.println("AEM has taken " + hours + " hours, " + minutes + " minutes, and " + seconds + " seconds.");
		
		System.out.print("\nPrinting results. ");
		PrintWriter pw = new PrintWriter(working_folder + "/p.all.results");
		
		pw.append("************************************************************\n");
		pw.append("\t\t\tPoolHap2.0\n");
		pw.append("De novo Haplotype Reconstruction from Pooled Sequencing Data\n\n");
		pw.append("************************************************************\n\n");
		
		pw.append("Input Files and Parameters: \n");
		pw.append("VEF List: " + file_veflist + "\n");
		pw.append("Variant Count File: " + file_data + "\n");
		pw.append("Variant Position File: " + file_pos + "\n");
		
		pw.append("num pools = " + num_pools + "\t");
		pw.append("est num ind per pool = " + N + "\t");
		pw.append("epsilon = " + epsilon + "\t");
		pw.append("rarity cutoff = " + rare_cutoff + "\n\n");

		pw.append("rjMCMC Algorithm Parameters: \n");		
		pw.append("num rounds rjmcmc = " + round + "\t");
		pw.append("num iter rjmcmc = " + max_iterations + "\t");
		pw.append("num burnin rjmcmc = " + burnin + "\t");
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
		pool2.best_Haps.print(pw);
		
		pw.append("\n");
		
		pw.append("+======BEST_INTRAPOOL=======+\n");
		pool2.print_intrapool(pw);
		pw.append("\n");
		
		// pw.append("+======All candidates======+\n");
		// pool2.print_buffer(rare_cutoff, pw);
		pw.close();
		startDate = new Date(System.currentTimeMillis());
		time = sdf.format(startDate);
		System.out.print("Finished. The current time is " + time + ".");
		pw.close();
	}	
}
