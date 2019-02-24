package PoolHap;

import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;
import org.apache.commons.math3.distribution.BetaDistribution;
import org.apache.commons.math3.linear.MatrixUtils;
import org.apache.commons.math3.linear.RealMatrix;
import org.apache.commons.math3.linear.SingularValueDecomposition;

import PoolHap.Parameters.*;
import PoolHap.HapConfig;

public class PoolSolver2_Testing {

	// Redundant variables for convenience
	public int num_snp; // number of SNPs in this region
	public int num_pools; // number of pools under study
	// Main structures
	public HapConfig initial_Haps;
	HapConfig[] best_Haps_buffer;
	HapConfig best_Haps;
	double logL_best;
	// Solver parameters
	McmcParameters mcmc_parameters;
	BetaDistribution new_old_haps_beta; // This is for suggesting proportions of haplotypes to section off for new
										// haplotypes.

	/*
	 * The constructor for the PoolSolver object i.e.: the core algorithm.
	 * 
	 * @param Estimated variant frequencies, positions of the primitive loci, list
	 * of each pool's VEFs.
	 * 
	 * @return To set up and adjust variant composition and inter/intra-pool
	 * frequencies of curr_Haps.
	 */
	public PoolSolver2_Testing(int num_snp, int num_pools, HapConfig hap_config, String parameter_file)	throws Exception {
		this.num_snp = num_snp;
		this.num_pools = num_pools;
		this.initial_Haps = hap_config;
		this.initial_Haps.solver = this; // TODOFIXED This doesn't pass down into the this.initial_Haps clones of
											// rjMCMC! So we're not calculating stuff correctly at all!
		analyze_a_region_aem(parameter_file);
		// analyze_a_region_mcmc(parameter_file);
	}

	/*
	 * The main method for running EM on haplotypes to estimate their global
	 * frequency. Inspired by the AEM function from Zhang et al., 2008.
	 * 
	 * @param this.initial_Haps, the set of possible sub-haplotypes in a specific
	 * region. aem_parameters include est_ind_pool, epsilon, rare_cutoff.
	 * 
	 * @return Updated guess of the overall (between-pool frequencies of the
	 * sub-haplotypes) frequencies i.e.: updated this.initial_Haps object.
	 */
	public void analyze_a_region_aem(String parameter_file) throws IOException{	
		AemParameters aem_parameters = new AemParameters(parameter_file);
// 		PrintWriter stdout = new PrintWriter(new FileWriter("/home/lmak/Documents/v0.7_test/aem_I_0.txt"));
		double[] p = this.initial_Haps.global_haps_freq.clone(); // The current estimate of global haplotype frequencies.
// 		stdout.write("Initial global frequencies:\n"); 
// 		for (double i : p) 	stdout.format("%.3f\t", i);
// 		stdout.write("\n"); 
		double[] Rh = p.clone();	// The importance factor of each haplotype is the previous iteration's estimate of its global frequency.
		this.initial_Haps.update_sigma_mu_logL();

// 		stdout.write("Initial mu:\n"); 
//		for(int v=0;v<this.initial_Haps.num_loci;v++)
// 			stdout.format("%.3f\t", this.initial_Haps.mu[v]);
// 		stdout.write("\nInitial sigma:\n");
//		for(int v=0;v<this.initial_Haps.num_loci;v++) {
//			for(int w=0;w<this.initial_Haps.num_loci;w++)
// 				stdout.format("%.3f\t", this.initial_Haps.sigma[v][w]);
// 			stdout.format("\n");
//		} 
		// Step 1) Find the (approximate) maximum likelihood global frequency of each haplotype by applying i) linear constraints on allele frequencies and ii) pairwise haplotype frequencies.
		for(int iter = 0; iter < aem_parameters.max_iteration; iter++) {	// For each ieration of AEM... 
			// Step 1a) Calculate ii) using the sigma matrix, which represents linkage between the alleles at different sites. 
// 			// stdout.write("\nIteration " + iter + " initial logL = " + this.initial_Haps.logL + "\n"); 	// REPORT
			SingularValueDecomposition svd = new SingularValueDecomposition(MatrixUtils.createRealMatrix(this.initial_Haps.sigma));	// Calculate the inverse singular value decomposition of the haplotype set's sigma (covariance) matrix.
			RealMatrix svd_inv = svd.getSolver().getInverse();
// 			stdout.write("\nInverse sigma matrix:\n");
// 			double[][] tmp = svd_inv.getData(); 
//			for (int r = 0; r < svd_inv.getRowDimension(); r++) {
//				for (int c = 0 ; c < svd_inv.getColumnDimension(); c++) {
//					if (Double.isNaN(tmp[r][c]))  {
// 						stdout.write(iter + ": Inverse sigma row " + r + " column " + c + " is NaN.");
//					}
// 					stdout.format("%.3f\t",tmp[r][c]); 
//				}
// 				stdout.write("\n");
//			}
			// Step 1b) Calculate i) by taking the distance of alternate alleles from the average haplotype of each pool. 
// 			stdout.write("\nDistance to the average haplotype matrix:\n");
			double[][] dist_to_mu = new double[this.num_pools][this.num_snp];			   //dist_to_mu=sweep(A/(2*N),2,omega, "-");
			for(int u = 0; u < this.num_pools; u++){
				for(int v = 0; v < this.num_snp; v++){
					dist_to_mu[u][v] = this.initial_Haps.inpool_site_freqs[v][u] - this.initial_Haps.mu[v];
// 					stdout.format("%.3f\t", dist_to_mu[u][v]); 
				}
// 				stdout.write("\n");
			}
			double[] IF = new double[this.initial_Haps.num_global_hap]; 
			for (int j = 0; j < this.initial_Haps.num_global_hap; j++){
				double[] hh = this.initial_Haps.global_haps[j];	 
				// rh1=exp(-1/(4*N)*t(omega-h)%*%svd_inv%*%(omega-h))
				double rh1 = Math.exp(-1/(2.0 * aem_parameters.est_ind_pool) * Algebra.quadratic_form(Algebra.minus(this.initial_Haps.mu, hh), svd_inv.getData()));	// The approximation of the importance factor of each haplotype, according to Appendix B of Zhang et al. 2008.
				// if (j % 50 == 0) {
// 					stdout.write("\nFor regional haplotype " + j + " of iteration " + iter + ":\n");
// 					stdout.format("rh1 = %.5f\nrh2 array:\t", rh1);
				// }
				// QUESTION Is this.aem_parameters.est_ind_pool making rh1 a lot smaller than it has to be? 
				// if(iter%100==0 && j==0) System.out.println("j:" + j + "\trh1=" + rh1);
				RealMatrix dist_mtx = MatrixUtils.createRealMatrix(dist_to_mu);// rh2=exp(-dist_to_mu%*%svd_inv%*%(omega-h)-diag( dist_to_mu%*%svd_inv%*%t(dist_to_mu)/2))
				double[] first_term = svd_inv.preMultiply(dist_mtx).operate(Algebra.minus(hh, this.initial_Haps.mu));
				double[] second_term = Algebra.diag(dist_mtx.multiply(svd_inv).multiply(dist_mtx.transpose()).scalarMultiply(0.5).getData());
				double[] rh2 = new double[first_term.length]; 
				for (int r = 0; r < rh2.length; r++) {
					rh2[r] = Math.exp(first_term[r] - second_term[r]);
				}				
				// double[] rh2 = Algebra.exp(Algebra.minus(svd_inv.preMultiply(dist_mtx).operate(Algebra.minus(hh, this.initial_Haps.mu)),
						// Algebra.diag(dist_mtx.multiply(svd_inv).multiply(dist_mtx.transpose()).scalarMultiply(0.5).getData())));
				// if (j % 50 == 0) {
// 					for (double tmp_rh2 : rh2) stdout.format("%.5f\t", tmp_rh2);
// 					stdout.write("\n");
//					if (j % 1 == 0) {
						// double[] first_term = svd_inv.preMultiply(dist_mtx).operate(Algebra.minus(hh, this.initial_Haps.mu));
						// double[] second_term = Algebra.diag(dist_mtx.multiply(svd_inv).multiply(dist_mtx.transpose()).scalarMultiply(0.5).getData());
//						for (int r = 0; r < rh2.length; r++) {
//							double exponent = Math.exp(first_term[r] - second_term[r]); 
							// System.out.println(first_term[r] + "\t" + second_term[r] + "\t" + exponent);						
//						}
						// System.out.println();
//					}
				// }
				double rh = rh1*Algebra.mean(rh2);
				IF[j] = rh; 
			}
			double delta = Algebra.sum(Algebra.abs(Algebra.minus(Rh, IF)));//sum(abs(Rh-IF))
			double delta1 = Algebra.sum(Algebra.abs(Algebra.add(Algebra.divide(Rh, IF),-1)));//sum(abs(IF/Rh-1)) 
			Rh = IF;
			double[] p_new= Algebra.times(Rh,p);
// 			stdout.format("\ndelta = %.5f\t\tdelta1 = %.5f\nImportance factor (rh * mean(rh2)):\n", delta, delta1);
// 			for (double impfac : IF) stdout.format("%.3f\t", impfac);
// 			stdout.write("\n");
			Algebra.normalize_ditribution(p_new);
			Algebra.rmlow_and_normalize(p_new, aem_parameters.rare_cutoff);	

// 			stdout.write("\nIteration " + iter + " updated global frequencies:\n");
// 			for (double freq : p_new) stdout.format("%.3f\t", freq);
// 			stdout.write("\n");
			
			if (delta < aem_parameters.epsilon || delta1 < aem_parameters.epsilon) {
// 				stdout.write(delta + " < " + aem_parameters.epsilon + " || " +  delta1 + " < " + aem_parameters.epsilon); 
				break;		
			}
			p=p_new.clone(); 
			this.initial_Haps.global_haps_freq = p.clone(); 
			this.initial_Haps.update_sigma_mu_logL();
			
// 			stdout.write("Updated mu:\n"); 
//			for(int v=0;v<this.initial_Haps.num_loci;v++)
// 				stdout.format("%.3f\t", this.initial_Haps.mu[v]);
// 			stdout.write("\nUpdated sigma:\n");
//			for(int v=0;v<this.initial_Haps.num_loci;v++) {
//				for(int w=0;w<this.initial_Haps.num_loci;w++)
// 					stdout.format("%.3f\t", this.initial_Haps.sigma[v][w]);
// 				stdout.write("\n");
//			}
		} 
// 		/* stdout.write("\nFinal AEM global frequencies:\n");
// 		for (double freq : p) stdout.format("%.3f\t", freq);
// 		stdout.write("\n"); */

		boolean[] list_rem_haps = new boolean[this.initial_Haps.num_global_hap];
		int num_rem_hap = 0; 
		for (int h = 0; h < this.initial_Haps.num_global_hap; h++)
			if (p[h] < aem_parameters.final_cutoff) {
				list_rem_haps[h] = true; 
				num_rem_hap++; 
			} else this.initial_Haps.global_haps_freq[h] = p[h]; 
		this.initial_Haps.remHaps(list_rem_haps, num_rem_hap);
// 		stdout.write("\n" + num_rem_hap + " haplotypes were removed due to zeroed frequency after AEM.");
		
		// System.out.println("\nThe pools are composed of " + this.initial_Haps.num_global_hap + " haplotypes.");
// 		stdout.close();
	}

/*	
	 * The main method for running multi-round rjMCMC on haplotypes to refine their variant compositions and global frequencies. 
	 * @param Various parameters for estimating parameters from a Hastings-within-Gibbs procedure. 
     * @return Refined guess of inter-pool haplotypes in best_Haps. 
	 
	public void analyze_a_region_mcmc(String parameter_file) throws IOException{
// 		PrintWriter stdout = new PrintWriter(new FileWriter("/home/lmak/Documents/v0.7_test/mcmc_I_0.txt"));
		this.mcmc_parameters = new McmcParameters(parameter_file);
		this.new_old_haps_beta=new BetaDistribution(this.mcmc_parameters.c_old, mcmc_parameters.c_new);	
		this.best_Haps_buffer = new HapConfig[this.mcmc_parameters.num_round]; 
// 		stdout.println("There are " + this.initial_Haps.num_global_hap + " initial haplotypes.\n");
		for (int round = 0; round < this.mcmc_parameters.num_round; round++){
			HapConfig curr_Haps = this.initial_Haps.clone(0);
			curr_Haps.solver = this;	// TODO make sure that the update functions that require the beta distribution actually access this.
			curr_Haps.checkrank_and_fullfill(mcmc_parameters.freqs_sum);
			System.out.println("Starting RJMCMC round " + round + ": logL = " + curr_Haps.logL + "\nThere are " + curr_Haps.num_global_hap + " initial haplotypes.");
// 			stdout.write("Starting RJMCMC round " + round + ": logL = " + curr_Haps.logL + "\nThere are " + curr_Haps.num_global_hap + " initial haplotypes.\n\n");
// 			this.best_Haps_buffer[round] = single_round_rjmcmc(curr_Haps, stdout);
		}
		int best_index = 0;
		if (mcmc_parameters.num_round != 0)
			for (int round = 1; round<mcmc_parameters.num_round;round++)
				if (this.best_Haps_buffer[round].logL > this.best_Haps_buffer[best_index].logL && this.best_Haps_buffer[round].logL < 0) best_index = round;
		// for (int h = 0; h < this.best_Haps_buffer[best_index].num_global_hap; h++)
			// System.out.println(this.best_Haps_buffer[best_index].hap_IDs[h]); 
		this.best_Haps = this.best_Haps_buffer[best_index].clone(mcmc_parameters.rare_cutoff);
		this.logL_best = best_Haps.logL;
// 		stdout.close();
	}
	
	
	 * The main method for running single-round rjMCMC. Inspired by the hippo program in Pirinen 2009.
	 * @param Various parameters for estimating parameters from a Hastings-within-Gibbs procedure. 
	 * @return Proposal for a refined guess of inter-pool haplotypes. 
	 
// 	public HapConfig single_round_rjmcmc(HapConfig curr_Haps, PrintWriter stdout){
		int accept_freq=0, accept_substit=0, accept_additions=0,accept_deletions=0; // ,accept_recomb=0,reject_recomb=0; also had rejection before
		try{
// 			stdout.write("Burning in for " + mcmc_parameters.burn_in + " iterations...\n");
			for (int iter = 0; iter < mcmc_parameters.burn_in; iter++) {	// These are the burn-in iterations. Nothing will be recorded here...
				curr_Haps.update_freqs(mcmc_parameters.beta_a, mcmc_parameters.beta_c, mcmc_parameters.alpha, iter); // !!!
				curr_Haps.mutate_a_hap();
				curr_Haps.add_mutant_or_coalesce(mcmc_parameters.alpha, mcmc_parameters.p_add, mcmc_parameters.gamma, mcmc_parameters.coalescing_mismatch);
			}
			System.out.println("Burnin done. Running for " + mcmc_parameters.max_iteration + " iterations...\n");
// 			stdout.write("Running for " + mcmc_parameters.max_iteration + " iterations...\n");
			for(int iter=0;iter< mcmc_parameters.max_iteration; iter++){
				accept_freq += curr_Haps.update_freqs(mcmc_parameters.beta_a, mcmc_parameters.beta_c, mcmc_parameters.alpha, iter); // !!!  		
				accept_substit += curr_Haps.mutate_a_hap();
				int add_or_del=curr_Haps.add_mutant_or_coalesce(mcmc_parameters.alpha, mcmc_parameters.p_add, mcmc_parameters.gamma, mcmc_parameters.coalescing_mismatch);
				if(add_or_del==1) accept_additions++;
				else if(add_or_del==2) accept_deletions++;
			}
		}catch(Exception e){e.printStackTrace();}	 
// 		stdout.write("Final number of haplotypes = " + curr_Haps.num_global_hap + "\t\tLog-likelihood = "+ curr_Haps.logL + "\n");
// 		stdout.write("Frequency Change: " + (double) accept_freq / mcmc_parameters.max_iteration + "\n"); 
// 		stdout.write("Single-locus Mutation: " + (double) accept_substit / mcmc_parameters.max_iteration + "\n"); 
// 		stdout.write("Haplotype Addition: " + (double) accept_additions / mcmc_parameters.max_iteration + "\n"); 
// 		stdout.write("Haplotype Deletion: " + (double) accept_deletions / mcmc_parameters.max_iteration + "\n\n"); 
		return curr_Haps;
	}

	public HapConfig report(){
		ArrayList<String> tmp_IDs = new ArrayList<String>(); 
		HashMap<String, ArrayList<Integer>> id2pools = new HashMap<String, ArrayList<Integer>>(); // Maps haplotypes to the (list of) pools they occur in. 
		for (int p = 0; p < this.num_pools; p++) {
			HapConfig haps_in = this.final_Haps[p]; 
			for (String curr_id : haps_in.hap_IDs) {
				if (!id2pools.containsKey(curr_id)) {
					id2pools.put(curr_id, new ArrayList<Integer>()); 
				}
				id2pools.get(curr_id).add(p); 
				tmp_IDs.add(curr_id);
			}
		}
		String[][] tmp_global_string = new String[tmp_IDs.size()][this.num_snp];
		double[][] tmp_in_pool_freq = new double[tmp_IDs.size()][this.num_pools]; 	// #global_hap x #pools
		double[] tmp_global_freq = new double[tmp_IDs.size()];	// #global_hap
		for (int h = 0; h < tmp_IDs.size(); h++) {
			String curr_ID = tmp_IDs.get(h);
			for (int p = 0; p < id2pools.get(curr_ID).size(); p++) {
				tmp_in_pool_freq[h][p] = id2pools.get(curr_ID).get(p);
				tmp_global_freq[h] += id2pools.get(curr_ID).get(p); 
			}
			double tmp = tmp_global_freq[h] / this.num_pools; 
			tmp_global_freq[h] = tmp; 
			for (int s = 0; s < this.num_snp; s++) 
				tmp_global_string[h][s] = tmp_IDs.get(h).charAt(s) + "";
		}
		String[] final_IDs = tmp_IDs.toArray(new String[tmp_IDs.size()]);
		// System.out.println(final_IDs);
		return new HapConfig(tmp_global_string, tmp_global_freq, tmp_in_pool_freq, this.initial_Haps.inpool_site_freqs, 
				this.initial_Haps.locusInfo, this.num_pools, final_IDs, this.initial_Haps.pool_IDs, this.initial_Haps.est_ind_pool);
	}*/
}