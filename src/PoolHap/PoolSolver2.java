package PoolHap;

import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import org.apache.commons.math3.distribution.BetaDistribution;
import org.apache.commons.math3.linear.MatrixUtils;
import org.apache.commons.math3.linear.RealMatrix;
import org.apache.commons.math3.linear.SingularValueDecomposition;

import PoolHap.Parameters.*;
import PoolHap.HapConfig;

public class PoolSolver2 {
	
	// Redundant variables for convenience
	public int num_snp;					// number of SNPs in this region
	public int num_pools;				// number of pools under study	
	// Msin structures
	HapConfig initial_Haps; 
	HapConfig[] final_Haps; 			// The best haplotype set for each pool.
	// Solver parameters
	McmcParameters mcmc_parameters;
	BetaDistribution new_old_haps_beta;	// This is for suggesting proportions of haplotypes to section off for new haplotypes.

	/*
	 * The constructor for the PoolSolver object i.e.: the core algorithm. 
	 * @param Estimated variant frequencies, positions of the primitive loci, list of each pool's VEFs. 
	 * @return To set up and adjust variant composition and inter/intra-pool frequencies of curr_Haps. 
	 */
	public PoolSolver2(int num_snp, int num_pools, HapConfig hap_config, String parameter_file) throws Exception{
		this.num_snp = num_snp;
		this.num_pools = num_pools;
		this.initial_Haps = hap_config;
		this.initial_Haps.solver = this;	// TODO This doesn't pass down into the this.initial_Haps clones of rjMCMC! So we're not calculating stuff correctly at all!
		this.final_Haps = new HapConfig[this.num_pools]; 
		analyze_a_region_aem(parameter_file); 
		analyze_a_region_mcmc(parameter_file);
	}

	/*
	 * The main method for running EM on haplotypes to estimate their global frequency. Inspired by the AEM function from Zhang et al., 2008.
	 * @param this.initial_Haps, the set of possible sub-haplotypes in a specific region. aem_parameters include est_ind_pool, epsilon, rare_cutoff.
	 * @return Updated guess of the overall (between-pool frequencies of the sub-haplotypes) frequencies i.e.: updated this.initial_Haps object. 
	 */
	public void analyze_a_region_aem(String parameter_file) throws IOException{	
		AemParameters aem_parameters = new AemParameters(parameter_file);
		double[] p = this.initial_Haps.global_haps_freq.clone(); // The current estimate of global haplotype frequencies.
		double[] Rh = p.clone();	// The importance factor of each haplotype is the previous iteration's estimate of its global frequency.
		/* double rh1_dist_to_mu = 0; 
		double rh2_0 = 0; 
		double rh2_mean = 0;
		double rh_dist_to_mu = 0; 
		double[] rh2array = new double[this.num_pools];
		double[] diffarray = new double[this.num_pools]; */ 
		// Step 1) Find the (approximate) maximum likelihood global frequency of each haplotype by applying i) linear constraints on allele frequencies and ii) pairwise haplotype frequencies.
		this.initial_Haps.update_sigma_mu_logL();
		AEM: for(int iter = 0; iter < aem_parameters.max_iteration; iter++) {	// For each ieration of AEM...
			// Step 1a) Calculate ii) using the sigma matrix, which represents linkage between the alleles at different sites. 
			if(iter % 10 == 0) System.out.println("\n" + iter+":"+this.initial_Haps.logL); 	// REPORT
			SingularValueDecomposition svd = new SingularValueDecomposition(MatrixUtils.createRealMatrix(this.initial_Haps.sigma));	// Calculate the inverse singular value decomposition of the haplotype set's sigma (covariance) matrix.
			RealMatrix svd_inv = svd.getSolver().getInverse();
			System.out.println("Inverse sigma matrix:");
			double[][] tmp = svd_inv.getData(); 
			for (int r = 0; r < svd_inv.getRowDimension(); r++) {
				for (int c = 0 ; c < svd_inv.getColumnDimension(); c++) {
					if (Double.isNaN(tmp[r][c]))  {
						System.out.println(iter + ": Inverse sigma row " + r + " column " + c + " is NaN.");
						break AEM;
					}
					System.out.printf("%.5f\t",tmp[r][c]); 
					if (Double.isNaN(tmp[r][c])) System.out.print("svd_inv[" + r + "][" + c + "] is 0. ");
				}
				System.out.println();
			}
			// Step 1b) Calculate i) by taking the distance of alternate alleles from the average haplotype of each pool. 
			System.out.println("\nDistance to the average haplotype matrix:");
			double[][] dist_to_mu = new double[this.num_pools][this.num_snp];			   //dist_to_mu=sweep(A/(2*N),2,omega, "-");
			for(int u = 0; u < this.num_pools; u++){
				// System.out.println(this.initial_Haps.inpool_site_freqs[u].length);
				for(int v = 0; v < this.num_snp; v++){
					dist_to_mu[u][v] = this.initial_Haps.inpool_site_freqs[v][u] - this.initial_Haps.mu[v];
					System.out.printf("%.5f\t", dist_to_mu[u][v]); 
					// if(iter%100==0 && j==0) System.out.println("dist_to_mu[" + u + "][" + v + "] = " + dist_to_mu[u][v] + "\t");
					// if (dist_to_mu[u][v] == 0)	System.out.print("dist_to_mu[" + u + "][" + v + "] is 0. ");
				}
				// if(iter%100==0 && j==0) System.out.println();
				System.out.println();
			}
			System.out.println("rh2 arrays:");
			double[] IF = new double[this.initial_Haps.num_global_hap]; 
			for (int j = 0; j < this.initial_Haps.num_global_hap; j++){
				double[] hh = this.initial_Haps.global_haps[j];	 
				// rh1=exp(-1/(4*N)*t(omega-h)%*%svd_inv%*%(omega-h))
				double rh1 = Math.exp(-1/(4.0*aem_parameters.est_ind_pool) * Algebra.quadratic_form(Algebra.minus(this.initial_Haps.mu, hh), svd_inv.getData()));	// The approximation of the importance factor of each haplotype, according to Appendix B of Zhang et al. 2008. 
				// QUESTION Is this.aem_parameters.est_ind_pool making rh1 a lot smaller than it has to be? 
				// if(iter%100==0 && j==0) System.out.println("j:" + j + "\trh1=" + rh1);
				RealMatrix dist_mtx = MatrixUtils.createRealMatrix(dist_to_mu);// rh2=exp(-dist_to_mu%*%svd_inv%*%(omega-h)-diag( dist_to_mu%*%svd_inv%*%t(dist_to_mu)/2))
				double[] rh2 = Algebra.exp(Algebra.minus(svd_inv.preMultiply(dist_mtx).operate(Algebra.minus(hh, this.initial_Haps.mu)),
						Algebra.diag(dist_mtx.multiply(svd_inv).multiply(dist_mtx.transpose()).scalarMultiply(0.5).getData())));
				for (double tmp_rh2 : rh2) System.out.printf("%.5f\t", tmp_rh2);
				System.out.println();
				double rh = rh1*Algebra.mean(rh2);
				/* if(j==0) {
					rh1_dist_to_mu = rh1; 
					rh2_mean = Algebra.mean(rh2); 
					for (int i = 0; i < rh2.length; i++) {
						rh2array[i] = rh2[i];
						diffarray[i] = (double) hh[i] - mu[i];
					}
					rh_dist_to_mu = rh; 
				} */
				IF[j] = rh; //   IF=c(IF, rh)
				/* for (int u = 0; u < this.num_pools; u++) {
					// System.out.println(j);
					this.initial_Haps.in_pool_haps_freq[j][u] = rh1 * rh2[u];
				} */
			}
			double delta = Algebra.sum(Algebra.abs(Algebra.minus(Rh, IF)));//sum(abs(Rh-IF))
			double delta1 = Algebra.sum(Algebra.abs(Algebra.add(Algebra.divide(Rh, IF),-1)));//sum(abs(IF/Rh-1)) 
			Rh = IF;
			double[] p_new= Algebra.times(Rh,p);
			if(iter%1==0) {
				System.out.printf(iter + "\tdelta = %.5f\tdelta1 = %.5f\n", delta, delta1);
				// System.out.printf("Hap_0\trh1 = %.5f\tmean_rh2 = %.5f\trh = %.5f\n", rh1_dist_to_mu, rh2_mean, rh_dist_to_mu);
				// for (double rh2dist_to_mu : rh2array) System.out.printf("%.5f\t", rh2dist_to_mu);
				// System.out.println();
				System.out.print("IF:\t");
				for (double impfac : IF) System.out.printf("%.3f\t", impfac);
				System.out.println();
				System.out.print("p_pre:\t");
				for (double freq : p_new) System.out.printf("%.3f\t", freq);
				System.out.println();
			}
			Algebra.normalize_ditribution(p_new);
			Algebra.rmlow_and_normalize(p_new, aem_parameters.rare_cutoff);	

			if(iter%1==0) System.out.print("p_new:\t");
			if(iter%1==0) for (double freq : p_new) System.out.printf("%.3f\t", freq);
			if(iter%1==0) System.out.println();
			
			if (delta < aem_parameters.epsilon || delta1 < aem_parameters.epsilon) {
				System.out.println(delta + " < " + aem_parameters.epsilon + " || " +  delta1 + " < " + aem_parameters.epsilon); 
				break;		
			}
			p=p_new.clone(); 
			this.initial_Haps.global_haps_freq = p.clone(); 
			this.initial_Haps.update_sigma_mu_logL();
		} 
		System.out.print("p:\t");
		for (double freq : p) System.out.printf("%.3f\t", freq);
		System.out.println();

		boolean[] list_rem_haps = new boolean[this.initial_Haps.num_global_hap];
		int num_rem_hap = 0; 
		for (int h = 0; h < this.initial_Haps.num_global_hap; h++)
			if (p[h] == 0) {
				list_rem_haps[h] = true; 
				num_rem_hap++; 
			} else this.initial_Haps.global_haps_freq[h] = p[h]; 
		this.initial_Haps.remHaps(list_rem_haps, num_rem_hap, false);
		System.out.println("\n" + num_rem_hap + " haplotypes were removed due to low frequency after AEM.");
		
		/* for (int u = 0; u < this.num_pools; u++) {
			double poolCount = 0; 
			for (int j = 0; j < this.initial_Haps.num_global_hap; j++) {
				this.initial_Haps.in_pool_haps_freq[j][u] = this.initial_Haps.in_pool_haps_freq[j][u] * p[j];
				poolCount += this.initial_Haps.in_pool_haps_freq[j][u]; 
			}
			for (int j = 0; j < this.initial_Haps.num_global_hap; j++) this.initial_Haps.in_pool_haps_freq[j][u] = this.initial_Haps.in_pool_haps_freq[j][u] / poolCount;
		} */
		System.out.println("\nThe pools are composed of " + this.initial_Haps.num_global_hap + " haplotypes.");
	}

	/*
	 * The main method for running multi-round rjMCMC on haplotypes to refine their variant compositions and within-pool frequencies. 
	 * @param hap_config, the set of possible sub-haplotypes in a specific region. mcmc_parameters make the running of a reverse-jump MCMC procedure possible. 
     * @return Refined guess of within-pool haplotype variant compositions and their within-pool frequencies. 
	 */
	public void analyze_a_region_mcmc(String parameter_file) throws IOException{
		this.mcmc_parameters = new McmcParameters(parameter_file);
		this.new_old_haps_beta=new BetaDistribution(mcmc_parameters.c_old, mcmc_parameters.c_new);	
		for (int pool = 0; pool < this.num_pools; pool++) {
			ArrayList<HapConfig> best_Haps_buffer=new ArrayList<HapConfig>();
			for(int round=0; round < mcmc_parameters.num_round;round++){
				HapConfig curr_Haps = this.initial_Haps.clone(0, false);
				curr_Haps.checkrank_and_fullfill(mcmc_parameters.freqs_sum, pool);
				System.out.println("Starting MCMC round " + round + " for pool " + pool + ": logL = " + curr_Haps.logL + "\n");
				best_Haps_buffer.add(population_freq_rjmcmc(curr_Haps, pool));
			}
			int best_index=0;
			if (mcmc_parameters.num_round == 0) best_index = 0;
			else for(int round=1; round<mcmc_parameters.num_round;round++)
				if (best_Haps_buffer.get(round).logL > best_Haps_buffer.get(best_index).logL)
					best_index=round;
			this.final_Haps[pool] = best_Haps_buffer.get(best_index).clone(mcmc_parameters.rare_cutoff, true);
		}
	}
	
	/*
	 * The main method for running single-round rjMCMC. Inspired by the hippo program in Pirinen 2009.
	 * @param Various parameters for estimating parameters from a Hastings-within-Gibbs procedure. 
	 * @return Prposal for a refined guess of inter-pool haplotypes. 
	 */
	public HapConfig population_freq_rjmcmc(HapConfig curr_Haps, int pool){
		int accept_freq=0, accept_substit=0, accept_additions=0,accept_deletions=0; // ,accept_recomb=0,reject_recomb=0; also had rejection before
		try{
			
			System.out.println("Burning in for " + mcmc_parameters.burn_in + " iterations...\n");
			for (int iter = 0; iter < mcmc_parameters.burn_in; iter++) {	// These are the burn-in iterations. Nothing will be recorded here...
				curr_Haps.update_freqs(mcmc_parameters.beta_a, mcmc_parameters.beta_c, mcmc_parameters.alpha, iter, pool); // !!!
				curr_Haps.mutate_a_hap(pool);
				curr_Haps.add_mutant_or_coalesce(mcmc_parameters.alpha, mcmc_parameters.p_add, mcmc_parameters.gamma, mcmc_parameters.coalescing_mismatch, pool);
			}
			System.out.println("Running for " + mcmc_parameters.max_iteration + " iterations...\n");
			for(int iter=0;iter< mcmc_parameters.max_iteration; iter++){
				accept_freq += curr_Haps.update_freqs(mcmc_parameters.beta_a, mcmc_parameters.beta_c, mcmc_parameters.alpha, iter, pool); // !!!  		
				accept_substit += curr_Haps.mutate_a_hap(pool);
				int add_or_del=curr_Haps.add_mutant_or_coalesce(mcmc_parameters.alpha, mcmc_parameters.p_add, mcmc_parameters.gamma, mcmc_parameters.coalescing_mismatch, pool);
				if(add_or_del==1) accept_additions++;
				else if(add_or_del==2) accept_deletions++;
			}
		}catch(Exception e){e.printStackTrace();}	 
		System.out.println("Number of haplotypes = " + curr_Haps.num_global_hap + "\tLog-likelihood = "+ curr_Haps.logL + "\nAcceptance Ratios:");
		System.out.println("Frequency Change: " + (double) accept_freq / mcmc_parameters.max_iteration); 
		System.out.println("Single-locus Mutation: " + (double) accept_substit / mcmc_parameters.max_iteration); 
		System.out.println("Haplotype Addition: " + (double) accept_additions / mcmc_parameters.max_iteration); 
		System.out.println("Haplotype Deletion: " + (double) accept_deletions / mcmc_parameters.max_iteration + "\n"); 
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
		return new HapConfig(tmp_global_string, tmp_global_freq, tmp_in_pool_freq, this.initial_Haps.inpool_site_freqs, 
				this.initial_Haps.locusInfo, this.num_pools, final_IDs, this.initial_Haps.pool_IDs, this.initial_Haps.est_ind_pool);
	}
}