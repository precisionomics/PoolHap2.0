package PoolHap;

import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.concurrent.ThreadLocalRandom;
import org.apache.commons.math3.linear.MatrixUtils;
import org.apache.commons.math3.linear.RealMatrix;
import org.apache.commons.math3.linear.SingularValueDecomposition;
import PoolHap.Parameters;

public class RegionEMSolver {
    // Redundant variables for convenience.
    public int num_loci; // number of SNPs in this region
    public int num_pools; // number of pools under study

    // Main structures.
    public HapConfig initial_Haps;
    public boolean failure = false;
    public HapConfig final_Haps;

    /**
     *  The constructor for the PoolSolver object i.e.: the core algorithm.
     *
     *  @param Estimated variant frequencies, positions of the primitive loci, list of each pool's
     *  VEFs.
     *
     *  @return To set up and adjust variant composition and inter/intra-pool frequencies of
     *  curr_Haps.
     */
    public RegionEMSolver(HapConfig hap_config, String parameter_file) throws Exception {
        this.initial_Haps = hap_config;
        this.num_loci = this.initial_Haps.num_loci;
        this.num_pools = this.initial_Haps.num_pools;
        this.initial_Haps.regional_EM_solver = this;
        analyze_a_region_aem(parameter_file);
        // analyze_a_region_mcmc(parameter_file);
    }

    /**
     *  The main method for running EM on haplotypes to estimate their global frequency. Inspired by
     *  the AEM function from Kuk et al., 2009.
     *
     *  @param this.initial_Haps, the set of possible sub-haplotypes in a specific region.
     *  aem_parameters include est_ind_pool, epsilon, rare_cutoff.
     *
     *  @return Updated guess of the overall (between-pool frequencies of the sub-haplotypes)
     *  frequencies i.e.: updated this.initial_Haps object.
     */
    
    public int num_of_alternate (String [] hap ) throws IOException {
    	int count= 0;
    	for (int i =0; i< hap.length;i++) {
    		if (hap[i].equals("1")) {
    			 count++;
    		}
    	}
    	return count;
    }
    public void analyze_a_region_aem(String parameter_file) throws IOException {
        Parameters aem_parameters = new Parameters(parameter_file);

        // The current estimate of global haplotype frequencies.
        double[] freq = this.initial_Haps.global_haps_freq.clone();
        
//        for (int i = 0; i  < freq.length; i++) {    		
//        	System.out.print(freq[i]  );
//        	System.out.print("\t" );
//        	for (int j = 0; j  < this.initial_Haps.global_haps_string[i].length; j++) { 
//        		System.out.print(this.initial_Haps.global_haps_string[i][j]  );
//        	}
//        	System.out.print("\n" );
//        }
        // The importance factor of each haplotype is the previous iteration's estimate of its
        // global frequency.
        double[] Rh = freq.clone();
        this.initial_Haps.update_sigma_mu_logL();
//        for (int i = 0; i  < this.initial_Haps.mu.length ;i++) {    	
//        	System.out.println( this.initial_Haps.mu[i]+"\t");
//        	
//        	for (int j = 0; j  < this.initial_Haps.sigma[i].length; j++) { 
//        		System.out.print( this.initial_Haps.sigma[i][j]+ "\t" );
//        	}
//        	System.out.println();
//        }
        /*
         *  Step 1: Find the (approximate) maximum likelihood global frequency of each haplotype by
         *  applying i) linear constraints on allele frequencies and ii) pairwise haplotype
         *  frequencies.
         */
        // For each AEM iteration...
        for (int iter = 0; iter < aem_parameters.aem_max_iteration; iter++) {
//        for (int iter = 0; iter < 2; iter++) {
        	
//    		String file_path= "/home/chencao/Desktop/test.txt";
//    		FileWriter mydata = new FileWriter(file_path,true);
//            PrintWriter pw = new PrintWriter(mydata);
//            pw.write(  Integer.toString(iter)+"\n");
//        	for (int i = 0; i  < freq.length; i++) {    		
//        		pw.write(Double.toString(freq[i])  );
//        		pw.write("\t" );
//            	for (int j = 0; j  < this.initial_Haps.global_haps_string[i].length; j++) { 
//            		pw.write(this.initial_Haps.global_haps_string[i][j]  );
//            	}
//            	pw.write("\n" );
//            }
//            pw.flush();
//            pw.close();
        	
            // Step 1a) Calculate ii) using the sigma matrix, which represents linkage between the
            // alleles at different sites.
            // TODO: [ReconEP]:: maybe each of these steps should be a helper.

            // Calculate the inverse singular value decomposition of the haplotype set's sigma
            // (covariance) matrix.
        	double total_sigma =0;
        	for (int j = 0; j < this.initial_Haps.sigma.length; j++) {
        		for (int k = 0; k < this.initial_Haps.sigma[j].length; k++) {
//        			System.out.print(this.initial_Haps.sigma[j][k]  );
        			total_sigma+= this.initial_Haps.sigma[j][k];
        		}
//        		System.out.println( );
        	}
        	
        	double total_sigma_cufoff= 0.02* (double)this.num_loci* (double)this.num_loci;
//        	total_sigma_cufoff =0;
//        	System.out.println(Double.toString(check)+"---" );
        	if (total_sigma< total_sigma_cufoff) {
        		double multiple = total_sigma_cufoff/ total_sigma;
        		for (int j = 0; j < this.initial_Haps.sigma.length; j++) {
            		for (int k = 0; k < this.initial_Haps.sigma[j].length; k++) {
            			this.initial_Haps.sigma[j][k] = this.initial_Haps.sigma[j][k]* multiple;
            		}
            	}
        	}
            SingularValueDecomposition svd = new SingularValueDecomposition(
                MatrixUtils.createRealMatrix(this.initial_Haps.sigma));
//            System.out.println( "**************************");
//            System.out.println( svd.getS());
            
            RealMatrix svd_inv = svd.getSolver().getInverse();

//            System.out.println(svd.getSolver().isNonSingular());            
//            System.out.println( svd_inv);
//            if (Double.isNaN(svd_inv.getColumn(0)[0])){
//            	System.exit(0);
//            }

            // Step 1b) Calculate i) by taking the distance of alternate alleles from the average
            // haplotype of each pool.
            // TODO: [ReconEP]:: extract to helper?
            double[][] dist_to_mu = new double[this.num_pools][this.num_loci];
            // dist_to_mu = sweep(A / (2*N), 2, omega, "-");
            for (int u = 0; u < this.num_pools; u++) {
                for (int v = 0; v < this.num_loci; v++) {
                    dist_to_mu[u][v] = this.initial_Haps.inpool_site_freqs[v][u]
                        - this.initial_Haps.mu[v];
                }
            }

            
            double[] IF = new double[this.initial_Haps.num_global_hap];
            for (int j = 0; j < this.initial_Haps.num_global_hap; j++) {
                double[] hh = this.initial_Haps.global_haps[j];
                // rh1 = exp(-1 / (4 * N) * t(omega - h) %*% svd_inv %*% (omega - h))
                // The approximation of the importance factor of each haplotype, according to
                // Appendix B of Zhang et al. 2008.
                double rh1 = Math.exp(-1 / (2.0 * aem_parameters.est_ind_pool)
                    * Algebra.quadratic_form(
                        Algebra.minus(this.initial_Haps.mu, hh),
                        svd_inv.getData()));
//                System.out.println("XXX");
                
                // QUESTION: is this.aem_parameters.est_ind_pool making rh1 a lot smaller than it
                // has to be?
                RealMatrix dist_mtx = MatrixUtils.createRealMatrix(dist_to_mu);
                // rh2 = exp(-dist_to_mu %*% svd_inv %*% (omega-h)
                //     - diag(dist_to_mu %*% svd_inv %*% t(dist_to_mu) / 2))

                
//                System.out.println(Double.toString( svd_inv.preMultiply(dist_mtx)
//                        .operate(Algebra.minus(hh, this.initial_Haps.mu))[0])+" @@"
//                		+Double.toString(Algebra.diag(dist_mtx.multiply(svd_inv)
//                        .multiply(dist_mtx.transpose())
//                        .scalarMultiply(0.5)
//                        .getData())[0])+" !!");
                
                double[] rh2 = Algebra.exp(
                    Algebra.minus(svd_inv.preMultiply(dist_mtx)
                        .operate(Algebra.minus(hh, this.initial_Haps.mu)),
                    Algebra.diag(dist_mtx.multiply(svd_inv)
                        .multiply(dist_mtx.transpose())
                        .scalarMultiply(0.5)
                        .getData())));

                double rh = rh1 * Algebra.mean(rh2);
                IF[j] = rh;
            }
            double delta = Algebra.sum(Algebra.abs(Algebra.minus(Rh, IF))); // sum(abs(Rh - IF)
            
            // sum(abs(IF / Rh - 1))
            double delta1 = Algebra.sum(Algebra.abs(Algebra.add(Algebra.divide(Rh, IF), -1)));
            Rh = IF;
            double[] p_new = Algebra.times(Rh, freq);

            Algebra.normalize_ditribution(p_new);
            Algebra.rmlow_and_normalize(p_new, aem_parameters.aem_zero_cutoff);
            if (delta < aem_parameters.aem_epsilon || delta1 < aem_parameters.aem_epsilon) {
                break;
            }
            freq = p_new.clone();
            this.initial_Haps.global_haps_freq = freq.clone();
            
// Chen: Ensure a non-singular matrix   
            double total_freq=1.0; 
            double min_freq_cufoff=0.0002; 
            
            for (int j = 0; j < this.initial_Haps.global_haps_freq.length; j++) {
            	if (num_of_alternate ( this.initial_Haps.global_haps_string[j]) ==1) {
            		if (this.initial_Haps.global_haps_freq[j]< min_freq_cufoff) {
            			this.initial_Haps.global_haps_freq[j]+= min_freq_cufoff;
            			total_freq+= min_freq_cufoff;
            		}
            	}
            }
            for (int j = 0; j < this.initial_Haps.global_haps_freq.length; j++) {
            	this.initial_Haps.global_haps_freq[j] = this.initial_Haps.global_haps_freq[j]
            			/total_freq;
            }
//            System.out.println("**");
//            for (int j = 0; j < this.initial_Haps.global_haps_freq.length; j++) {
//            	System.out.print(this.initial_Haps.global_haps_freq[j]);
//            }
//            System.out.println("");
            this.initial_Haps.update_sigma_mu_logL();
            
//            for (int j = 0; j < this.initial_Haps.sigma.length; j++) {
//        		for (int k = 0; k < this.initial_Haps.sigma[j].length; k++) {
//        			
//        			System.out.print(this.initial_Haps.sigma[j][k]);
//        		}
//        		System.out.println();
//        	}
        }

        for (double f : freq) {
            if (Double.isNaN(f)) {
                failure = true;
//                System.exit(0);
            }
        }

        // TODO: [ReconEP]:: the below shoudl be split into multiple helpers?
        if (!failure) {
            boolean[] list_remove_haps = new boolean[freq.length];
            int num_remove_hap = 0;
            for (int h = 0; h < freq.length; h++) {
                if (freq[h] < aem_parameters.aem_regional_cross_pool_freq_cutoff) {
                    list_remove_haps[h] = true;
                    num_remove_hap++;
                } else {
                    this.initial_Haps.global_haps_freq[h] = freq[h];
                }
            }
            double actual_cutoff = aem_parameters.aem_regional_cross_pool_freq_cutoff;
            // If too many of them are below the regional frequency minimum,
            // i.e., too few haps are remained 
            if (aem_parameters.aem_hapset_size_min > freq.length - num_remove_hap) {
                if(freq.length >= aem_parameters.aem_hapset_size_min) {
                    list_remove_haps = Algebra.permute_sort_and_remove(freq.clone(), 
                        aem_parameters.aem_hapset_size_min);
                }
            }
            // Or, if too many haps are remained...
            else if (aem_parameters.aem_hapset_size_max < freq.length - num_remove_hap) {
                list_remove_haps = Algebra.permute_sort_and_remove(freq.clone(), 
                    aem_parameters.aem_hapset_size_max);
            }
            num_remove_hap = 0;
            actual_cutoff=1;
            for(int h=0; h < freq.length; h++) {
                if(list_remove_haps[h]) {
                    num_remove_hap++;
                } else {
                    this.initial_Haps.global_haps_freq[h] = freq[h];
                    if(actual_cutoff > freq[h]) actual_cutoff = freq[h];
                }
            }
            this.initial_Haps.remHaps(list_remove_haps, num_remove_hap);
            System.out.println("Of the "
                + freq.length
                + " regional haplotypes, "
                + num_remove_hap
                + " were removed. The frequency cutoff was "
                + actual_cutoff
                + ".");           
            this.final_Haps = this.initial_Haps;
        }

    }
}
