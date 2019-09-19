package PoolHap;

import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Random;
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
    
    
    
    public double freq_dist (double [] x, double [] y) {
    	double n =(double) x.length;
    	double value=0;
    	for (int i=0;i< x.length;i++) {
    		value += Math.abs(x[i]-y[i]);
    	}
    	return value/ n;
    }
    
    public double [] var_freq_cal (String [][] hap, double [] hap_freq ) throws IOException {
    	double [] var_freq =new double [hap[0].length ];
    	for (int i=0;i< var_freq.length; i++) {
    		double value =0;
    		for (int j=0; j< hap.length; j++) {
    			if (hap[j][i].equals("1")) {
    				value =value+ hap_freq[j];
    			}
    		}
    		var_freq[i]=value;
    	}
    	return var_freq;
    }
    
    public double [] var_freq_cal2 ( double [][] site_freq ) throws IOException {
    	double [] var_freq =new double [site_freq.length ];
    	for (int i=0;i< var_freq.length; i++) {
    		double value =0;
    		for (int j=0; j< site_freq[i].length; j++) {
    			value =value+ site_freq[i][j];

    		}
    		var_freq[i]=value/ (double)site_freq[i].length ;
    	}
    	return var_freq;
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

        int iter=0;
        int diagonal_matrix_freq_index =0;
        double[] diagonal_matrix_freq =new double [] {0.0001, 0.001, 0.01, 0.1, 1.0, 2.0, 5.0};
//        for (iter = 0; iter < this.num_loci*100; iter++) {
        while ((iter< this.num_loci*100)  &&
        		(diagonal_matrix_freq_index < diagonal_matrix_freq.length ) )  {

        	double total_sigma =0;
        	for (int j = 0; j < this.initial_Haps.sigma.length; j++) {
        		for (int k = 0; k < this.initial_Haps.sigma[j].length; k++) {
        			total_sigma+= this.initial_Haps.sigma[j][k];
        		}
        	}
        	double total_sigma_cufoff= 0.02* (double)this.num_loci* (double)this.num_loci;
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

            
            RealMatrix svd_inv = svd.getSolver().getInverse();
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
                
                // QUESTION: is this.aem_parameters.est_ind_pool making rh1 a lot smaller than it
                // has to be?
                RealMatrix dist_mtx = MatrixUtils.createRealMatrix(dist_to_mu);
                // rh2 = exp(-dist_to_mu %*% svd_inv %*% (omega-h)
                //     - diag(dist_to_mu %*% svd_inv %*% t(dist_to_mu) / 2))

                
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
            
            freq = p_new.clone();
            
            for (double f : freq) {
                if (Double.isNaN(f)) {
                	double[] r_f=new double [freq.length] ;
                	double total_r= 0;
                	for (int j = 0; j < freq.length; j++) {
    	        	    double s = 1.0;
    	        	    r_f[j]=s;
    	        	    total_r +=s;
                	}
                	for (int j = 0; j < freq.length; j++) {
                		freq[j]=(double) r_f[j]/ (double) total_r;
                	}
                }
            }
            if (delta < aem_parameters.aem_epsilon || delta1 < aem_parameters.aem_epsilon) {
            	double [] predicted_var_freq =new double [ this.initial_Haps.num_loci];
            	double [] real_var_freq =new double [ this.initial_Haps.num_loci];
            	predicted_var_freq= var_freq_cal(this.initial_Haps.global_haps_string,freq);
            	real_var_freq =var_freq_cal2(this.initial_Haps.inpool_site_freqs);
            	if (freq_dist( predicted_var_freq,real_var_freq )< 0.12) {
            		break;
            	}
            }
            
            
            iter++;
            double total_freq=0;
//            if (iter< -10) {
            if (iter== this.num_loci*100) {
            	iter=0;
            	double [] predicted_var_freq =new double [ this.initial_Haps.num_loci];
            	double [] real_var_freq =new double [ this.initial_Haps.num_loci];
            	predicted_var_freq= var_freq_cal(this.initial_Haps.global_haps_string,freq);
            	real_var_freq =var_freq_cal2(this.initial_Haps.inpool_site_freqs);
//            	for (int j=0;j<predicted_var_freq.length;j++ ) {
//            		System.out.print(predicted_var_freq[j]+" " );
//            	}
//            	System.out.println();
//            	System.out.println(freq_dist( predicted_var_freq,real_var_freq )+"***");
            	if (freq_dist( predicted_var_freq,real_var_freq ) < 0.12) {
            		this.initial_Haps.global_haps_freq = freq.clone();
            		break;
            	}
            	else {
            		for (int j=0;j<freq.length;j++ ) {
            			if (num_of_alternate(this.initial_Haps.global_haps_string[j])==1) {
            				freq[j]= diagonal_matrix_freq[diagonal_matrix_freq_index];
            				total_freq+= diagonal_matrix_freq[diagonal_matrix_freq_index];
            			}else {
            				freq[j]= 1.0;
            				total_freq+=1.0;
            			}
            		}
            		diagonal_matrix_freq_index++;
            		for (int j = 0; j < freq.length; j++) {
                    	freq[j] = freq[j]/ total_freq;
                    }
            	}
            }
            
            this.initial_Haps.global_haps_freq = freq.clone();
            
// Chen: Ensure a non-singular matrix   
            total_freq=1.0; 
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
            this.initial_Haps.update_sigma_mu_logL();
            
        }
        
        
        

        
        
        for (double f : freq) {
            if (Double.isNaN(f)) {
                failure = true;
            }
        }
        // TODO: [ReconEP]:: the below shoudl be split into multiple helpers?
        if (!failure) {
//        	for (int j = 0; j < this.initial_Haps.global_haps_string.length; j++) {
//        		String tmp="";
//        		for (int k = 0; k < this.initial_Haps.global_haps_string[j].length; k++) {
//        			tmp=tmp+ this.initial_Haps.global_haps_string[j][k];
//        		}
//        		System.out.println( tmp+" "+ this.initial_Haps.global_haps_freq[j]);
//        	}
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
