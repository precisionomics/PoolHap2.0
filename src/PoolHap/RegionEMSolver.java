package PoolHap;

import java.io.IOException;
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
    public void analyze_a_region_aem(String parameter_file) throws IOException {
        Parameters aem_parameters = new Parameters(parameter_file);

        // The current estimate of global haplotype frequencies.
        double[] freq = this.initial_Haps.global_haps_freq.clone();
        // The importance factor of each haplotype is the previous iteration's estimate of its
        // global frequency.
        double[] Rh = freq.clone();
        this.initial_Haps.update_sigma_mu_logL();
        /*
         *  Step 1: Find the (approximate) maximum likelihood global frequency of each haplotype by
         *  applying i) linear constraints on allele frequencies and ii) pairwise haplotype
         *  frequencies.
         */
        // For each AEM iteration...
        for (int iter = 0; iter < aem_parameters.aem_max_iteration; iter++) {
            // Step 1a) Calculate ii) using the sigma matrix, which represents linkage between the
            // alleles at different sites.
            // TODO: [ReconEP]:: maybe each of these steps should be a helper.

            // Calculate the inverse singular value decomposition of the haplotype set's sigma
            // (covariance) matrix.
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
            if (delta < aem_parameters.aem_epsilon || delta1 < aem_parameters.aem_epsilon) {
                break;
            }
            freq = p_new.clone();
            this.initial_Haps.global_haps_freq = freq.clone();
            this.initial_Haps.update_sigma_mu_logL();
        }

        for (double f : freq) {
            if (Double.isNaN(f)) {
                failure = true;
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
