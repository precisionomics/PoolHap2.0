package PoolHap;

import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.concurrent.ThreadLocalRandom;
import org.apache.commons.math3.linear.MatrixUtils;
import org.apache.commons.math3.linear.RealMatrix;
import org.apache.commons.math3.linear.SingularValueDecomposition;
import PoolHap.Parameters.AemParameters;

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
        this.initial_Haps.solver = this;
        analyze_a_region_aem(parameter_file);
        // analyze_a_region_mcmc(parameter_file);
    }

    /**
     *  The main method for running EM on haplotypes to estimate their global frequency. Inspired by
     *  the AEM function from Zhang et al., 2008.
     *
     *  @param this.initial_Haps, the set of possible sub-haplotypes in a specific region.
     *  aem_parameters include est_ind_pool, epsilon, rare_cutoff.
     *
     *  @return Updated guess of the overall (between-pool frequencies of the sub-haplotypes)
     *  frequencies i.e.: updated this.initial_Haps object.
     */
    public void analyze_a_region_aem(String parameter_file) throws IOException {
        AemParameters aem_parameters = new AemParameters(parameter_file);

        // TODO: [LEFTOVER]
        // PrintWriter stdout = new PrintWriter(
        //     new FileWriter("/home/lmak/Documents/v0.7_test/aem_I_0.txt"));

        // The current estimate of global haplotype frequencies.
        double[] freq = this.initial_Haps.global_haps_freq.clone();

        // TODO: [LEFTOVER]
        // stdout.write("Initial global frequencies:\n");
        // for (double i : p) {
        //     stdout.format("%.3f\t", i);
        // }
        // stdout.write("\n");

        // The importance factor of each haplotype is the previous iteration's estimate of its
        // global frequency.
        double[] Rh = freq.clone();
        this.initial_Haps.update_sigma_mu_logL();

        // TODO: [LEFTOVER]
        // stdout.write("Initial mu:\n");
        // for (int v = 0; v < this.initial_Haps.num_loci; v++) {
        //     stdout.format("%.3f\t", this.initial_Haps.mu[v]);
        // }
        // stdout.write("\nInitial sigma:\n");
        // for (int v = 0; v < this.initial_Haps.num_loci; v++) {
        //     for (int w = 0; w < this.initial_Haps.num_loci; w++) {
        // 	       stdout.format("%.3f\t", this.initial_Haps.sigma[v][w]);
        //     }
        //     stdout.format("\n");
        // }


        /*
         *  Step 1: Find the (approximate) maximum likelihood global frequency of each haplotype by
         *  applying i) linear constraints on allele frequencies and ii) pairwise haplotype
         *  frequencies.
         */
        // For each AEM iteration...
        for (int iter = 0; iter < aem_parameters.max_iteration; iter++) {
            // Step 1a) Calculate ii) using the sigma matrix, which represents linkage between the
            // alleles at different sites.
            // TODO: [ReconEP]:: maybe each of these steps should be a helper.

            // TODO: [LEFTOVER]
            // REPORT
            // stdout.write("\nIteration "
            //     + iter
            //     + " initial logL = "
            //     + this.initial_Haps.logL
            //     + "\n");

            // Calculate the inverse singular value decomposition of the haplotype set's sigma
            // (covariance) matrix.
            SingularValueDecomposition svd = new SingularValueDecomposition(
                MatrixUtils.createRealMatrix(this.initial_Haps.sigma));

            RealMatrix svd_inv = svd.getSolver().getInverse();

            // TODO: [LEFTOVER]
            // stdout.write("\nInverse sigma matrix:\n");
            // double[][] tmp = svd_inv.getData();
            // for (int r = 0; r < svd_inv.getRowDimension(); r++) {
            //     for (int c = 0; c < svd_inv.getColumnDimension(); c++) {
            //         if (Double.isNaN(tmp[r][c])) {
            //            stdout.write(iter
            //                + ": Inverse sigma row "
            //                + r
            //                + " column "
            //                + c
            //                + " is NaN.");
            //
            // 		   }
            //         stdout.format("%.3f\t",tmp[r][c]);
            //     }
            //     stdout.write("\n");
            // }

            // TODO: [LEFTOVER]
            // stdout.write("\nDistance to the average haplotype matrix:\n");

            // Step 1b) Calculate i) by taking the distance of alternate alleles from the average
            // haplotype of each pool.
            // TODO: [ReconEP]:: extract to helper?
            double[][] dist_to_mu = new double[this.num_pools][this.num_loci];

            // TODO: [LEFTOVER]
            // dist_to_mu = sweep(A / (2*N), 2, omega, "-");

            for (int u = 0; u < this.num_pools; u++) {
                for (int v = 0; v < this.num_loci; v++) {
                    dist_to_mu[u][v] = this.initial_Haps.inpool_site_freqs[v][u]
                        - this.initial_Haps.mu[v];

                    // TODO: [LEFTOVER]
                    // stdout.format("%.3f\t", dist_to_mu[u][v]);

                }

                // TODO: [LEFTOVER]
                // stdout.write("\n");

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

                // TODO: [LEFTOVER]
                // if (j % 50 == 0) {
                //     stdout.write("\nFor regional haplotype "
                //        + j
                //        + " of iteration "
                //        + iter
                //        + ":\n");
                //
                //     stdout.format("rh1 = %.5f\nrh2 array:\t", rh1);
                // }
                //
                // QUESTION: is this.aem_parameters.est_ind_pool making rh1 a lot smaller than it
                // has to be?
                // if (iter % 100 == 0 && j == 0) {
                //     System.out.println("j:" + j + "\trh1=" + rh1);
                // }

                RealMatrix dist_mtx = MatrixUtils.createRealMatrix(dist_to_mu);

                // TODO: [LEFTOVER]
                // rh2 = exp(-dist_to_mu %*% svd_inv %*% (omega-h)
                //     - diag(dist_to_mu %*% svd_inv %*% t(dist_to_mu) / 2))
                //
                // double[] first_term = svd_inv.preMultiply(dist_mtx)
                //     .operate(Algebra.minus(hh, this.initial_Haps.mu));
                //
                // double[] second_term = Algebra.diag(dist_mtx.multiply(svd_inv)
                //         .multiply(dist_mtx.transpose())
                //         .scalarMultiply(0.5)
                //         .getData());
                //
                // double[] rh2 = new double[first_term.length];
                // for (int r = 0; r < rh2.length; r++) {
                //     rh2[r] = Math.exp(first_term[r] - second_term[r]);
                // }

                double[] rh2 = Algebra.exp(
                    Algebra.minus(svd_inv.preMultiply(dist_mtx)
                        .operate(Algebra.minus(hh, this.initial_Haps.mu)),
                    Algebra.diag(dist_mtx.multiply(svd_inv)
                        .multiply(dist_mtx.transpose())
                        .scalarMultiply(0.5)
                        .getData())));

                // TODO: [LEFTOVER]
                // if (j % 50 == 0) {
                //     for (double tmp_rh2 : rh2) {
                //         stdout.format("%.5f\t", tmp_rh2);
                //     }
                //     stdout.write("\n");
                //     if (j % 1 == 0) {
                //         double[] first_term = svd_inv.preMultiply(dist_mtx)
                //                 .operate(Algebra.minus(hh, this.initial_Haps.mu));
                //
                //         double[] second_term = Algebra.diag(dist_mtx.multiply(svd_inv)
                //                 .multiply(dist_mtx.transpose())
                //                 .scalarMultiply(0.5)
                //                 .getData());
                //
                //         for (int r = 0; r < rh2.length; r++) {
                //             double exponent = Math.exp(first_term[r] - second_term[r]);
                //             System.out.println(first_term[r]
                //                 + "\t"
                //                 + second_term[r]
                //                 + "\t"
                //                 + exponent);
                //
                //         }
                //         System.out.println();
                //     }
                // }

                double rh = rh1 * Algebra.mean(rh2);
                IF[j] = rh;
            }
            double delta = Algebra.sum(Algebra.abs(Algebra.minus(Rh, IF))); // sum(abs(Rh - IF)

            // sum(abs(IF / Rh - 1))
            double delta1 = Algebra.sum(Algebra.abs(Algebra.add(Algebra.divide(Rh, IF), -1)));
            Rh = IF;
            double[] p_new = Algebra.times(Rh, freq);

            // TODO: [LEFTOVER]
            // stdout.format(
            //     "\ndelta = %.5f\t\tdelta1 = %.5f\nImportance factor (rh * mean(rh2)):\n",
            //     delta,
            //     delta1);
            //
            // for (double impfac : IF) {
            //    stdout.format("%.3f\t", impfac);
            // }
            // stdout.write("\n");

            Algebra.normalize_ditribution(p_new);
            Algebra.rmlow_and_normalize(p_new, aem_parameters.rare_cutoff);

            // TODO: [LEFTOVER]
            // stdout.write("\nIteration " + iter + " updated global frequencies:\n");
            // for (double freq : p_new) {
            //     stdout.format("%.3f\t", freq);
            // }
            // stdout.write("\n");

            if (delta < aem_parameters.epsilon || delta1 < aem_parameters.epsilon) {

                // TODO: [LEFTOVER]
                // stdout.write(delta
                //     + " < "
                //     + aem_parameters.epsilon
                //     + " || "
                //     + delta1
                //     + " < "
                //     + aem_parameters.epsilon);

                break;
            }

            freq = p_new.clone();
            this.initial_Haps.global_haps_freq = freq.clone();
            this.initial_Haps.update_sigma_mu_logL();

            // TODO: [LEFTOVER]
            // stdout.write("Updated mu:\n");
            // for (int v = 0; v < this.initial_Haps.num_loci; v++) {
            //     stdout.format("%.3f\t", this.initial_Haps.mu[v]);
            // }
            // stdout.write("\nUpdated sigma:\n");
            // for (int v = 0; v < this.initial_Haps.num_loci; v++) {
            //     for (int w=0; w < this.initial_Haps.num_loci; w++) {
            //         stdout.format("%.3f\t", this.initial_Haps.sigma[v][w]);
            //     }
            //     stdout.write("\n");
            // }

        }

        // TODO: [LEFTOVER]
        // stdout.write("\nFinal AEM global frequencies:\n");
        // for (double freq : p) {
        //     stdout.format("%.3f\t", freq);
        // }
        // stdout.write("\n");

        for (double f : freq) {
            if (Double.isNaN(f)) {
                failure = true;
            }
        }

        // TODO: [ReconEP]:: the below shoudl be split into multiple helpers?
        if (!failure) {
            boolean[] list_rem_haps = new boolean[freq.length];
            int num_rem_hap = 0;
            double actual_cutoff = aem_parameters.final_cutoff;
            for (int h = 0; h < freq.length; h++) {
                if (freq[h] < actual_cutoff) {
                    list_rem_haps[h] = true;
                    num_rem_hap++;

                } else {
                    this.initial_Haps.global_haps_freq[h] = freq[h];
                }
            }

            // If too many of them are below the regional frequency minimum...
            if (num_rem_hap > freq.length - aem_parameters.hapset_size_max) {
                list_rem_haps = new boolean[freq.length];
                num_rem_hap = 0;
                double[] haps_freq_copy = new double[freq.length];	// deep instead of shallow copy
                for (int f = 0; f < freq.length; f++) {
                    haps_freq_copy[f] = freq[f];
                }

                Arrays.sort(haps_freq_copy);

                // Take the haplotypes that are the top 5 frequencies or more common.
                actual_cutoff = haps_freq_copy[haps_freq_copy.length
                    - aem_parameters.adhoc_freq_cutoff];

                // Don't want too many because GC gets confused.
                for (int h = 0; h < freq.length; h++) {
                    if (freq[h] < actual_cutoff) {
                        list_rem_haps[h] = true;
                        num_rem_hap++;
                    } else {
                        this.initial_Haps.global_haps_freq[h] = freq[h];
                    }
                }
            }

            // ...and then if all of them are about the same low frequency...
            if (num_rem_hap < freq.length - aem_parameters.hapset_size_min) {
                list_rem_haps = new boolean[freq.length];
                num_rem_hap = 0;

                // AEM starts with 2^l haplotypes, so we need to keep much less than 20%.
                ArrayList<Integer> haps_to_keep = new ArrayList<Integer>();
                for (int h = 0; h < aem_parameters.hapset_size_rand; h++) {
                    haps_to_keep.add((int) ThreadLocalRandom.current().nextDouble() * freq.length);
                }
                for (int h = 0; h < freq.length; h++) {
                    if (!haps_to_keep.contains(h)) {
                        list_rem_haps[h] = true;
                        num_rem_hap++;
                    } else {
                        this.initial_Haps.global_haps_freq[h] = freq[h];
                    }
                }
                actual_cutoff = Double.NaN;
            }

            System.out.println("Of the "
                + freq.length
                + " regional haplotypes, "
                + num_rem_hap
                + " were removed. The frequency cutoff was "
                + actual_cutoff
                + ".");

            this.initial_Haps.remHaps(list_rem_haps, num_rem_hap);
            this.final_Haps = this.initial_Haps;
        }

        // TODO: [LEFTOVER]
        // stdout.write("\n"
        //     + num_rem_hap
        //     + " haplotypes were removed due to zeroed frequency after AEM.");
        //
        // System.out.println("\nThe pools are composed of "
        //     + this.initial_Haps.num_global_hap
        //     + " haplotypes.");
        //
        // stdout.close();

    }
}
