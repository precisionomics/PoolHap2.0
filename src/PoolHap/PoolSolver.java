package PoolHap;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.util.ArrayList;
import java.util.Collection;
import java.util.HashMap;
import java.util.Set;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;

import org.apache.commons.math3.distribution.BetaDistribution;
import org.apache.commons.math3.linear.MatrixUtils;
import org.apache.commons.math3.linear.RealMatrix;
import org.apache.commons.math3.linear.SingularValueDecomposition;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

import PoolHap.GraphUV.SubGraph;

public class PoolSolver {
	
	static double _MAX_ = 1E20;
	double[][] data;	// Estimated variant frequencies from BAMFormatter.jar.
	int num_snp;		// Number of primitive loci (SNPs and short indels).
	int num_pool;
	int est_ind_pool;	// Estimated number of individuals in each pool.
	
	HapConfig initial_Haps; 
	HapConfig curr_Haps; 
	HapConfig best_Haps; 
	ArrayList<HapConfig> best_Haps_buffer=new ArrayList<HapConfig>();

	public double[] mu;
	public double[][] sigma;
	double logL_best;
	BetaDistribution new_old_haps_beta;	// This is for suggesting proportions of haplotypes to section off for new haplotypes.

	/*
	 * The constructor for the PoolSolver object i.e.: the core algorithm. 
	 * @param Estimated variant frequencies, positions of the primitive loci, list of each pool's VEFs. 
     * @return To set up and adjust variant composition and inter/intra-pool frequencies of curr_Haps. 
	 */
	public PoolSolver(double[][] data, int[] posArray, int est_ind, String vef_list, String minisat_path) throws Exception{
		this.est_ind_pool=est_ind;
		this.data=data;
		this.num_snp=data[0].length;
		System.out.println(this.num_snp);
		this.num_pool=data.length;
		this.mu=new double[this.num_snp];
		this.sigma=new double[this.num_snp][this.num_snp];
		this.initial_Haps = gc(vef_list, posArray, this.num_snp, this.num_pool, minisat_path);
		this.initial_Haps.solver = this;
	}

	/*
	 * The main method for graph colouring. Originally written by Jia Wang. 
	 * @param list of each pool's VEF file, positions of the primitive loci, path to the miniSAT program.  
     * @return Initial guess of inter-pool haplotypes in the curr_Haps HapConfig object. 
	 */
	public static HapConfig gc(String file_veflist, int[] posArray, int snps, int num_pool, String minisat_path) throws Exception	{
		HashMap<String,Double> knownHapSet = new HashMap<String,Double>();
		BufferedReader br = new BufferedReader(new FileReader(file_veflist)); 
		String currLine = br.readLine();
		int poolNum = 0; 
		while (currLine != null) {
			String[] currArray = currLine.split("\t");
			ExecutorService workers = getWorkers(Integer.parseInt(System.getProperty("reads.vc.workers", "0")));
			
			double sample = (currArray.length > 1)?
				Double.parseDouble(currArray[1]): 2;

			String method = (currArray.length > 2)? currArray[2]: "md";
			try
			{
				System.out.println("\nGuessing for pool number " + poolNum + "...");
				// knownHapSet.addAll(
				HashMap<String,Double> tmpHaps = run(new File(currArray[0]), posArray, sample, method, workers, minisat_path);
				for (String hapVars : tmpHaps.keySet()) {
					if (knownHapSet.containsKey(hapVars)) {
						double tmpFreq = knownHapSet.get(hapVars) + tmpHaps.get(hapVars);
						knownHapSet.put(hapVars, tmpFreq); 
					} else knownHapSet.put(hapVars, tmpHaps.get(hapVars)); 
				}
			} catch (Exception e){
				e.printStackTrace();
				System.out.println("The graph-colouring problem for pool number " + poolNum + " is unsatisfiable.");	
			} finally {
				workers.shutdownNow();
				poolNum++;
			}
			currLine = br.readLine();
		}
		br.close();
		String[][] initialHaps = new String[knownHapSet.size()][snps];
		double[] initialFreqs = new double[knownHapSet.size()]; 
		Set<String> knownHaps = knownHapSet.keySet();
		
		int currHap = 0;
	    for (String tmpKey : knownHaps) {
	    	String[] tmpHap = tmpKey.split(""); 
	    	for (int i = 0; i < tmpHap.length; i++) {
	    		initialHaps[currHap][i] = tmpHap[i];
	    	}
	    	initialFreqs[currHap] = (knownHapSet.get(tmpKey) / num_pool); 
	    	currHap++; 
	    }
	    
	    LocusAnnotation[] locusInfo = new LocusAnnotation[posArray.length]; 
	    for (int p = 0; p < posArray.length; p++) {
	    	locusInfo[p] = new LocusAnnotation(false, 0, posArray[p], posArray[p], new String[]{"0","1"}); 
	    }
	    HapConfig initialConfig = new HapConfig(initialHaps, initialFreqs, null, locusInfo, num_pool, null, null); 
		
	    return initialConfig;
	}
	
	/*
	 * The module selector for graph colouring. Can choose to colour with or without knowledge of existing haplotypes. 
	 * Originally written by Jia Wang. 
	 * @param each pool's VEF file, positions of the primitive loci, path to the miniSAT program.  
     * @return Initial guess of inter-pool haplotypes as a HashMap of strings (composition) and doubles (frequencies). 
	 */
	private static HashMap<String,Double> run(File input, int[] posArray, 
			double sample, String method, ExecutorService workers, String minisat_path) throws Exception
		{
			String name = (sample > 1)? input.getName():
				String.format("%s.%.0f", input.getName(), sample*100);
			
			Timer t = new Timer();
			
			ReadsDB db = new ReadsDB();
			db.load(input, sample);
			
			Collection<Subtype> goldens = db.getGoldens();
			if (!goldens.isEmpty())
			{
				db.showSubtypes("goldens", goldens, false);
				Subtype.saveSubtypes(
					new File(name+".goldens"), goldens);
			}
			
			ReadsGraph rG = new ReadsGraph(db, workers);
			GraphUV G = rG.getGraph();
			//writeGraph(new File(name+".graph"), G);
		
			int[] colors = null;

			if (method.equals("assist"))
			{
				GoldenDB gdb = new GoldenDB();
				gdb.load(new File(System.getProperty(
					"reads.vc.goldenDB")));
				
				MaxClique4 mc4 = new MaxClique4(G);
				int[] clique = mc4.solve();

				HyperGraphV ccs = G.bfsCC();
				System.out.printf("@ccs %d, clique %d%n",
					ccs.numEdges(), clique.length);
				logger.info("CCs: {}", ccs.numEdges());
				
				AssistedColoring assist
					= new AssistedColoring(name, G, clique, gdb,
						v -> rG.getRead(v));
				
				colors = assist.run();
			}
			else
			{
				SubsetReduction ssr = new SubsetReduction(G);
				SubGraph sG = ssr.apply();
				sG.saveDIMACS(new File(name+".ssr.graph"));

				MaxClique4 mc4 = new MaxClique4(sG);
				int[] clique = mc4.solve();
				
				HyperGraphV subCCs = sG.bfsCC();
				System.out.printf("@ccs %d, clique %d%n",
					subCCs.numEdges(), clique.length);
				logger.info("CCs: {}", subCCs.numEdges());
				
				SATVertexColoring sat
					= new SATVertexColoring(name, sG, clique);
				if (!sat.run(clique.length, method, minisat_path))
					return null;
				
				colors = ssr.resolveColors(sat.getColors());
			}
			
			if (colors != null)
			{
				ArrayList<Subtype> subtypes = rG.createSubtypes(colors);
				
				int numColors = subtypes.size();

				System.out.printf("@num_colors %d, real %.3f%n",
					numColors, t.now());

				logger.info("{}: colors {}, time {}",
					name, numColors, t.now());

				db.showSubtypes("vc", subtypes, true);
				Subtype.saveSubtypes(new File(name+".vc"), subtypes);
				return Subtype.reportSubtypes(subtypes, db, posArray);
			}
			return null; 
		}
		
	private static ExecutorService getWorkers(int n) {
		int m = (n != 0)? n:
			Math.max(2, Runtime.getRuntime().availableProcessors());
		
		logger.info("create {} workers", m);
		
		return Executors.newFixedThreadPool(m);
	}
		
	private static final Logger logger
		= LoggerFactory.getLogger(Main2.class);

	/*
	 * The main method for running multi-round rjMCMC. 
	 * @param Various parameters for estimating parameters from a Hastings-within-Gibbs procedure. 
     * @return Refined guess of inter-pool haplotypes in best_Haps. 
	 */
	public void population_freq_rjmcmc_multirun(int num_round, int max_iterations, int burn_in, 
			double alpha, double beta_a, double beta_c, double gamma, double c_old, double c_new, 
			double p_add, int coalescing_mismatch, double rare_cutoff){

		this.new_old_haps_beta=new BetaDistribution(c_old, c_new);	
		for(int round=0; round<num_round;round++){
			double freqs_sum=0.001;	// TODO should be an args. 
			this.curr_Haps = this.initial_Haps.clone(0);
			this.curr_Haps.checkrank_and_fullfill(freqs_sum);
		    System.out.println("Starting MCMC round " + round + ": logL = " + this.curr_Haps.logL + "\n");
			population_freq_rjmcmc(max_iterations, burn_in, alpha, beta_a, beta_c, gamma, p_add, coalescing_mismatch);
		}
		int best_index=0;
		if (num_round == 0) best_index = 0;
		else for(int round=1; round<num_round;round++)
			if(this.best_Haps_buffer.get(round).logL > this.best_Haps_buffer.get(best_index).logL)
				best_index=round;
		this.best_Haps=this.best_Haps_buffer.get(best_index).clone(rare_cutoff);
		this.logL_best = this.best_Haps_buffer.get(best_index).logL;
	}
	
	/*
	 * The main method for running single-round rjMCMC. Inspired by the hippo program (Pirinen 2009).
	 * @param Various parameters for estimating parameters from a Hastings-within-Gibbs procedure. 
     * @return Prposal for a refined guess of inter-pool haplotypes. 
	 */
	public void population_freq_rjmcmc(int max_iterations, int burn_in, double alpha, double beta_a, double beta_c, 
			double gamma, double p_add, int coalescing_mismatch){			
    	int accept_freq=0, accept_substit=0, accept_additions=0,accept_deletions=0; // ,accept_recomb=0,reject_recomb=0; also had rejection before
		try{
			System.out.println("Burning in for " + burn_in + " iterations...\n");
	    	for (int iter = 0; iter<burn_in;iter++) {	// These are the burn-in iterations. Nothing will be recorded here...
				this.curr_Haps.update_freqs(beta_a, beta_c, alpha, iter); // !!!
				this.curr_Haps.mutate_a_hap();
				this.curr_Haps.add_mutant_or_coalesce(alpha, p_add, gamma, coalescing_mismatch);
	    	}
			System.out.println("Running for " + max_iterations + " iterations...\n");
	    	for(int iter=0;iter<max_iterations;iter++){
	    		accept_freq += this.curr_Haps.update_freqs(beta_a, beta_c, alpha, iter); // !!!  		
	    		accept_substit += this.curr_Haps.mutate_a_hap();
	    		int add_or_del=this.curr_Haps.add_mutant_or_coalesce(alpha, p_add, gamma, coalescing_mismatch);
	    		if(add_or_del==1) accept_additions++;
	    		else if(add_or_del==2) accept_deletions++;
	    	}
	    }catch(Exception e){e.printStackTrace();}	 
	    this.best_Haps_buffer.add(this.curr_Haps);
		System.out.println("Number of haplotypes = " + this.curr_Haps.num_global_hap + "\tLog-likelihood = "+ this.curr_Haps.logL + "\nAcceptance Ratios:");
		System.out.println("Frequency Change: " + (double) accept_freq / max_iterations); 
		System.out.println("Single-locus Mutation: " + (double) accept_substit / max_iterations); 
		System.out.println("Haplotype Addition: " + (double) accept_additions / max_iterations); 
		System.out.println("Haplotype Deletion: " + (double) accept_deletions / max_iterations + "\n"); 
	}

	public void population_freq_aem(double epsilon, double rare_cutoff, int max_iteration){	
		double[] p= this.best_Haps.global_haps_freq.clone(); // Changed from num_all_hap because input from rjMCMC.
		double[] Rh=p.clone();
		this.best_Haps.in_pool_haps_freq = new double[this.best_Haps.num_global_hap][this.best_Haps.num_loci]; 
    	/* double rh1_tmp = 0; 
    	double rh2_0 = 0; 
    	double rh2_mean = 0;
    	double rh_tmp = 0; 
    	double[] rh2array = new double[this.num_pool];
    	double[] diffarray = new double[this.num_pool]; */ 
    	AEM: for(int ii=0;ii<max_iteration;ii++){
    		if(ii%100==0) System.out.println("\n" + ii+":"+this.best_Haps.logL);	// !!! Formerly, logL_normal.
		    SingularValueDecomposition svd=new SingularValueDecomposition(MatrixUtils.createRealMatrix(this.best_Haps.sigma));
		    RealMatrix si= svd.getSolver().getInverse();
	        System.out.println("Inverse sigma matrix:");
		    double[][] si_tmp = si.getData(); 
	        for (int r = 0; r < si.getRowDimension(); r++) {
	        	for (int c = 0 ; c < si.getColumnDimension(); c++) {
	        		if (Double.isNaN(si_tmp[r][c]))  {
	        			System.out.println(ii + ": Inverse sigma row " + r + " column " + c + " is NaN.");
	        			break AEM;
	        		}
	        		System.out.printf("%.5f\t",si_tmp[r][c]); 
	        		if (Double.isNaN(si_tmp[r][c])) System.out.print("si[" + r + "][" + c + "] is 0. ");
	        	}
	        	System.out.println();
	        }
	        System.out.println("TMP matrix:");
	        double[][] tmp=new double[this.num_pool][this.num_snp];               //tmp=sweep(A/(2*N),2,omega, "-");
	        for(int k=0;k<this.num_pool;k++){
	        	for(int q=0;q<this.num_snp;q++){
	        		tmp[k][q]=this.data[k][q]/this.est_ind_pool-this.best_Haps.mu[q];
	        		System.out.printf("%.5f\t",tmp[k][q]); 
		        	// if(ii%100==0 && j==0) System.out.println("tmp[" + k + "][" + q + "] = " + tmp[k][q] + "\t");
		        	// if (tmp[k][q] == 0)	System.out.print("tmp[" + k + "][" + q + "] is 0. ");
	        	}
	        	// if(ii%100==0 && j==0) System.out.println();
	        	System.out.println();
	        }
	        System.out.println("rh2 arrays:");
		    double[] IF=new double[this.best_Haps.num_global_hap]; 
		    for (int j=0;j<this.best_Haps.num_global_hap;j++){
		        float[] hh=this.best_Haps.global_haps[j];	 
		        // rh1=exp(-1/(4*N)*t(omega-h)%*%si%*%(omega-h))
		        double rh1=Math.exp(-1/(4.0*this.est_ind_pool)*Algebra.quadratic_form(Algebra.minus(this.best_Haps.mu, hh), si.getData()));	
		        // if(ii%100==0 && j==0) System.out.println("j:" + j + "\trh1=" + rh1);
		        RealMatrix TMP=MatrixUtils.createRealMatrix(tmp);// rh2=exp(-tmp%*%si%*%(omega-h)-diag( tmp%*%si%*%t(tmp)/2))
		        double[] rh2=Algebra.exp(Algebra.minus(si.preMultiply(TMP).operate(Algebra.minus(hh, this.best_Haps.mu)),
		        		Algebra.diag(TMP.multiply(si).multiply(TMP.transpose()).scalarMultiply(0.5).getData())));
		    	for (double rh2tmp : rh2) System.out.printf("%.5f\t", rh2tmp);
		    	System.out.println();
		        double rh=rh1*Algebra.mean(rh2);
		        /* if(j==0) {
		        	rh1_tmp = rh1; 
		        	rh2_mean = Algebra.mean(rh2); 
		        	for (int i = 0; i < rh2.length; i++) {
		        		rh2array[i] = rh2[i];
		        		diffarray[i] = (double) hh[i] - mu[i];
		        	}
		        	rh_tmp = rh; 
		        } */
		        IF[j]=rh; //   IF=c(IF, rh)
		        for (int k = 0; k < this.num_pool; k++) {
		        	// System.out.println(j);
		        	this.best_Haps.in_pool_haps_freq[j][k] = rh1 * rh2[k];
		        }
		    }
		    double delta = Algebra.sum(Algebra.abs(Algebra.minus(Rh, IF)));//sum(abs(Rh-IF))
		    double delta1 =Algebra.sum(Algebra.abs(Algebra.add(Algebra.divide(Rh, IF),-1)));//sum(abs(IF/Rh-1)) 
		    Rh=IF;
		    double[] p_new= Algebra.times(Rh,p);
		    if(ii%1==0) {
		    	System.out.printf(ii + "\tdelta = %.5f\tdelta1 = %.5f\n", delta, delta1);
		    	// System.out.printf("Hap_0\trh1 = %.5f\tmean_rh2 = %.5f\trh = %.5f\n", rh1_tmp, rh2_mean, rh_tmp);
		    	// for (double rh2tmp : rh2array) System.out.printf("%.5f\t", rh2tmp);
		    	// System.out.println();
		    	System.out.print("IF:\t");
		    	for (double impfac : IF) System.out.printf("%.3f\t", impfac);
		    	System.out.println();
		    	System.out.print("p_pre:\t");
		    	for (double freq : p_new) System.out.printf("%.3f\t", freq);
		    	System.out.println();
		    }
		    Algebra.normalize_ditribution(p_new);
		    Algebra.rmlow_and_normalize(p_new, 0.001);	// TODO should be an args

		    if(ii%1==0) System.out.print("p_new:\t");
		    if(ii%1==0) for (double freq : p_new) System.out.printf("%.3f\t", freq);
		    if(ii%1==0) System.out.println();
		    
		    if (delta < epsilon || delta1 < epsilon) {
		    	System.out.println(delta + " < " + epsilon + " || " +  delta1 + " < " + epsilon); 
		    	break;		
		    }
		    p=p_new.clone(); 
		    this.best_Haps.global_haps_freq = p.clone(); 
			this.best_Haps.update_sigma_mu_logL();;
		} 
		System.out.print("p:\t");
	    for (double freq : p) System.out.printf("%.3f\t", freq);
	    System.out.println();

	    boolean[] list_rem_haps = new boolean[this.best_Haps.num_global_hap];
	    int num_rem_hap = 0; 
		for (int h = 0; h < this.best_Haps.num_global_hap; h++)
			if (p[h] == 0) {
				list_rem_haps[h] = true; 
				num_rem_hap++; 
			} else this.best_Haps.global_haps_freq[h] = p[h]; 
		this.best_Haps.remHaps(list_rem_haps, num_rem_hap);
		
		for (int k = 0; k < this.num_pool; k++) {
			double poolCount = 0; 
			for (int j = 0; j < this.best_Haps.num_global_hap; j++) {
				this.best_Haps.in_pool_haps_freq[j][k] = this.best_Haps.in_pool_haps_freq[j][k] * p[j];
				poolCount += this.best_Haps.in_pool_haps_freq[j][k]; 
			}
			for (int j = 0; j < this.best_Haps.num_global_hap; j++) this.best_Haps.in_pool_haps_freq[j][k] = this.best_Haps.in_pool_haps_freq[j][k] / poolCount;
		}
		System.out.println("\nThe pools are composed of " + this.best_Haps.num_global_hap + " haplotypes.");
	}
}