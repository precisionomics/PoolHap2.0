package PoolHap; 

import PoolHap.Timer;
import PoolHap.GraphUV.SubGraph;

import java.io.BufferedReader;
import java.util.ArrayList;
import java.util.HashSet;
import java.util.Iterator;
import java.io.File;
import java.io.FileReader;
import java.io.PrintWriter;
import java.util.Collection;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;

import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

import org.apache.commons.math3.distribution.BetaDistribution;
import org.apache.commons.math3.linear.MatrixUtils;
import org.apache.commons.math3.linear.RealMatrix;
import org.apache.commons.math3.linear.SingularValueDecomposition;

public class PoolAnalyzer {

	/*
	How to get the graph-colouring algorithm to work:
	1) sudo apt install minisat; NOTE- Can only use Minisat on the command-line, so to use on Windows have to...
	2) Export from Eclipse as a runnable JAR file so that the Manifest contains the link to the proper main class. 
	- The results are in the file input_name.vc.
	*/
	
	static double _MAX_ = 1E20;
	double[] pop_freq;
	double[][] data;
	int[][] known_haps; 
	int num_snp;
	int num_pool;
	int num_all_hap;
	double[][] inpool_freq;	
	int num_hap_inpool;
	double freq_proportion_unknown; 

	CurrentHap curr_Haps;
	CurrentHap best_Haps;
	ArrayList<CurrentHap> best_Haps_buffer=new ArrayList<CurrentHap>();
	
	public double[] mu;
	public double[][] sigma;
	double logl_best;
	BetaDistribution new_old_haps_beta;
	
	public PoolAnalyzer(double[][] data, int[] posArray, int M, String vef_list, double freq_proportion_unknown, String minisat_path) throws Exception{
		this.freq_proportion_unknown=freq_proportion_unknown;
		this.num_hap_inpool=M;
		this.data=data;
		this.num_snp=data[0].length;
		this.num_pool=data.length;
		this.mu=new double[this.num_snp];
		this.sigma=new double[this.num_snp][this.num_snp];
		this.known_haps = graph_colouring(vef_list, posArray, data[0].length, minisat_path);
		// From trial p0 in PoolHap_Testing.
		// new int[][]{{0,1,1,0,0,0,0,1,1,1,1,1,0,1,0,1,1,0,1,1,1,1,0,0,0,0,1,0,1,1,1,0,1,1,1,1,0,0,0,0,1,0,0,1,1,0,1,1,1,0,0,0,1,1,0,0,0,1,0,1,1,1,1,1,0,1,0,1,0,1,1,1,1,1,1,1,0,1,1,0,1,1,0,0,1,1,0,1,1,1,1,1},
		// {1,0,0,0,0,0,1,0,0,0,0,0,1,1,0,1,0,1,0,1,0,0,0,0,1,1,0,1,0,0,1,1,0,0,1,1,1,1,1,0,1,1,0,0,0,0,0,1,0,0,1,1,1,0,1,0,0,0,1,0,0,0,0,0,0,1,0,0,1,0,1,0,0,0,0,0,1,0,0,0,0,0,0,1,1,0,0,0,0,0,0,0},
		// {0,0,0,1,1,1,0,0,0,0,1,0,0,0,1,0,0,1,1,0,0,0,1,1,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,1,1,1,1,0,0,1,0,1,0,1,0,0,0,0,0,1,1,0,0,0,0,0,0,1,1,1,1,1,1,1,0,0,0,0,0,0,0,0,0,1,0,0,1,0,1,0,1,0,0,0,0,0}};				
		/*
		for (int i = 0; i < this.known_haps.length; i++) {
			for (int j = 0; j < this.known_haps[0].length; j++) System.out.print(this.known_haps[i][j] + "\t");
			System.out.println();
		}
		*/
		this.setup_currhaps(this.known_haps, freq_proportion_unknown);
	}

	public static int[][] graph_colouring(String file_veflist, int[] posArray, int snps, String minisat_path) throws Exception	{
		HashSet<String> knownHapSet = new HashSet<String>();
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
				knownHapSet.addAll(run(new File(currArray[0]), posArray, sample, method, workers, minisat_path));
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
		int[][] initialHaps = new int[knownHapSet.size()][snps];
		Iterator<String> knownHapIt = knownHapSet.iterator();
		int currHap = 0;
	    while (knownHapIt.hasNext()) {
	    	String[] tmpHap = knownHapIt.next().split(""); 
	    	for (int i = 0; i < tmpHap.length; i++) {
	    		initialHaps[currHap][i] = Integer.parseInt(tmpHap[i]);
	    	}
	    	currHap++; 
	    }
	    return initialHaps;
	}

	private static HashSet<String> run(File input, int[] posArray, 
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
			return Subtype.reportSubtypes(subtypes, posArray);
		}
		return null; 
	}
	
	private static ExecutorService getWorkers(int n)
	{
		int m = (n != 0)? n:
			Math.max(2, Runtime.getRuntime().availableProcessors());
		
		logger.info("create {} workers", m);
		
		return Executors.newFixedThreadPool(m);
	}
	
	private static final Logger logger
		= LoggerFactory.getLogger(Main.class);
	
	public void setup_currhaps(int[][] known_haps, double freq_proportion_unknown){
		this.curr_Haps=new CurrentHap(known_haps, freq_proportion_unknown, this.num_snp);
		this.curr_Haps.analyzer=this;
	}
	
	public void update_sigma_mu_fixed(double[] p){
		this.mu=new double[num_snp];
	    for(int q=0;q<num_snp;q++){ // assign omega based on p
	    	for(int h=0;h<best_Haps.num_curr_H;h++){	// Changed from num_all_hap because input from rjMCMC.
	    		this.mu[q]=this.mu[q]+p[h]*best_Haps.nodes.get(h).hap[q];
	    	}
	    }//  apply(H*p,2, sum); 
	    double[][] eta=new double[num_snp][num_snp];
	    for(int q1=0;q1<num_snp;q1++){
	    	 for(int q2=q1;q2<num_snp;q2++){
	    		 for(int h=0;h<best_Haps.num_curr_H;h++){
	 		    	eta[q1][q2]+=(best_Haps.nodes.get(h).hap[q1]*best_Haps.nodes.get(h).hap[q2]*p[h]);
	 		     }eta[q2][q1]=eta[q1][q2];	
			 }
	    }  // eta=t(H)%*%diag(p)%*%H //gsl_blas_dsyr(CblasLower, p[i], H[i], sigma); //sigma+=p[i]*H[i]*H[i]^T
	    this.sigma=new double[num_snp][num_snp];
	    for(int q1=0;q1<num_snp;q1++){
	    	for(int q2=0;q2<num_snp;q2++){
	    		sigma[q1][q2]=eta[q1][q2]-mu[q1]*mu[q2];
	    	}
	    } //gsl_blas_dsyr(CblasLower,-1.0,mu,sigma);//sigma-=mu*mu^T 
	}
			
	public void population_freq_rjmcmc_multirun(double alpha, double beta_a, double beta_c, double gamma, double c_old, 
			double c_new, double p_add, int max_iterations, int burn_in,int coalescing_mismatch,
			int num_round, String working_folder, double rare_cutoff){

		this.new_old_haps_beta=new BetaDistribution(c_old, c_new);

		for(int round=0; round<num_round;round++){
			double freqs_sum=0.0001;
			this.setup_currhaps(this.known_haps, freq_proportion_unknown);
			this.curr_Haps.checkrank_and_fullfill(freqs_sum);	// NOTE: This is only called once in the Hippo program.
		    this.curr_Haps.update_sigma_mu_current();
		    this.curr_Haps.loglikelihood=Algebra.logL_aems(Algebra.times(sigma, num_hap_inpool),
					Algebra.times(mu, num_hap_inpool), data);	// !!! Formerly, logL_normal. logL_rjmcmc(this.sigma, this.mu, this.data, this.num_hap_inpool)
		    this.logl_best=this.curr_Haps.loglikelihood;
		    System.out.println("Starting MCMC round " + round + ": logL = "+this.curr_Haps.loglikelihood);
			population_freq_rjmcmc(alpha, beta_a, beta_c, gamma, c_old, c_new, p_add, 
					max_iterations, burn_in, coalescing_mismatch, round, working_folder);
		}
		int best_index=0;
		if (num_round == 0) {
			this.best_Haps = curr_Haps.clone_highfreq(rare_cutoff);
		} else {
			for(int round=1; round<num_round;round++){
				if(this.best_Haps_buffer.get(round).loglikelihood>this.best_Haps_buffer.get(best_index).loglikelihood)
					best_index=round;
			}
			this.best_Haps=this.best_Haps_buffer.get(best_index).clone_highfreq(rare_cutoff); // clone_highfreq(rare_cutoff)
		}
	}
		
	public void population_freq_rjmcmc(double alpha, double beta_a, double beta_c, double gamma, double c_old, 
			double c_new, double p_add, int max_iterations, int burn_in,int coalescing_mismatch, int round_index,
			String working_folder){	
		/* Removed to see if this helps prioritize the original haplotypes.
		if (round_index == 0) {	// !!To check accuracy of rjMCMC!!! Sets all of the haplotype frequencies to the same in the first round.
			for (int h = 0; h < this.curr_Haps.num_curr_H; h++) {
				this.curr_Haps.nodes.get(h).change_my_freq(1.0 / (double)this.curr_Haps.num_curr_H); 
			}
		}
		*/
	    try{
			System.out.println("Burning in for " + burn_in + " rounds...\n");
	    	for (int iter=0;iter<burn_in;iter++) {	// These are the burn-in iterations. Nothing will be recorded here...
				@SuppressWarnings("unused")
				int freq_change=this.curr_Haps.update_freqs(beta_a, beta_c, alpha, iter, num_hap_inpool); // !!!
	    		@SuppressWarnings("unused")
				int add_or_del=this.curr_Haps.add_mutant_or_coalesce(alpha, p_add, c_old, c_new, gamma, coalescing_mismatch, num_hap_inpool);
	    		@SuppressWarnings("unused")
				int substitute=this.curr_Haps.mutate_a_hap(num_hap_inpool);
	    	}
	    	int accept_freq=0, reject_freq=0, accept_additions=0,reject_additions=0,accept_deletions=0,
	    			reject_deletions=0,accept_substit=0,reject_substit=0; // ,accept_recomb=0,reject_recomb=0;
	    	for(int iter=0;iter<max_iterations;iter++){
	    		// if (iter%1000==0) System.out.println("Currently on iteration " + iter + "...\n");
	    		int report = max_iterations / 5; 
	    		if(iter%report==0) {
	    			System.out.println(round_index+":"+iter+" n_haps="+this.curr_Haps.num_curr_H+", logL="+this.logl_best+"\n");
	    			System.out.println(" acc/rej_frequency="+accept_freq+"/"+reject_freq+",  "
		    	  		+ "acc/rej_additions="+accept_additions+"/"+reject_additions+", \n "
		    	  		+ "acc/rej_deletions="+accept_deletions+"/"+reject_deletions+",   "
		    	  		+ "acc/rej_substitutions="+accept_substit+"/"+reject_substit+"\n ");
	    				// + "acc/rej_recombinations="+accept_recomb+"/"+reject_recomb+"\n");
	    		}
				int freq_change=this.curr_Haps.update_freqs(beta_a, beta_c, alpha, iter, num_hap_inpool); // !!!
	    		if(freq_change==1)accept_freq++;
	    		else if(freq_change==-1)reject_freq++;
	    		
	    		int add_or_del=this.curr_Haps.add_mutant_or_coalesce(alpha, p_add, c_old, c_new, gamma, coalescing_mismatch, num_hap_inpool);
	    		if(add_or_del==1)accept_additions++;
	    		else if(add_or_del==-1)reject_additions++;
	    		else if(add_or_del==2)accept_deletions++;
	    		else if(add_or_del==-2)reject_deletions++;
	    		
	    		int substitute=this.curr_Haps.mutate_a_hap(num_hap_inpool);
	    		if(substitute==1)accept_substit++;
	    		else if(substitute==-1)reject_substit++;
	    		
	    		/*
	    		int recomb_or_not=this.curr_Haps.recombine_two_haps(alpha, gamma, iter);
	    		if(recomb_or_not==1)accept_recomb++;
	    		else if(recomb_or_not==-1)reject_recomb++;
	    		*/
	    		if(iter%2000==0) System.out.println("iter " + iter + " done...");
	    	}
	    }catch(Exception e){e.printStackTrace();}	 
	    this.best_Haps_buffer.add(curr_Haps.clone());
	    System.out.println();
	    // System.out.println("curr_Haps likelihood = " + curr_Haps.loglikelihood + "\tSize of best_Haps_buffer = " + best_Haps_buffer.size());
	}

	public void population_freq_aem(double[][] data, double epsilon, double rare_cutoff, int max_iteration){	
		double[] p=new double[best_Haps.num_curr_H]; // Changed from num_all_hap because input from rjMCMC.
		for (int f = 0; f < p.length; f++) p[f] = best_Haps.nodes.get(f).freq; 
		// System.out.print("freqs_init: ");
		// for (double freq : p) System.out.print(freq + "\t");
		// System.out.println();
		double[] p_record = new double[best_Haps.num_curr_H]; 
		double[] Rh=p.clone();
		this.inpool_freq = new double[best_Haps.num_curr_H][num_pool];
    	double rh1_tmp = 0; 
    	double rh2_0 = 0; 
    	double rh2_mean = 0;
    	double rh_tmp = 0; 
    	double[] rh2array = new double[this.num_pool];
    	double[] diffarray = new double[this.num_pool]; 
    	AEM: for(int ii=0;ii<max_iteration;ii++){
			this.update_sigma_mu_fixed(p);
    		if(ii%100==0) System.out.println("\n" + ii+":"+Algebra.logL_aems(Algebra.times(sigma, num_hap_inpool),
					Algebra.times(mu, num_hap_inpool), data));	// !!! Formerly, logL_normal.
		    SingularValueDecomposition svd=new SingularValueDecomposition(MatrixUtils.createRealMatrix(sigma));
		    RealMatrix si= svd.getSolver().getInverse();
	        System.out.println("Inverse sigma matrix:");
		    double[][] si_tmp = si.getData(); 
	        for (int r = 0; r < si.getRowDimension(); r++) {
	        	for (int c = 0 ; c < si.getColumnDimension(); c++) {
	        		if (Double.isNaN(si_tmp[r][c]))  {
	        			System.out.println(ii + ": Inverse sigma row " + r + " column " + c + " is NaN.");
	        		    p=p_record.clone();      
	        			break AEM;
	        		}
	        		System.out.printf("%.5f\t",si_tmp[r][c]); 
	        		if (Double.isNaN(si_tmp[r][c])) System.out.print("si[" + r + "][" + c + "] is 0. ");
	        	}
	        	System.out.println();
	        }
	        System.out.println("TMP matrix:");
	        double[][] tmp=new double[num_pool][num_snp];               //tmp=sweep(A/(2*N),2,omega, "-");
	        for(int k=0;k<num_pool;k++){
	        	for(int q=0;q<num_snp;q++){
	        		tmp[k][q]=data[k][q]/this.num_hap_inpool-mu[q];
	        		System.out.printf("%.5f\t",tmp[k][q]); 
		        	// if(ii%100==0 && j==0) System.out.println("tmp[" + k + "][" + q + "] = " + tmp[k][q] + "\t");
		        	// if (tmp[k][q] == 0)	System.out.print("tmp[" + k + "][" + q + "] is 0. ");
	        	}
	        	// if(ii%100==0 && j==0) System.out.println();
	        	System.out.println();
	        }
	        System.out.println("rh2 arrays:");
		    double[] IF=new double[best_Haps.num_curr_H]; 
		    for (int j=0;j<best_Haps.num_curr_H;j++){
		        int[] hh=best_Haps.nodes.get(j).hap;	 
		        // rh1=exp(-1/(4*N)*t(omega-h)%*%si%*%(omega-h))
		        double rh1=Math.exp(-1/(4.0*num_hap_inpool)*Algebra.quadratic_form(Algebra.minus(mu, hh), si.getData()));	
		        // if(ii%100==0 && j==0) System.out.println("j:" + j + "\trh1=" + rh1);
		        RealMatrix TMP=MatrixUtils.createRealMatrix(tmp);// rh2=exp(-tmp%*%si%*%(omega-h)-diag( tmp%*%si%*%t(tmp)/2))
		        double[] rh2=Algebra.exp(Algebra.minus(si.preMultiply(TMP).operate(Algebra.minus(hh, mu)),
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
		        for (int k = 0; k < num_pool; k++) this.inpool_freq[j][k] = rh1 * rh2[k];
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
		    Algebra.rmlow_and_normalize(p_new, rare_cutoff/100);

		    if(ii%1==0) System.out.print("p_new:\t");
		    if(ii%1==0) for (double freq : p_new) System.out.printf("%.3f\t", freq);
		    if(ii%1==0) System.out.println();
		    
		    if (delta < epsilon || delta1 < epsilon) {
		    	System.out.println(delta + " < " + epsilon + " || " +  delta1 + " < " + epsilon); 
		    	break;		
		    }
		    p_record = p.clone();    
		    p=p_new.clone();      
		} this.pop_freq= p.clone();
		System.out.print("p:\t");
	    for (double freq : p) System.out.printf("%.3f\t", freq);
	    System.out.println();
		for (int k = 0; k < num_pool; k++) {
			double poolCount = 0; 
			for (int j = 0; j < best_Haps.num_curr_H; j++) {
				this.inpool_freq[j][k] = this.inpool_freq[j][k] * p[j];
				poolCount += this.inpool_freq[j][k]; 
			}
			for (int j = 0; j < best_Haps.num_curr_H; j++) this.inpool_freq[j][k] = this.inpool_freq[j][k] / poolCount;
		}
		for (int f = 0; f < p.length; f++) best_Haps.nodes.get(f).change_my_freq(p[f]);
		int maxIndex = best_Haps.nodes.size() - 1; 
		for (int i = maxIndex; i >= 0; i--) {
			if (best_Haps.nodes.get(i).freq == 0) {
				// System.out.println(i + "\t" + best_Haps.num_curr_H + "\t" + best_Haps.nodes.get(i).freq);
				best_Haps.delete(i);
			}			
		}
		System.out.println("\nThere are " + best_Haps.num_curr_H + " haplotypes left.");
		best_Haps.loglikelihood = Algebra.logL_aems(Algebra.times(sigma, num_hap_inpool),
				Algebra.times(mu, num_hap_inpool), data); 
	}

	public void print_intrapool(PrintWriter pw) {
		for(int h=0;h<best_Haps.num_curr_H;h++){
			if (best_Haps.nodes.get(h).freq == 0) continue; 
			pw.append("hap_"+h+"\t"+best_Haps.nodes.get(h).output_nofreq("") + "\t");
			for (int k = 0; k < num_pool; k++) pw.format("%.3f\t",this.inpool_freq[h][k]); 
			pw.append("\n");
		}
	}

	public void print_buffer(double rare_cutoff, PrintWriter pw){
		for(int k=0;k<this.best_Haps_buffer.size();k++){
			pw.append("No."+k+" round" + "\n");
			this.best_Haps_buffer.get(k).print(pw);
		}	
	}
}