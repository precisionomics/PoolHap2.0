package PoolHap; 

import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.concurrent.ThreadLocalRandom;

import org.apache.commons.math3.distribution.BetaDistribution;
import org.apache.commons.math3.exception.NumberIsTooSmallException; 
import org.apache.commons.math3.linear.MatrixUtils;
import org.apache.commons.math3.linear.SingularValueDecomposition;

/*
 * An ArrayList maintaining current haplotypes during the MCMC sampling
 * It generates the mutation likelihoods for acceptance probability in rjMCMC.
 */
public class CurrentHap {
	ArrayList<OneHap> nodes;
	int num_curr_H;
	int num_snps;
	double loglikelihood;
	
	PoolAnalyzer analyzer;
	
	public CurrentHap clone(){
		return (new CurrentHap(nodes, num_curr_H, num_snps, loglikelihood));
	}
	
	public CurrentHap clone_highfreq(double freq_cutoff){
		ArrayList<Integer> common_hap_index = new ArrayList<Integer>();
		ArrayList<OneHap> common_haps = new ArrayList<OneHap>();
		for (int h = 0; h < num_curr_H; h++) {
			if (this.nodes.get(h).freq >= freq_cutoff) common_hap_index.add(h); 
		}
		/*
		HashMap<Double,Integer> freq2index = new HashMap<Double,Integer>();
		for (int h = 0; h < num_curr_H; h++) freq2index.put(this.nodes.get(h).freq, h);
		ArrayList<Double> sortedFreq = new ArrayList<Double>(freq2index.keySet());
		Collections.sort(sortedFreq, Collections.reverseOrder());
		for (int f = 0; f < 20; f++) common_haps.add(this.nodes.get(freq2index.get((sortedFreq.get(f)))));
		*/
		for (Integer h : common_hap_index) common_haps.add(this.nodes.get(h)); 
		return (new CurrentHap(common_haps, common_haps.size(), num_snps, loglikelihood));
	}
	
	public CurrentHap(ArrayList<OneHap> nodes, int num_curr_H, int num_snps, double loglikelihood){
		this.loglikelihood=loglikelihood;
		this.nodes=new ArrayList<OneHap>();
		this.num_curr_H=num_curr_H;
		this.num_snps=num_snps;
		for(int k=0;k<num_curr_H;k++){
			OneHap the_node=new OneHap(nodes.get(k).hap.clone(), nodes.get(k).freq);
			this.nodes.add(the_node);
		}
	}
	
	public CurrentHap(int[][] known_haps, double freq_proportion_unknown, int num_snp){
		this.num_snps=num_snp;
		ArrayList<Double> freqs=new ArrayList<Double>();
		this.nodes=new ArrayList<OneHap>();
		double the_freq=Double.NaN;
		for (int h = 0; h < known_haps.length; h++) {
			int[] the_hap = new int[num_snp];
			if(this.num_snps!=known_haps[h].length)	System.out.println("Length of known haps not equals to #SNP.");	
			for (int k = 0; k < this.num_snps; k++) the_hap[k] = known_haps[h][k];
			this.nodes.add(new OneHap(the_hap, the_freq));
			freqs.add(the_freq);
			this.num_curr_H++;
		}
		/*
		for (int h = 0; h < known_haps.length; h++) {
			for (int k = 0; k < this.num_snps; k++) System.out.print(this.nodes.get(h).hap[k] + "\t");
			System.out.println();
		}
		*/
		Algebra.normalize_ditribution(freqs, freq_proportion_unknown);
		for(int k=0;k<this.num_curr_H;k++) {
			this.nodes.get(k).freq=freqs.get(k);
		}
	}			
	
	public CurrentHap(int[][] current_H, double[] current_freq, int num_snp){
		this.num_curr_H=current_H.length;
		this.num_snps=num_snp;
		this.nodes=new ArrayList<OneHap>();
		for(int k=0;k<this.num_curr_H;k++){
			this.nodes.add(new OneHap(current_H[k].clone(), current_freq[k]));
		}
	}
	
	public void checkrank_and_fullfill(double freqs_sum){
		if(this.num_curr_H<this.num_snps)this.fullfill(freqs_sum);
		else{
			double[][] the_H_array=new double[num_curr_H][num_snps];
			for(int k=0;k<num_curr_H;k++){
				for(int i=0;i<num_snps;i++)the_H_array[k][i]=this.nodes.get(k).hap[i];
			}
			SingularValueDecomposition svd=new SingularValueDecomposition(MatrixUtils.createRealMatrix(the_H_array));
			double[] sv=svd.getSingularValues();
			int rank=0;
			for(int i=0;i<sv.length;i++){
				if(sv[i]!=0)rank++;
			}if(rank<num_snps)this.fullfill(freqs_sum);
		}
	}
	
	public void fullfill(double freqs_sum){
		 
		for(int k=0;k<this.num_curr_H;k++)this.nodes.get(k).freq*=(1-freqs_sum);
		int fillRank = this.num_snps - this.num_curr_H;
		for (int p = 0; p < fillRank; p++) {
			int[] the_new_hap=new int[num_snps];
			the_new_hap[ThreadLocalRandom.current().nextInt(0,num_snps)]=1;
			OneHap the_new_node=new OneHap(the_new_hap, freqs_sum/num_snps);
			this.append(the_new_node);
		}
		/* The original fulfill function, which generates as many haplotypes as there are SNPs. The new one generates as many haplotypes as are missing. 
		 * 		for(int i=0;i<num_snps;i++){
			int[] the_new_hap=new int[num_snps];
			the_new_hap[i]=1;
			OneHap the_new_node=new OneHap(the_new_hap, freqs_sum/num_snps);
			this.append(the_new_node);
		}
		 */
	}
	
	void delete(int index){
		if(this.num_curr_H==1){
			System.out.println("Can't delete the only hap in the current distribution. Continued...");
			return;
		}
		if(index>this.num_curr_H-1){
			System.out.println("Trying to delete the i_th hap, while i>=this.num_curr_H: Wrong");
			return;
		}
		this.nodes.set(index, this.nodes.get(this.num_curr_H-1));
		this.nodes.remove(num_curr_H-1);
		
		this.num_curr_H--;
	}
	
	void append(OneHap the_node){
		this.nodes.add(the_node);
		this.num_curr_H++;
	}
	
	/*
	 * Update \sigma and \mu based on the current hap frequencies, ready for Log-likelihood calculation.
	 * Its difference to the above one is the multiplier of number of haps in each pool
	 */
	public void update_sigma_mu_current(){
		analyzer.mu=new double[num_snps];
	    for(int q=0;q<this.num_snps;q++){ // assign omega based on p
	    	for(int h=0;h<num_curr_H;h++){		    
	    		analyzer.mu[q]=analyzer.mu[q]+this.nodes.get(h).freq*this.nodes.get(h).hap[q];
	    		// System.out.println(this.nodes.get(h).freq + "\t" + this.nodes.get(h).hap[q]);
	    	}
	    } 
	    double[][] eta=new double[num_snps][num_snps];
	    for(int q1=0;q1<num_snps;q1++){
	    	 for(int q2=0;q2<num_snps;q2++){
	    		 for(int h=0;h<num_curr_H;h++){
	 		    	eta[q1][q2]+=(this.nodes.get(h).hap[q1]*this.nodes.get(h).hap[q2]*this.nodes.get(h).freq);
	 		     }
			 }
	    }
	    analyzer.sigma=new double[num_snps][num_snps];
	    for(int q1=0;q1<num_snps;q1++){
	    	 for(int q2=0;q2<num_snps;q2++){
	    		 analyzer.sigma[q1][q2]=eta[q1][q2]-analyzer.mu[q1]*analyzer.mu[q2];
	    	 }
	    } //gsl_blas_dsyr(CblasLower,-1.0,mu,sigma);//sigma-=mu*mu^T 
	}
	
	/*
	 * return two indexes of haps, ensuring i1 < i2 
	 */
	int[] sample_two_haps(){
		int i1=(int)(ThreadLocalRandom.current().nextDouble()*this.num_curr_H);
	    int i2=(int)(ThreadLocalRandom.current().nextDouble()*this.num_curr_H);
	    if(i1==this.num_curr_H)i1--;
	    if(i2==this.num_curr_H)i2--;
	    if(i1>i2){int tmp=i1;i1=i2;i2=tmp;}
	    if(i1==i2){
	    	if(i2!=this.num_curr_H-1)i2++;
	    	else i1--;
	    }
	    int[] indexes=new int[2];
	    indexes[0]=i1; indexes[1]=i2;
	    return indexes;
	}

	int[] sample_two_haps_prob(){	// !!!
		double cumul_freq = 0;
		double cumul_lim1 = ThreadLocalRandom.current().nextDouble();	// Real number between 0 and 1; 'floor' of maximum cumulative frequency.
		int i1 = 0, i2 = 0; 
		while ((cumul_lim1 > cumul_freq) && (i1 < (this.nodes.size()-1))) {
			cumul_freq += this.nodes.get(i1).freq;		// pop_freq is only updated in the EM algorithm. Need a different one for MCMC testing.
			i1++; 
		}
		double cumul_lim2 = ThreadLocalRandom.current().nextDouble() * (1 - this.nodes.get(i1).freq);
		// System.out.println(cumul_lim2);
		if (i1 == 0) i2 = 1;
		cumul_freq = 0; 
		while ((cumul_lim2 > cumul_freq) && (i2 < (this.nodes.size()-1))) {
			if (i2 == i1) i2++; 
			cumul_freq += this.nodes.get(i2).freq; 
			i2++; 
		}
		int[] indexes = new int[2];
		if (i2 < i1) {indexes[0] = i2; indexes[1] = i1;}
		if (i2 == this.nodes.size()) i2--;
		else {indexes[0] = i1; indexes[1] = i2;}
		// System.out.println(this.nodes.size() + "\t" + indexes[0] + "\t" + indexes[1]);
		return indexes; 
	}
	
	int sample_one_hap(){
		int i1=(int)(ThreadLocalRandom.current().nextDouble()*this.num_curr_H);
	    if(i1==this.num_curr_H)i1--;
	    return i1;
	}
	
	double coalescing_probability(int i1, int i2){
		double sum=0;
		for(int i=0;i<this.num_curr_H;i++){
			int diff=0;
			for(int l=0;l<this.num_snps;l++){
				if(nodes.get(i1).hap[l]!=nodes.get(i).hap[l])diff++;
				if(diff>1) break;
			}
			if(diff==1)sum+=nodes.get(i).freq;
			if(i==i2 && diff!=1)System.out.println("ERROR: in coalescing probability, d(h(i1),h(12))>1\n");
	    }
		return nodes.get(i2).freq/sum;
	}
		
	void print(PrintWriter pw){
		pw.append("haps:\t"+this.num_curr_H+"\tsnps:\t"+this.num_snps+"\n");
		pw.append("logL\t" + this.loglikelihood + "\n");
		for(int h=0;h<this.num_curr_H;h++){
			pw.append("hap_"+h+"\t"+this.nodes.get(h).output("") + "\n");
		}
	}
	
	void print_stdout(){
		for(int h=0;h<this.num_curr_H;h++)System.out.println("hap_"+h+"\t"+this.nodes.get(h).output(""));
	}

	void print_stdout(double freq_cutoff){
		for(int h=0;h<this.num_curr_H;h++)
			if(this.nodes.get(h).freq>freq_cutoff) System.out.println("hap_"+h+"\t"+this.nodes.get(h).output(""));
	}
	
	/*
	 * search a hap with at most "coalescing_mismatch".  
	 */
	int seek_coalescence(int index, int coalescing_mismatch, double[] return_proba){
		int unacceptable=coalescing_mismatch+1;
		// double sum_freq_closest=0;
		int[] distance=new int[this.num_curr_H];
		int smallest_found=unacceptable; // initialize
		distance[index]=unacceptable;  // itself will be excluded.
		int[] the_target=this.nodes.get(index).hap;
		for(int h=0;h<this.num_curr_H;h++){
			if(h==index)continue;
			int[] to_compare=this.nodes.get(h).hap;
			for(int k=0;k<this.num_snps;k++){				
				if(the_target[k]!=to_compare[k]){
					distance[h]++;
					if(distance[h]==unacceptable)break;
				}
			}if(distance[h]<smallest_found)smallest_found=distance[h];
		}if(smallest_found==unacceptable)return -1; // no hap found
		ArrayList<Integer> the_closest_ones=new ArrayList<>();		
		for(int h=0;h<this.num_curr_H;h++){
			if(smallest_found==distance[h]){
				the_closest_ones.add(h);
				// sum_freq_closest+=nodes.get(h).freq;
			}
		}// randomly choose an index from the closest ones:
		int the_chosen=(int)(ThreadLocalRandom.current().nextDouble()*the_closest_ones.size());
		if(the_chosen==the_closest_ones.size())the_chosen--;
		return_proba[0]=1.0/the_closest_ones.size();
		return the_closest_ones.get(the_chosen);
	}
	
	public int search_a_hap(int[] hap){
		for(int i=0;i<this.num_curr_H;i++){
			boolean match=true;
			int[] the_candidate=this.nodes.get(i).hap;
			for(int k=0;k<this.num_snps;k++){
				if(hap[k]!=the_candidate[k]){
					match=false;break;
				}
			}if(match)return i;
		}return -1;
	}
	
	/*
	 * Update frequencies of two haplotypes 
	 */
	int update_freqs(double beta_a, double beta_c, double alpha, int iter, int N){	// !!!
		int[] index = new int[2]; 
		if ((iter % 3) == 0) {	// !!! 
			index = this.sample_two_haps();	
		} else {
			index = this.sample_two_haps_prob();
		}
		OneHap n1=nodes.get(index[0]);
		OneHap n2=nodes.get(index[1]);
		double ori_freq_1=n1.freq;
		double ori_freq_2=n2.freq; // backup for returning to the original of rejected. 
	    BetaDistribution beta_dist=new BetaDistribution(beta_a+n1.freq*beta_c, beta_a+ n2.freq*beta_c);
	    double css= beta_dist.sample(); //gsl_ran_beta(rng,beta_a+p[i1]*beta_c,beta_a+p[i2]*beta_c);
	    if(css>0.99999) css=0.99999;
	    if(css<0.00001) css=0.00001;
	    n1.freq=(n1.freq+n2.freq)*css;
	    n2.freq=(ori_freq_1+ori_freq_2)-n1.freq;
	    this.update_sigma_mu_current();
	    double logl2=Algebra.logL_aems(Algebra.times(analyzer.sigma, N), Algebra.times(analyzer.mu,N), analyzer.data);	// logL_rjmcmc(analyzer.sigma, analyzer.mu, analyzer.data, analyzer.num_hap_inpool)
	    double logHR=logl2-this.loglikelihood;	
	    BetaDistribution beta_dist2=new BetaDistribution(beta_a+n1.freq*beta_c,beta_a+n2.freq*beta_c);  
	    double prop_dens = 0;	// !!!
	    try {	// !!!
	    	double below_min = 0; 
	    	if (ori_freq_2 > Math.pow(10, -17)) below_min = ori_freq_2; 
	    	prop_dens = beta_dist2.logDensity(ori_freq_1/(ori_freq_1+below_min))-beta_dist.logDensity(n1.freq/(n1.freq+n2.freq));
	    } catch (NumberIsTooSmallException e) {
	    }
	    logHR=logHR+prop_dens;    
	    logHR+=(alpha-1)*(Math.log(n1.freq)+Math.log(n2.freq)-Math.log(ori_freq_1)-Math.log(ori_freq_2));
	    if ((iter % 3) != 0) { // !!!
	    	logHR+=Math.log(n1.freq)+Math.log(n2.freq)+Math.log(1-ori_freq_2)+Math.log(1-ori_freq_1)-(Math.log(1-n1.freq)+Math.log(1-n2.freq)+Math.log(ori_freq_2)+Math.log(ori_freq_1));
	    }
		if( (!(Double.isNaN(logl2))) && Math.log(ThreadLocalRandom.current().nextDouble())<=logHR){//proposal accepted
			// System.out.println("====Freq Accepted. Logl===: "+logl2+"/"+this.loglikelihood);
			this.loglikelihood=logl2;
			if(this.loglikelihood>analyzer.logl_best)analyzer.logl_best=this.loglikelihood;
			return 1;
		}
		else{//proposal rejected
			n1.freq=ori_freq_1;
			n2.freq=ori_freq_2;
			return -1;
		}		
	}
	
	/*
	 * Add one hap by mutating an existing hap, or remove one hap by coalescing it with a close one.
	 * return 0 if nothing done, 1 if mutated, -1 if mutant rejected, 2 if coalesced, -2 if coalescence rejected.  
	 */
	int add_mutant_or_coalesce(double alpha, double p_add, double c_old, double c_new, double gamma, int coalescing_mismatch, int N){
		int index=this.sample_one_hap();
		boolean proceed_mut;
		if(this.num_curr_H<=2)proceed_mut=true;
		else proceed_mut=ThreadLocalRandom.current().nextDouble()<p_add;
		if(proceed_mut){
			double proportion=analyzer.new_old_haps_beta.sample();
			// move some frequency from the old hap to the new one. 
			double ori_freq=nodes.get(index).freq;
			OneHap n1=nodes.get(index).generate_mutant((1.0-proportion)*nodes.get(index).freq);
			if(search_a_hap(n1.hap)!=-1){return 0;}
			nodes.get(index).freq=proportion*nodes.get(index).freq;
			this.append(n1);		
			// calculating new logL:
		    this.update_sigma_mu_current();
		    double logl2=Algebra.logL_aems(Algebra.times(analyzer.sigma, N), Algebra.times(analyzer.mu,N), analyzer.data);	// logL_rjmcmc(analyzer.sigma, analyzer.mu, analyzer.data, analyzer.num_hap_inpool)
		    double logHR=logl2-this.loglikelihood;	
		    logHR-=gamma;
		    logHR+=(alpha-1)*(Math.log(nodes.get(index).freq)+Math.log(n1.freq)-Math.log(ori_freq)); //prior for p
		    logHR+=-Math.log(this.num_curr_H+1.0)+Math.log(1-p_add)+Math.log(coalescing_probability(index, num_curr_H-1));//inverse transition probability
		    logHR+=Math.log(this.num_curr_H*this.num_snps*1.0)-Math.log(analyzer.new_old_haps_beta.density(proportion)); //transition probability
		    logHR+=Math.log(ori_freq); //|Jacobian|    
		    
			if((!(Double.isNaN(logl2))) && Math.log(ThreadLocalRandom.current().nextDouble())<=logHR){//proposal accepted				
				// System.out.println("=====Add Accepted. Logl====: "+logl2+"/"+this.loglikelihood);
				this.loglikelihood=logl2;
				if(this.loglikelihood>analyzer.logl_best)analyzer.logl_best=this.loglikelihood;
				return 1;
			}
			else{//proposal rejected
				nodes.get(index).freq=ori_freq;
				this.delete(num_curr_H-1);
				return -1;
			}	
		}else{ //coalesce index with the closest hap.
			double[] sampling_prob=new double[1];
			int coal_index=seek_coalescence(index, coalescing_mismatch, sampling_prob);
			if(coal_index==-1)return 0; // no hap selected, return 0.
			double ori_freq=nodes.get(index).freq;
			double ori_freq2=nodes.get(coal_index).freq;
			nodes.get(index).freq=0;
			nodes.get(coal_index).freq+=ori_freq;
			// calculating new logL:
		    this.update_sigma_mu_current();
		    double logl2=Algebra.logL_aems(Algebra.times(analyzer.sigma, N), Algebra.times(analyzer.mu,N), analyzer.data);	// logL_rjmcmc(analyzer.sigma, analyzer.mu, analyzer.data, analyzer.num_hap_inpool)
		    //System.out.println("ori_p:"+index+":"+ori_freq);
		    //System.out.println("tested_p:"+index+":"+nodes.get(index).freq+"/combined to:"+nodes.get(coal_index).freq);
		 // calculating acceptance prob: 
		    double logHR=logl2-this.loglikelihood;	
		    logHR+=gamma;
		    logHR+=(alpha-1)*(Math.log(nodes.get(coal_index).freq)-Math.log(ori_freq)-Math.log(ori_freq2)); //prior for p
		    logHR-=-Math.log(this.num_curr_H)+Math.log(1-p_add)+Math.log(sampling_prob[0]);//transition probability
		    logHR+=-Math.log((this.num_curr_H-1)*this.num_snps*1.0)+
		    		Math.log(analyzer.new_old_haps_beta.density(ori_freq2/nodes.get(coal_index).freq)); //inverse transition probability
			logHR+=-Math.log(ori_freq+ori_freq2); //|Jacobian|
			
			if((!(Double.isNaN(logl2))) && Math.log(ThreadLocalRandom.current().nextDouble())<=logHR){//proposal accepted				
				if(this.loglikelihood>analyzer.logl_best)analyzer.logl_best=this.loglikelihood;
				// System.out.println("=====Del Accepted. Logl====: "+logl2+"/"+this.loglikelihood);
				this.loglikelihood=logl2;
				this.delete(index);
				return 2;
			}
			else{//proposal rejected
				nodes.get(index).freq=ori_freq;
				nodes.get(coal_index).freq=ori_freq2;
				return -2;
			}	
		}					
	}

	/*
	 * try to mutate a hap (with frequency unchanged), return 0 if the mutant is in the current haps
	 */
	int mutate_a_hap(int N){
		int index=this.sample_one_hap(); 
		OneHap ori_hap=nodes.get(index);			
		OneHap n1=nodes.get(index).generate_mutant(nodes.get(index).freq);
		if(search_a_hap(n1.hap)!=-1)return 0;
		this.nodes.set(index, n1);		
		// calculating new logL:
		this.update_sigma_mu_current();
	    double logl2=Algebra.logL_aems(Algebra.times(analyzer.sigma, N), Algebra.times(analyzer.mu,N), analyzer.data);	// logL_rjmcmc(analyzer.sigma, analyzer.mu, analyzer.data, analyzer.num_hap_inpool)
		// calculating acceptance prob: 
		double logHR=logl2-this.loglikelihood;	
		if((!(Double.isNaN(logl2))) && Math.log(ThreadLocalRandom.current().nextDouble())<=logHR){//proposal accepted			
			// System.out.println("====Mut Accepted. Logl====: "+logl2+"/"+this.loglikelihood);
			this.loglikelihood=logl2;
			if(this.loglikelihood>analyzer.logl_best)analyzer.logl_best=this.loglikelihood;
			return 1;
		}
		else{//proposal rejected
			nodes.set(index,ori_hap);
			return -1;
		}	
	}
	
	// NOTE: Hippo doesn't actually have this function so the calculation of the probability density ratio may not be correct! 
	int recombine_two_haps(double alpha, double gamma, int iter, int N){ 
		int[] index = new int[2]; 
		if ((iter % 3) == 0) {	// !!! 
			index = this.sample_two_haps();	
		} else {
			index = this.sample_two_haps_prob();
		}
		OneHap n1=nodes.get(index[0]);
		OneHap n2=nodes.get(index[1]);
		double recomb_freq = Math.min(n1.freq, n2.freq) * analyzer.new_old_haps_beta.sample(); 
		double ori_freq_1=n1.freq;
		double ori_freq_2=n2.freq; // backup for returning to the original of rejected.
		int recomb_loc = ThreadLocalRandom.current().nextInt(0, this.num_snps);	// Needs the ThreadLocalRandom import.

		OneHap[] recomb_haps = OneHap.generate_recombinants(n1, n2, recomb_loc, recomb_freq);
		if(search_a_hap(recomb_haps[0].hap)!=-1)return 0;
		if(search_a_hap(recomb_haps[1].hap)!=-1)return 0;
		this.append(recomb_haps[0]);
		this.append(recomb_haps[1]);	
	    n1.change_my_freq(n1.freq - recomb_freq);
	    n2.change_my_freq(n2.freq - recomb_freq);

	    this.update_sigma_mu_current();
	    double logl2=Algebra.logL_aems(Algebra.times(analyzer.sigma, N), Algebra.times(analyzer.mu,N), analyzer.data);	// logL_rjmcmc(analyzer.sigma, analyzer.mu, analyzer.data, analyzer.num_hap_inpool)
	    double logHR=logl2-this.loglikelihood;	
	    // Treating this like the mutation part of add_mutant_or_coalesce . Why? Technically just rearranging the linkage, not introducing new mutations 
	    // or increasing concentration of a particular pattern of linkage.  
	    logHR-=gamma;
	    logHR+=(alpha-1)*(Math.log(recomb_haps[0].freq)+Math.log(n1.freq)-Math.log(ori_freq_1)); //prior for p, first recombinant
	    logHR+=(alpha-1)*(Math.log(recomb_haps[1].freq)+Math.log(n2.freq)-Math.log(ori_freq_2)); //prior for p, second recombinant	    
	    logHR+=-Math.log(this.num_curr_H+2.0)+Math.log(coalescing_probability(index[0], num_curr_H-2));//inverse transition probability, first recombinant
	    logHR+=-Math.log(this.num_curr_H+2.0)+Math.log(coalescing_probability(index[1], num_curr_H-2));//inverse transition probability, second recombinant
	    logHR+=2.0 * (Math.log(this.num_curr_H*this.num_snps*1.0)-(1-Math.log(analyzer.new_old_haps_beta.density(0)))); //transition probability, both recombinants
	    logHR+=Math.log(ori_freq_1); //|Jacobian|, first recombinant
	    logHR+=Math.log(ori_freq_2); //|Jacobian|, second recombinant

		if((!(Double.isNaN(logl2))) && Math.log(ThreadLocalRandom.current().nextDouble())<=logHR){//proposal accepted				
			this.loglikelihood=logl2;
			if(this.loglikelihood>analyzer.logl_best)analyzer.logl_best=this.loglikelihood;
			return 1;
		}
		else{//proposal rejected
		    n1.change_my_freq(ori_freq_1);
		    n2.change_my_freq(ori_freq_2);
			this.delete(num_curr_H-1);
			this.delete(num_curr_H-1);
			return -1;
		}
	}
}
