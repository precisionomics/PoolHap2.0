package PoolHap;

import java.util.concurrent.ThreadLocalRandom;

public class OneHap {
	int[] hap;
	double freq;
	double[] posteriori;
	double psum;
	double p2sum;
	
	OneHap next;
	
	public OneHap(int[] hap, double freq){
		this.hap=hap;
		this.freq=freq;
	}	
	
	public OneHap(int[] hap, double freq, double psum, double p2sum){
		this.hap=hap;
		this.freq=freq;
		this.psum=psum;
		this.p2sum=p2sum;
	}
	
	public void setNext(OneHap next){
		this.next=next;
	}
	
	public void mutate_myself(int loc_index, int allele){
		if(this.hap[loc_index]==allele){
			System.out.println("The mutated allele is the same as the original: Error?");
		}
		this.hap[loc_index]=allele;
	}
	
	public void change_my_freq(double new_fereq){
		this.freq=new_fereq;
	}
	
	public OneHap generate_mutant(int loc_index, int allele, double new_fereq){
		if(this.hap[loc_index]==allele){
			System.out.println("The mutated allele is the same as the original: Error?");
			return null;
		}
		OneHap new_hap=new OneHap(this.hap.clone(), new_fereq);
		new_hap.hap[loc_index]=allele;
		return new_hap;
	}	
	
	public OneHap generate_mutant(double new_fereq){
		int loc_index=(int) (ThreadLocalRandom.current().nextDouble()*this.hap.length);
		if(loc_index==this.hap.length)loc_index--;
		int allele=0;
		if(this.hap[loc_index]==0)allele=1; //TODO if more than two alleles!! 
		return generate_mutant(loc_index, allele, new_fereq);
	}
	
	/*
	 *  recombine two haps between loc and loc+1
	 */
	public static OneHap[] generate_recombinants(OneHap h1, OneHap h2, int loc, double recomb_freq){
		if(h1.hap.length!=h2.hap.length){
			System.out.println("generate_recombinants: h1.hap.length!=h2.hap.length");
			return null;
		}
		int num_snps=h1.hap.length;
		int[] hap1=new int[num_snps];
		int[] hap2=new int[num_snps];
		for(int k=0;k<=loc;k++){
			hap1[k]=h1.hap[k];
			hap2[k]=h2.hap[k];
		}for(int k=loc+1;k<num_snps;k++){
			hap1[k]=h2.hap[k];
			hap2[k]=h1.hap[k];
		}
		OneHap[] recombinants=new OneHap[2];
		recombinants[0]=new OneHap(hap1, recomb_freq);
		recombinants[1]=new OneHap(hap2, recomb_freq);
		return recombinants;
	}
	
	String output(String sep){
		String out=""+this.hap[0];
		for(int k=1;k<this.hap.length;k++)out=out+sep+this.hap[k];
		out=out+"\t"+this.freq;
		return out;
	}

	String output_nofreq(String sep){
		String out=""+this.hap[0];
		for(int k=1;k<this.hap.length;k++)out=out+sep+this.hap[k];
		out=out+"; ";
		return out;
	}
}
