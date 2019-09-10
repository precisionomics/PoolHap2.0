package MiscFunctions;

import java.io.*;
import java.util.*;

import PoolHap.Entrance;

/**
 * @author Chen Cao 2019-08
 * Random Generate all files which could be used for PoolHapX testing.
 * Including the gold standard files, vef files, sample name file,
 * Each pool variant frequency file (In paratice, this file could be 
 * called by GATK/ SamTools). 
 * To simulate the evolution, each round n (last parameter) snps changes. 
 */  

public class RandomSimulator {
	String project_name;
	String folder_path;
	int num_pools;
	int num_total_haps;
	double ave_haps_inpool;
	int genome_len;
	int num_vars;
	int read_len;
	int outer_dist;
	int mut_each_generation;
	double ave_coverage;
	/*	1	0	1	0	0	0	1
	 * 	1	0	0	0	0	0	1
	 * 	0	0	1	0	0	1	0
	 * 3 genomes, the 1st one is 1010001
	 */	
	int[][] sim_genome;   // num_total_haps* num_vars
	double [][] hap_freq_inpool; // num_pools* num_total_haps 
	int[] var_pos; 
	
	public int  mut(int  base)  throws IOException {
		if (base==0 ) {
			return 1;
		}else {
			return 0;
		}
	}
	public void GenomeSimulator() throws IOException {
		this.sim_genome = new int [this.num_total_haps][this.num_vars];
		
		for (int i = 0; i < this.sim_genome[0].length; i++) {
			Random random = new Random(); 
			int genotype  = random.nextInt(2);
			this.sim_genome[0][i] = genotype ;
		}
		
		
		for (int i = 1; i < this.sim_genome.length; i++) {
			HashSet<Integer> pos_set = new HashSet<Integer>();
			for (int j = 0; j < this.sim_genome[i].length; j++) {
				this.sim_genome[i][j] = this.sim_genome[i-1][j];
			}
			int count =0;
			while (count < this.mut_each_generation) { 
				
				Random random = new Random(); 
				int pos   = random.nextInt(this.num_vars);
				if  (!pos_set.contains(pos)){
					count++;
					sim_genome[i][pos]=mut (sim_genome[i-1][pos] );
					pos_set.add(pos);
				}
				
			}
		}	
		
		
		Boolean[] array = new Boolean[this.genome_len];
		this.var_pos = new int [this.num_vars];
		Arrays.fill(array, Boolean.FALSE);
		int count =0;
		while (count < this.num_vars) {
			Random random = new Random(); 
			int value   = random.nextInt(this.genome_len- 3*this.read_len-this.outer_dist)+ 
					this.read_len;
			if (array[value]==false ) {
				count++ ;
				array[value]= true ;
			}
		}
		count=0;
		for (int i = 0; i < array.length; i++) {
			if (array[i]== true) {
				this.var_pos[count]=i+1;
				count=count+1;
			}
		}
		
		String file_path= this.folder_path+"/"+ this.project_name+"/gold_standard/"+
				this.project_name+ "_haps.txt";
		FileWriter mydata = new FileWriter(file_path,false);
        PrintWriter pw = new PrintWriter(mydata);
        String pw_line ="";
        for (int i = 0; i < this.sim_genome.length; i++) {
        	pw_line="@Haplotype_"+ Integer.toString(i);
        	pw.write(pw_line+"\n");
        	pw_line="";
        	for (int j = 0; j < this.sim_genome[i].length; j++) {
        		pw_line=  pw_line+ Integer.toString(this.sim_genome[i][j]) ; 
        	}
        	pw.write(pw_line+"\n");
        }
        pw.flush();
        pw.close();
		return;	
	}
	
	public void HapFreqSimulator() throws IOException {
		double [] total_value_arr= new double [this.num_pools];
		this.hap_freq_inpool = new double [this.num_pools][this.num_total_haps];
		for (int i = 0; i < this.hap_freq_inpool.length; i++) {
			double total_value = 0.0;
			for (int j = 0; j < this.hap_freq_inpool[i].length; j++) {
				Random random = new Random(); 
				int value   = random.nextInt(10*this.num_total_haps);
				if (value< 10*(this.num_total_haps- this.ave_haps_inpool )) {
					value =0;
				}
				this.hap_freq_inpool[i][j]= value;
				total_value+= Double.valueOf(value) ;
			}
			total_value_arr[i]= total_value;
		}
		
		for (int i = 0; i < this.hap_freq_inpool.length; i++) {
			for (int j = 0; j < this.hap_freq_inpool[i].length; j++) {
				if (this.hap_freq_inpool[i][j] / total_value_arr[i] > 0.0005) {
					this.hap_freq_inpool[i][j] =  (double)Math.round(this.hap_freq_inpool[i][j] / 
							total_value_arr[i]*1000)/1000 ;
				}else {
					this.hap_freq_inpool[i][j]= 0.0;
				}
			}
		}
		
		String file_path= this.folder_path+"/"+ this.project_name+"/gold_standard/"+
				this.project_name+ "_haps.intra_freq.txt";
		FileWriter mydata = new FileWriter(file_path,false);
        PrintWriter pw = new PrintWriter(mydata); 
        String pw_line = "Hap_ID";
        for (int i = 0; i < this.num_total_haps; i++) { 
        	pw_line = pw_line+  "\t"+ "h"+Integer.toString(i);
        }
        pw.write( pw_line + "\n"); 
        for (int i = 0; i < this.hap_freq_inpool.length; i++) {
        	pw_line = this.project_name+"p"+Integer.toString(i);
        	for (int j = 0; j < this.hap_freq_inpool[i].length; j++) {
        		pw_line =pw_line +"\t"+  Double.toString(this.hap_freq_inpool[i][j]);
        	}
        	pw.write( pw_line + "\n"); 
        }
        pw.flush();
        pw.close();
        
       
        file_path= this.folder_path+"/"+ this.project_name+"/gold_standard/"+
				this.project_name+ "_haps.inter_freq_haps.txt";
        FileWriter mydata2 = new FileWriter(file_path,false);
        PrintWriter pw2 = new PrintWriter(mydata2); 
        pw_line = "Hap_ID";
        for (int i = 0; i < this.num_total_haps; i++) { 
        	pw_line = pw_line+  "\t"+ "h"+Integer.toString(i);
        }
        pw2.write( pw_line + "\n"); 
        pw_line= "Freq";
        for (int j = 0; j < this.hap_freq_inpool[0].length; j++) {
        	double total_value = 0.0;
        	for (int i = 0; i < this.hap_freq_inpool.length; i++) {
        		total_value+=  this.hap_freq_inpool[i][j];
			}
        	pw_line=pw_line+"\t" + Double.toString(total_value/ Double.valueOf(this.num_pools)); 
		}
        pw2.write( pw_line + "\n"); 
        for (int j = 0; j < this.sim_genome[0].length; j++) {
        	pw_line = "0:"+ Integer.toString(this.var_pos[j])+";"+  Integer.toString(this.var_pos[j])
    		+";0:1";
        	for (int i = 0; i < this.sim_genome.length; i++) {
        		pw_line=pw_line + "\t"+ Integer.toString(sim_genome[i][j]) ; 
        	}
        	pw2.write( pw_line + "\n"); 
        }
        pw2.flush();
        pw2.close();
		return;
	}
	
	public void VefSimulator() throws IOException {
		new File(this.folder_path + "/"+this.project_name + "/intermediate/vef").mkdir();
		for (int i = 0; i < this.num_pools; i++) { 
			String file_path= this.folder_path+"/"+ this.project_name+"/"+
					"intermediate/vef/"+ this.project_name+ "_p"+Integer.toString(i)+".vef" ;
			FileWriter mydata = new FileWriter(file_path,false);
			PrintWriter pw = new PrintWriter(mydata); 
			String pw_line ="";
//			int[][] sim_genome;   // num_total_haps* num_vars
//			double [][] hap_freq_inpool; // num_pools* num_total_haps
			int [] hap_count = new int [this.num_total_haps ];
			int roulette_total =0;
			HashMap<Integer, Integer> hap_dict = new HashMap<Integer, Integer>(); 
			for (int j = 0; j < this.num_total_haps; j++) { 
				int count = (int) (hap_freq_inpool[i][j]*1000);
				hap_count[j]= count ;	
				for (int k = 0; k < count ; k++) {
					hap_dict.put(roulette_total, j);
					roulette_total++;
				}
			}
			
			int total_read = (int) ((this.ave_coverage* this.genome_len)/ (this.read_len*2 ));
			for (int t = 0; t < total_read; t++) { 
				Random random = new Random(); 
//				System.out.println(roulette_total);
				int value   = random.nextInt( roulette_total ) ;
				int hap = hap_dict.get(value);
				Random random_pos = new Random(); 
				int pos =  random_pos.nextInt(this.genome_len ) ;
//				@Haplotype_0_4618_4869_0_1_0_0_0:0:0_0:0:0_1c/1
				int start_1=pos;
				int end_1=pos+this.read_len;
				int start_2=pos+ this.read_len+ this.outer_dist;
				int end_2=pos+ 2*this.read_len+ this.outer_dist;
				pw_line="@Haplotype_"+ Integer.toString(hap)+"_"+ Integer.toString(t)+"\t";
				boolean flag= false;
				HashMap<Integer, Integer> var_dict = new HashMap<Integer, Integer>();
				for (int j = 0; j < this.num_vars; j++) { 
					var_dict.put(var_pos[j],sim_genome[hap][j] );
				}
				for (int j = start_1; j < end_1 ; j++) {
					if (var_dict.containsKey(j)) {
						pw_line=pw_line+ Integer.toString(j)+"="+var_dict.get(j)+";";
						flag=true;
					}
				}
				for (int j = start_2; j < end_2 ; j++) {
					if (var_dict.containsKey(j)) {
						pw_line=pw_line+ Integer.toString(j)+"="+var_dict.get(j)+";";
						flag=true;
					}
				}
				if (flag==true ) {
					pw_line=pw_line+ "\t"+ "//\t"+ Integer.toString(start_1)+"\t"+ 
							Integer.toString(end_1) +"\t"+
							Integer.toString(start_2)+"\t"+Integer.toString(end_2);
					pw.write( pw_line + "\n"); 
				}
			}
			pw.flush();
	        pw.close();
			
		}
		
	}
	
	public void VarFreqSimulator() throws IOException {
		String file_path= this.folder_path+"/"+ this.project_name+"/"+"intermediate/"
				+ this.project_name+ "_sample_names.txt";
		FileWriter mydata0 = new FileWriter(file_path,false);
        PrintWriter pw0 = new PrintWriter(mydata0); 
        for (int i = 0; i < this.num_pools; i++) { 
        	pw0.write( this.project_name +"_p"+ Integer.toString(i ) + "\n");  
        }
        pw0.flush();
        pw0.close();
        
        file_path= this.folder_path+"/"+ this.project_name+"/"+"intermediate/" +
				this.project_name+ "_vars.intra_freq.txt";
		FileWriter mydata = new FileWriter(file_path,false);
		PrintWriter pw = new PrintWriter(mydata); 
		String pw_line="Pool_ID";
		for (int i = 0; i < this.num_pools; i++) { 
			pw_line=pw_line+ "\t"+ this.project_name +"p"+ Integer.toString(i ) ;
		}
		pw.write(pw_line+"\n");
		for (int i = 0; i < this.num_vars; i++) {
			pw_line = "0;"+ Integer.toString(this.var_pos[i])+";"+  Integer.toString(this.var_pos[i])
    		+";0:1";
			for (int j = 0; j < this.num_pools; j++) {
				double freq = 0.0;
//				int[][] sim_genome;   // num_total_haps* num_vars
//				double [][] hap_freq_inpool; // num_pools* num_total_haps
				for (int k = 0; k < this.num_total_haps; k++) {
					freq+= hap_freq_inpool[j][k]* sim_genome [k][i];
				}
				freq= (double)Math.round(freq *1000)/1000;
				pw_line=pw_line + "\t"+ Double.toString(freq )  ;
			}
			pw.write(pw_line+"\n");
		}
        
		pw.flush();
        pw.close();
	}
	
	public RandomSimulator(String[] args) throws IOException {
		this.folder_path= args[0];
		this.project_name= args[1];
		this.num_pools = Integer.parseInt(args[2]);
		this.num_total_haps = Integer.parseInt(args[3]);
		this.ave_haps_inpool= Double.parseDouble( args[4]);
		this.genome_len= Integer.parseInt(args[5]);
		this.num_vars = Integer.parseInt(args[6]);
		this.read_len = Integer.parseInt(args[7]);
		this.outer_dist = Integer.parseInt(args[8]);
		this.ave_coverage= Double.parseDouble( args[9]);
		this.mut_each_generation = Integer.parseInt(args[10]);
		new File(this.folder_path + "/"+this.project_name ).mkdir();
		new File(this.folder_path + "/"+this.project_name + "/"+ "/gold_standard/").mkdir();
		new File(this.folder_path + "/"+this.project_name + "/"+ "/input/").mkdir();
		new File(this.folder_path + "/"+this.project_name + "/"+ "/output/").mkdir();
		new File(this.folder_path + "/"+this.project_name + "/"+ "/intermediate/").mkdir();
		GenomeSimulator (); 
		System.out.println(" Genome Simulation Finished!");
		HapFreqSimulator();
		System.out.println(" Haplotype Frequency for Each Pool Simulation Finished!");
		VarFreqSimulator();
		System.out.println(" Variants Frequency for Each Haplotype Simulation Finished!");
		VefSimulator();
		System.out.println(" Vef File Simulation Finished!");
	}
	

	
	public static  void main(String[] args) throws IOException {
		// /home/chencao/Desktop  sim001	10	20	15	5000	50	150	50	100		5
		//parameter 0: folder path
		//parameter 1: project name
		//parameter 2: number of pools (or generations)
		//parameter 3: number of total haplotypes
		//parameter 4: average number of haplotypes in each pool
		//parameter 5: genome length 
		//parameter 6: number of variants
		//parameter	7: read length
		//parameter 8: outer distance between the two ends for pairs
		//parameter 9: average coverage
		//parameter 10: number of mutations each generation
		
		RandomSimulator rs = new RandomSimulator(args);
	}
	

}
