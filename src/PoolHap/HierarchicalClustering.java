package PoolHap;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileOutputStream;
import java.io.FileReader;
import java.io.IOException;
import java.io.OutputStreamWriter;
import java.io.UnsupportedEncodingException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;

/**
 * @author Chen Cao: 2019-09
 * Hierarchical Clustering
 */  

public class HierarchicalClustering {
	
	public HashMap<String, Integer> output_ref_arr ;
    public HashMap<String, String> conf_ref_arr;
    public int num_loci;
    public int num_pools;
    public LocusAnnotation[] locusInfo; // # of loci, note that a locus can be a SNP or a region
    public double[][] inpool_site_freqs; // # of loci x # of pools; Added by Quan Dec. 2018.
    public HashMap<Integer, Integer> loci_index_dict;
    
	
	public HierarchicalClustering(String[] haps_ID, String[][] global_haps_string,
			HashMap<String,Integer> hapID2index , String outputfile,  double disThreshold)
					throws IOException{
		
		ArrayList<String> haplist = new ArrayList<>();
		ArrayList<Integer> hapindex = new ArrayList<>();
		int count=0;
		HashMap<Integer, String> index2hapID = new HashMap<Integer, String>();
		for (int i=0; i<global_haps_string.length;i++ ) {
			index2hapID.put(count, haps_ID[i]);
			String tmp ="";
			for (int j=0; j<global_haps_string[i].length;j++ ) {
				tmp =tmp+ global_haps_string[i][j];
			}
			haplist.add(tmp ); 
			hapindex.add(count);
			count++;
		}
		
		ArrayList<String[]> hapvec =  new ArrayList<String[]> ();
		for (int i=0;i< global_haps_string.length; i++) {
			hapvec.add(global_haps_string[i].clone() );
		}
		
		double[][] simmatrix = this.CalSimMatrix(hapvec);
		List<Cluster> cluste_res = this.starAnalysis(hapindex, simmatrix, disThreshold);
		this.writeClusterToFile(outputfile, cluste_res, haplist, index2hapID);
	}
	
	
	public  double[][] CalSimMatrix(ArrayList<String[]> vectorlist) {
		double[][] sim = new double[vectorlist.size()][vectorlist.size()];
		for (int i = 0; i < vectorlist.size(); i++) {
			String[] vec1 = vectorlist.get(i);
			for (int j = i + 1; j < vectorlist.size(); j++) {
				String[] vec2 = vectorlist.get(j);
				sim[i][j] = getSimilarity(vec1, vec2);
				sim[j][i] = sim[i][j];
			}
		}
		return sim;
	}
	
	public  double getSimilarity(String[] x, 
			String[] y) {
		if (x.length!= y.length) {
			System.out
            .println("Error: The haplotypes have different length for clustering.\n");
			System.exit(0);
		}
		double len = (double )x.length;
		double iden= 0;
		for (int i=0;i < x.length;i++) {
			if (x[i].equals(y[i] )){
				iden=iden+ 1;
			}
		}
		return iden/ len;
	}

	//Initialization
	public  List<Cluster> initialCluster(ArrayList<Integer> CatalogList){
		List<Cluster> originalClusters = new ArrayList<Cluster>();
		for(int i = 0; i < CatalogList.size(); i ++){
			int index = CatalogList.get(i);
			List<Integer> tempDataPoints = new ArrayList<Integer>();
			tempDataPoints.add(index);  
            Cluster tempCluster = new Cluster();  
            tempCluster.setClusterName("Cluster " + String.valueOf(i));  
            tempCluster.setDataPoints(tempDataPoints);  
            originalClusters.add(tempCluster);  
		}
		
		return originalClusters;
	}
	
	
	 /** 
     * Merge two clusters
     * @param clusters 
     * @param mergeIndexA 
     * @param mergeIndexB 
     * @return 
     */  
	 private List<Cluster> mergeCluster(List<Cluster> clusters, int mergeIndexA,  int mergeIndexB) {  
	        if (mergeIndexA != mergeIndexB) {  
	           
	            Cluster clusterA = clusters.get(mergeIndexA);  
	            Cluster clusterB = clusters.get(mergeIndexB);  
	            List<Integer> dpA = clusterA.getDataPoints();  
	            List<Integer> dpB = clusterB.getDataPoints(); 
	            for (int dp : dpB) {
	                dpA.add(dp); 
	            }	            
	            clusterA.setDataPoints(dpA); 
	            clusters.remove(mergeIndexB);  
	        }  
	        return clusters;  
	    }  
	
	/**
	 * Clustering the haplotypes
	 * @param CatalogList
	 * @param sim
	 * @param threshold
	 * @return
	 * @throws IOException
	 */
	public List<Cluster> starAnalysis(ArrayList<Integer> CatalogList, double[][] sim,double threshold)
			throws IOException{
		List<Cluster> finalClusters = new ArrayList<Cluster>();
		List<Cluster> originalClusters = initialCluster(CatalogList);
		finalClusters = originalClusters;
		
		while(finalClusters.size() >  1){
			double max = Double.MIN_VALUE; 
			int mergeIndexA = 0;  
            int mergeIndexB = 0;
            for (int i = 0; i < finalClusters.size(); i++) {  
                for (int j = i + 1; j < finalClusters.size(); j++) {   
                    Cluster clusterA = finalClusters.get(i);//
                    Cluster clusterB = finalClusters.get(j);//
                    List<Integer> dataPointsA = clusterA.getDataPoints();  
                    List<Integer> dataPointsB = clusterB.getDataPoints();  
                          
                    //
                    double minTempWeight = Double.MAX_VALUE;
                    boolean flag = false;
                    for (int m = 0; m < dataPointsA.size(); m++) {  
                    	for (int n = 0; n < dataPointsB.size(); n++) {//
                            double tempWeight = sim[dataPointsA.get(m)][dataPointsB.get(n)];     
                            if(tempWeight < minTempWeight){
                                minTempWeight = tempWeight;
                                flag = true;
                            }
                                                               
                        }//end_for  
                    }//end_for
                                              
                        
                    if (minTempWeight > max && flag == true) {  
                        max = minTempWeight;  
                        mergeIndexA = i;// --
                        mergeIndexB = j;// --
                    }//end_if  
                } // end for j  
            }// end for i  
            //Merge cluster[mergeIndexA] and cluster[mergeIndexB]
            if(max < threshold) {  
//                System.out.println("The similarity is le！");  
                break;  
            }  
            finalClusters = mergeCluster(finalClusters, mergeIndexA, mergeIndexB);  
        }// end while  
        return finalClusters;  
           		
	}
	
	
/** 
 * Output clustering results
 * @param clusterFilePath 
 * @param clusters 
 * @throws FileNotFoundException 
 * @throws UnsupportedEncodingException 
 */
	public  void writeClusterToFile(String clusterFilePath,List<Cluster>clusters, 
			ArrayList<String> sentencelist,HashMap<Integer, String> index2hapID ) throws IOException{
		File f = new File(clusterFilePath);
		this.output_ref_arr = new HashMap<String,Integer>();
		BufferedWriter bw;
		int count0 = 0;
		try {
			bw = new BufferedWriter(new OutputStreamWriter(new FileOutputStream(f),"UTF-8"));
			for (Cluster cl : clusters){
			
				count0++;
				try {
					bw.write("Cluster:\t" + count0);
					bw.newLine();
					bw.flush();	
				} catch (IOException e1) {
					// TODO Auto-generated catch block
					e1.printStackTrace();
				}
				
				List<Integer> tempDps = cl.getDataPoints();
				ArrayList<String>  hap_list  = new  ArrayList<String>();
				for(int tempdp:tempDps){				
					try { 				
						bw.write(index2hapID.get(tempdp) + ":\t" + sentencelist.get(tempdp));
						hap_list.add(sentencelist.get(tempdp));
						bw.newLine();
						bw.flush();
					} catch (IOException e) {
					// TODO Auto-generated catch block
						e.printStackTrace();
					}			
				}
				if (hap_list.size()>0){
					double max_sim= -1.0;
					int index =0 ;
					for (int i=0;i<hap_list.size();i++ ) {
						double i_sim =0.0;
						for (int j=0;j<hap_list.size();j++ ) {
							String[] x = hap_list.get(i).split("");
							String[] y = hap_list.get(j).split("");
							i_sim+=this.getSimilarity(x,y);
			
						}
						if (i_sim> max_sim) {
							max_sim= i_sim;
							index=i;
						}
					}
//					System.out.println(hap_list.get(index));
//					this.output_ref_arr.put("1111111110000", 1);
					this.output_ref_arr.put(hap_list.get(index), 1);	
				}
			}
		} catch (UnsupportedEncodingException e1) {
		// TODO Auto-generated catch block
			e1.printStackTrace();
		} catch (FileNotFoundException e1) {
		// TODO Auto-generated catch block
			e1.printStackTrace();
		}
		System.out.println("Number of Clusters "
				+ ":" + clusters.size());
	}
	
	public void gc_solver(String gs_var_pos ) throws IOException {
        /*
         *  Initialize variables.
         */
        // HashMap<Seg_Site, Index>
        HashMap<Integer, Integer> pos_dict = new HashMap<Integer, Integer>();
        Integer pos_index = 0;

        /*
         *  Read through gold standard variant positions file to load into dictionary and count
         *  total number of loci.
         *
         *  also this needs to be done in SiteInPoolFreqAnno as well, refactor to a static method?
         *  retain number of loci as a global variable?
         */
        // Read into file.
        BufferedReader br = new BufferedReader(new FileReader(gs_var_pos));
        String currLine = br.readLine(); // skip header
        currLine = br.readLine();

        // Read each line into a dictionary.
        while (currLine != null) {
            pos_dict.put(Integer.parseInt(currLine.split(";")[1]), pos_index);

            // TODO: [LEFTOVER]
            // System.out.println(currLine.split(";")[1] + "\t" +  pos_index);

            pos_index++; // move to next variant position index
            currLine = br.readLine(); // read next line
        }
        br.close();
        this.loci_index_dict = pos_dict;
        this.num_loci = pos_dict.size(); // number of loci in gold standard variants file
        

        /*
         *  Read through gold standard variant positions file again to load loci info into a matrix.
         */
        // Read into file.
        br = new BufferedReader(new FileReader(gs_var_pos));
        currLine = br.readLine(); // skip header
        currLine = br.readLine();

        // Initialize variables.
        int loci_index = 0;
        this.num_pools = currLine.split("\t").length - 1;
        this.locusInfo = new LocusAnnotation[this.num_loci];
        this.inpool_site_freqs = new double[this.num_loci][this.num_pools];
        
        while (currLine != null) {
            String[] tmp = currLine.split("\t");

            // The first column is the locus-info.
            this.locusInfo[loci_index] = new LocusAnnotation(tmp[0]);

            // For each pool, load locus frequency into matrix.
            for (int p = 0; p < this.num_pools; p++) {
                this.inpool_site_freqs[loci_index][p] = Double.parseDouble(tmp[p + 1]);
            }

            loci_index++; // move to next locus index
            currLine = br.readLine(); // read next line
        }

	}
	
	
	public HapConfig hapOut(String[] pool_IDs) {
        //
        int num_global_hap = this.output_ref_arr.size();
        String[][] global_haps_string = new String[num_global_hap][num_loci];
        int[] global_haps_ct = new int[num_global_hap];
        int tot_hap_ct = 0;
        int hap_index = 0;

        //
        for (String entry : this.output_ref_arr.keySet()) {
            String[] var_comp = entry.split("");
            for (int v = 0; v < num_loci; v++) {
                global_haps_string[hap_index][v] = var_comp[v];
            }

            global_haps_ct[hap_index] = this.output_ref_arr.get(entry);
            tot_hap_ct += this.output_ref_arr.get(entry);
            hap_index++;
        }

        //
        double[] global_haps_freq = new double[num_global_hap];
        System.out.println("There are " + num_global_hap + " center haplotypes "
        		+ "generated by HierarchicalClustering.");
        for (int h = 0; h < num_global_hap; h++) {
            global_haps_freq[h] = (double) global_haps_ct[h] / (double) tot_hap_ct;
        }

        return new HapConfig(
            global_haps_string,
            global_haps_freq,
            null,
            this.inpool_site_freqs,
            this.locusInfo,
            this.num_pools,
            null,
            pool_IDs,
            0);

    }
}
