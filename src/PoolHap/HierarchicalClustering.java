package PoolHap;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.OutputStreamWriter;
import java.io.UnsupportedEncodingException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;

/**
 * @author Chen Cao 2019-09
 * Hierarchical Clustering
 */  

public class HierarchicalClustering {
	
	public HierarchicalClustering(String[] haps_ID, String[][] global_haps_string,
			HashMap<String,Integer> hapID2index , String outputfile, double disThreshold)
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
		HierarchicalClustering.writeClusterToFile(outputfile, cluste_res, haplist, index2hapID);
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
	
	public static double getSimilarity(String[] x, 
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
	public static void writeClusterToFile(String clusterFilePath,List<Cluster>clusters, 
			ArrayList<String> sentencelist,HashMap<Integer, String> index2hapID ) throws IOException{
		File f = new File(clusterFilePath);
		
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
				for(int tempdp:tempDps){				
					try { 				
						bw.write(index2hapID.get(tempdp) + ":\t" + sentencelist.get(tempdp));
						bw.newLine();
						bw.flush();
					} catch (IOException e) {
					// TODO Auto-generated catch block
						e.printStackTrace();
					}			
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
}
