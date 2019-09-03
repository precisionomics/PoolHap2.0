package PoolHap;

import java.io.*;
import java.util.*;

public class Evaluate {
	public ArrayList<String> gold_haps= new ArrayList<String>();
	public ArrayList<String> compare_haps= new ArrayList<String>();
	
	public double mcc(String A_hap, String B_hap) {
//		MCC=(TP×TN−FP×FN)/ sqrt((TP+FP)*(TP+FN)*(TN+FP)* (TN+FN))
		double tp = 0.0;
		double tn = 0.0;
		double fp = 0.0;
		double fn = 0.0;
		for (int i = 0; i < A_hap.length(); i++) {
			if ((A_hap.substring(i, i+1).equals("1")) &&  
					(B_hap.substring(i, i+1).equals("1"))) {
				tp= tp +1 ;
			}else if ((A_hap.substring(i, i+1).equals("0")) &&  
					(B_hap.substring(i, i+1).equals("0"))){
				tn=tn+1 ;
			}else if ((A_hap.substring(i, i+1).equals("0")) && 
					(B_hap.substring(i, i+1).equals("1"))){
				fp=fp+1;
			}else if ((A_hap.substring(i, i+1).equals("1")) &&  
					(B_hap.substring(i, i+1).equals("0"))){
				fn=fn+1;
			}
		}
		double numerator= tp*tn-fp*fn;
		double denominator= Math.sqrt((tp+fp)*(tp+fn)*(tn+fp)* (tn+fn)) ;		
		return numerator/ denominator;
		
	}
	public void AemEvaluate(String gold_file, String aem_file) throws IOException {
		BufferedReader bufferedreader = new BufferedReader(new FileReader(gold_file));
        String line = "";
        while ((line = bufferedreader.readLine()) != null) {
        	line =line.replace("\n", "").replace("\r", "");
//        	System.out.println(line);
        	if (!line.startsWith("@")){
        		this.gold_haps.add(line );
        	}
        }
        bufferedreader.close();
        ArrayList<ArrayList<String >>  geno_2D = new ArrayList<ArrayList<String>>();
        BufferedReader bufferedreader2 = new BufferedReader(new FileReader(aem_file));
        while ((line = bufferedreader2.readLine()) != null) {
        	line =line.replace("\n", "").replace("\r", "");
        	if ((!line.startsWith("Hap_ID")) &&  (!line.startsWith("Freq") )){
        		String[] tmp = line.split("\t");
        		ArrayList<String> tmp_arr = new ArrayList<String>();
        		for (int i = 1; i < tmp.length; i++) {
        			tmp_arr.add(tmp[i]);
        		}
        		geno_2D.add(tmp_arr);
        	}
        }
        for (int j = 0; j < geno_2D.get(0).size(); j++) {
        	String tmp_str="";
        	for (int i = 0; i < geno_2D.size(); i++) {
        		tmp_str=tmp_str+ geno_2D.get(i).get(j);
        	}
        	this.compare_haps.add(tmp_str);
        }
        bufferedreader2.close();
        
//        for (int i = 0; i < this.compare_haps.size(); i++) {
//        	System.out.println(this.compare_haps.get(i));
//        }
        double total_mcc=0;
        for (int i = 0; i < this.gold_haps.size(); i++) {
        	double max_mcc= -1;
        	String nearest_hap="";
        	for (int j = 0; j < this.compare_haps.size(); j++) {
        		if( mcc(this.gold_haps.get(i), this.compare_haps.get(j)) >  max_mcc) {
        			nearest_hap= this.compare_haps.get(j);
        			max_mcc= mcc(this.gold_haps.get(i), this.compare_haps.get(j));
        		}
        	}
        	total_mcc += max_mcc;
        	System.out.println("Haplotype_"+ Integer.toString(i) );
        	System.out.println( this.gold_haps.get(i) );
        	System.out.println( "MaxMcc:\t"+ Double.toString(max_mcc)  );
        	System.out.println(  nearest_hap );
        	System.out.println( "-----------------------------------------");
        }
        System.out.println( "Average MCC:\t" + total_mcc/ Double.valueOf(this.gold_haps.size())  );
		return;
	}
	
	public void LassoEvaluate(String gold_file, String aem_file) throws IOException {
		return;
	}
	
	public Evaluate() throws IOException {
		
	}
}
