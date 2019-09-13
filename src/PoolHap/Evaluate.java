package PoolHap;

import java.io.*;
import java.util.*;

public class Evaluate {
	public ArrayList<String> gold_haps= new ArrayList<String>();
	public ArrayList<String> compare_haps= new ArrayList<String>();
	public String proj_name; 
	
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
	
	public int NumofMismatch (String x, String y) throws IOException {
		int mismatch= 0;
		for (int i=0;i < x.length();i++){
			if (!x.substring(i, i+1).equals( y.substring(i, i+1))) {
				mismatch++ ;
			}
		}
		return mismatch;
	}
			
	
	public void AemEvaluate(String pj_name, String dc_file, String gold_file, 
			String aem_folder)  throws IOException {
		this.proj_name= pj_name;
		this.gold_haps.clear();
        this.compare_haps.clear();
//		int num_match = 0;
		int num_mismatch = 0;
		ArrayList<Integer > level_I_start = new ArrayList<Integer>();
		ArrayList<Integer > level_I_end = new ArrayList<Integer>();
		ArrayList<Integer > level_II_start = new ArrayList<Integer>();
		ArrayList<Integer > level_II_end = new ArrayList<Integer>();
		int count = 0;
		BufferedReader bufferedreader = new BufferedReader(new FileReader(dc_file));
        String line = "";
        while ((line = bufferedreader.readLine()) != null) {
        	line =line.replace("\n", "").replace("\r", "");
//        	System.out.println(line);
        	if (count== 1) {
        		String[] tmp = line.split("\t");
        		for (int i = 0; i < tmp.length; i++) {
        			String[] tmp2= tmp[i].split(":");
        			level_I_start.add(Integer.parseInt(tmp2[0] )) ;
        			level_I_end.add(Integer.parseInt(tmp2[1] )) ;
        		}
        	}else if (count== 3) {
        		String[] tmp = line.split("\t");
        		for (int i = 0; i < tmp.length; i++) {
        			String[] tmp2= tmp[i].split(":");
        			level_II_start.add(Integer.parseInt(tmp2[0] )) ;
        			level_II_end.add(Integer.parseInt(tmp2[1] )) ;
        		}
        	}
        	count ++;
        }
        bufferedreader.close();
        
        BufferedReader bufferedreader2 = new BufferedReader(new FileReader(gold_file));
        line = "";
        while ((line = bufferedreader2.readLine()) != null) {
        	line =line.replace("\n", "").replace("\r", "");
//        	System.out.println(line);
        	if (!line.startsWith("@")){
        		this.gold_haps.add(line );
        	}
        }
        bufferedreader2.close();
        
        
        
        
        for (int i = 0; i < level_I_start.size(); i++) {
        	String aem_file= aem_folder+ this.proj_name+"_level_1_region_"+Integer.toString(i)
        	+".inter_freq_haps.txt";
        	ArrayList<ArrayList<String >>  geno_2D = new ArrayList<ArrayList<String>>();
        	BufferedReader br = new BufferedReader(new FileReader(aem_file));
        	while ((line = br.readLine()) != null) {
            	line =line.replace("\n", "").replace("\r", "");
            	if ((!line.startsWith("Hap_ID")) &&  (!line.startsWith("Freq") )){
            		String[] tmp = line.split("\t");
            		ArrayList<String> tmp_arr = new ArrayList<String>();
            		for (int j = 1; j < tmp.length; j++) {
            			tmp_arr.add(tmp[j]);
            		}
            		geno_2D.add(tmp_arr);
            	}
            }
        	br.close();
        	
        	for (int j = 0; j < geno_2D.get(0).size(); j++) {
            	String tmp_str="";
            	for (int k = 0; k < geno_2D.size(); k++) {
            		tmp_str=tmp_str+ geno_2D.get(k).get(j);
            	}
            	this.compare_haps.add(tmp_str);
            }
        	
//        	for (int j = 0; j < this.compare_haps.size(); j++) {
//        		System.out.println(this.compare_haps.get(j));
//        		
//        	}
        	for (int j = 0; j < this.gold_haps.size(); j++) {
            	int  max_mismatch=  level_I_end.get(i)- level_I_start.get(i)+1;
            	for (int k = 0; k < this.compare_haps.size(); k++) {
            		if( NumofMismatch (this.gold_haps.get(j).substring(level_I_start.get(i), level_I_end.get(i)+1)
            				, this.compare_haps.get(k)) <  max_mismatch) {
            			max_mismatch= NumofMismatch(this.gold_haps.get(j).substring(level_I_start.get(i), 
            					level_I_end.get(i)+1), this.compare_haps.get(k));
            		}
            	}
            	System.out.println(this.gold_haps.get(j).substring(level_I_start.get(i), level_I_end.get(i)+1)
            			+"\tFor Level I: Region "+ Integer.toString(i)+": Mismatch: "
            			+Integer.toString(max_mismatch));
            	if (max_mismatch!=0 ) {
            		num_mismatch++;
            	}
            }
        	
        	this.compare_haps.clear();	
        }
        
        for (int i = 0; i < level_II_start.size(); i++) {
        	String aem_file= aem_folder+ this.proj_name+"_level_2_region_"+Integer.toString(i)
        	+".inter_freq_haps.txt";
        	ArrayList<ArrayList<String >>  geno_2D = new ArrayList<ArrayList<String>>();
        	BufferedReader br = new BufferedReader(new FileReader(aem_file));
        	while ((line = br.readLine()) != null) {
            	line =line.replace("\n", "").replace("\r", "");
            	if ((!line.startsWith("Hap_ID")) &&  (!line.startsWith("Freq") )){
            		String[] tmp = line.split("\t");
            		ArrayList<String> tmp_arr = new ArrayList<String>();
            		for (int j = 1; j < tmp.length; j++) {
            			tmp_arr.add(tmp[j]);
            		}
            		geno_2D.add(tmp_arr);
            	}
            }
        	br.close();
        	
        	for (int j = 0; j < geno_2D.get(0).size(); j++) {
            	String tmp_str="";
            	for (int k = 0; k < geno_2D.size(); k++) {
            		tmp_str=tmp_str+ geno_2D.get(k).get(j);
            	}
            	this.compare_haps.add(tmp_str);
            }
        	
//        	for (int j = 0; j < this.compare_haps.size(); j++) {
//        		System.out.println(this.compare_haps.get(j));
//        		
//        	}
        	for (int j = 0; j < this.gold_haps.size(); j++) {
            	int  max_mismatch=  level_II_end.get(i)- level_II_start.get(i)+1;
            	for (int k = 0; k < this.compare_haps.size(); k++) {
            		if( NumofMismatch (this.gold_haps.get(j).substring(level_II_start.get(i), level_II_end.get(i)+1)
            				, this.compare_haps.get(k)) <  max_mismatch) {
            			max_mismatch= NumofMismatch(this.gold_haps.get(j).substring(level_II_start.get(i), 
            					level_II_end.get(i)+1), this.compare_haps.get(k));
            		}
            	}
            	System.out.println(this.gold_haps.get(j).substring(level_II_start.get(i), level_II_end.get(i)+1)
            			+"\tFor Level II: Region "+ Integer.toString(i)+": Mismatch: "
            			+Integer.toString(max_mismatch));
            	if (max_mismatch!=0 ) {
            		num_mismatch++;
            	}
            }
        	
        	this.compare_haps.clear();	
        }
        

        System.out.println("****The total number of disagree haplotypes is:\t"+ 
        		Integer.toString(num_mismatch)) ;
        
        
       
		return;
		
	}
	
	
	public void Level_III_AemEvaluate(String pj_name, String dc_file, String gold_file, 
			String aem_folder)  throws IOException {
		this.proj_name= pj_name;
		this.gold_haps.clear();
        this.compare_haps.clear();
//		int num_match = 0;
		int num_mismatch = 0;
		ArrayList<Integer > level_III_start = new ArrayList<Integer>();
		ArrayList<Integer > level_III_end = new ArrayList<Integer>();
		ArrayList<Integer > level_II_start = new ArrayList<Integer>();
		ArrayList<Integer > level_II_end = new ArrayList<Integer>();
		int count = 0;
		BufferedReader bufferedreader = new BufferedReader(new FileReader(dc_file));
        String line = "";
        while ((line = bufferedreader.readLine()) != null) {
        	line =line.replace("\n", "").replace("\r", "");
//        	System.out.println(line);
        	if (count== 5) {
        		String[] tmp = line.split("\t");
        		for (int i = 0; i < tmp.length; i++) {
        			String[] tmp2= tmp[i].split(":");
        			level_III_start.add(Integer.parseInt(tmp2[0] )) ;
        			level_III_end.add(Integer.parseInt(tmp2[1] )) ;
        		}
        	}else if (count== 3) {
        		String[] tmp = line.split("\t");
        		for (int i = 0; i < tmp.length; i++) {
        			String[] tmp2= tmp[i].split(":");
        			level_II_start.add(Integer.parseInt(tmp2[0] )) ;
        			level_II_end.add(Integer.parseInt(tmp2[1] )) ;
        		}
        	}
        	count ++;
        }
        bufferedreader.close();
        
        BufferedReader bufferedreader2 = new BufferedReader(new FileReader(gold_file));
        line = "";
        while ((line = bufferedreader2.readLine()) != null) {
        	line =line.replace("\n", "").replace("\r", "");
//        	System.out.println(line);
        	if (!line.startsWith("@")){
        		this.gold_haps.add(line );
        	}
        }
        bufferedreader2.close();
        
        
        
        
        for (int i = 0; i < level_III_start.size(); i++) {
        	String aem_file= aem_folder+ this.proj_name+"_level_3_region_"+Integer.toString(i)
        	+".inter_freq_haps.txt";
        	ArrayList<ArrayList<String >>  geno_2D = new ArrayList<ArrayList<String>>();
        	BufferedReader br = new BufferedReader(new FileReader(aem_file));
        	while ((line = br.readLine()) != null) {
            	line =line.replace("\n", "").replace("\r", "");
            	if ((!line.startsWith("Hap_ID")) &&  (!line.startsWith("Freq") )){
            		String[] tmp = line.split("\t");
            		ArrayList<String> tmp_arr = new ArrayList<String>();
            		for (int j = 1; j < tmp.length; j++) {
            			tmp_arr.add(tmp[j]);
            		}
            		geno_2D.add(tmp_arr);
            	}
            }
        	br.close();
        	
        	for (int j = 0; j < geno_2D.get(0).size(); j++) {
            	String tmp_str="";
            	for (int k = 0; k < geno_2D.size(); k++) {
            		tmp_str=tmp_str+ geno_2D.get(k).get(j);
            	}
            	this.compare_haps.add(tmp_str);
            }
        	
//        	for (int j = 0; j < this.compare_haps.size(); j++) {
//        		System.out.println(this.compare_haps.get(j));
//        		
//        	}
        	for (int j = 0; j < this.gold_haps.size(); j++) {
            	int  max_mismatch=  level_III_end.get(i)- level_III_start.get(i)+1;
            	for (int k = 0; k < this.compare_haps.size(); k++) {
            		if( NumofMismatch (this.gold_haps.get(j).substring(level_III_start.get(i), level_III_end.get(i)+1)
            				, this.compare_haps.get(k)) <  max_mismatch) {
            			max_mismatch= NumofMismatch(this.gold_haps.get(j).substring(level_III_start.get(i), 
            					level_III_end.get(i)+1), this.compare_haps.get(k));
            		}
            	}
            	System.out.println(this.gold_haps.get(j).substring(level_III_start.get(i), level_III_end.get(i)+1)
            			+"\tFor Level III: Region "+ Integer.toString(i)+": Mismatch: "
            			+Integer.toString(max_mismatch));
            	if (max_mismatch!=0 ) {
            		num_mismatch++;
            	}
            }
        	
        	this.compare_haps.clear();	
        }
        
       
        

        System.out.println("****The total number of disagree haplotypes is:\t"+ 
        		Integer.toString(num_mismatch)) ;
        
        
       
		return;
		
	}
	
	public void Level_IV_AemEvaluate(String pj_name, String dc_file, String gold_file, 
			String aem_folder)  throws IOException {
		this.proj_name= pj_name;
		this.gold_haps.clear();
        this.compare_haps.clear();
//		int num_match = 0;
		int num_mismatch = 0;
		ArrayList<Integer > level_III_start = new ArrayList<Integer>();
		ArrayList<Integer > level_III_end = new ArrayList<Integer>();
		ArrayList<Integer > level_II_start = new ArrayList<Integer>();
		ArrayList<Integer > level_II_end = new ArrayList<Integer>();
		int count = 0;
		BufferedReader bufferedreader = new BufferedReader(new FileReader(dc_file));
        String line = "";
        while ((line = bufferedreader.readLine()) != null) {
        	line =line.replace("\n", "").replace("\r", "");
//        	System.out.println(line);
        	if (count== 7) {
        		String[] tmp = line.split("\t");
        		for (int i = 0; i < tmp.length; i++) {
        			String[] tmp2= tmp[i].split(":");
        			level_III_start.add(Integer.parseInt(tmp2[0] )) ;
        			level_III_end.add(Integer.parseInt(tmp2[1] )) ;
        		}
        	}else if (count== 3) {
        		String[] tmp = line.split("\t");
        		for (int i = 0; i < tmp.length; i++) {
        			String[] tmp2= tmp[i].split(":");
        			level_II_start.add(Integer.parseInt(tmp2[0] )) ;
        			level_II_end.add(Integer.parseInt(tmp2[1] )) ;
        		}
        	}
        	count ++;
        }
        bufferedreader.close();
        
        BufferedReader bufferedreader2 = new BufferedReader(new FileReader(gold_file));
        line = "";
        while ((line = bufferedreader2.readLine()) != null) {
        	line =line.replace("\n", "").replace("\r", "");
//        	System.out.println(line);
        	if (!line.startsWith("@")){
        		this.gold_haps.add(line );
        	}
        }
        bufferedreader2.close();
        
        
        
        
        for (int i = 0; i < level_III_start.size(); i++) {
        	String aem_file= aem_folder+ this.proj_name+"_level_4_region_"+Integer.toString(i)
        	+".inter_freq_haps.txt";
        	ArrayList<ArrayList<String >>  geno_2D = new ArrayList<ArrayList<String>>();
        	BufferedReader br = new BufferedReader(new FileReader(aem_file));
        	while ((line = br.readLine()) != null) {
            	line =line.replace("\n", "").replace("\r", "");
            	if ((!line.startsWith("Hap_ID")) &&  (!line.startsWith("Freq") )){
            		String[] tmp = line.split("\t");
            		ArrayList<String> tmp_arr = new ArrayList<String>();
            		for (int j = 1; j < tmp.length; j++) {
            			tmp_arr.add(tmp[j]);
            		}
            		geno_2D.add(tmp_arr);
            	}
            }
        	br.close();
        	
        	for (int j = 0; j < geno_2D.get(0).size(); j++) {
            	String tmp_str="";
            	for (int k = 0; k < geno_2D.size(); k++) {
            		tmp_str=tmp_str+ geno_2D.get(k).get(j);
            	}
            	this.compare_haps.add(tmp_str);
            }
        	
//        	for (int j = 0; j < this.compare_haps.size(); j++) {
//        		System.out.println(this.compare_haps.get(j));
//        		
//        	}
        	for (int j = 0; j < this.gold_haps.size(); j++) {
            	int  max_mismatch=  level_III_end.get(i)- level_III_start.get(i)+1;
            	for (int k = 0; k < this.compare_haps.size(); k++) {
            		if( NumofMismatch (this.gold_haps.get(j).substring(level_III_start.get(i), level_III_end.get(i)+1)
            				, this.compare_haps.get(k)) <  max_mismatch) {
            			max_mismatch= NumofMismatch(this.gold_haps.get(j).substring(level_III_start.get(i), 
            					level_III_end.get(i)+1), this.compare_haps.get(k));
            		}
            	}
            	System.out.println(this.gold_haps.get(j).substring(level_III_start.get(i), level_III_end.get(i)+1)
            			+"\tFor Level IV: Region "+ Integer.toString(i)+": Mismatch: "
            			+Integer.toString(max_mismatch));
            	if (max_mismatch!=0 ) {
            		num_mismatch++;
            	}
            }
        	
        	this.compare_haps.clear();	
        }
        
       
        

        System.out.println("****The total number of disagree haplotypes is:\t"+ 
        		Integer.toString(num_mismatch)) ;
        
        
       
		return;
		
	}
	
	public void LassoEvaluate(String gold_file, String aem_file) throws IOException {
		this.gold_haps.clear();
        this.compare_haps.clear();
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
	
	public void GcAemEvaluate(String gold_file, String aem_file) throws IOException {
		BufferedReader bufferedreader = new BufferedReader(new FileReader(gold_file));
        String line = "";
        this.gold_haps.clear();
        this.compare_haps.clear();
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
//        	this.compare_haps.add(tmp_str);
        	if ((j%1) ==0) {
        		this.compare_haps.add(tmp_str);
        	}
        }
        bufferedreader2.close();
        
		String file_path= "/home/chencao/Desktop/test2.txt";
		FileWriter mydata = new FileWriter(file_path,true);
        PrintWriter pw = new PrintWriter(mydata);
        
    	for (int i = 0; i  < this.compare_haps.size(); i++) {    		
    		pw.write(compare_haps.get(i)+"\n");
        
        }
        pw.flush();
        pw.close();
        
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
	
	
	public void StringEvaluate(String gold_file, String str_file) throws IOException {
		this.gold_haps.clear();
        this.compare_haps.clear();
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
        
        
		BufferedReader bufferedreader2 = new BufferedReader(new FileReader(str_file));
        line = "";
        while ((line = bufferedreader2.readLine()) != null) {
        	line =line.replace("\n", "").replace("\r", "");
//        	System.out.println(line);
        	if (!line.startsWith("@")){
        		this.compare_haps.add(line );
        	}
        }
        
        bufferedreader2.close();
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

	
	public Evaluate() throws IOException {
		
	}
}
