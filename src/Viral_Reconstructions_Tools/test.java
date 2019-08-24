package Viral_Reconstructions_Tools;

import java.util.ArrayList;
import java.util.HashMap;

import org.netlib.util.doubleW;

import shapeless.newtype;
import spire.optional.intervalGeometricPartialOrder;

public class test {

	public static void main(String[] args) {
		HashMap<String, double[]> hap2poolfre = new HashMap<String, double[]>();
		
//		String key_0 = "0_key";
//		String key_1 = "1_key";
//		String key_2 = "1_key";
//		String key_3 = "2_key";
		hap2poolfre.put("key_0", new double[4]);
		hap2poolfre.put("key_1", new double[4]);
		hap2poolfre.put("key_2", new double[4]);
		hap2poolfre.get("key_0")[0] = 0.1;
		hap2poolfre.get("key_0")[3] = 0.2;
		System.out.println(hap2poolfre.get("key_0"));
		for (int i=0; i< 4; i++) {
			System.out.println(hap2poolfre.get("key_0")[i]);
		}
		
//		System.out.println(hap2poolfre.size());
//		System.out.println(hap2poolfre.containsKey("0_key"));
//		System.out.println(hap2poolfre.containsKey("3_key"));
//		System.out.println(hap2poolfre.get(key_0).contains(0.1));
		
		ArrayList<String> hap_string_list=new ArrayList<String>();
		
		System.out.println(hap2poolfre.keySet());
		for ( String key : hap2poolfre.keySet() ) {
		    System.out.println( key );
		    hap_string_list.add(key);
		}
		
		System.out.print(hap_string_list);
		
		
		
		
		
		
//
//		
//		
//		ArrayList<ArrayList<String>> ref_seq_listlist=new ArrayList<ArrayList<String>>();
//		ArrayList<String> hap_string_list=new ArrayList<String>();
//		ArrayList<String> hap_string_list_2=new ArrayList<String>();
//		ArrayList<String> hap_string_list_3=new ArrayList<String>();
//		hap_string_list.add("1");
//		hap_string_list.add("2");
//		hap_string_list.add("3");
//		hap_string_list_3.add("1");
//		hap_string_list_3.add("2");
//		hap_string_list_3.add("3");
//		hap_string_list_2.add("4");
//		hap_string_list_2.add("5");
//		ref_seq_listlist.add(hap_string_list);
//		ref_seq_listlist.add(hap_string_list_2);
//		ref_seq_listlist.add(hap_string_list_3);
//		System.out.println(ref_seq_listlist);	
//		ref_seq_listlist.get(0).set(1,"0");
//		
//		System.out.println(ref_seq_listlist);
		
		
	}

}
