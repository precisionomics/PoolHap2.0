package Viral_Reconstructions_Tools;

import java.util.ArrayList;
import java.util.HashMap;

import shapeless.newtype;
import spire.optional.intervalGeometricPartialOrder;

public class test {

	public static void main(String[] args) {
		String d0String= "9.0E-4";
		double d=Double.parseDouble(d0String);
		if(d>0.01) {
			System.out.println(d);
		}
		
		ArrayList<HashMap<String, Integer>> pool=new ArrayList<HashMap<String, Integer>>();
		HashMap<String, Integer> has1HashMap= new HashMap<String, Integer>();
		HashMap<String, Integer> has2HashMap= new HashMap<String, Integer>();
//		String[] str1Strings= {"0","0","0","0"};
//		String[] str2Strings= {"0","0","0","0"};
//		System.out.println(str1Strings);
//		System.out.println(str2Strings);
		has1HashMap.put("1000", 3);
		pool.add(has1HashMap);
		has2HashMap.put("000", 2);
		has2HashMap.put("100", 3);
//		int update = has2HashMap.get("000")+2;
		System.out.println(has2HashMap);
		pool.add(has2HashMap);
		pool.get(0).put("1111", 4);
		System.out.println(pool);
		System.out.println(pool.get(0).containsKey("1111"));
		System.out.println(pool.get(0).get("1000"));
//		System.out.println(pool.get(0).size());
		
		ArrayList<String> actualhaplist = new ArrayList<String>();
		
		for(int i=0;i<pool.size();i++) {
			int num_hap=0;
			for(String h:pool.get(i).keySet()) {
				num_hap=num_hap+pool.get(i).get(h);
				actualhaplist.add(h);
			}
			System.out.println(num_hap);
		}
		System.out.println(actualhaplist);
		int[] array = {15,7,10,5,2,9,22};
		int[] array2 = {1,2,3,4,5,6,7};
		for(int i = 0; i < array.length - 1;i++) {
			for(int j = 0; j < array.length - 1 - i; j++) {
				if(array[j] > array[j+1]) {
					int temp = array[j];
					array[j] = array[j+1];
					array[j+1] = temp;
					int temp2 = array2[j];
					array2[j] = array2[j+1];
					array2[j+1] = temp2;
				}
			}
		}
		for(int i=0;i<array.length;i++) {
			System.out.println(array[i]);

		}
		for(int i=0;i<array.length;i++) {
			System.out.println(array2[i]);

		}
	}
}
