package Viral_Reconstructions_Tools;

import java.util.ArrayList;
import java.util.HashMap;

import shapeless.newtype;

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
		String[] str1Strings= {"0","0","0","0"};
		String[] str2Strings= {"0","0","0","0"};
		System.out.println(str1Strings);
		System.out.println(str2Strings);
		has1HashMap.put("1000", 3);
		pool.add(has1HashMap);
		has2HashMap.put("000", 2);
		has2HashMap.put("100", 3);
		pool.add(has2HashMap);
		pool.get(0).put("1111", 4);
		System.out.println(pool);
	}
}
