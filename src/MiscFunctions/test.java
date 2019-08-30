package MiscFunctions;

import java.util.ArrayList;

import org.netlib.util.doubleW;

import PoolHap.HapConfig;

public class test {

	public static void main(String[] args) {
		
//		 ArrayList<String> ori = new ArrayList<>();
//		 ori.add("1");
//		 ori.add("3");
//		 ori.add("4");
//		 ori.remove(0);
//		 System.out.println(ori);
		 
		 
		 ArrayList<String> new_hapIDs_list = new ArrayList<>();
		 ArrayList<String> final_hapIDs_list = new ArrayList<>();
		 double[] hap_freq_array = {0.01,0.1,0.01,0.2,0.3,0.2,0.4,0.2,0.01,0.4,0.01};
		 ArrayList<Double> new_freq_list = new ArrayList<>();
		 ArrayList<Double> final_freq_list = new ArrayList<>();
		 ArrayList<Integer> identical_hap= new ArrayList<>();
		 new_hapIDs_list.add("h0");
		 new_freq_list.add(0.01);
		 new_hapIDs_list.add("h1");
		 new_freq_list.add(0.1);
		 new_hapIDs_list.add("h0");
		 new_freq_list.add(0.01);
		 new_hapIDs_list.add("h2");
		 new_freq_list.add(0.2);
		 new_hapIDs_list.add("h3");
		 new_freq_list.add(0.3);
		 new_hapIDs_list.add("h2");
		 new_freq_list.add(0.2);
		 new_hapIDs_list.add("h4");
		 new_freq_list.add(0.4);
		 new_hapIDs_list.add("h2");
		 new_freq_list.add(0.2);
		 new_hapIDs_list.add("h0");
		 new_freq_list.add(0.01);
		 new_hapIDs_list.add("h4");
		 new_freq_list.add(0.4);
		 new_hapIDs_list.add("h0");
		 new_freq_list.add(0.01);
		 System.out.println(new_hapIDs_list);
		 System.out.println(new_hapIDs_list.size());
		 for(int h1=0;h1<new_hapIDs_list.size()-1;h1++) {
			 System.out.println("h1+"+h1);
				for(int h2=h1+1;h2<new_hapIDs_list.size();h2++) {
					System.out.println("h2+"+h2+".."+new_hapIDs_list.get(h2));
					if(new_hapIDs_list.get(h2).equals(new_hapIDs_list.get(h1))) {
						//new_hapIDs_list.remove(h2);
						identical_hap.add(h2);
						new_freq_list.set(h1,new_freq_list.get(h1)+new_freq_list.get(h2));
						hap_freq_array[h1]=hap_freq_array[h1]+hap_freq_array[h2];
						//new_freq_list.remove(h2);
					}
				}
				System.out.println(new_hapIDs_list);
				System.out.println(identical_hap);
			}
		 System.out.println(identical_hap);
		 
//		 System.out.println(new_hapIDs_list);
		 System.out.println(new_freq_list);
		 
		 for(int i=0;i<new_hapIDs_list.size();i++) {
			 if(!identical_hap.contains(i)) {
				 final_hapIDs_list.add(new_hapIDs_list.get(i));
				 //final_freq_list.add(new_freq_list.get(i));
				 final_freq_list.add(hap_freq_array[i]);
			 }
		 }
		 System.out.println(final_hapIDs_list);
		 System.out.println(final_freq_list);

	}

}
