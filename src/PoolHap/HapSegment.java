package PoolHap;

import java.util.HashMap;

/*
 * Created by Quan Long Jan 2019.
 * 
 * This class records one segment of a haplotype. 
 * 
 * This is used in the DivideConquer.java that form haplotype configurations based on GC outcome. 
 * It puts indexes and sequences in the same object to simplify the the structure of the method
 * 
 * public HapConfig generate_hapconfig_from_gc(int[] the_region, String[][] gc_outcome2, SiteInPoolFreqAnno in_pool_sites2)
 * 
 */
public class HapSegment {
	
	public int start_index;
	public int end_index;
	public String sequence;
	
	public HapSegment(int start_index, int end_index, String sequence){
		if(start_index>end_index){
			System.out.println("ERROR HapSegment: start_index("+start_index+")>end_index("+end_index+")");
		}
		this.sequence=sequence;
		this.start_index=start_index;
		this.end_index=end_index;
	}
	
	/*
	 * Combine two segments into one. 
	 * It will fail if these two are not adjacent.  
	 */
	
	public static HapSegment combine(HapSegment hs1, HapSegment hs2){
		if(hs1.end_index+1==hs2.start_index){ //hs1 comes first, the preferred pattern
			return new HapSegment(hs1.start_index, hs2.end_index, hs1.sequence+hs2.sequence);
		}else if(hs2.end_index+1==hs1.start_index){ // hs2 comes first, which is not encouraged, but OK. 
			return new HapSegment(hs2.start_index, hs1.end_index, hs2.sequence+hs1.sequence);
		}else{
			System.out.println("ERROR: two segment haplotypes to be combined are not adjacent!");
			return null;
		}
	}
	
	/*
	 * Override equals
	 * @see java.lang.Object#equals(java.lang.Object)
	 */
	public boolean equals(Object seg2){
		if(seg2==this)return true;
		if (!(seg2 instanceof HapSegment)) { 
            return false; 
        } 
		else {
			HapSegment casted=(HapSegment) seg2;
			return (this.start_index==casted.start_index 
				&& this.end_index==casted.end_index 
				&& this.sequence.equals(casted.sequence));
		}
	}
	
	/*
	 * override hashCode by using the hashCode of the unique string of this object. 
	 * @see java.lang.Object#hashCode()
	 */
	public int hashCode(){
		String unique_string_of_this_object=this.start_index+"_"+this.end_index+"_"+this.sequence;
		return unique_string_of_this_object.hashCode();
	}
	
	/*
	 * Output the String representation that is easier to interpret, instead of [PoolHap.HapSegment@Integer.toHexString(hashCode())]
	 * @see java.lang.Object#toString()
	 */
	public String toString(){
		return this.start_index+"_"+this.end_index+"_"+this.sequence;
	}
	
	public static void add2hashmap(HashMap<HapSegment, Integer> hapSeg_set, HapSegment seg) {
		if(hapSeg_set.containsKey(seg)) {
			int the_increased_num=hapSeg_set.get(seg)+1;
			hapSeg_set.put(seg, the_increased_num);
		}else {
			hapSeg_set.put(seg, 1);
		}
	}
//	public static void main(String[] args){
//		HashSet<HapSegment> set=new HashSet<HapSegment>();
//		HapSegment s1=new HapSegment(1,2,"0101");
//		HapSegment s2=new HapSegment(1,3,"0101");
//		HapSegment s3=new HapSegment(1,2,"0111");
//		HapSegment s4=new HapSegment(1,2,"0101");
//		
//		set.add(s1);
//		System.out.println(set);
//		System.out.println(Integer.toHexString("1_2_0101".hashCode()));
//		System.out.println(Integer.toHexString(s1.hashCode()));
//		set.add(s2);
//		System.out.println(set);
//		set.add(s3);
//		System.out.println(set);
////		set.add(s4);
////		System.out.println(set);
//		System.out.println(set.contains(s4));
//	}
	
	
}
