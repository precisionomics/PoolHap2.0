package MiscFunctions;

import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;

import spire.optional.intervalGeometricPartialOrder;


public class PoolSimulatorSliM {
	
	public static void processing_standard_outcome() throws IOException{
//		BufferedReader br = new BufferedReader(new FileReader(gs_dir + project_name + "_slim.txt")); 
//		String file = "D:\\PhD-Studying\\Informatics\\Project\\HIV_project\\PoolHapX_testing\\input\\panmictic_diploid.out";
		String file = "D:\\PhD-Studying\\Informatics\\Project\\HIV_project\\PoolHapX_testing\\input\\island_haploid_p9.out";
		ArrayList<String> index2varpos = new ArrayList<String>();
		ArrayList<HashMap<String, Integer>> pool2allhapList=new ArrayList<HashMap<String, Integer>>();
		BufferedReader br = new BufferedReader(new FileReader(file)); 
		String currLine = br.readLine(); // header
		String[] tmpcurrpos = currLine.split(" ");
		Boolean is_single_population = false;
		if(is_single_population==true) {
			HashMap<String, Integer> hapforpool= new HashMap<String, Integer>();
			pool2allhapList.add(hapforpool);
			while(!tmpcurrpos[0].equals("Mutations:")) { // Read those lines before "Mutations:"
				currLine = br.readLine();
				tmpcurrpos = currLine.split(" ");
			}
			currLine = br.readLine();
			tmpcurrpos = currLine.split(" ");
		}else {
			while(!tmpcurrpos[0].equals("Populations:")) { // Read those lines before "Populations:"
				currLine = br.readLine();
				tmpcurrpos = currLine.split(" ");
			}
			currLine = br.readLine();
			while(!tmpcurrpos[0].equals("Mutations:")) { // Read those lines before "Mutations:"
				HashMap<String, Integer> hapforpool= new HashMap<String, Integer>();
				pool2allhapList.add(hapforpool);
				currLine = br.readLine();
				tmpcurrpos = currLine.split(" ");
			}
			currLine = br.readLine();
			tmpcurrpos = currLine.split(" ");
		}
		System.out.println(pool2allhapList);
		System.out.println(currLine);
		System.out.println(tmpcurrpos[0]);
		while(!tmpcurrpos[0].equals("Genomes:")&&!tmpcurrpos[0].equals("Individuals:")) {
			 index2varpos.add(tmpcurrpos[0]+"_"+tmpcurrpos[3]);
			 currLine = br.readLine();
			 tmpcurrpos = currLine.split(" ");
		}
		System.out.println(index2varpos);
		int num_var_pos=index2varpos.size();
		int[] sim_var_pos=new int[num_var_pos];
		for(int p=0;p<index2varpos.size();p++) {
			String[] index2varposarray = index2varpos.get(p).split("_");
			sim_var_pos[Integer.parseInt(index2varposarray[0])]=Integer.parseInt(index2varposarray[1]);
		}
		for(int i=0;i<sim_var_pos.length;i++) {
			//System.out.println(sim_var_pos[i]);
		}
		if(!tmpcurrpos[0].equals("Genomes:")) {
			while(!tmpcurrpos[0].equals("Genomes:")) {
				 currLine = br.readLine();
				 tmpcurrpos = currLine.split(" ");
			}
		}
		//System.out.println(currLine);
		currLine = br.readLine();
		HashMap<String, Integer> hapsHS = new HashMap<String, Integer>();
		while(currLine!=null) {
			tmpcurrpos = currLine.split(" ");
			int curr_pool_index= Integer.parseInt(tmpcurrpos[0].split("")[1])-1;
			if(tmpcurrpos.length>2) {
				ArrayList<Integer> currvarpos = new ArrayList<Integer>();
				String curr_hap="";
				for (int p=2;p<tmpcurrpos.length;p++){
					currvarpos.add(Integer.parseInt(tmpcurrpos[p]));
				}
				for(int p=0;p<num_var_pos;p++) {
					if(currvarpos.contains(p)) {
						curr_hap=curr_hap+1;
					}else {
						curr_hap=curr_hap+0;
					}
				}
				if (!hapsHS.containsKey(curr_hap)) {
					hapsHS.put(curr_hap, 1);
				}else {
					int tmpCt =  hapsHS.get(curr_hap) + 1;
					hapsHS.put(curr_hap, tmpCt);
				}
				if(!pool2allhapList.get(curr_pool_index).containsKey(curr_hap)) {
					pool2allhapList.get(curr_pool_index).put(curr_hap, 1);
				}else if(pool2allhapList.get(curr_pool_index).containsKey(curr_hap)) {
				    int tmpCt = pool2allhapList.get(curr_pool_index).get(curr_hap)+1;
					pool2allhapList.get(curr_pool_index).put(curr_hap, tmpCt);
				}
			}
			currLine = br.readLine();
		}
		System.out.println(pool2allhapList);
		br.close();
	}

	public static void main(String[] args) throws IOException {
		processing_standard_outcome();

	}

}
