package MiscFunctions;

import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;


public class PoolSimulatorSliM {
	
	public static void processing_standard_outcome() throws IOException{
//		BufferedReader br = new BufferedReader(new FileReader(gs_dir + project_name + "_slim.txt")); 
		String file = "D:\\PhD-Studying\\Informatics\\Project\\HIV_project\\PoolHapX_testing\\input\\panmictic_diploid.out";
		ArrayList<String> index2varpos = new ArrayList<String>();
		ArrayList<HashMap<String, Integer>> pool2allhapList=new ArrayList<HashMap<String, Integer>>();
		BufferedReader br = new BufferedReader(new FileReader(file)); 
		String currLine = br.readLine(); // header
		String[] tmpVarPos = currLine.split(" ");
		Boolean is_single_population = true;
		if(is_single_population==true) {
			HashMap<String, Integer> hapforpool= new HashMap<String, Integer>();
			pool2allhapList.add(hapforpool);
			while(!tmpVarPos[0].equals("Mutations:")) { // Read those lines before "Mutations:"
				currLine = br.readLine();
				tmpVarPos = currLine.split(" ");
			}
			currLine = br.readLine();
			tmpVarPos = currLine.split(" ");
		}else {
			while(!tmpVarPos[0].equals("Populations:")) { // Read those lines before "Populations:"
				currLine = br.readLine();
				tmpVarPos = currLine.split(" ");
			}
			currLine = br.readLine();
			while(!tmpVarPos[0].equals("Mutations:")) { // Read those lines before "Mutations:"
				HashMap<String, Integer> hapforpool= new HashMap<String, Integer>();
				pool2allhapList.add(hapforpool);
				currLine = br.readLine();
				tmpVarPos = currLine.split(" ");
			}
			currLine = br.readLine();
			tmpVarPos = currLine.split(" ");
		}
		System.out.println(pool2allhapList);
		System.out.println(currLine);
		System.out.println(tmpVarPos[0]);
		while(!tmpVarPos[0].equals("Genomes:")&&!tmpVarPos[0].equals("Individuals:")) {
			 index2varpos.add(tmpVarPos[0]+"_"+tmpVarPos[3]);
			 currLine = br.readLine();
			 tmpVarPos = currLine.split(" ");
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
		if(!tmpVarPos[0].equals("Genomes:")) {
			while(!tmpVarPos[0].equals("Genomes:")) {
				 currLine = br.readLine();
				 tmpVarPos = currLine.split(" ");
			}
		}
		//System.out.println(currLine);
		currLine = br.readLine();
		while(currLine!=null) {
			tmpVarPos = currLine.split(" ");
			
			currLine = br.readLine();
		}
		
		br.close();
	}

	public static void main(String[] args) throws IOException {
		processing_standard_outcome();

	}

}
