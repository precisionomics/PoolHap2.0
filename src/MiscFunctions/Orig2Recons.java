package MiscFunctions;

import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;

public class Orig2Recons {

	// args[0] = simhaps.mutations.txt
	// args[1] = p.all.pos
	// args[2] = p.all.results
	// args[3] = num_pts
	
	public static void main(String[] args) throws IOException {
		
		System.out.println("Results file:\t" + args[2] + "\n");
		
		System.out.println("Profiling all of the original haplotypes...\n");
		BufferedReader br1 = new BufferedReader(new FileReader(args[0]));
		HashMap<Integer,Boolean> origPosMap = new HashMap<Integer,Boolean>();
		HashMap<Integer,HashMap<Integer,Integer>> origVarMap = new HashMap<Integer,HashMap<Integer,Integer>>();
		String currLine = br1.readLine(); 
		while (currLine != null) {
			int hapID = Integer.parseInt(currLine.split("[\t_]")[1]); 
			int varPos = Integer.parseInt(currLine.split("[\t_]")[2]);
			origPosMap.put(varPos, false); 
			if (!origVarMap.containsKey(hapID)) origVarMap.put(hapID, new HashMap<Integer,Integer>()); 
			origVarMap.get(hapID).put(varPos, 1);
			currLine = br1.readLine(); 
		}
		br1.close();
		
		System.out.println("Comparing all of the called variant positions...");
		BufferedReader br2 = new BufferedReader(new FileReader(args[1]));
		String[] calledPos = br2.readLine().split("\t"); 
		ArrayList<Integer> verifiedPos = new ArrayList<Integer>();
		ArrayList<Integer> vcfPos = new ArrayList<Integer>();
		br2.close(); 
		int falseVars = 0; 
		for (String p : calledPos) {
			if (origPosMap.containsKey(Integer.parseInt(p))) {
				verifiedPos.add(Integer.parseInt(p));
				vcfPos.add(Integer.parseInt(p)); 
				origPosMap.put(Integer.parseInt(p), true); 
			} else {
				falseVars++;
				vcfPos.add(Integer.parseInt(p)); 
			}
		}
		int missedVars = 0; 
		ArrayList<Integer> missedPos = new ArrayList<Integer>();
		for (int p : origPosMap.keySet()) {
			if (!origPosMap.get(p)) {
				missedVars++; 
				missedPos.add(p); 
			}
		}
		System.out.println("Number of correctly called variants\t" + verifiedPos.size());
		System.out.println("Number of incorrectly called variants\t" + falseVars);
		System.out.println("Number of missed original variants\t" + missedVars); 

		System.out.println("\nComparing the original and reconstructed haplotypes...\n");
		ArrayList<Integer> origPosSorted = new ArrayList<>(origPosMap.keySet());
		Collections.sort(origPosSorted);
		for (int p : missedPos) {
			for (int h : origVarMap.keySet()){
				if (origVarMap.get(h).containsKey(p)) origVarMap.get(h).remove(p); 
			}
		}
		ArrayList<Boolean> includeOrNot = new ArrayList<Boolean>();  
		for (int p : vcfPos) {
			if (origPosMap.containsKey(p)) includeOrNot.add(origPosMap.get(p)); 
			else includeOrNot.add(false); 
			// System.out.println(p + "\t" + includeOrNot.get(vcfPos.indexOf(p)));
		}
		// System.out.println(vcfPos.size() + "\t" + includeOrNot.size()); 
		for (int p : verifiedPos) {
			for (int h : origVarMap.keySet()) origVarMap.get(h).putIfAbsent(p, 0);
			System.out.print(p + "\t");
		}
		System.out.println("\nThese are the original haplotypes...");
		for (int h : origVarMap.keySet()) {
			System.out.print("orig_" + h + "\t");
			for (int p : verifiedPos) {
				System.out.print(origVarMap.get(h).get(p) + "\t"); 
			}
			System.out.println();
		}
		BufferedReader br3 = new BufferedReader(new FileReader(args[2]));
		for (int l = 0; l < 20; l++) currLine = br3.readLine(); 
		int numHaps = Integer.parseInt(currLine.split("\t")[1]);
		int numReconVars = Integer.parseInt(currLine.split("\t")[3]);
		HashMap<Integer,HashMap<Integer,Integer>> reconVarMap = new HashMap<Integer,HashMap<Integer,Integer>>();
		HashMap<Integer,String> reconNames = new HashMap<Integer,String>(); 
		currLine = br3.readLine();
		for (int h = 0; h < numHaps; h++) {
			currLine = br3.readLine();
			reconNames.put(h,currLine.split("\t")[0]);
			String[] currHapVars= currLine.split("\t")[1].split("");
			// System.out.println(numReconVars + "\t" + currHapVars.length + "\t" + currLine.split("\t")[1]);
			reconVarMap.put(h, new HashMap<Integer,Integer>()); 
			int tmp = 0; 
			for (int v = 0; v < numReconVars; v++)  {
				// System.out.println(v + "\t" + includeOrNot.get(v));
				if (includeOrNot.get(v)) {
					reconVarMap.get(h).put(verifiedPos.get(tmp),Integer.parseInt(currHapVars[v])); 
					tmp++; 
				}
				// System.out.println(tmp + "\t" + includeOrNot.size() + "\tdone");
			}
		}
		/*
		currLine = br3.readLine();
		currLine = br3.readLine();
		int numPts = Integer.parseInt(args[3]); 
		double[][] intrapoolFreqs = new double[numHaps][numPts];
		for (int h = 0; h < numHaps; h++) {
			currLine = br3.readLine();
			String[] currHapFreqs= currLine.split("\t"); 
			for (int p = 0; p < numPts; p++) intrapoolFreqs[h][p] = Double.parseDouble(currHapFreqs[p + 2]); 
		}
		*/
		br3.close();
		System.out.println("\nThese are the reconstructed haplotypes...");
		for (int h : reconVarMap.keySet()) {
			System.out.print("recon_" + h + "\t");
			for (int p : verifiedPos) {
				System.out.print(reconVarMap.get(h).get(p) + "\t"); 
			}
			System.out.println();
		}
		// 		HashMap<Integer,HashMap<Integer,Integer>> origVarMap = new HashMap<Integer,HashMap<Integer,Integer>>();
		// 		HashMap<Integer,HashMap<Integer,Integer>> reconVarMap = new HashMap<Integer,HashMap<Integer,Integer>>();
		// 		HashMap<Integer,String> reconNames = new HashMap<Integer,String>(); 
		HashMap<Integer,String> minDiffHap = new HashMap<Integer,String>();
		HashMap<Integer,Integer> minDiffNum = new HashMap<Integer,Integer>();
		HashMap<Integer, ArrayList<Integer>> missedPosLists = new HashMap<Integer, ArrayList<Integer>>(); 
		for (int orig : origVarMap.keySet()) {
			int minHap = 0;
			int minNum = verifiedPos.size(); 
			for (int recon : reconVarMap.keySet()) {
				int diff = 0; 
				missedPosLists.put(orig, new ArrayList<Integer>()); 
				for (int pos : verifiedPos) {
					if (origVarMap.get(orig).get(pos) != reconVarMap.get(recon).get(pos)) {
						diff++;
						missedPosLists.get(orig).add(0);
					} else {
						missedPosLists.get(orig).add(1);						
					}
				}
				if (diff < minNum) {
					minHap = recon; 
					minNum = diff; 
				}
			}
			minDiffHap.put(orig,reconNames.get(minHap));
			minDiffNum.put(orig,minNum); 
		}
		System.out.println("\nOriginal_Hap\tClosest_Recons\tNum_Diff");
		for (int orig : minDiffHap.keySet()) {
			System.out.print(orig + "\t" + minDiffHap.get(orig) + "\t" + minDiffNum.get(orig) + "\t"); 
			for (int p : missedPosLists.get(orig)) System.out.print(p + "\t");
			System.out.println();
		}
	}
}