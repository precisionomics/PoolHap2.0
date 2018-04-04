package MiscFunctions;

import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.Set;

public class Orig2Recons {

	// args[0] = simhaps.mutations.txt
	// args[1] = p.all.pos
	// args[2] = num_pts
	// args[3] = folder
	// args[4] = maximum threshold of differing positions above the minimum for OH == RH 
		// (i.e.: considered part of the quasispecies 'cloud') 
	// args[5] = maximum threshold of differing frequencies for OH == RH
	// args[6] = output.file 
	// args[7] = p.all.results

	public static void main(String[] args) throws IOException {

		PrintWriter pw = new PrintWriter(args[3] + "/full_analysis.txt");
		pw.append("Folder name: " + args[3] + "\n\n"); 
		
		BufferedReader br1 = new BufferedReader(new FileReader(args[0]));
		HashMap<Integer,Boolean> origPosMap = new HashMap<Integer,Boolean>();
		HashMap<Integer,HashMap<Integer,Integer>> origVarMap = new HashMap<Integer,HashMap<Integer,Integer>>();
		String currLine = br1.readLine(); 
		int prevID = 0; 
		int numHapsSkipped = 0; 
		while (currLine != null) {
			int hapID = Integer.parseInt(currLine.split("[\t_]")[1]); 
			int varPos = Integer.parseInt(currLine.split("[\t_]")[2]);
			origPosMap.put(varPos, false); 
			if (!origVarMap.containsKey(hapID)) {
				if (hapID > prevID + 1) {
					numHapsSkipped = hapID - prevID - 1; 
					for (int h = 1; h <= numHapsSkipped; h++) {
						origVarMap.put(prevID + h, new HashMap<Integer,Integer>());	
						// These will be empty i.e.: completely reference-base haplotypes.  
					}
				}
				origVarMap.put(hapID, new HashMap<Integer,Integer>()); 
			}
			origVarMap.get(hapID).put(varPos, 1);
			currLine = br1.readLine(); 
			prevID = hapID; 
		}
		br1.close();
		
		pw.append("Comparing all of the called variant positions...\n");
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
		pw.append("Number of correctly called variants\t" + verifiedPos.size() + "\n");
		pw.append("Number of incorrectly called variants\t" + falseVars + "\n");
		pw.append("Number of missed original variants\t" + missedVars + "\n"); 
		System.out.println(origVarMap.keySet().size()); // Number of OH.
		double mutRate = Double.parseDouble(args[6].split("_")[1]);
		System.out.println(mutRate + "\n" + mutRate); // The mutation and recombination rates (identical).
		System.out.println(verifiedPos.size()); // Number of verified variant positions.
		double propAccCall = (double) verifiedPos.size() / (verifiedPos.size() + missedVars); 
		System.out.printf("%.3f\n",propAccCall);	// Fraction of actual variants that are also called (i.e.: verified).
		double propFalseCall = (double) falseVars / (verifiedPos.size() + falseVars);
		System.out.printf("%.3f\n", propFalseCall);	// 	Fraction of called variant positions that are also inaccurate (i.e.: not in the actual set). 

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
		}
		// Purpose: To print the intersection of actual and called (i.e.: verified) variant positions.
		// System.out.print("\nPositions: \t");
		for (int p : verifiedPos) {
			for (int h : origVarMap.keySet()) origVarMap.get(h).putIfAbsent(p, 0);
			// System.out.print(p + "\t");
		}
		int numPts = Integer.parseInt(args[2]); 
		double[][] origCts = new double[origVarMap.keySet().size()][numPts];
		double[][] origFreqs = new double[origVarMap.keySet().size()][numPts];
		int interSum = 0; 
		for (int p = 0; p < numPts; p++) {
			int ptSum = 0;
			BufferedReader brp = new BufferedReader(new FileReader(args[3] + "/p" + p + ".fa"));
			currLine = brp.readLine(); 
			while (currLine != null) {
				if (currLine.contains(">")) {
					int hapNum = Integer.parseInt(currLine.split("_")[1].replace(" ",""));
					origCts[hapNum][p]++; 
					ptSum++;
				}
				currLine = brp.readLine(); 
			}
			brp.close();
			for (int h = 0; h < origVarMap.keySet().size(); h++) origFreqs[h][p] = origCts[h][p] / ptSum;
			interSum += ptSum; 
		}
		
		double[] origInterFreqs = new double[origVarMap.keySet().size()]; 
		for (int h = 0; h < origVarMap.keySet().size(); h++) {
			int hapCt = 0; 
			for (int p = 0; p < numPts; p++) {
				hapCt += origCts[h][p]; 
			}
			origInterFreqs[h] = (double) hapCt/interSum; 
		}
		int maxPosDiff = Integer.parseInt(args[4]);

		pw.append("\nOutput file:\t" + args[7] + "\n");
		pw.append("\nComparing the original and reconstructed haplotypes...\n");
		BufferedReader brOut = new BufferedReader(new FileReader(args[6]));
		String gcResults = "initial haplotypes spanning";
		currLine = brOut.readLine(); 
		while (currLine.indexOf(gcResults) == -1) currLine = brOut.readLine();
		int gcHaps = Integer.parseInt(currLine.split(" ")[2]);
		System.out.println(gcHaps); // Number of GC-reconstructed haplotypes. 
		int numGCVars = Integer.parseInt(currLine.split(" ")[6]);
		double[] gcInterFreqs = new double[gcHaps];
		HashMap<Integer,HashMap<Integer,Integer>> gcVarMap = new HashMap<Integer,HashMap<Integer,Integer>>();
		HashMap<Integer,String> gcNames = new HashMap<Integer,String>(); 
		for (int h = 0; h < gcHaps; h++) {
			currLine = brOut.readLine();
			gcNames.put(h,currLine.split("\t")[0]);
			String[] currHapVars= currLine.split("\t")[1].split("");
			gcVarMap.put(h, new HashMap<Integer,Integer>()); 
			int tmp = 0; 
			for (int v = 0; v < numGCVars; v++)  {
				if (includeOrNot.get(v)) {
					gcVarMap.get(h).put(verifiedPos.get(tmp),Integer.parseInt(currHapVars[v])); 
					tmp++; 
				}
			}
			gcInterFreqs[h] = Double.parseDouble(currLine.split("\t")[2]); 
		}
		brOut.close();
		OutputReporter(pw, gcVarMap, origVarMap, verifiedPos, gcNames, maxPosDiff, gcHaps, gcInterFreqs, origInterFreqs); 
			// GC results.

		pw.append("\nResults file:\t" + args[7] + "\n");
		pw.append("\nComparing the original and reconstructed haplotypes...\n");
		BufferedReader brRes = new BufferedReader(new FileReader(args[7]));
		for (int l = 0; l < 20; l++) currLine = brRes.readLine(); 
		int numHaps = Integer.parseInt(currLine.split("\t")[1]);
		System.out.println(numHaps); // Number of PHX-reconstructed haplotypes. 
		int numReconVars = Integer.parseInt(currLine.split("\t")[3]);
		double[] reconInterFreqs = new double[numHaps];
		HashMap<Integer,HashMap<Integer,Integer>> reconVarMap = new HashMap<Integer,HashMap<Integer,Integer>>();
		HashMap<Integer,String> reconNames = new HashMap<Integer,String>(); 
		currLine = brRes.readLine();
		for (int h = 0; h < numHaps; h++) {
			currLine = brRes.readLine();
			reconNames.put(h,currLine.split("\t")[0]);
			String[] currHapVars= currLine.split("\t")[1].split("");
			reconVarMap.put(h, new HashMap<Integer,Integer>()); 
			int tmp = 0; 
			for (int v = 0; v < numReconVars; v++)  {
				if (includeOrNot.get(v)) {
					reconVarMap.get(h).put(verifiedPos.get(tmp),Integer.parseInt(currHapVars[v])); 
					tmp++; 
				}
			}
			reconInterFreqs[h] = Double.parseDouble(currLine.split("\t")[2]); 
		}
		currLine = brRes.readLine();
		currLine = brRes.readLine();
		double[][] intrapoolFreqs = new double[numHaps][numPts];			
		for (int h = 0; h < numHaps; h++) {
			currLine = brRes.readLine();
			String[] currHapFreqs= currLine.split("\t"); 
			for (int p = 0; p < numPts; p++) intrapoolFreqs[h][p] = Double.parseDouble(currHapFreqs[p + 2]); 
		}
		brRes.close();
		double maxFreqDiff = Double.parseDouble(args[5]);
		for (int f = 0; f <= maxPosDiff; f++) {	// GC + rjMCMC + EM results.
			// System.out.println(f + " position differences...");
			ResultsReporter(pw, reconVarMap, origVarMap, verifiedPos, reconNames, f, numHaps, reconInterFreqs, 
					origInterFreqs, numPts, intrapoolFreqs, origFreqs, maxFreqDiff);
		}
		pw.close();
	}
	
	public static void OutputReporter(PrintWriter pw, HashMap<Integer,HashMap<Integer,Integer>> reconVarMap, 
			HashMap<Integer,HashMap<Integer,Integer>> origVarMap, ArrayList<Integer> verifiedPos, 
			HashMap<Integer,String> reconNames, int maxPosDiff, int numHaps, double[] reconInterFreqs, 
			double[] origInterFreqs) throws IOException {
				
		/* Purpose: To print the variant composition of the original haplotypes at (only) the verified positions.
		System.out.println("\nThese are the original haplotypes...");
		for (int h : origVarMap.keySet()) {
			System.out.print("orig_" + h + "\t");
			for (int p : verifiedPos) {
				System.out.print(origVarMap.get(h).get(p) + "\t"); 
			}
			System.out.println();
		}
		*/
		/* Purpose: To print the variant composition of the reconstructed haplotypes at (only) the verified positions.
		System.out.println("\nThese are the reconstructed haplotypes...");
		for (int h : reconVarMap.keySet()) {
			System.out.print("recon_" + h + "\t");
			for (int p : verifiedPos) {
				System.out.print(reconVarMap.get(h).get(p) + "\t"); 
			}
			System.out.println();
		}
		*/
		HashMap<Integer,String> minDiffHap = new HashMap<Integer,String>();
		HashMap<Integer, Integer> minDiffHapInd = new HashMap<Integer,Integer>();
		HashMap<Integer,Integer> minDiffNum = new HashMap<Integer,Integer>();
		Set<Integer> tmpRecons = reconVarMap.keySet(); 
		HashMap<Integer,ArrayList<Integer>> allDiffNum = new HashMap<Integer,ArrayList<Integer>>(); 
		HashMap<Integer, ArrayList<Integer>> maxDiffOrLess = new HashMap<Integer, ArrayList<Integer>>(); 
		HashMap<Integer, ArrayList<Integer>> diffAnalysis = new HashMap<Integer, ArrayList<Integer>>(); 
		for (int orig : origVarMap.keySet()) {
			int minHap = 0;
			int minNum = verifiedPos.size(); 
			allDiffNum.put(orig, new ArrayList<Integer>()); 
			for (int recon : tmpRecons) {
				int diff = 0; 
				for (int pos : verifiedPos) {
					if (origVarMap.get(orig).get(pos) != reconVarMap.get(recon).get(pos)) diff++;
				}
				allDiffNum.get(orig).add(diff); 
				if (diff < minNum) {
					minHap = recon; 
					minNum = diff; 
				}
			}
			minDiffHap.put(orig,reconNames.get(minHap));
			minDiffHapInd.put(orig,minHap);
			minDiffNum.put(orig,minNum); 
			maxDiffOrLess.put(orig, new ArrayList<Integer>()); 
			diffAnalysis.put(orig, new ArrayList<Integer>()); 
			for (int i = 0; i <= maxPosDiff; i++) diffAnalysis.get(orig).add(0); 
			for (int recon = 0; recon < allDiffNum.get(orig).size(); recon++) { 
				if (allDiffNum.get(orig).get(recon) <= minDiffNum.get(orig) + maxPosDiff) maxDiffOrLess.get(orig).add(recon);
				if (allDiffNum.get(orig).get(recon) <= (minNum + maxPosDiff)) {
					int tmp = diffAnalysis.get(orig).get(allDiffNum.get(orig).get(recon) - minNum) + 1;
					diffAnalysis.get(orig).set(allDiffNum.get(orig).get(recon) - minNum, tmp);
				}
			}
		}
		pw.append("Original_Hap\tClosest_Recons\tNum_Diff\tProp_Diff\t");
		for (int i = 0; i <= maxPosDiff; i++) {
			pw.append("+" + i + "\t");
		}
		pw.append("Prop_Quasi\n");
		double propQuasiSum = 0; 
		double propErrorSum = 0; 
		for (int orig : minDiffHap.keySet()) {
			double propDiff = (double) minDiffNum.get(orig) / verifiedPos.size();  
			pw.append(orig + "\t" + minDiffHap.get(orig) + "\t" + minDiffNum.get(orig) + "\t");
			pw.format("%.3f\t", propDiff);
			// Purpose: To print the number of differences between each reconstructed haplotype and the current original haplotype.
			// for (int d : allDiffNum.get(orig)) System.out.print(d + "\t");
			double ctQuasi =  0; 
			for (int d : diffAnalysis.get(orig)) { 
				pw.append(d + "\t");
				ctQuasi += d;
			}
			propQuasiSum += ctQuasi / numHaps;
			propErrorSum += propDiff; 
			pw.format("%.3f\n", ctQuasi / numHaps);
		}
		double propQuasi = propQuasiSum / minDiffHap.keySet().size(); 
		double propError = propErrorSum / minDiffHap.keySet().size(); 
		int tmpProp = (int) ((double)maxPosDiff / verifiedPos.size() * 100); 
		pw.append("Average proportion of minimum reconstruction error = ");
		pw.format("%.3f\n", propError);
		System.out.printf("%.3f\n", propError);	// The proportion of the variant composition of the most similar reconstructed haplotype that is incorrect.
		pw.append("Average proportion of reconstructed at " + tmpProp + "% similarity = ");
		pw.format("%.3f\n", propQuasi);
		System.out.printf("%.3f\n", propQuasi);	// The proportion of all reconstructed haplotypes that are more similar than the threshold.	

		/* Purpose: To print the original and reconstructed haplotype interpool frequencies respectively.
		System.out.println("\nThese are the interpool frequencies of the original haplotype in each patient...");
		for (int h = 0; h < origVarMap.keySet().size(); h++) {
			System.out.print("orig_" + h + "\t");
 			for (int p = 0; p < numPts; p++) System.out.prinf("%.3f\t", origFreqs[h][p]);
			System.out.println();
		}
		System.out.println("\nThese are the reconstructed haplotype frequencies in each patient...");
		for (int h = 0; h < numHaps; h++) {
			System.out.print("recon_" + h + "\t");
			System.out.printf("%.3f\t", reconInterFreqs[h]);
			System.out.println();
		}
		*/
		
		pw.append("\nThe difference between the reconstructed and original interpool frequencies and the reconstructed ones...\n");
		double sumFreqDiff = 0; 
		for (int h = 0; h < origVarMap.keySet().size(); h++) {
			pw.append("orig_" + h + "\t");
			double maxDiffOne = 0;
			for (int r : maxDiffOrLess.get(h)) maxDiffOne += reconInterFreqs[r];
			double ipDiff = maxDiffOne - origInterFreqs[h]; 
			pw.format("%.3f\n", ipDiff);
			sumFreqDiff += Math.abs(ipDiff); 
		}
		double avgFreqDiff = sumFreqDiff / origVarMap.keySet().size();
		pw.format("Average difference between reconstructed and original interpool frequencies = %.3f\n", avgFreqDiff);		
		System.out.printf("%.3f\n", avgFreqDiff); 	// The average difference between the interpool reconstructed and original frequencies.
	}
	
	public static void ResultsReporter(PrintWriter pw, HashMap<Integer,HashMap<Integer,Integer>> reconVarMap, 
			HashMap<Integer,HashMap<Integer,Integer>> origVarMap, ArrayList<Integer> verifiedPos, 
			HashMap<Integer,String> reconNames, int maxPosDiff, int numHaps, double[] reconInterFreqs, 
			double[] origInterFreqs, int numPts, double[][] intrapoolFreqs, double[][] origFreqs, double maxFreqDiff) 
					throws IOException {
				
		/* Purpose: To print the variant composition of the original haplotypes at (only) the verified positions.
		System.out.println("\nThese are the original haplotypes...");
		for (int h : origVarMap.keySet()) {
			System.out.print("orig_" + h + "\t");
			for (int p : verifiedPos) {
				System.out.print(origVarMap.get(h).get(p) + "\t"); 
			}
			System.out.println();
		}
		*/
		/* Purpose: To print the variant composition of the reconstructed haplotypes at (only) the verified positions.
		System.out.println("\nThese are the reconstructed haplotypes...");
		for (int h : reconVarMap.keySet()) {
			System.out.print("recon_" + h + "\t");
			for (int p : verifiedPos) {
				System.out.print(reconVarMap.get(h).get(p) + "\t"); 
			}
			System.out.println();
		}
		*/
		HashMap<Integer,String> minDiffHap = new HashMap<Integer,String>();
		HashMap<Integer, Integer> minDiffHapInd = new HashMap<Integer,Integer>();
		HashMap<Integer,Integer> minDiffNum = new HashMap<Integer,Integer>();
		Set<Integer> tmpRecons = reconVarMap.keySet(); 
		HashMap<Integer,ArrayList<Integer>> allDiffNum = new HashMap<Integer,ArrayList<Integer>>(); 
		HashMap<Integer, ArrayList<Integer>> maxDiffOrLess = new HashMap<Integer, ArrayList<Integer>>(); 
		HashMap<Integer, ArrayList<Integer>> diffAnalysis = new HashMap<Integer, ArrayList<Integer>>(); 
		for (int orig : origVarMap.keySet()) {
			int minHap = 0;
			int minNum = verifiedPos.size(); 
			allDiffNum.put(orig, new ArrayList<Integer>()); 
			for (int recon : tmpRecons) {
				int diff = 0; 
				for (int pos : verifiedPos) {
					if (origVarMap.get(orig).get(pos) != reconVarMap.get(recon).get(pos)) diff++;
				}
				allDiffNum.get(orig).add(diff); 
				if (diff < minNum) {
					minHap = recon; 
					minNum = diff; 
				}
			}
			minDiffHap.put(orig,reconNames.get(minHap));
			minDiffHapInd.put(orig,minHap);
			minDiffNum.put(orig,minNum); 
			maxDiffOrLess.put(orig, new ArrayList<Integer>()); 
			diffAnalysis.put(orig, new ArrayList<Integer>()); 
			for (int i = 0; i <= maxPosDiff; i++) diffAnalysis.get(orig).add(0); 
			for (int recon = 0; recon < allDiffNum.get(orig).size(); recon++) { 
				if (allDiffNum.get(orig).get(recon) <= minDiffNum.get(orig) + maxPosDiff) maxDiffOrLess.get(orig).add(recon);
				if (allDiffNum.get(orig).get(recon) <= (minNum + maxPosDiff)) {
					int tmp = diffAnalysis.get(orig).get(allDiffNum.get(orig).get(recon) - minNum) + 1;
					diffAnalysis.get(orig).set(allDiffNum.get(orig).get(recon) - minNum, tmp);
				}
			}
		}
		pw.append("Original_Hap\tClosest_Recons\tNum_Diff\tProp_Diff\t");
		for (int i = 0; i <= maxPosDiff; i++) {
			pw.append("+" + i + "\t");
		}
		pw.append("Prop_Quasi\n");
		double propQuasiSum = 0; 
		double propErrorSum = 0; 
		for (int orig : minDiffHap.keySet()) {
			double propDiff = (double) minDiffNum.get(orig) / verifiedPos.size();  
			pw.append(orig + "\t" + minDiffHap.get(orig) + "\t" + minDiffNum.get(orig) + "\t");
			pw.format("%.3f\t", propDiff);
			// Purpose: To print the number of differences between each reconstructed haplotype and the current original haplotype.
			// for (int d : allDiffNum.get(orig)) System.out.print(d + "\t");
			double ctQuasi =  0; 
			for (int d : diffAnalysis.get(orig)) { 
				pw.append(d + "\t");
				ctQuasi += d;
			}
			propQuasiSum += ctQuasi / numHaps;
			propErrorSum += propDiff; 
			pw.format("%.3f\n", ctQuasi / numHaps);
		}
		double propQuasi = propQuasiSum / minDiffHap.keySet().size(); 
		double propError = propErrorSum / minDiffHap.keySet().size(); 
		int tmpProp = (int) ((double)maxPosDiff / verifiedPos.size() * 100); 
		pw.append("Average proportion of minimum reconstruction error = ");
		pw.format("%.3f\n", propError);
		System.out.printf("%.3f\n", propError);	// The proportion of the variant composition of the most similar reconstructed haplotype that is incorrect.
		pw.append("Average proportion of reconstructed at " + tmpProp + "% similarity = ");
		pw.format("%.3f\n", propQuasi);
		System.out.printf("%.3f\n", propQuasi);	// The proportion of all reconstructed haplotypes that are more similar than the threshold.	

		/* Purpose: To print the original and reconstructed haplotype interpool frequencies respectively.
		System.out.println("\nThese are the interpool frequencies of the original haplotype in each patient...");
		for (int h = 0; h < origVarMap.keySet().size(); h++) {
			System.out.print("orig_" + h + "\t");
 			for (int p = 0; p < numPts; p++) System.out.prinf("%.3f\t", origFreqs[h][p]);
			System.out.println();
		}
		System.out.println("\nThese are the reconstructed haplotype frequencies in each patient...");
		for (int h = 0; h < numHaps; h++) {
			System.out.print("recon_" + h + "\t");
			System.out.printf("%.3f\t", reconInterFreqs[h]);
			System.out.println();
		}
		*/
		
		pw.append("\nThe difference between the reconstructed and original interpool frequencies and the reconstructed ones...\n");
		double sumFreqDiff = 0; 
		for (int h = 0; h < origVarMap.keySet().size(); h++) {
			pw.append("orig_" + h + "\t");
			double maxDiffOne = 0;
			for (int r : maxDiffOrLess.get(h)) maxDiffOne += reconInterFreqs[r];
			double ipDiff = maxDiffOne - origInterFreqs[h]; 
			pw.format("%.3f\n", ipDiff);
			sumFreqDiff += Math.abs(ipDiff); 
		}
		double avgFreqDiff = sumFreqDiff / origVarMap.keySet().size();
		pw.format("Average difference between reconstructed and original interpool frequencies = %.3f\n", avgFreqDiff);		
		System.out.printf("%.3f\n", avgFreqDiff); 	// The average difference between the interpool reconstructed and original frequencies.
		
		// Purpose: To print the difference between the intrapool frequencies of the original and their closest reconstructed haploypes.
		pw.append("\nThe difference between the reconstructed and original interpool frequencies and the reconstructed ones...\n");
		int numOver = 0;
		int numUnder = 0;
		int numNP = 0; 
		int numAcc = 0; 
		for (int h = 0; h < origVarMap.keySet().size(); h++) {
			pw.append("orig_" + h + "\t");
			int minHap = minDiffHapInd.get(h); 
			for (int p = 0; p < numPts; p++) {
				double ipDiff = intrapoolFreqs[minHap][p] - origFreqs[h][p]; 
				pw.format("%.3f\t", ipDiff);
			}
			for (int p = 0; p < numPts; p++) {
				double maxDiffOne = 0;
				for (int r : maxDiffOrLess.get(h)) maxDiffOne += intrapoolFreqs[r][p];
				double ipDiff = maxDiffOne - origFreqs[h][p];
				if (origFreqs[h][p] > 0) {
					if (ipDiff > (origFreqs[h][p] * maxFreqDiff)) numOver++;
					else if (ipDiff < -(origFreqs[h][p] * maxFreqDiff)) numUnder++;
					else numAcc++; 
				} else {
					if (ipDiff < maxFreqDiff) numAcc++;
					else numNP++;
				}
			}
			pw.append("\n");
		}

		double freqOver = (double) numOver / (origVarMap.keySet().size() * numPts);
		double freqUnder = (double) numUnder / (origVarMap.keySet().size() * numPts);
		double freqNP = (double) numNP / (origVarMap.keySet().size() * numPts); 
		double freqAcc = (double) numAcc / (origVarMap.keySet().size() * numPts); 
		pw.append("\nProportion of reconstructed intrapool frequencies that were (in)accurate...\n");
		pw.append("Over\tUnder\tNotPres\tAccurate\n");
		pw.format("%.3f\t%.3f\t%.3f\t%.3f\n\n", freqOver, freqUnder, freqNP, freqAcc);
		System.out.printf("%.3f\n%.3f\n%.3f\n%.3f\n", freqOver, freqUnder, freqNP, freqAcc);
			// The proportions of intrapool frequency estimates that are (in)accurate. 
	}
}