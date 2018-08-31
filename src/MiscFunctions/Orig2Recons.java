package MiscFunctions;

import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.Set;

public class Orig2Recons {

	// args[0] = number of iterations
	// args[1] = folder prefix for each iteration
	// args[2] = num_pts
	// args[3] = maximum threshold of differing positions above the minimum for OH == RH 
		// (i.e.: considered part of the quasispecies 'cloud') 
	// args[4] = maximum threshold of differing frequencies for OH == RH
	
	public static void main(String[] args) throws IOException {
		ArrayList<ArrayList<Double>> resultTable = new ArrayList<ArrayList<Double>>();
		int iter = Integer.parseInt(args[0]);
		String folderRoot = args[1];
		int numPts = Integer.parseInt(args[2]);
		int maxPosDiff = Integer.parseInt(args[3]);
		double maxFreqDiff = Double.parseDouble(args[4]);
		for (int i = 0; i < iter; i++) {
			System.out.println("Processing iteration folder " + i + "...");
			String folderFull = folderRoot + "/"; // i + TODO return
			resultTable.add(AnalysisManager(folderFull, numPts, maxPosDiff, maxFreqDiff)); 
		}
		ArrayList<String> headings = new ArrayList<String>(Arrays.asList("General","General","General","General","GC","GC","GC","GC","MCMC","MCMC","MCMC","MCMC"));
		ArrayList<String> subheadings = new ArrayList<String>(Arrays.asList("orig_haps","ver_var_pos","prop_ver","prop_false","num_haps","prop_error","prop_quasi","inter_freq_diff","num_haps","prop_error","prop_quasi","inter_freq_diff"));
		ArrayList<String> subhForPHX = new ArrayList<String>(Arrays.asList("num_haps","prop_error","prop_quasi","inter_freq_diff","freq_over","freq_under","freq_np","freq_acc")); 
		for (int tmp = 0; tmp <= maxPosDiff; tmp++) {
			subheadings.addAll(subhForPHX);
			for (int tmp2 = 0; tmp2 < ((maxPosDiff + 1) * 8); tmp2++) {
				headings.add(Integer.toString(tmp));		
			}
		}
		PrintWriter pw = new PrintWriter(folderRoot + "full_analysis.txt");
		for (int row = 0; row < headings.size(); row++) {
			pw.append(headings.get(row) + "\t" + subheadings.get(row) + "\t");
			// System.out.print(headings.get(row) + "\t" + subheadings.get(row) + "\t"); 
			for (int col = 0; col < iter; col++) {
				pw.format("%.3f", resultTable.get(col).get(row));
				// System.out.printf("%.3f\t", resultTable.get(col).get(row));
				pw.append("\t");
			}
			pw.append("\n");
			// System.out.println();
		}
		pw.close(); 
	}
	
	public static ArrayList<Double> AnalysisManager(String folder, int numPts, int maxPosDiff, double maxFreqDiff) throws IOException {

		ArrayList<Double> analysisReporter = new ArrayList<Double>(); 
		// 1) Two records made from the ms-generated list of true variants for each haplotype.
		//		i) origPosMap = (variant position, called or not)
		//		ii) origVarMap = (each haplotype, (variant position, alternate allele i.e.: 1)
		BufferedReader br1 = new BufferedReader(new FileReader(folder + "simhaps.mutations.txt"));
		HashMap<Integer,Boolean> origPosMap = new HashMap<Integer,Boolean>();
		HashMap<Integer,HashMap<Integer,Integer>> origVarMap = new HashMap<Integer,HashMap<Integer,Integer>>();
		String currLine = br1.readLine(); 
		int prevID = 0; 
		int numHapsSkipped = 0; 
		int hapID = Integer.parseInt(currLine.split("[\t_]")[1]); 
		if (hapID != 0) origVarMap.put(0, new HashMap<Integer,Integer>());
		while (currLine != null) {
			hapID = Integer.parseInt(currLine.split("[\t_]")[1]); 
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
		
		// 2a) Three records made from the BAM2VEF-generated list of called variants for each haplotype.
		//		i) verifiedPos = (called true variant positions)
		//		ii) vcfPos = (all called variant positions)
		//		iii) missedPos = (missed true variant positions)
		BufferedReader br2 = new BufferedReader(new FileReader(folder + "p.all.pos"));
		String[] calledPos = br2.readLine().split("\t"); 
		ArrayList<Integer> verifiedPos = new ArrayList<Integer>();
		ArrayList<Integer> vcfPos = new ArrayList<Integer>();
		br2.close(); 
		int falseVars = 0; 
		for (String p : calledPos) {
			if (origPosMap.containsKey(Integer.parseInt(p))) {
				verifiedPos.add(Integer.parseInt(p));
				origPosMap.put(Integer.parseInt(p), true); 
			} else
				falseVars++;
			vcfPos.add(Integer.parseInt(p)); 
		}
		int missedVars = 0; 
		ArrayList<Integer> missedPos = new ArrayList<Integer>();
		for (int p : origPosMap.keySet()) {
			if (!origPosMap.get(p)) {
				missedVars++; 
				missedPos.add(p); 
			}
		}
		// 2b) Stats to do with variant calling accuracy printed to summary file and STDOUT.  
		analysisReporter.add((double) origVarMap.keySet().size()); // Number of OH.
		analysisReporter.add((double) verifiedPos.size()); // Number of verified variant positions.
		double propAccCall = (double) verifiedPos.size() / (verifiedPos.size() + missedVars); 
		analysisReporter.add(propAccCall);	// Fraction of actual variants that are also called (i.e.: verified).
		double propFalseCall = (double) falseVars / (verifiedPos.size() + falseVars);
		analysisReporter.add(propFalseCall);	// 	Fraction of called variant positions that are also inaccurate (i.e.: not in the actual set). 

		// 3) Updates to origVarMap: Remove missed true variant positions from each haplotype carrying it so that 
		//    only the called true variant positions are compared to each other. For every haplotype, at every
		//    called true variant position, add the variant position and corresponding reference allele code i.e. 0. 
		//    i) includeOrNot = (trues at indexes corresponding to called true variant positions out of all called
		//       variant positions)
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
		for (int p : verifiedPos) {
			for (int h : origVarMap.keySet()) origVarMap.get(h).putIfAbsent(p, 0);
		}
		
		// 4) Three records made from the FullSimulator-generated FastA of each patient's haplotype composition.
		//	  i and ii) origCts and origFreqs = (x = haplotypes, y = patients, cells = quantity of haplotypes)
		//    iii) origInterFreqs = (frequency of each haplotype in total)
		double[][] origCts = new double[origVarMap.keySet().size()][numPts];
		double[][] origFreqs = new double[origVarMap.keySet().size()][numPts];
		int interSum = 0; 
		for (int p = 0; p < numPts; p++) {
			int ptSum = 0;
			BufferedReader brp = new BufferedReader(new FileReader(folder + "p" + p + ".fa"));
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
		
		// 5) Get the GC-predicted haplotype variant compositions and inter-patient frequencies. The following 
		//    variables are associated with filtering for called true variant positions:
		//    tmp (index for CTVPs only), numGCVars, v (VCF position), includeOrNot (checking for truth of v). 
		//    Three records made from the output.file which contains the STDOUT of the GC part of PHX.
		//    i) gcInterFreqs = (frequency of each haplotype overall)
		//    ii) gcVarMap = (each haplotype, (called true variant position, reference 0 or alternate 1)
		//    iii) gcNames = (each haplotype, its name in String form)
		BufferedReader brOut = new BufferedReader(new FileReader(folder + "stdout.txt_1"));
		String gcResults = "initial haplotypes spanning";
		currLine = brOut.readLine(); 
		while (currLine.indexOf(gcResults) == -1) currLine = brOut.readLine();
		int gcHaps = Integer.parseInt(currLine.split(" ")[2]);
		analysisReporter.add((double) gcHaps); // Number of GC-reconstructed haplotypes. 
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
		analysisReporter.addAll(OutputReporter(gcVarMap, origVarMap, verifiedPos, gcNames, maxPosDiff, gcHaps, 
				gcInterFreqs, origInterFreqs)); // GC results.

		String rjResults = "refined haplotypes spanning";
		currLine = brOut.readLine(); 
		while (currLine.indexOf(rjResults) == -1) currLine = brOut.readLine();
		int rjHaps = Integer.parseInt(currLine.split(" ")[2]);
		analysisReporter.add((double) rjHaps); // Number of rjMCMC-reconstructed haplotypes.
		int numRJVars = Integer.parseInt(currLine.split(" ")[6]);
		double[] rjInterFreqs = new double[rjHaps];
		HashMap<Integer,HashMap<Integer,Integer>> rjVarMap = new HashMap<Integer,HashMap<Integer,Integer>>();
		HashMap<Integer,String> rjNames = new HashMap<Integer,String>(); 
		for (int h = 0; h < rjHaps; h++) {
			currLine = brOut.readLine();
			rjNames.put(h,currLine.split("\t")[0]);
			String[] currHapVars= currLine.split("\t")[1].split("");
			rjVarMap.put(h, new HashMap<Integer,Integer>()); 
			int tmp = 0; 
			for (int v = 0; v < numRJVars; v++)  {
				if (includeOrNot.get(v)) {
					rjVarMap.get(h).put(verifiedPos.get(tmp),Integer.parseInt(currHapVars[v])); 
					tmp++; 
				}
			}
			rjInterFreqs[h] = Double.parseDouble(currLine.split("\t")[2]); 
		}
		brOut.close();
		analysisReporter.addAll(OutputReporter(rjVarMap, origVarMap, verifiedPos, rjNames, maxPosDiff, rjHaps, rjInterFreqs, origInterFreqs)); // rjMCMC results.

		// 6) Get the full PHX-predicted haplotype variant compositions and intra/inter-patient frequencies. The 
		//    following variables are associated with filtering for called true variant positions:
		//    tmp (index for CTVPs only), numReconVars, v (VCF position), includeOrNot (checking for truth of v). 
		//    Two records made from the output.file which contains the STDOUT of the GC part of PHX.
		//    i) reconInterFreqs = (frequency of each haplotype overall)
		//    ii) reconVarMap = (each haplotype, (called true variant position, reference 0 or alternate 1)
		//    iii) reconNames = (each haplotype, its name in String form)
		//    iv) intrapoolFreqs = (x = haplotypes, y = patients, cell = frequency of haplotype in that patient)
		BufferedReader brRes = new BufferedReader(new FileReader(folder + "p.all.results_1"));
		for (int l = 0; l < 20; l++) currLine = brRes.readLine(); 
		int numHaps = Integer.parseInt(currLine.split("\t")[1]);
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
		for (int f = 0; f <= maxPosDiff; f++) {	// GC + rjMCMC + EM results.
			// System.out.println(f + " position differences...");
			analysisReporter.add((double) numHaps); // Number of PHX-reconstructed haplotypes. 
			// analysisReporter.add((double) f); // Number of differences such that reconstruction is a quasispecies. 
			analysisReporter.addAll(ResultsReporter(reconVarMap, origVarMap, verifiedPos, reconNames, f, numHaps, 
					reconInterFreqs, origInterFreqs, numPts, intrapoolFreqs, origFreqs, maxFreqDiff));
		}
		return analysisReporter;
	}
	
	public static ArrayList<Double> OutputReporter(HashMap<Integer,HashMap<Integer,Integer>> reconVarMap, 
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
		
		// 5a) Six records comparing the original and GC-reconstructed haplotypes.
		//     i, ii, and iii) minDiffHap, minDiffHapInd, minDiffNum = (Closest reconstruction name, index, or the  
		//     number of variant mismatches)
		//     iv) allDiffNum = (original haplotype, (distance from each reconstruction))
		//     v) maxDiffOrLess = (original haplotype, (to keep track of haplotype indices that qualify as quasispecies))
		//     vi) diffAnalysis = (original haplotype, (zero to max differences, ))
		ArrayList<Double> outputData = new ArrayList<Double>(); 
		HashMap<Integer,String> minDiffHap = new HashMap<Integer,String>();
		HashMap<Integer, Integer> minDiffHapInd = new HashMap<Integer,Integer>();
		HashMap<Integer,Integer> minDiffNum = new HashMap<Integer,Integer>();
		Set<Integer> tmpRecons = reconVarMap.keySet(); 
		HashMap<Integer,ArrayList<Integer>> allDiffNum = new HashMap<Integer,ArrayList<Integer>>(); 
		HashMap<Integer, ArrayList<Integer>> maxDiffOrLess = new HashMap<Integer, ArrayList<Integer>>(); 
		HashMap<Integer, ArrayList<Integer>> diffAnalysis = new HashMap<Integer, ArrayList<Integer>>(); 
		for (int orig : origVarMap.keySet()) { // For each original haplotype (index)...
			int minHap = 0;
			int minNum = verifiedPos.size(); 
			allDiffNum.put(orig, new ArrayList<Integer>()); 
			for (int recon : tmpRecons) { // For each reconstructed haplotype index...
				int diff = 0; 
				for (int pos : verifiedPos) {
					if (origVarMap.get(orig).get(pos) != reconVarMap.get(recon).get(pos)) diff++;
				}
				allDiffNum.get(orig).add(diff); 
				if (diff < minNum) { // ...figure out the closest reconstructed haplotype.   
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
			// Count up the number of reconstructions at each distance less than the maximum to the original. 
			for (int tempDiff = 0; tempDiff < allDiffNum.get(orig).size(); tempDiff++) { 
				if (allDiffNum.get(orig).get(tempDiff) <= minDiffNum.get(orig) + maxPosDiff) 
					maxDiffOrLess.get(orig).add(tempDiff);
				if (allDiffNum.get(orig).get(tempDiff) <= (minNum + maxPosDiff)) {
					int tmp = diffAnalysis.get(orig).get(allDiffNum.get(orig).get(tempDiff) - minNum) + 1;
					diffAnalysis.get(orig).set(allDiffNum.get(orig).get(tempDiff) - minNum, tmp);
				}
			}
		}
		
		// 5b) Report reconstruction success as the closest distance and how many other reconstructions are as similar.  
		double propQuasiSum = 0; 
		double propErrorSum = 0; 
		for (int orig : minDiffHap.keySet()) {
			double propDiff = (double) minDiffNum.get(orig) / verifiedPos.size();  
			// Purpose: To print the number of differences between each reconstructed haplotype and the current original haplotype.
			// for (int d : allDiffNum.get(orig)) System.out.print(d + "\t");
			double ctQuasi =  0; 
			for (int d : diffAnalysis.get(orig)) ctQuasi += d;
			propQuasiSum += ctQuasi / numHaps;
			propErrorSum += propDiff; 
		}
		double propQuasi = propQuasiSum / minDiffHap.keySet().size(); 
		double propError = propErrorSum / minDiffHap.keySet().size(); 
		outputData.add(propError);	// The proportion of the variant composition of the most similar reconstructed haplotype that is incorrect.
		outputData.add(propQuasi);	// The proportion of all reconstructed haplotypes that are more similar than the threshold.	

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
		
		double sumFreqDiff = 0; 
		for (int h = 0; h < origVarMap.keySet().size(); h++) {
			double maxDiffOne = 0;
			for (int r : maxDiffOrLess.get(h)) maxDiffOne += reconInterFreqs[r];
			double ipDiff = maxDiffOne - origInterFreqs[h]; 
			sumFreqDiff += Math.abs(ipDiff); 
		}
		double avgFreqDiff = sumFreqDiff / origVarMap.keySet().size();
		outputData.add(avgFreqDiff); 	// The average difference between the interpool reconstructed and original frequencies.
		return outputData; 
	}
	
	public static ArrayList<Double> ResultsReporter(HashMap<Integer,HashMap<Integer,Integer>> reconVarMap, 
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
		ArrayList<Double> resultsData = new ArrayList<Double>(); 
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
		double propQuasiSum = 0; 
		double propErrorSum = 0; 
		for (int orig : minDiffHap.keySet()) {
			double propDiff = (double) minDiffNum.get(orig) / verifiedPos.size();
			// Purpose: To print the number of differences between each reconstructed haplotype and the current original haplotype.
			double ctQuasi =  0; 
			for (int d : diffAnalysis.get(orig)) ctQuasi += d;
			propQuasiSum += ctQuasi / numHaps;
			propErrorSum += propDiff; 
			// System.out.println(minDiffNum.get(orig) + " \t" + propDiff + " \t" + propErrorSum);
		}
		double propQuasi = propQuasiSum / (double) minDiffHap.keySet().size(); 
		double propError = propErrorSum / (double) minDiffHap.keySet().size(); 
		// System.out.println(propError);
		resultsData.add(propError);	// The proportion of the variant composition of the most similar reconstructed haplotype that is incorrect.
		resultsData.add(propQuasi);	// The proportion of all reconstructed haplotypes that are more similar than the threshold.	

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
		
		double sumFreqDiff = 0; 
		for (int h = 0; h < origVarMap.keySet().size(); h++) {
			double maxDiffOne = 0;
			for (int r : maxDiffOrLess.get(h)) maxDiffOne += reconInterFreqs[r];
			double ipDiff = maxDiffOne - origInterFreqs[h]; 
			sumFreqDiff += Math.abs(ipDiff); 
		}
		double avgFreqDiff = sumFreqDiff / origVarMap.keySet().size();
		resultsData.add(avgFreqDiff); 	// The average difference between the interpool reconstructed and original frequencies.
		
		// Purpose: To print the difference between the intrapool frequencies of the original and their closest reconstructed haploypes.
		int numOver = 0;
		int numUnder = 0;
		int numNP = 0; 
		int numAcc = 0; 
		for (int h = 0; h < origVarMap.keySet().size(); h++) {
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
		}

		double freqOver = (double) numOver / (origVarMap.keySet().size() * numPts);
		double freqUnder = (double) numUnder / (origVarMap.keySet().size() * numPts);
		double freqNP = (double) numNP / (origVarMap.keySet().size() * numPts); 
		double freqAcc = (double) numAcc / (origVarMap.keySet().size() * numPts); 
		resultsData.addAll(Arrays.asList(freqOver, freqUnder, freqNP, freqAcc)); // The proportions of intrapool frequency estimates that are (in)accurate. 
		return resultsData;
	}
}