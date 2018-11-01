package MiscFunctions;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.concurrent.ThreadLocalRandom;
import org.apache.commons.lang3.ArrayUtils;

public class FullSimulator {
	
	// args[0] is the reference sequence.
	// args[1] is the variant composition of the simulated haplotypes.
	// args[2] is the rate of indels given the rate of mutations.
	// args[3] is the probability that an insertion is extended.
	// args[4] is the output directory.
	// args[5] is the number of patients.
	// args[6] is the number of haplotypes assigned to each patient.

	public static void main(String[] args) throws IOException {
		
		// 1) HashMap<totPos,ref> and ArrayList<availPos>. 
		// Put all reference bases in an indexed object.
		BufferedReader br1 = new BufferedReader(new FileReader(args[0])); 
		HashMap<Integer,String> refSequence = new HashMap<Integer,String>();
		String currLine = br1.readLine();
		int tmpPos = 0; 
		while (currLine != null) {
			if (currLine.contains(">")) {
				currLine = br1.readLine();
				continue; 
			}
			for (String b : currLine.split("")) {
				refSequence.put(tmpPos,b); 
				tmpPos++;
			}
			currLine = br1.readLine();
		}
		br1.close();
		ArrayList<Integer> availPos = new ArrayList<Integer>();
		availPos.addAll(refSequence.keySet());
		int refSeqLen = availPos.size(); 

		// 2) int[varPos] and int[hap][varPos] = refOrAlt.
		// Put all of the simulated haplotype information (from the ms simulation)
		// program into into accessible objects.
		BufferedReader br2 = new BufferedReader(new FileReader(args[1])); 
		currLine = br2.readLine();
		int rawNumHaps = Integer.parseInt(currLine.split(" ")[1]); 
		for (int i = 0; i < 6; i++) currLine = br2.readLine();
		int numVarPos = Integer.parseInt(currLine.split(" ")[1]);
		currLine = br2.readLine();
		String[] tmpVarPos = currLine.split(" "); 
		ArrayList<Integer> allVarPos = new ArrayList<Integer>(); 
		for (int p = 1; p < numVarPos + 1; p++) {
			allVarPos.add((int) Math.floor(Double.parseDouble(tmpVarPos[p]) * refSeqLen));
			// System.out.println(tmpVarPos.length + "\t" + p);			
		}
		currLine = br2.readLine();
		ArrayList<String> allOrigHaps = new ArrayList<String>();
		ArrayList<Integer> allOrigCts = new ArrayList<Integer>();  
		for (int h = 0; h < rawNumHaps; h++) {
			if (allOrigHaps.contains(currLine)) {
				int tmpCt =  allOrigCts.get(allOrigHaps.indexOf(currLine)) + 1;
				allOrigCts.set(allOrigHaps.indexOf(currLine), tmpCt);
			} else {
				allOrigHaps.add(currLine);
				allOrigCts.add(1);
			}
			currLine = br2.readLine();
		}
		br2.close();
		int numHaps = allOrigHaps.size(); 
		// PrintWriter tmp = new PrintWriter(args[4] + "/inpool_vars.txt");
		int[][] simHapVC = new int[numHaps][numVarPos]; 
		for (int h = 0; h < numHaps; h++) {
			String[] tmpHap = allOrigHaps.get(h).split(""); 
			for (int p = 0; p < numVarPos; p++) {
				simHapVC[h][p] = Integer.parseInt(tmpHap[p]); 
				// tmp.append(simHapVC[h][p] + "\t");
			}
			// tmp.append("\n");
		}
		// tmp.close();

		// 3) HashMap<varPos,altAllele>. Simulate mutations on the reference.
		HashMap<Integer,String> allAltAlleles = new HashMap<Integer,String>();
		double rateIndels = Double.parseDouble(args[2]); 
		double extendInsert = Double.parseDouble(args[3]); 
		for (int p = 0; p < numVarPos; p++)  {
			allAltAlleles.put(allVarPos.get(p),
				simVariant(refSequence.get(allVarPos.get(p)),rateIndels, extendInsert));
		}

		// 4) HashMap<hap,ArrayList<base>>. 
		// Add simulated mutations to finish the simulated haplotypes.
		HashMap<Integer,ArrayList<String>> allSimHaps = new HashMap<Integer,ArrayList<String>>();
		PrintWriter pw0 = new PrintWriter(args[4] + "/simhaps.mutations.txt"); 
		for (int h = 0; h < numHaps; h++) {
			allSimHaps.put(h, new ArrayList<String>()); 
			for (int p = 0; p < refSeqLen; p++) {
				if (allVarPos.indexOf(p)!=-1) {
					int varPosIndex = allVarPos.indexOf(p);
					if (simHapVC[h][varPosIndex]==1) {
						allSimHaps.get(h).add(allAltAlleles.get(p));
						int oneIndexed = p + 1;
						pw0.append("Haplotype_" + h + "\t" + oneIndexed + "\t" + refSequence.get(p) +"\t");
						if (allAltAlleles.get(p).isEmpty()) pw0.append("*\n");
						else pw0.append(allAltAlleles.get(p) + "\n");
					} else 
						allSimHaps.get(h).add(refSequence.get(p)); 
				} else allSimHaps.get(h).add(refSequence.get(p));
			}
		}
		pw0.close();

		// 5) Simulate the patients and give a random number of haplotypes.
		int numPts = Integer.parseInt(args[5]);	
		int ptHaps = Integer.parseInt(args[6]); 
		double[][] origCts = new double[numHaps][numPts];
		double[][] origFreqs = new double[numHaps][numPts];
		double[][] poolCts = new double[numVarPos][numPts];
		double[][] poolFreqs = new double[numVarPos][numPts];
		double[] globalCt = new double[numHaps]; 
		double[] globalFreq = new double[numHaps]; 
		for (int p = 0; p < numPts; p++) {
			PrintWriter pw = new PrintWriter(args[4] + "/p" + p + ".fa");
			for (int h = 0; h < ptHaps; h++) {
				int currHap = ThreadLocalRandom.current().nextInt(0, numHaps); // Randomly pick an index from the list of haplotypes.
				while (allOrigCts.get(currHap) == 0) {
					currHap = ThreadLocalRandom.current().nextInt(0, numHaps);
				}
				origCts[currHap][p]++; 
				pw.append(">Haplotype_" + currHap + " \n");
				for (String s : allSimHaps.get(currHap)) pw.append(s);
				pw.append("\n\n"); 
				for (int v = 0; v < numVarPos; v++) {
					poolCts[v][p] += simHapVC[h][v]; 
				}
				globalCt[currHap]++; 
				int newCt = allOrigCts.get(currHap) - 1; 
				allOrigCts.set(currHap, newCt);
			}
			pw.close();
			for (int h = 0; h < numHaps; h++) origFreqs[h][p] = origCts[h][p] / ptHaps;
			for (int v = 0; v < numVarPos; v++) poolFreqs[v][p] = poolCts[v][p] / ptHaps;
		}
		for (int h = 0; h < numHaps; h++) globalFreq[h] = globalCt[h] / rawNumHaps;
		
		BufferedWriter bv_glob= new BufferedWriter(new FileWriter(args[4] + "simhaps.inter_freq_vars.txt"));
		bv_glob.write("Hap_ID");
		for(int h=0;h<numHaps;h++) bv_glob.write("\t"+h);
		bv_glob.write("\nFreq");
		for(int h=0;h<numHaps;h++) bv_glob.write("\t"+globalFreq[h]); // Global freq
		for(int v=0;v<numVarPos;v++){
			int oneIndexed = allVarPos.get(v) + 1; 
			bv_glob.write("\n0;" + oneIndexed + ";" + oneIndexed + ";0:1\t"); // TODO This only allows for biallelic simple loci (single alternate allele) for now. (int h=0;h<numHaps;h++)
			for(int h=0;h<numHaps;h++) bv_glob.write("\t"+simHapVC[h][v]);
		}
		bv_glob.close();
		
		BufferedWriter bv_loc= new BufferedWriter(new FileWriter(args[4] + "simhaps.intra_freq.txt"));
		bv_loc.write("Hap_ID");
		for(int h=0;h<numHaps;h++) bv_loc.write("\t"+h);
		bv_loc.write("\n");
		for(int p=0;p<numPts;p++){
			bv_loc.write(p + "");
			for(int h=0;h<numHaps;h++)
				bv_loc.write("\t"+origFreqs[h][p]);
			bv_loc.write("\n");
		}
		bv_loc.close();
		
		BufferedWriter bvp= new BufferedWriter(new FileWriter(args[4] + "simvars.intra_freq.txt"));
		bvp.write("Var_ID\t");
		for(int p=0;p<numPts;p++) bvp.write(p + "\t");
		for(int v=0;v<numVarPos;v++) {
			bvp.write("\n0;" + allVarPos.get(v) + ";" + allVarPos.get(v) + ";0:1\t"); // TODO This only allows for biallelic simple loci (single alternate allele) for now. 
			for(int p=0;p<numPts;p++) bvp.write("\t"+poolFreqs[v][p]);
		}
		bvp.close();
	}

	static String simVariant(String refBase, double rateIndels, double extendInsert) {
		String[] bases = new String[] {"A","C","G","T"};
		if (ThreadLocalRandom.current().nextDouble() < rateIndels) {
			if (ThreadLocalRandom.current().nextDouble() > 0.5) {
				return "";
			} else { 
				StringBuilder insertBases = new StringBuilder(); 
				insertBases.append(refBase); 
				insertBases.append(bases[ThreadLocalRandom.current().nextInt(0, 4)]); 
				double extension = ThreadLocalRandom.current().nextDouble(); 
				extension = ThreadLocalRandom.current().nextDouble(); 
				while (extension < extendInsert) {
					insertBases.append(bases[ThreadLocalRandom.current().nextInt(0, 4)]); 
					extension = ThreadLocalRandom.current().nextDouble(); 
				}
				return insertBases.toString(); 
				
			}
		} else {
			bases = ArrayUtils.removeElement(bases, refBase);
			// for (String b : bases) System.out.print(b);
			// System.out.println();
			return bases[ThreadLocalRandom.current().nextInt(0, 3)]; 
		}
	}
}