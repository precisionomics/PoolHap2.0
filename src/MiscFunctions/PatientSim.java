package MiscFunctions;

import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.Random;

public class PatientSim {

	public static void main(String[] args) throws IOException {
		BufferedReader br1 = new BufferedReader(new FileReader(args[0])); // args[0] needs to be $prefix.fa
		String currLine1 = br1.readLine();
		String currKey = null; 
		HashMap<String,ArrayList<String>> head2seq = new HashMap<String,ArrayList<String>>(); 
		while (currLine1 != null) {
			if (currLine1.contains(">")) {
				head2seq.put(currLine1.replace(">",""), new ArrayList<String>()); 
				currKey = currLine1.replace(">",""); 
			} else { 
				String[] tmpSeq = currLine1.split("");
				// System.out.println(tmpSeq.length + "\t" + currKey);
				for (String s : tmpSeq)  {
					// System.out.print(s + " ");
					head2seq.get(currKey).add(s);
				}
				// System.out.println("Done with " + currKey);
			}
			currLine1 = br1.readLine();
		}
		br1.close();
		BufferedReader br2 = new BufferedReader(new FileReader(args[1])); // args[1] needs to be $prefix.mutations.txt
		String currLine2 = br2.readLine();
		int refLength; 
		while (currLine2 != null) {
			String[] currRow = currLine2.split("\t");
			if (currRow[3].equals("-")) {	// Deletion
				refLength = currRow[2].length();
				for (int pos = Integer.parseInt(currRow[1]) - 1; pos < Integer.parseInt(currRow[1])- 1 + refLength; pos++)
					head2seq.get(currRow[0]).set(pos,"0"); 
			}
			else if (currRow[2].equals("-")) {	// Insertion
				String newIns = head2seq.get(currRow[0]).get(Integer.parseInt(currRow[1]) - 1) + currRow[3];
				head2seq.get(currRow[0]).set(Integer.parseInt(currRow[1]) - 1,newIns);
			}
			else head2seq.get(currRow[0]).set(Integer.parseInt(currRow[1]) - 1,currRow[3]);	// Substitution
			currLine2 = br2.readLine();
		}
		br2.close();
		int numPts = Integer.parseInt(args[2]);	// args[2] needs to be the number of patients.
		int numHaps = Integer.parseInt(args[3]);	// args[3] needs to be the number of haplotypes.
		Random rand = new Random();
		for (int p = 0; p < numPts; p++) {
			int ptHaps = rand.nextInt(13) + 3; 	// Between 3 and 15 individual viruses in each patient.
			PrintWriter pw = new PrintWriter(args[4] + "/p" + p + ".fa"); // args[4] needs to be $outdir
			for (int h = 0; h < ptHaps; h++) {
				int currHap = rand.nextInt(numHaps); // Randomly pick an index from the list of haplotypes.
				pw.append(">Haplotype_" + currHap + " \n");
				int seqLength = head2seq.get("Haplotype_" + currHap).size(); 
				for (int s = 0; s < seqLength; s++) {
					if (head2seq.get("Haplotype_" + currHap).get(s).equals("0")) continue;	// If this position is supposed to be deleted, don't add it.
					pw.append(head2seq.get("Haplotype_" + currHap).get(s));
				}
				pw.append("\n"); 
			}
			pw.close(); 
		}
	}
}
