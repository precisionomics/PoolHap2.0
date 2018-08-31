package MiscFunctions;

import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;
import java.io.PrintWriter;

public class HaploWriter {
	
	public static void main(String[] args) throws IOException {
		BufferedReader br = new BufferedReader(new FileReader(args[0])); // args[0] needs to be param[3]
		String currLine = br.readLine();
		PrintWriter pw = new PrintWriter(args[1] + "/simhaps.fa"); // args[1] needs to be $outdir
		StringBuilder sequence = new StringBuilder();
		while (currLine != null) {
			if (currLine.contains(">")) {
				currLine = br.readLine();
				continue; 
			}
			else sequence.append(currLine.replace("\n", ""));
			currLine = br.readLine();
		}
		br.close();
		int numHaps = Integer.parseInt(args[2]); // args[2] needs to be the number of haplotypes.
		for (int h = 0; h < numHaps; h++) {
			pw.append(">Haplotype_" + h + "\n" + sequence.toString() + "\n");
		}
		pw.close(); 
	}
}