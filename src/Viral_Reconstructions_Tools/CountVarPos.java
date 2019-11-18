package Viral_Reconstructions_Tools;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;

public class CountVarPos {
	
	public static void countvars(String model, int num_file) throws IOException{
		String outfile = "varcount.out";
		BufferedWriter bw = new BufferedWriter(new FileWriter(outfile));
		for(int f=1;f<=num_file;f++) {
			String infile = "0_"+f+"_"+model+"_haploid.out";
			BufferedReader br = new BufferedReader(new FileReader(infile)); 
			String currLine = br.readLine(); // header
			String[] tmpcurrpos = currLine.split(" "); 
			if(model.equals("panmictic")) {
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
					currLine = br.readLine();
					tmpcurrpos = currLine.split(" ");
				}
				currLine = br.readLine();
				tmpcurrpos = currLine.split(" ");
			}
			int num_vars =0;
			while(!tmpcurrpos[0].equals("Genomes:")&&!tmpcurrpos[0].equals("Individuals:")) {
				//Read those lines under the "Mutations" section
				num_vars++;
				currLine = br.readLine();
				tmpcurrpos = currLine.split(" ");
			}
			bw.write("0_"+f+"\t"+num_vars+"\n");
			br.close();
		}
		bw.close();
		
	}

	public static void main(String[] args)throws IOException {
		String model = args[0];
		int num_file = Integer.parseInt(args[1]);
		countvars(model, num_file);
	}

}
