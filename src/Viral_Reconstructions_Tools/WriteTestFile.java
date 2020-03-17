package Viral_Reconstructions_Tools;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;

public class WriteTestFile {

	public static void main(String[] args) throws IOException{
		BufferedReader br = new BufferedReader(new FileReader("D:\\PhD-Studying"
				+ "\\Informatics\\Project\\Meta_Learning\\Prediction\\dataset\\"+
				"phg000010.MajorDepression.genotype-calls.Perlegen500K.v1.p1.c1.PRSC.filtered.matrixfmt.genotype.num.csv"));
		PrintWriter pw = new PrintWriter(new FileWriter("D:\\PhD-Studying"
				+ "\\Informatics\\Project\\Meta_Learning\\Prediction\\dataset\\"+ "phg000010.MajorDepression.genotype_test.csv", false)); 
		String currline=br.readLine();
		int line_count=0;
		while (line_count<1000) {
			pw.append(currline+"\n");
			currline=br.readLine();
			line_count++;
		}
		br.close();
		pw.close();
	}

}
