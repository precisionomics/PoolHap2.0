package Viral_Reconstructions_Tools;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.FileInputStream;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.io.InputStream;
import java.util.Properties;

public class Transfer_TO_Fasta {
	String output_dir;
	String outfile_name;
    String fasta_file_name;
    String final_output_dir;
	String referece_seqence;
	
	public Transfer_TO_Fasta(String parameter_file) throws IOException {
		 InputStream is = new FileInputStream(parameter_file);
	     Properties prop = new Properties();
	     prop.load(is);
	     this.output_dir = prop.getProperty("Output_Dir");
	     this.outfile_name = prop.getProperty("Out_File_Name");
	     this.final_output_dir=prop.getProperty("Final_Output_Dir");
	     this.fasta_file_name=prop.getProperty("Fasta_File_Name");
	     this.referece_seqence = prop.getProperty("Ref_Seq_Path");
	     is.close();
	     
	}

	public void output_to_fasta() throws IOException, InterruptedException{
		BufferedReader br_ref_file = new BufferedReader(new FileReader(
				referece_seqence ));
		BufferedWriter bw_fasta_file = new BufferedWriter(new FileWriter(
					final_output_dir + fasta_file_name + ".fasta"));
		bw_fasta_file.write(">Reference"+"\n");
		String ref_line=br_ref_file.readLine();
		ref_line=br_ref_file.readLine();
		while(ref_line!=null) {
			String [] each_position=ref_line.split("");
			for(int i =0 ;i < each_position.length; i++) {
				bw_fasta_file.write(each_position[i]);
			}
			ref_line=br_ref_file.readLine();
		}
		br_ref_file.close();
		int num_haps_per_pool = 0;
		BufferedReader br = new BufferedReader(new FileReader(output_dir 
				+ outfile_name));
		String currline = br.readLine(); // read the first line
		while(currline!=null) {
			String[] each_position = currline.split("");
			if(each_position[0].equals(">")) {
				bw_fasta_file.write("\n");
				String[] strain_fre_line=currline.split("_");
				bw_fasta_file.write(">"+ "Hap" + num_haps_per_pool +"_" 
					+ strain_fre_line[strain_fre_line.length-1] + "\n" );
				num_haps_per_pool++;
			}else {
				for(int i =0; i < each_position.length; i++) {
					bw_fasta_file.write(each_position[i]);
				}
					
			}
		currline = br.readLine();
		}
		br.close();
		bw_fasta_file.close();
	}
	public static void main(String[] args)throws IOException, InterruptedException  {
		String parameter = args[0];//"D:\\PhD-Studying\\Informatics\\Project\\HIV project\\Viral_reconstruction\\Other_tools_results\\CliqueSNV\\fasta\\FS3.properties";
		Transfer_TO_Fasta gf = new Transfer_TO_Fasta(parameter);
		gf.output_to_fasta();
	}

}
