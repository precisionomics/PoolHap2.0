package Viral_Reconstructions_Tools;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.FileInputStream;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.io.InputStream;
import java.util.ArrayList;
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
	
	public void output_to_fasta_predicthaplo() throws IOException, InterruptedException{
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
		//First read the fasta_file_name.out to get the range of largest 
		//reconstruct global haplotype
		BufferedReader br_out = new BufferedReader(new FileReader(output_dir 
				+ fasta_file_name + ".out"));
		int line_count=0;
		String line = br_out.readLine();;
		while(line_count!=44) {
			line = br_out.readLine();
			line_count++;
		}
		String[] curr_sentence = line.split(". Last");
		String first_loci="";
		String final_loci="";
		String[] curr_word = line.split(" ");
		if(curr_word[0].equals("First")&&curr_word[1].equals("read")) {
			first_loci=curr_sentence[0].split(" ")[curr_sentence[0].
			                                              split(" ").length-1];
			final_loci=curr_word[curr_word.length-1];

		}else if(first_loci.equals("")||final_loci.equals("")){
			System.out.println("Didn't find the first loci or the final loci location");
		}
//			System.out.println(first_loci);
//			System.out.println(final_loci);
		br_out.close();
		
		//Read the global haplotype output file
		BufferedReader br = new BufferedReader(new FileReader(output_dir 
				+ fasta_file_name + "_global_"+ first_loci+"_"+final_loci+".fas"));
		int num_haps_per_pool = 0;
		String currline = br.readLine(); // read the first line
		while(currline!=null) {
			String[] each_position = currline.split("");
			if(each_position[0].equals(">")) {
				bw_fasta_file.write("\n");
				currline = br.readLine();
				String[] strain_fre_line=currline.split(":");
				bw_fasta_file.write(">"+ "Hap" + num_haps_per_pool +"_" 
					+ strain_fre_line[strain_fre_line.length-1] + "\n" );
				num_haps_per_pool++;
				while(!currline.equals(";EndOfComments")) {
					currline = br.readLine();
				}
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
	
	public void output_to_fasta_TenSQR() throws IOException, InterruptedException{
		BufferedReader br_ref_file = new BufferedReader(new FileReader(
				referece_seqence ));
		BufferedWriter bw_fasta_file = new BufferedWriter(new FileWriter(
					final_output_dir + fasta_file_name + ".fasta"));
		bw_fasta_file.write(">Reference"+"\n");
		ArrayList<String> ref_seq_list =new ArrayList<String>();
		String ref_line=br_ref_file.readLine();
		ref_line=br_ref_file.readLine();
		while(ref_line!=null) {
			String [] each_position=ref_line.split("");
			for(int i =0 ;i < each_position.length; i++) {
				bw_fasta_file.write(each_position[i]);
				ref_seq_list.add(each_position[i]);
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
			String[] each_word = currline.split(" ");
			if(each_word[0].equals("Viral")) {
				bw_fasta_file.write("\n");
				bw_fasta_file.write(">"+ "Hap" + num_haps_per_pool +"_" 
					+ each_word[each_word.length-1] + "\n" );
				num_haps_per_pool++;
			}else {
				for(int i =0; i < each_position.length; i++) {
					if(each_position[i].equals("*")) {
						bw_fasta_file.write(ref_seq_list.get(i));
					}else {
						bw_fasta_file.write(each_position[i]);
					}
				}					
			}
		currline = br.readLine();
		}
		br.close();
		bw_fasta_file.close();
	}
	public static void main(String[] args)throws IOException, InterruptedException  {
		String parameter = args[0];//"D:\\PhD-Studying\\Informatics\\Project\\HIV_project\\Viral_reconstruction\\PredictHaplo\\O2R.properties";//"D:\\PhD-Studying\\Informatics\\Project\\HIV project\\Viral_reconstruction\\Other_tools_results\\CliqueSNV\\fasta\\FS3.properties";
		InputStream is = new FileInputStream(parameter);
        Properties prop = new Properties();
        prop.load(is);
	    String tool_name = prop.getProperty("Tool_Name");
	    is.close(); 
		Transfer_TO_Fasta gf = new Transfer_TO_Fasta(parameter);
		if(tool_name.equals("CliqueSNV")||tool_name.equals("QuasiRecomb")
				||tool_name.equals("RegressHaplo")) {
			gf.output_to_fasta();
		}else if(tool_name.equals("PredictHaplo")) {
			gf.output_to_fasta_predicthaplo();
		}else if(tool_name.equals("TenSQR")) {	
			gf.output_to_fasta_TenSQR();
		}else {
			System.out.println("This tools are not in the list. Please choose "
					+ "from CliqueSNV,QuasiRecomb,RegressHaplo,PredictHaplo "
					+ "and TenSQR.");
		}

	}

}
