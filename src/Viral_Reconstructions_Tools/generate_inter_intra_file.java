package Viral_Reconstructions_Tools;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.FileInputStream;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.io.InputStream;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Properties;

import spire.optional.intervalGeometricPartialOrder;

public class generate_inter_intra_file {
    String output_dir;
    String fasta_folder;
    String gs_dir;
	String project_name;
	int num_pools;
	ArrayList<ArrayList<String>> ref_seq_listlist=new ArrayList<ArrayList<String>>();
	ArrayList<ArrayList<String>> hap_seq_listlist=new ArrayList<ArrayList<String>>();
	ArrayList<String> hap_string_list=new ArrayList<String>();
	ArrayList<Integer> gs_variant_position = new ArrayList<Integer>(); 
	//ArrayList<Integer> rc_variant_position = new ArrayList<Integer>(); 
	HashSet<Integer> rc_variant_position = new HashSet<>();
	HashMap<String, ArrayList<Double>> hap2poolfre = new HashMap<String, ArrayList<Double>>();
	
	public generate_inter_intra_file(String parameter_file) throws IOException {
		 InputStream is = new FileInputStream(parameter_file);
	     Properties prop = new Properties();
	     prop.load(is);
	     this.output_dir = prop.getProperty("Output_Dir");
	     this.fasta_folder=this.output_dir+"fasta/";
	     this.gs_dir = prop.getProperty("Gold-Standard_Dir");
	     this.project_name=prop.getProperty("Proj_Name");
	     this.num_pools= Integer.parseInt(prop.getProperty("Num_Pools"));
	     
	     
	}
	
	public void true_var_position() throws IOException, InterruptedException{
		BufferedReader br = new BufferedReader(new FileReader(
				gs_dir+project_name+".inter_freq_vars.txt"));
		String currline = br.readLine();
		currline = br.readLine();
		currline = br.readLine();// read the third line
		while(currline != null) {
			String[] curr_positionStrings = currline.split("\t");
			String[] var_position = curr_positionStrings[0].split(";");
			// var_position_int is the position in the reference sequence
			this.gs_variant_position.add(Integer.parseInt(var_position[2]));
			currline = br.readLine();
		}
		br.close();
	}
	
	public void generate_hap2poolfre_hashmap() throws IOException, InterruptedException {
		BufferedReader br_fa = new BufferedReader(new FileReader(
				fasta_folder + project_name +".fa"));
		BufferedWriter bw_file=new BufferedWriter(new FileWriter(
				output_dir + project_name + "_first_output.txt"));
		String fa_line = br_fa.readLine();
		fa_line = br_fa.readLine();  // read the second sequence line
		String[] each_base=fa_line.split(""); 
		// This step is to add reference sequence into ref_seq_listlist
		// "-" that occurred in the reference sequence means insertion 
		while(!each_base[0].equals(">")) {
			ArrayList<String> sequence_row=new ArrayList<String>();
			for(int i=0;i<each_base.length;i++) {	
				sequence_row.add(each_base[i]);
				//this.seg_site_number.add(true_seg_site);	// New. 
				//if (!each_base[i+1].equals("-"))
				//	true_seg_site++;
				// New. If the reference sequence does not align with an insertion 
				//the segregating site index is incremented.
			}
			this.ref_seq_listlist.add(sequence_row);
			fa_line = br_fa.readLine();
			each_base=fa_line.split("");
		}
		// TO DO: When there is an insertion in the reference
		// The variant position is changed
		while(fa_line!=null) {
			if(each_base[0].equals(">")) {
				// write the header for each haplotype
				//bw_file.write(fa_line + "\n"); 
				String[] fre_position = fa_line.split("_");
				ArrayList<String> hap_sequence=new ArrayList<String>();
				for(int i=0; i < fre_position.length; i++) {
					hap_sequence.add(fre_position[i]);
				}
				fa_line=br_fa.readLine();//read the next line, the sequence line
				each_base=fa_line.split("");
				int sequence_line_count=0;
				while(!each_base[0].equals(">")) {
				// For diG testing to differentiate between a non-reference 
			    // insertion and a deletion here.
					for(int i=0;i<each_base.length;i++) {
						if(each_base[i].equals("-")) { 
							if(each_base[i].equals(ref_seq_listlist.get(
									sequence_line_count).get(i))) {  
						// 1st Situation:ref and each_base[i] both equals "-"
						//	bw_file.write("-");		
							} else { // 2nd Situation: Deletion
								hap_sequence.add("-");
								//bw_file.write("-"); 
							}
						}else if(!each_base[i].equals("-")){
							if( // 3rd Situation: It's the same as the ref
								each_base[i].equals(ref_seq_listlist.get(
										sequence_line_count).get(i))) {
							//bw_file.write("0");
							 hap_sequence.add("0");
							}else if(each_base[i].equals("*")) { 
						// 4th Situation: missing from the reconstruct sequence
						// bw_file.write("*");
							hap_sequence.add("*");
							}else if(ref_seq_listlist.get(
									sequence_line_count).get(i).equals("-")) {
						// 5th Situation: There is an insertion
						//bw_file.write("+"); 
							}else {
								hap_sequence.add("1");
							}
						}// end of if
//							}
					}//end of for
					fa_line = br_fa.readLine();
					if(fa_line!=null) {  
					// !!!this step is important, otherwise, 
					// when read till the last line, the line can not be splitted
						each_base=fa_line.split("");
						sequence_line_count++;
					}else {
						break;
					}
				} // end of while
				this.hap_seq_listlist.add(hap_sequence);
			}// end of if
		} // end of while
		
		// if the variant positions that shown as "-" 
		// or "1" are not in the gs_varian_position,
		// we need to get rid of this variant position,
		// and set it as "0"
		for (int i=0; i< hap_seq_listlist.size(); i++) {
			for (int j =3;j< hap_seq_listlist.get(i).size();j++ ) {
				if(hap_seq_listlist.get(i).get(j).equals("-")||hap_seq_listlist.get(i).get(j).equals("1")) {
					for(int pos = 0; pos < gs_variant_position.size(); pos++) {
						if(j == (gs_variant_position.get(pos)-1+3)) {
							this.rc_variant_position.add(gs_variant_position.get(pos));
						}else {
							hap_seq_listlist.get(i).set(j,"0");
						}						
					}
				}
			}
		}
		
		// Convert listlist to String
		for (int i=0; i< hap_seq_listlist.size(); i++) {
			String tmp_hap = "";
			for (int j =3;j< hap_seq_listlist.get(i).size();j++ ) {
				tmp_hap = tmp_hap + hap_seq_listlist.get(i).get(j);
			}
			this.hap_string_list.add(tmp_hap);
		}
		
		// Generate hap2poolfre_hashmap
//		for(int h =0; h < hap_string_list.size(); h++) {
//			if(hap2poolfre)
//		}

		
		
		System.out.println(gs_variant_position);
		System.out.println(gs_variant_position.size());
		System.out.println(hap_seq_listlist.size());
//		System.out.println(hap_seq_listlist.get(1));
//		System.out.println(hap_seq_listlist.get(hap_seq_listlist.size()-2));
//		System.out.println(hap_seq_listlist.get(0).size());
//		System.out.println(hap_seq_listlist.get(hap_seq_listlist.size()-1).size());
		System.out.println(rc_variant_position);
		System.out.println(rc_variant_position.size());
//		System.out.println(hap_string_list);
		System.out.println(hap_string_list.size());
	}
	
	

	public static void main(String[] args) throws IOException, InterruptedException{
		String parameter = "D:\\PhD-Studying\\Informatics\\Project\\HIV project\\Viral_reconstruction\\Other_tools_results\\CliqueSNV\\fasta\\FS3.properties";
		generate_inter_intra_file gf = new generate_inter_intra_file(parameter);
		gf.true_var_position();
		gf.generate_hap2poolfre_hashmap();

		
	}

}
