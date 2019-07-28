package MiscFunctions;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.FileReader;
import java.io.FileWriter;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashSet;

public class TSR2SD {
	public static void main(String[] args) {
		//String folder=args[0]; 
		String output_fa_file=args[0]+".fa";//"XXX.fa";
		String first_output_file=args[0]+"_first_outcome.txt"; //"XXX_first_outcome.txt"
		String final_output_file=args[0]+"_final_outcome.txt";//"XXx_final_outcome.txt"
		transfer_file(output_fa_file,first_output_file,final_output_file);
	}
	
	/*
	 * transfer_file():
	 * The fa_file is generated by multiple sequence alignment tools (Clustal Omega)
	 * 1. From the fa_file, we first read the reference sequence into ArrayList<ArrayList<String>> (reference_listlist)
	 * 2. Then we read the following lines of each strains, the line started with ">S" are the name of the each strain and its frequency
	 * we save the strain frequency into the stran_frequency ArrayList, and calculate the number of strains(num_of_hap);
	 * 3. Next, the following lines after this title line till the next line stated with ">S" are sequence data for one strain;
	 * We compare the value at each base with the value at the same location in reference_listlist, and write the first_output_file,
	 * where each row is the sequence data for each strain. If they are the same, write "0",
	 * if they are different, write "-" and save the position in the variant_position HashSet, 
	 * if the value is "-", write "-", if the value is "*", write "*";
	 * 4.To be noticed, when it's "1", it means that this position have variants. 
	 * Therefore, first we save all positions that appear "1" into a variant_position Hashset(therefore, there will be no repetition)
	 * Then, transfer this variant_position Hashset into an Arraylist variant_position_list in oder to sort
	 * Finally, sort the variant_position_list
	 * Now, the variant_position_list contains all the positions that could have variants from the lowest position to the highest
	 * 5.For each strain, we calculate the recombinate rate and save it into recombinate_rate_list
	 */
	public static void transfer_file(String output_fa_file, String first_output_file, String final_output_file) {
		ArrayList<String> strain_frequency= new ArrayList<String>();
		ArrayList<ArrayList<String>> reference_listlist=new ArrayList<ArrayList<String>>();
		ArrayList<Integer> seg_site_number = new ArrayList<Integer>();	// For diG testing.
		HashSet<Integer> variant_position=new HashSet<Integer>(); 
		ArrayList<Double> recombinate_rate_list= new ArrayList<Double>();
		try {
			BufferedReader br_fa_file=new BufferedReader(new FileReader(output_fa_file));
			BufferedWriter bw_file=new BufferedWriter(new FileWriter(first_output_file));
			String fa_line=br_fa_file.readLine();//read the first line
			fa_line=br_fa_file.readLine();// read the second line, the reference sequence
			String[] each_base=fa_line.split("");
			int true_seg_site = 0; // For diG testing.
			while(!each_base[0].equals(">")) {
				ArrayList<String> sequence_row=new ArrayList<String>();
				for(int i=0;i<each_base.length;i++) {	
					sequence_row.add(each_base[i]);
					seg_site_number.add(true_seg_site); // For diG testing.
					if (!each_base[i+1].equals("-")) true_seg_site++; // If the reference sequence does not align with an insertion, the segregating site index is incremented.
				}
				reference_listlist.add(sequence_row);
				fa_line=br_fa_file.readLine();
				each_base=fa_line.split("");
			}
			// for (int i: seg_site_number) System.out.print(i + "\t");
			int num_of_hap=0;
			while(fa_line!=null) {
				if(each_base[0].equals(">")) {
					String[] each_position=fa_line.split(" ");
					strain_frequency.add(each_position[each_position.length-1]);
					num_of_hap++;
					fa_line=br_fa_file.readLine();//read the next line, the sequence line
					each_base=fa_line.split("");
					int sequence_line_count=0;
					int num_of_one=0;
					int num_of_zero=0;
					int num_of_empty=0;
					double recombinant_rate=0.0;
					while(!each_base[0].equals(">")) {
						for(int i=0;i<each_base.length;i++) {
							if(each_base[i].equals("-")) {
								if(each_base[i].equals(reference_listlist.get(sequence_line_count).get(i))) {  	// For diG testing to differentiate between a non-reference insertion and a deletion here.
									bw_file.write("0");
									num_of_zero++;			
								} else {
									bw_file.write("-");
									if(sequence_line_count!=0) {
										variant_position.add(sequence_line_count*reference_listlist.get(sequence_line_count-1).size()+i);
									}else if (sequence_line_count==0) {
										variant_position.add(i);
									}
									num_of_one++;
								}
								// bw_file.write("-");		// For TenSQR comparisons.
								// num_of_empty++;			// For TenSQR comparisons.
							}else if(each_base[i].equals(reference_listlist.get(sequence_line_count).get(i))) {
								bw_file.write("0");
								num_of_zero++;
							}else if(each_base[i].equals("*")) {
								bw_file.write("*");
								num_of_empty++;
							}else {
								if (reference_listlist.get(sequence_line_count).get(i).equals("-")) bw_file.write("+"); // For diG testing to identify a non-reference insertion.
								else bw_file.write("1");	// For TenSQR comparisons, delete 'else'.
								if(sequence_line_count!=0) {
									variant_position.add(sequence_line_count*reference_listlist.get(sequence_line_count-1).size()+i);
								}else if (sequence_line_count==0) {
									variant_position.add(i);
								}
								num_of_one++;
							}
						}
						fa_line=br_fa_file.readLine(); 
						if(fa_line!=null) {  // !!!this step is important, otherwise, when read till the last line, the line can not be splitted
							each_base=fa_line.split("");
							sequence_line_count++;
						}else {
							break;
						}
					}// end of while
					recombinant_rate=(double)(num_of_zero+num_of_one)/(double)(num_of_empty+num_of_one+num_of_zero);
					recombinate_rate_list.add(recombinant_rate);
					bw_file.write("\n");
				}// end of if
			}
			System.out.println(recombinate_rate_list);
			bw_file.close();
			br_fa_file.close();
			ArrayList<Integer> variant_position_list = new ArrayList<Integer>(variant_position); // hashset convert to arraylist
			Collections.sort(variant_position_list);
			final_output_file(strain_frequency, variant_position_list, seg_site_number, first_output_file,final_output_file,num_of_hap,recombinate_rate_list);
		}catch(Exception e) {e.printStackTrace();}
	}
	
	/*
	 * final_output_file: 
	 * 1. read the first_output_file,but this time, we only need to see what the value(1,0 or *) in those variant_position, the value will be saved in the variant_in_hap array
	 * 2. The variant_in_hap: Each column refers each Haplotype, and each rows refers to a variant_position and the values on the variant_position for each Haplotype
	 * 3. write the final_output_file: 
	 *   the first line: the Haplotype ID; 
	 *   the second line: frequency for each haplotype; 
	 *   the third line and line after: the value on each variant_position for each haplotype
	 */
	public static void final_output_file(ArrayList<String> strain_frequency, ArrayList<Integer> variant_position_list,ArrayList<Integer> seg_site_number, String first_output_file,String final_output_file,int num_of_hap,ArrayList<Double> recombinate_rate_list) {
		try {
			BufferedReader br_file=new BufferedReader(new FileReader(first_output_file));
			String[][] variant_in_hap=new String[variant_position_list.size()][num_of_hap];
			int line_count=0;
			String line=br_file.readLine();
			while(line!=null) {
				String [] each_position=line.split("");
				for(int i=0;i<variant_position_list.size();i++) {
					/*if(each_position[variant_position_list.get(i)].equals("1")) { //the base at the position of variant_position
						variant_in_hap[i][line_count]="1";
					}else if(each_position[variant_position_list.get(i)].equals("0")) {
						variant_in_hap[i][line_count]="0";
					}else {
							variant_in_hap[i][line_count]="*";
					}*/	// For TenSQR comparisons
					variant_in_hap[i][line_count]=each_position[variant_position_list.get(i)]; // For diG testing.
				}
				line=br_file.readLine();
				line_count++;
			}
			BufferedWriter bw_file=new BufferedWriter(new FileWriter(final_output_file));
			bw_file.write("Hap_ID"+"\t"); // write the first line
			for(int i=0;i<num_of_hap-1;i++) {
				bw_file.write(i+"\t");
			}
			bw_file.write((num_of_hap-1)+"\n");
			//write the second line
			bw_file.write("freq"+"\t");
			for(int i=0;i<strain_frequency.size()-1;i++) {
				bw_file.write(strain_frequency.get(i)+"\t");
			}
			bw_file.write(strain_frequency.get(strain_frequency.size()-1)+"\n");
			//write the following line
			for(int i=0;i<variant_position_list.size();i++) {
				// The diG testing version.
				bw_file.write("0;"+Integer.toString(seg_site_number.get(variant_position_list.get(i)) + 1)+";"+Integer.toString(seg_site_number.get(variant_position_list.get(i)) + 1)+";"+"0:1"+"\t");
				for(int j=0;j<num_of_hap-1;j++) {
					bw_file.write(variant_in_hap[i][j]+"\t");
				}
				bw_file.write(variant_in_hap[i][num_of_hap-1]+"\n");
			}
			bw_file.write("Recombinate Rate"+"\t");
			for(int i=0;i<recombinate_rate_list.size()-1;i++) {
				bw_file.write(String.format("%.3f",recombinate_rate_list.get(i))+"\t"); //double transfer to string
			}
			bw_file.write(String.format("%.3f",recombinate_rate_list.get(recombinate_rate_list.size()-1)));
			br_file.close();
			bw_file.close();
		}catch(Exception e) {e.printStackTrace();}	
	}
	
}
