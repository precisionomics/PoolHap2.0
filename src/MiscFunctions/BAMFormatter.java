package MiscFunctions;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.PrintWriter;
import java.util.HashMap;
import java.util.Scanner;
import java.util.Set;
import java.util.SortedSet;
import java.util.TreeSet;

public class BAMFormatter {
	
	public static void main(String[] args) throws FileNotFoundException {
		
		String BAMPrefix = args[0]; 
		File BAMPath = new File(BAMPrefix + ".bqsr.sam"); 	// args[0] needs to be $pfile. 
		Scanner BAMScanner = new Scanner(BAMPath);
		BAMScanner.useDelimiter("\t");
		String currLine = BAMScanner.nextLine(); 
		// System.out.println(currLine);
		while (currLine.matches("@(.*)")) {
			currLine = BAMScanner.nextLine(); 
		}
		// System.out.println(currLine);
		String QNAME, CIGAR, SEQ, tempPos = "",currAltAllele;
		Integer startPOS = 1, posToAdd, currentPos, matchToIndel, endPOS, skipSBases = 0;
		Character tempChar; 
		VarObj currVarObj; 
		HashMap<Integer, HashMap<Integer,String>> hmMATCH = new HashMap<Integer, HashMap<Integer,String>>();
		HashMap<Integer, HashMap<Integer,String>> hmINS = new HashMap<Integer, HashMap<Integer,String>>();
		HashMap<Integer, HashMap<Integer,String>> hmDEL = new HashMap<Integer, HashMap<Integer,String>>();
		HashMap<Integer,HashMap<String,VarObj>> variantEncoder = VariantMapping(args[1], BAMPrefix, true, args[4]);	// args[1] needs to be $outdir/p.all.vcf.
		Set<Integer> keyMATCH, keyINS, keyDEL, variantKS = variantEncoder.keySet(), tempSet, indelKS; 
		SortedSet<Integer> variantSS = new TreeSet<Integer>();
		variantSS.addAll(variantKS); 
		HashMap<Integer,String> defaultHM, tempFinderHM;
		HashMap<String,VarObj> tempReporterHM; 
		HashMap<Integer,Integer> hmVarCount = new HashMap<Integer,Integer>(); 
		for (Integer i : variantSS) {
			hmVarCount.put(i,0);
		}
		int readCount = 0;  
		
		PrintWriter VEFFile = new PrintWriter(BAMPrefix + ".vef"); 

		while (BAMScanner.hasNextLine()) {
			readCount++; 
			QNAME = BAMScanner.next().replace(":", "");
			for (int s = 0; s < 2; s++) {
				BAMScanner.next(); 	// Skip the FLAG, RNAME.
			}
			startPOS = BAMScanner.nextInt();
			currentPos = startPOS; 
			BAMScanner.next(); 		// Skip the MAPQ.			
			CIGAR = BAMScanner.next();
			
			if (CIGAR.charAt(0) == '*') {
				BAMScanner.nextLine(); 	// * refers to the fact that no CIGAR information is available. 
				continue;	 			// The above line skips the rest of the information in this read entirely.
			}
			for (int c = 0; c < CIGAR.length(); c++) {
				tempChar = CIGAR.charAt(c); 
				if (tempChar.compareTo('M') == 0) {
					posToAdd = Integer.parseInt(tempPos);
					for (int addPos = currentPos; addPos < currentPos + posToAdd; addPos++) {
						defaultHM = new HashMap<Integer,String>(); 
						defaultHM.put(1, "N"); 
						hmMATCH.put(addPos,defaultHM);
					}
					currentPos += posToAdd;
					tempPos = ""; 
				} else if (tempChar.compareTo('I') == 0) { 
					// Command to check for insertions: samtools view HIV_p1_sample_1.procd.bam | awk '{if ($6 ~ /I/) print $6;}' -
					posToAdd = Integer.parseInt(tempPos);
					matchToIndel = currentPos - 1; 	// This is the first position of the insertion, which it is recognized by. 
					hmMATCH.remove(matchToIndel); 	// Remove the first position of the insertion from the match HM.
					defaultHM = new HashMap<Integer,String>(); 
					defaultHM.put(posToAdd + 1, "N");	// This accounts for the entire length of the insertion.
					hmINS.put(matchToIndel, defaultHM); 
					currentPos = matchToIndel + 1;		// Since the insertion starts at the previous currentPos - 1, the next match along starts at matchToIndel + 1 = previous currentPos.
					tempPos = ""; 
				} else if (tempChar.compareTo('D') == 0) {
					posToAdd = Integer.parseInt(tempPos);
					matchToIndel = currentPos - 1; 	// This is the first position of the deletion, which it is recognized by. 
					hmMATCH.remove(matchToIndel); 	// Remove the first position of the deletion from the match HM.
					defaultHM = new HashMap<Integer,String>(); 
					defaultHM.put(posToAdd + 1, "N");	// This accounts for the entire length of the deletion.
					hmDEL.put(matchToIndel, defaultHM); 
					currentPos += posToAdd; 	// Since the deletion starts at the previous currentPos - 1, the next match along starts at currentPos + matchToIndel - 1.
					tempPos = ""; 
				} else if (tempChar.compareTo('S') == 0) {
					skipSBases = Integer.parseInt(tempPos);	// These bases do appear in SEQ BUT they are not aligned to the reference genome. Therefore these bases do not contribute to the currentPOS. 
					tempPos = ""; 				
				} else if (tempChar.compareTo('H') == 0) {
					tempPos = ""; // These bases DO NOT appear in SEQ and they are not aligned to the reference genome. Therefore these bases do not contribute to the currentPOS. 				
				} else { 
					tempPos += tempChar;
				}
			}

			keyMATCH = hmMATCH.keySet(); 
			keyINS = hmINS.keySet();
			keyDEL = hmDEL.keySet(); 
			for (int s = 0; s < 3; s++) {
				BAMScanner.next(); 	// Skip the RNEXT, PNEXT, TLEN.
			}
			SEQ = BAMScanner.next();
			endPOS = startPOS + SEQ.length();
			currentPos = startPOS; 
			// 4. Test the following to see if the SEQ-parsing code works properly, and the proper bases are assigned to the correct Hash Maps. 
			for (int s = skipSBases; s < SEQ.length(); s++) {	// Skip all of the bases that were soft-clipped but still reported at the beginning of SEQ. In the event there is no 'S' in the CIGAR string, this is 0 and we start reporting from the end of SEQ.  
				if (keyMATCH.contains(currentPos)) {
					hmMATCH.get(currentPos).put(1,SEQ.substring(s,s+1));
					currentPos++; 
				} else if (keyINS.contains(currentPos)) {
					tempSet = hmINS.get(currentPos).keySet();
					for (Integer length : tempSet) {
						hmINS.get(currentPos).put(length,SEQ.substring(s,s+length)); 
					}
					currentPos++; 
				} else if (keyDEL.contains(currentPos)) {
					tempSet = hmDEL.get(currentPos).keySet();
					for (Integer length : tempSet) {
						hmDEL.get(currentPos).put(length,SEQ.substring(s,s+1));
						currentPos += length; 
					}
				}
			}
			skipSBases = 0; // Reset the 'start' position of SEQ-processing to 0 in case there aren't any S-es in the next CIGAR string. 
			
			StringBuilder readInfo = new StringBuilder(); 
			for (Integer VCFVarPos : variantSS.subSet(startPOS,endPOS+1)) {	// subSet method is (inclusive,exclusive)
				tempReporterHM = variantEncoder.get(VCFVarPos);
				if (keyMATCH.contains(VCFVarPos)) {
					currAltAllele = hmMATCH.get(VCFVarPos).get(1);
					if (tempReporterHM.containsKey(currAltAllele)) {
						currVarObj = tempReporterHM.get(currAltAllele);
						// VEFFile.print(currVarObj.intM + "=" + currVarObj.varCode + ";");
						readInfo.append(currVarObj.intM + "=" + currVarObj.varCode + ";");
						hmVarCount.put(VCFVarPos, hmVarCount.get(VCFVarPos) + 1);
					}
				} else if (keyINS.contains(VCFVarPos)) {
					tempFinderHM = hmINS.get(VCFVarPos);
					indelKS = tempFinderHM.keySet(); 
					for (Integer indel : indelKS) {
						currAltAllele = tempFinderHM.get(indel);
						if (tempReporterHM.containsKey(currAltAllele)) {
							currVarObj = tempReporterHM.get(currAltAllele);
							// VEFFile.print(currVarObj.intM + "=" + currVarObj.varCode + ";");
							readInfo.append(currVarObj.intM + "=" + currVarObj.varCode + ";");
							hmVarCount.put(VCFVarPos, hmVarCount.get(VCFVarPos) + 1);
						}
					}
				} else if (keyDEL.contains(VCFVarPos)) {
					tempFinderHM = hmDEL.get(VCFVarPos);
					indelKS = tempFinderHM.keySet(); 
					for (Integer indel : indelKS) {
						currAltAllele = tempFinderHM.get(indel);
						if (tempReporterHM.containsKey(currAltAllele)) {
							currVarObj = tempReporterHM.get(currAltAllele);
							// VEFFile.print(currVarObj.intM + "=" + currVarObj.varCode + ";");
							readInfo.append(currVarObj.intM + "=" + currVarObj.varCode + ";");
							hmVarCount.put(VCFVarPos, hmVarCount.get(VCFVarPos) + 1);
						}
					}
				}
				else readInfo.append(VCFVarPos + "=0;");			
			}
			// System.out.println(readInfo.toString());
			if (readInfo.toString().isEmpty()) {
				BAMScanner.nextLine(); 	// Skip the QUAL.
				hmMATCH.clear();
				hmINS.clear();
				hmDEL.clear(); 
				continue;
			}
			VEFFile.println(QNAME + ":\t" + readInfo + "\t//\t" + startPOS + "\t" + endPOS); 
			BAMScanner.nextLine(); 	// Skip the QUAL.
			hmMATCH.clear();
			hmINS.clear();
			hmDEL.clear(); 
		}
		VEFFile.close();
		BAMScanner.close();
		double readLength = Double.parseDouble(args[2]); 
		double fullSeqLength = Double.parseDouble(args[3]);
		double estNumInd = readCount * readLength / fullSeqLength; 
		// System.out.println(readCount + " \t" + readLength + " \t" + fullSeqLength+ " \t" + estNumInd);
		PrintWriter CTFile = new PrintWriter(BAMPrefix + ".ct"); 
		// System.out.println(hmVarCount.size());
		for (Integer i : variantSS) {
			double normalized_varcount = hmVarCount.get(i) / estNumInd; 
			// System.out.println(hmVarCount.get(i));
			CTFile.append(normalized_varcount + "\t"); 
		}
		CTFile.append("\n");
		CTFile.close();
	}
	
	private static HashMap<Integer,HashMap<String,VarObj>> VariantMapping(String VCFFile, String BAMPrefix, boolean biallelic, String outdir) throws FileNotFoundException { 
		
		File VCFPath = new File(VCFFile); 
		Scanner VCFScanner = new Scanner(VCFPath); 	 
		VCFScanner.useDelimiter("\t"); // In Linux, the separator is a '\t', while in Windows, it is a ','.
		String currLine = VCFScanner.nextLine();
		// System.out.println(currLine);
		while (!currLine.contains("#CHROM")) {
			currLine = VCFScanner.nextLine();
		}
		// System.out.println(currLine);
		PrintWriter PosFile = new PrintWriter(outdir + "/p.all.pos"); 
		Integer pos;
		Scanner altAlleleScanner = new Scanner("");
		String variantNow;
		Integer variantCode = 1;  
		HashMap<String,VarObj> varEncAtPos; 
		HashMap<Integer,HashMap<String,VarObj>> variantEncoder = new HashMap<Integer,HashMap<String,VarObj>>();
		VarObj altAlleleAtPos; 
		while (VCFScanner.hasNextLine()) {
			VCFScanner.next();
			pos = VCFScanner.nextInt();
			PosFile.append(pos + "\t");
			varEncAtPos = new HashMap<String,VarObj>(); // This HashMap is specific to each variant position, and will be 'cleared' at each line of the VCF. 
			variantEncoder.put(pos,varEncAtPos);
			for (int s = 0; s < 2; s++) {
				VCFScanner.next(); 
			} 
			altAlleleScanner = new Scanner(VCFScanner.next());	// This Scanner runs through the String containing all of the alternate allele(s). There is at least one. 
			altAlleleScanner.useDelimiter(",");
			while(altAlleleScanner.hasNext()) {					// In case of highly polymorphic (>=2 alternate alleles) positions.
				variantNow = altAlleleScanner.next().replace("\"","");
				altAlleleAtPos = new VarObj(pos,variantCode);	// This is the information that will be reported when the alternate allele at this position in the read is accessed.
				varEncAtPos.put(variantNow,altAlleleAtPos);			// The information is mapped to the alternate allele. 
				variantCode++;										// The code corresponding to the alternate allele is incremented so no two alternate alleles in the same location have the same code.
				if (biallelic == true) break; 	// If PoolHap has not been adjusted for 3+ alleles. 
			}  
			variantCode = 1;	// Reset for the next patient-specific variant position.
			VCFScanner.nextLine(); 
		}
		VCFScanner.close();
		altAlleleScanner.close();
		PosFile.close();
		return variantEncoder;
	}
}

class VarObj {

	public Integer intM; 	// The position of the variant according to the metareference.
	public Integer varCode; // The numerical code corresponding to the desired alternative allele.  
	
	public VarObj(Integer i, Integer c){	// This is the constructor. 
		intM = i;
		varCode = c; 
	}
}

