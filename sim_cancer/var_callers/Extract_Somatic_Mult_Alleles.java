import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.text.DecimalFormat;
import java.util.zip.GZIPInputStream;
import java.io.FileInputStream;
import java.io.FileOutputStream;
import java.util.zip.GZIPOutputStream;
import java.io.Reader;
import java.io.Writer;
import java.io.FileWriter;
import java.io.BufferedWriter;
import java.io.InputStreamReader;
import java.util.HashMap;
import java.util.Map;


public class Extract_Somatic_Mult_Alleles {

	// concatenates two strings arrays
	public static String[] ConcatenateTwoStringArrays( String[] firstString, String[] secondString) {
		String[] combinedString = new String[firstString.length + secondString.length];
		// copy first half
		System.arraycopy(firstString, 0, combinedString, 0, firstString.length);
		// copy second half
		System.arraycopy(secondString, 0, combinedString, firstString.length, secondString.length);
		if(combinedString.length != (firstString.length + secondString.length)){
			System.err.println("Error! Srtings have different lengths.");
			return(null);
		}

		return(combinedString);
	}

	// concatenate string and string array
	public static String[] ConcatenateStringAndStringArrays(String firstString, String[] secondString){
		//String[] combinedString = new String[1 + secondString.length];
		StringBuilder builder = new StringBuilder();
		builder.append(firstString + "\t");
		for(int i=0; i < secondString.length; i++){
			String s = secondString[i];
			if(i==secondString.length-1){
				builder.append(s);
			} else {
				builder.append(s + "\t");
			}
		}
		//return(combinedString);
		//String[] combinedStringPlusEmpty = builder.toString().split("\t");
		//String[] combinedString = new String[combinedStringPlusEmpty.length-1];
		//int i=0;
		//for(String s : combinedStringPlusEmpty){
		//	if(i==combinedStringPlusEmpty.length-1){
		//		break;
		//	}
		//	combinedString[i] = s;
		//	System.out.println("Curr string: " + s);
		//	i++;
		//}
		String[] combinedString = builder.toString().split("\t");
		//System.out.println("This is the concatenated string array:");
		//for(String s : combinedString){
		//	System.out.println(s);
		//}
		return(combinedString);
	}

	public static void main(String[] args) throws Exception {
			if(args.length < 1){
				System.err.println("Usage: java Extract_Somatic_Mult_Alleles <VCF> <filteredVCF>");
				System.err.println("VCF\t\t\t=\tthe vcf file to filter");
				System.err.println("filteredVCF\t\t=\tthe filtered vcf file, which will be generated");
				System.err.println("Output:");
				System.err.println("filteredVCF");
				System.err.println("Description: This script takes a VCF file from somaticSniper, and filters out germline mutations, which might be still in the VCF. These are the loci, where there are several alternate alleles reported, e.g. A -> C,G, but only one of them is actually somatic, and the other one is a germline mutation. The variant might still get a high score, but is of course counted as false positive because it is germline. Here, we split alternate alleles that are reported with a comma, and only write out the ones which are somatic. The status (germline/somatic), we do get from the genotype field.");
				return;
			}
			
			String VCFList=args[0];
			String filteredVCFList=args[1];

			String aLine = null;
			String altAllele = null;
			String genotypeNO = null;
			String genotypeTU = null;
			int columnNO = 9;
			int columnTU = 10;

			File outFile = new File(filteredVCFList);
	
			// continue only if the output files do not yet exist
			if(outFile.exists()){
				System.out.println( args[1] + " already exists.");
				return;
			}

			// prepare the output file
			BufferedWriter outWriter = new BufferedWriter(new FileWriter(outFile));

			// prepare reading in the vcf
			GZIPInputStream inVCFList;
			Reader decoderVCFList;
			BufferedReader VCFListReader;

			if (VCFList.endsWith(".gz")){ // gzipped input
				inVCFList = new GZIPInputStream(new FileInputStream(VCFList));
				decoderVCFList = new InputStreamReader(inVCFList);
				VCFListReader = new BufferedReader(decoderVCFList);
			} else { // not gzipped input
				VCFListReader = new BufferedReader(new FileReader(new File(VCFList)));
			}


			// go through the vcf and check for each case, where we have several alternate alleles
			int TotalCnt=0;
			int GermlineCnt=0;
			int SomaticCnt=0;
			while((aLine = VCFListReader.readLine())!=null){
				String[] splits = aLine.split("\t");
				if (splits[0].startsWith("#")){ 	// a header line will be written to the output file
					outWriter.write(aLine + "\n");
				} else {
					TotalCnt++;
					altAllele = splits[4];
					if (altAllele.contains(",")){ // in the case we have several alternate alleles (they are always comma separated)
						//split the alt allele, add the ref in the front, and then we can with the genotype easily extract those alleles that are new in the tumor!
						String refAllele = splits[3];
						String[] allAltAlleles = altAllele.split(",");
						//System.out.println("Curr ref allele: " + refAllele);
						//System.out.println("Curr alternate allele: " + altAllele);
						String[] allAlleles = ConcatenateStringAndStringArrays(refAllele,allAltAlleles);
						//System.out.println("These are all alleles:");
						//for(String s : allAlleles){
						//	System.out.println(s);
						//}
						genotypeNO = splits[columnNO].split(":")[0]; // e.g. 0/1
						genotypeTU = splits[columnTU].split(":")[0]; // e.g. 1/2
						//System.out.println("Normal genotype: " + genotypeNO);
						//System.out.println("Tumor genotype: " + genotypeTU);
						int[] genotypeNOnums = new int[2];
						int[] genotypeTUnums = new int[2];
						

						String[] splittedGT_NO = genotypeNO.split("");
						String[] splittedGT_TU = genotypeTU.split("");
						//for(String s : splittedGT_NO){
						//	System.out.println(s);
						//}

						int j=0;
						for(String currGT : genotypeNO.split("")){
							if(currGT.isEmpty() || currGT.equals("/")){
								continue;
							}
							if(currGT.equals(".")){
								genotypeNOnums[j]=0; // if there is no coverage at that position, we assume the genotype is the reference allele
							} else {
								genotypeNOnums[j]=Integer.parseInt(currGT);
							}
							//System.out.println("Current normal genotpe is: " + genotypeNOnums[j] );
							j++;
						}
						j=0;
						for(String currGT : genotypeTU.split("")){
							if(currGT.isEmpty() || currGT.equals("/")){
								continue;
							}
							if(currGT.equals(".")){
								genotypeTUnums[j]=0; // if there is no coverage at that position, we assume the genotype is the reference allele
							} else {
								genotypeTUnums[j]=Integer.parseInt(currGT);
							}
							//System.out.println("Current tumor genotpe is: " + genotypeTUnums[j] );
							j++;
						}

						// Now find which alleles are somatic (only in tumor, not in normal)
						for (int altAlleleNum=0; altAlleleNum < allAlleles.length; altAlleleNum++){ // loop over all alleles at the loci
							String currAlt = allAlleles[altAlleleNum];
							boolean inNormal = false;
							boolean inTumor = false;
							// check whether it occurs in the normal genotype
							for (j=0; j < genotypeNOnums.length; j++ ){
								if(altAlleleNum==genotypeNOnums[j]){
									inNormal=true;
								}
							}
							// check whether it occurs in the tumor genotype
							for(j=0; j<genotypeTUnums.length;j++){
								if(altAlleleNum == genotypeTUnums[j]){
									inTumor=true;
								}
							}
							if(inTumor && !inNormal){ // current allele is not in the normal genotype, but it is in the tumor genotype --> somatic!
								SomaticCnt++;
								// #CHROM  POS     ID      REF     ALT     QUAL    FILTER  INFO    FORMAT  NO      TU
								// write it out:
								for (j=0; j<splits.length-1; j++){
									if(j==4){
										outWriter.write(currAlt+"\t");
									} else {
										outWriter.write(splits[j]+"\t");
									}
								}
								outWriter.write(splits[splits.length-1]+"\n");
							} else {
								GermlineCnt++; // do not write out the allele
								//System.out.println(currAlt + " is germline:");
								//System.out.println(aLine);
							}
						}
					} else {
						SomaticCnt++;
						outWriter.write(aLine + "\n"); // write out line as is, because there is only one alternate allele
					}
				}
			}
			VCFListReader.close();
			outWriter.close();

			System.out.println("Found in total " + TotalCnt + " mutations (=loci; each can have several variants). Of these, we determined " +GermlineCnt+ " to be germline mutations, and " +SomaticCnt+ " to be somatic mutations.");
			if(SomaticCnt+GermlineCnt < TotalCnt){
				System.err.println("Error! SomaticCnt+GermlineCnt < TotalCnt");
				return;
			}
	}
}

