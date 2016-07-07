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


public class WriteDifference_TU_NO_asSomaticFilter {

	// concatenates two strings
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


	public static void main(String[] args) throws Exception {
			// description:
			// 	this function filters the vcf from SAMtools, GATK Unified Genotyper and GATK Haplotype Caller. It is necessary because these tools also report many
			// 	germline mutations and give them very high scores -> therefore they are not comparable to the other callers deepSNV, JointSNVMix2, somaticSniper,
			//	VarScan2 which are sorted by some sort of somatic score --> no germline mutations or they have a very low score
			// inputs: 	VCF Tumor, VCF Normal, output_germline.vcf, output_somatic.vcf
			// outputs:	output_germline.vcf, output_somatic.vcf will be written full with the corresponding germline and somatic mutations
			// assumptions: a variant is counted as overlap, if the tumor sample does not have additional alternate alleles compared to the normal sample, at a given
			//	 position. the alternate alleles are comma-separated. the following column order is assumed:
			//	chr	pos	reference_allele	alternate_allele(s)
			
			if(args.length != 4){
				System.err.println("Usage: java WriteDifference_TU_NO_asSomaticFilter <VCF_tumor> <VCF_normal> <VCF_TUMORgermline> <VCF_TUMORsomatic>");
				return;
			}			

			String VCFListTumor=args[0];
			String VCFListNormal=args[1];
			String aLine = null;
			String chrTU, chrNO, refTU, refNO = null;
			String varNO, varTU = null;
			String posTU, posNO = null;
			Map<String,String> map1 = new HashMap();
			File VCF_TUMORgermlineFile = new File(args[2]);
			File VCF_TUMORsomaticFile = new File(args[3]);
			String delimiter="__";
	
			// continue only if the output files do not yet exist
			if(VCF_TUMORgermlineFile.exists() || VCF_TUMORsomaticFile.exists()){
				System.out.println( args[2] + " and/or " + args[3] + " already exist(s).");
				return;
			}

			// prepare the output files
			BufferedWriter VCF_TUMORgermlineWriter = new BufferedWriter(new FileWriter(VCF_TUMORgermlineFile));
			BufferedWriter VCF_TUMORsomaticWriter = new BufferedWriter(new FileWriter(VCF_TUMORsomaticFile));

			// prepare reading in the tumor and normal vcfs
			GZIPInputStream inVCFListTumor, inVCFListNormal;
			Reader decoderVCFListTumor, decoderVCFListNormal;
			BufferedReader VCFListTumorreader, VCFListNormalreader;

			if (VCFListTumor.endsWith(".gz")){ // gzipped input, as for SAMtools
				inVCFListTumor = new GZIPInputStream(new FileInputStream(VCFListTumor));
				decoderVCFListTumor = new InputStreamReader(inVCFListTumor);
				VCFListTumorreader = new BufferedReader(decoderVCFListTumor);

				// also the normal vcf:
				inVCFListNormal = new GZIPInputStream(new FileInputStream(VCFListNormal));
				decoderVCFListNormal = new InputStreamReader(inVCFListNormal);
				VCFListNormalreader = new BufferedReader(decoderVCFListNormal);

			} else { // not gzipped input, as is the case for GATK
				VCFListTumorreader = new BufferedReader(new FileReader(new File(VCFListTumor)));
				VCFListNormalreader = new BufferedReader(new FileReader(new File(VCFListNormal)));
			}


			// read in the genes from the normal mut list and store them in map1
			while((aLine = VCFListNormalreader.readLine())!=null){
				String[] splits = aLine.split("\t");
				if (!splits[0].startsWith("#")){ 	// a header line will be skipped; each other line will be stored in the hash map:
					chrNO = splits[0];
					posNO = splits[1];
					refNO = splits[3];
					varNO = splits[4];
					
					String hashID = chrNO + delimiter + posNO + delimiter + refNO;

					if(map1.containsKey(hashID)){
						System.out.println("The same position occurs twice in the normal list: " + hashID);
						System.out.println("Call it <hashID>"+delimiter+"Number2");
						String newhashID = hashID + delimiter + "Number2";
						if(map1.containsKey(newhashID)){
							System.out.println("The same position occurs three times in the normal list: " + hashID);
							System.out.println("Call it <hashID>"+delimiter+"Number3");
							newhashID = hashID +delimiter + "Number3";
							if(map1.containsKey(newhashID)){
								System.out.println("Error! The same position occurs more than three times in the list: "+ hashID);
								return;
							} else {
								map1.put(newhashID, hashID + delimiter + varNO);
							}
						} else {
							map1.put(newhashID, hashID + delimiter + varNO);
						}
					} else {
						map1.put(hashID,hashID + delimiter + varNO);
					}
				}
			}
			VCFListNormalreader.close();

			int cntGermline = 0;
			int cntSomatic = 0;
			int cntTumorMuts = 0;

			// read in the tumor vcf and compare to what is in the hash map
			while((aLine = VCFListTumorreader.readLine())!=null){
				String[] splits = aLine.split("\t");
				if (splits[0].startsWith("#")){ 	// a header line will be written to both output files
					VCF_TUMORgermlineWriter.write(aLine + "\n");
					VCF_TUMORsomaticWriter.write(aLine + "\n");
				} else {
					cntTumorMuts++;
					chrTU = splits[0];
					posTU = splits[1];
					refTU = splits[3];
					varTU = splits[4];
						
					String hashID = chrTU + delimiter + posTU + delimiter + refTU;
					if(map1.containsKey(hashID)){ // compare the var cases
						String[] NormalMut = map1.get(hashID).split(delimiter);
						// chr, pos, and ref have to be identical here in tumor and normal vcf
						if(!NormalMut[0].equals(chrTU) || ! NormalMut[1].equals(posTU) || !NormalMut[2].equals(refTU)){
							System.err.println("Error! chr, pos, and/or ref are not the same!");
							System.err.println("NormalMut[0] = " + NormalMut[0]);
							System.err.println("chrTU = " + chrTU);
							System.err.println("NormalMut[1] = "+NormalMut[1]);
							System.err.println("posTU = " + posTU);
							System.err.println("NormalMut[2] = " + NormalMut[2]);
							System.err.println("refTU = " + refTU);
							return;
						}
						
						// there could be several alternate alleles at the same position (separated by comma)
						String[] NormalVars = NormalMut[3].split(",");
						String[] TumorVars = varTU.split(",");
						
						// check if there are more mutations listed at the same position
						// if yes, read in the additional alternate normal alleles 
						if(map1.containsKey(hashID +delimiter +"Number2")){
							NormalMut = map1.get(hashID +delimiter + "Number2").split(delimiter);
							String[] addNormalVars = NormalMut[3].split(",");
							NormalVars = ConcatenateTwoStringArrays(NormalVars, addNormalVars);
						}
						if(map1.containsKey(hashID +delimiter + "Number3")){
							NormalMut = map1.get(hashID +delimiter +"Number3").split(delimiter);
							String[] addNormalVars = NormalMut[3].split(",");
							NormalVars = ConcatenateTwoStringArrays(NormalVars, addNormalVars);
						}

						// check for each tumor alternate allele whether it is also present in the normal
						boolean Germline = true;
						String alternateSOMATICalleles = null;
						for (String currTUvar : TumorVars ){
							boolean currVarIsNO = false;
							for (String currNOvar : NormalVars ){
								// if at least one alternate tumor allele is not present in the normal, set Germline to false
								if (currTUvar.equals(currNOvar)){
									currVarIsNO = true;
								} else if (currTUvar.contains(currNOvar) || currNOvar.contains(currTUvar)){
									currVarIsNO = true;
								}
							}
							if(!currVarIsNO){
								Germline=false;
								if(alternateSOMATICalleles==null){
									alternateSOMATICalleles=currTUvar;
								} else {
									alternateSOMATICalleles=alternateSOMATICalleles+","+currTUvar;
								}
							}
						}
						if (Germline){
							VCF_TUMORgermlineWriter.write(aLine + "\n");
							cntGermline++;
						} else {
							for(int i=0; i<4; i++){
								VCF_TUMORsomaticWriter.write(splits[i] + "\t");
							}
							VCF_TUMORsomaticWriter.write(alternateSOMATICalleles + "\t");
							for(int i=5; i<9; i++){
								VCF_TUMORsomaticWriter.write(splits[i] + "\t");
							}
							VCF_TUMORsomaticWriter.write(splits[9] + "\n");
							cntSomatic++;
						}
					} else {
						VCF_TUMORsomaticWriter.write(aLine + "\n");
						cntSomatic++;
					}
				}
			}
			VCFListTumorreader.close();
			VCF_TUMORgermlineWriter.close();
			VCF_TUMORsomaticWriter.close();
		
			System.out.println("Found " + cntTumorMuts + " mutations in the tumor list.");
			System.out.println("Counted " + cntGermline + " germline mutations, and "+ cntSomatic +" somatic mutations.");
			if(cntTumorMuts != (cntGermline+cntSomatic)){
				System.err.println("Error! cntTumorMuts != (cntGermline+cntSomatic)");
				return;
			}
	}
}

