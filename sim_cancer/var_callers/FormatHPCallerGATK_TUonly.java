import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;


public class FormatHPCallerGATK_TUonly {

	public static void main(String[] args) throws Exception {
		// description:
		// 	this function rewrites the vcf from the GATK HaplotypeCaller. It is necessary
		// 	because some of the GT:AD:DP:GQ:PL fields are malformed
		// inputs: 	GATK_Haplotype_Caller_VCF
		// outputs:	the rewritten vcf is given out on the command line
		
		if(args.length < 1){
			System.err.println("Usage: java -jar FormatHPCallerGATK.jar <HPCallerGATK.vcf>");
			return;
		}
		
		String HPCallerList=args[0];
		
		BufferedReader HPCallerListreader = new BufferedReader(new FileReader(new File(HPCallerList)));
		String aLine = null;
		String newSplits9;
		String newSplits8 = "GT:AD:DP:GQ:PL";
		
		int cnt_dot = 0;
		while((aLine = HPCallerListreader.readLine())!=null){
			String[] splits = aLine.split("\t");
			if (splits[0].startsWith("#")){ 	// a header line
				System.out.println(aLine);  	// will be just printed out
			} else {							// each other line will be checked:
				String[] splits_GT_info9 = splits[9].split(":");
								
				if (splits_GT_info9.length == 5){
					if (splits_GT_info9[1].equals(".")){
						splits_GT_info9[1] = "0,0";
						cnt_dot++;
					}
					if (splits_GT_info9[2].equals(".")){
						splits_GT_info9[2] = "0";
						cnt_dot++;
					}
					newSplits9 = splits_GT_info9[0] + ":" + splits_GT_info9[1] + ":" + splits_GT_info9[2] + ":" + splits_GT_info9[3] + ":" + splits_GT_info9[4];
				} else {
					if (splits[9].equals("./.")){
						newSplits9 = ".:0,0:0:0:0,0,0";
					} else {
						newSplits9 = splits_GT_info9[0] + ":0,0:0:" + splits_GT_info9[1] + ":" + splits_GT_info9[2];
					}
				}
				System.out.println(splits[0] + "\t" + splits[1] + "\t" + splits[2] + "\t" + splits[3] + "\t" + splits[4] + "\t" + splits[5] + "\t" + splits[6] + "\t" + splits[7] + "\t" + newSplits8 + "\t" + newSplits9);	
			}
		}
		//System.out.println("The number of times the dot was replaced: " + cnt_dot);
		HPCallerListreader.close();
	}
}
