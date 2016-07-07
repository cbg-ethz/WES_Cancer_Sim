import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.text.DecimalFormat;

public class RewriteSomaticSniperQualField {

	public static void main(String[] args) throws Exception {
			// description:
			// 	this function rewrites the vcf from the somaticSniper. It is necessary
			// 	because there is no information in the QUAL field
			// inputs: 	somaticSniper_VCF
			// outputs:	the rewritten VCF is given out on the command line
			
			if(args.length < 1){
				System.err.println("Usage: java -jar RewriteSomaticSniperQualField.jar <SomaticSniper.vcf>");
				return;
			}
			
			String SomaticSniperList=args[0];
			BufferedReader SomaticSniperListreader = new BufferedReader(new FileReader(new File(SomaticSniperList)));
			String aLine = null;
				
			while((aLine = SomaticSniperListreader.readLine())!=null){
				String[] splits = aLine.split("\t");
				if (splits[0].startsWith("#")){ 	// a header line
					System.out.println(aLine);  	// will be just printed out
				} else {							// each other line will be checked:
					String[] splits_info = splits[splits.length-1].split(":");
					String ssc = splits_info[splits_info.length-1];
					int somScore = Integer.parseInt(ssc);
				
					System.out.println(splits[0] + "\t" + splits[1] + "\t" + splits[2] + "\t" + splits[3] + "\t" + splits[4] + "\t" + somScore + "\t" + splits[6] + "\t" + splits[7]+ "\t" + splits[8] + "\t" + splits[9] + "\t" + splits[10]);	
				}
			}
			SomaticSniperListreader.close();
		}
}




