import java.io.FileWriter;
import java.io.BufferedWriter;
import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.text.DecimalFormat;
import java.util.zip.GZIPInputStream;
import java.io.FileInputStream;
import java.io.FileOutputStream;
import java.util.zip.GZIPOutputStream;
import java.io.Reader;
import java.io.InputStreamReader;
import java.util.HashMap;
import java.util.Map;


public class createOriginalVCF {

	public static void main(String[] args) throws Exception {
			
			if(args.length < 3){
				System.err.println("Usage: java createOriginalVCF <posList> <originalVCF> <outVCF> [< --delimiter <char> >]");
				System.err.println("posList\t\t\t=\tthe file with the positions to extract; one line is e.g. chr1_123123");
				System.err.println("originalVCF\t\t=\tthe original vcf file, from which rows will be etracted");
				System.err.println("outVCF\t\t\t=\tthe generated vcf file, with only the positions from <posList>");
				System.err.println("Optional parameters:");
				//System.err.println("--delimiter <char>\t=\tthe delimiter which is between chr and the base position in the posList. [_]");
				System.err.println("Output:");
				System.err.println("outVCF");
				return;
			}
			
			String PosList=args[0];
			String VCFList=args[1];
			String outList=args[2];

			// Read in additional arguments if exist
			String delimiter="_"; 
			if(args.length>3){
				for(int j=3; j < args.length; j++){
					if(args[j].equals("--delimiter")){
						//delimiter=args[j+1].charAt(0);
						delimiter=args[j+1];
					}
				}
			}

			BufferedReader VCFListreader = new BufferedReader(new FileReader(new File(VCFList)));
			BufferedReader PosListreader = new BufferedReader(new FileReader(new File(PosList)));
			Map<String,String> map1 = new HashMap();
			String aLine = null;

			File outFile = new File(outList);
			// continue only if the output files do not yet exist
			if(outFile.exists()){
				System.out.println( outList + " already exists.");
				return;
			}
	
			// create the output file
			BufferedWriter outFileWriter = new BufferedWriter(new FileWriter(outFile));
	
			// read in the positions from the first list and store them in map1
			while((aLine = PosListreader.readLine())!=null){
				//String[] splits = aLine.split("_");
				map1.put(aLine, aLine);
			}
			PosListreader.close();
		
			// now read in the SNVs from the VCF list and check if they are in map1
			// if yes, write them out 
			while((aLine = VCFListreader.readLine())!=null){
				String[] splits = aLine.split("\t");
				if (splits[0].startsWith("#")){ // header line
					outFileWriter.write(aLine + "\n");
				} else {
					String chr_pos=splits[0] + delimiter + splits[1];
					if(map1.containsKey(chr_pos)){
						outFileWriter.write(aLine+"\n");
					}
				}
			}
			VCFListreader.close();
			outFileWriter.close();
	}
}
