import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;


public class SNV_Overlapper_Ordered_PlusAlleles {

	public static void main(String[] args) throws Exception{
		String list1=args[0];
		String list2=args[1];
		String overlap=args[2];
		String private_list=args[3];
		
		BufferedReader list1reader = new BufferedReader(new FileReader(new File(list1)));
		String aLine = null;
		String bLine = null;
		String lastBLine = null;
		BufferedWriter overlapWriter = new BufferedWriter(new FileWriter(new File(overlap)));
		BufferedWriter privateWriter = new BufferedWriter(new FileWriter(new File(private_list)));
		BufferedReader list2reader = new BufferedReader(new FileReader(new File(list2)));
		int cntPrivate=0;
		int cntOverlap=0;
		int cntTotal=0;

		lastBLine = list2reader.readLine();	
		while((aLine = list1reader.readLine())!=null){
			String[] splitsA = aLine.split("\t");
			String ACurrChrString = splitsA[0];
			int ACurrPos = Integer.parseInt(splitsA[1]);
			String ACurrRefAllele = splitsA[3];
			String ACurrAltAllele = splitsA[4];
			//int ACurrChrNum = Integer.parseInt(ACurrChrString.substring(3)); // String chr123123 -> int 123123

			//System.out.println(ACurrChrString + "\t" + ACurrPos);
			boolean isOverlap = false;

			bLine = lastBLine;
			if(bLine!=null){
				
				String[] splitsB = bLine.split("\t");
				String BCurrChrString = splitsB[0];
				int BCurrPos = Integer.parseInt(splitsB[1]);
				String BCurrRefAllele = splitsB[3];
				String BCurrAltAllele = splitsB[4];
				//int BCurrChrNum = Integer.parseInt(BCurrChrString.substring(3)); // String chr123123 -> int 123123
		
				// the B file is lower and needs to be incremented
				// example: int result = str1.compareTo(str2); if result < 0 => str1 is before str2
				while(BCurrChrString.compareTo(ACurrChrString)<0  || ( BCurrChrString.compareTo(ACurrChrString)==0 && BCurrPos < ACurrPos) ){ 
					bLine = list2reader.readLine();
					if(bLine==null){
						break; // should only break the inner while loop
					}
					splitsB = bLine.split("\t");
					BCurrChrString = splitsB[0];
					BCurrPos = Integer.parseInt(splitsB[1]);
					BCurrRefAllele = splitsB[3];
					BCurrAltAllele = splitsB[4];
					//BCurrChrNum = Integer.parseInt(BCurrChrString.substring(3)); // String chr123123 -> int 123123
				}
				lastBLine = bLine;

				//// the A file is lower and needs to be incremented; EDIT: This happens automatically
				//if(BCurrChrNum > ACurrChrNum || BCurrPos > ACurrPos ){
				//	continue;
				//}

				if(bLine!=null){
					// they are the same; hence overlap!
					if(BCurrChrString.equals(ACurrChrString) && BCurrPos == ACurrPos && ACurrRefAllele.equals(BCurrRefAllele) && ACurrAltAllele.equals(BCurrAltAllele)){
						isOverlap = true;
					}
				}
			}
			
			cntTotal++;
			if(isOverlap){
				overlapWriter.write(aLine + "\n");
				cntOverlap++;
			} else {
				privateWriter.write(aLine + "\n");
				cntPrivate++;
			}
		}
		if((cntOverlap+cntPrivate)!=cntTotal){
			System.err.println("Error! (cntOverlap+cntPrivate)!=cntTotal");
		}
		list2reader.close();
		list1reader.close();
		overlapWriter.close();
		privateWriter.close();	
	}
}



