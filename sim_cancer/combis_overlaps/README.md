###
# scripts to run for intersecting/overlapping vcfs:
overlaps_main.sh: creates the 10%, 5%, and 1% overlap files and extracts the numbers of the vennDiagram, also computes the entire overlap files of selected tools for evaluation; calls the script extract_nums_forVennDiagram.sh and the programs SNV_Overlapper_Ordered.java, createOriginalVCF.java

extractOverlaps.sh: generates the vcf file of the overlap, in order to be able to evaluate the overlap of all tools; calls the script getSets.R and the porgram createOriginalVCF.java

###
# script for running all combinations of all nine tools:
start_all_combinations.sh: starts all combis; calls the script combininglistsbetter_Automated.R

###
# Now, the vcfs of the combis and overlaps can be evaluated. see ../eval_vcf/
