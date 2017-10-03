# rank-combination

Important:
[The latest version of the rank-combination.R script can be found in our NGS-Pipe scripts folder](https://github.com/cbg-ethz/NGS-pipe/blob/master/scripts/rank_combination.R "rank-combination.R")


Usage:

```
/pathTo/bin/Rscript  rank_combination.R  combined_out_file.txt  tool_1.vcf  tool_2.vcf ... tool_n.vcf
```

Command line arguments:

combined_out_file.txt 			= the output file with the combined and ranked variants to be generated

tool_1.vcf tool_2.vcf ... tool_n.vcf	= the vcf files of the tools to be combined. Should be at least two files.

Assumptions:
- The vcf files have a sixth column with the confidence score from the variant callers for ranking. For the
  confidence score, the assumption is that higher scores represent a higher confidence in the variant call.
- In the case of MuTect, the seventh column is considered, which has the "ACCEPT" or "REJECT" label. 
  And the script recognizes that the vcf is from MuTect, by finding the substring "mutect" in its name.



Recommended filtering procedure afterwards:
The rank-combination.R provides a re-ranking of the variants from the different vcf files which were provided. The higher the
rank-combination score, the higher the confidence. We recommend prioritizing according to this new rank-combination score. Also, for each
of these variants, we still know the confidence score (i.e. the p-value/score/label) that they got from the original variant caller. 
So in order to get a final filtered list of high confidence variants, we recommend the following:
Starting with the variant that got the highest rank-combination score, then the next highest, etc. going down in the descending order of the
rank-combination score taking into the final list those variants that are still within acceptable limits of at least one of the original 
variant caller's confidence score. Stop, when none of the original variant callers had an acceptable confidence score any more.
A less stringent but more sensitive variation of this would be, to allow a certain percentage of lower confidence scores in the final
lists (e.g. 5%). If there are questions about the rank-combination.R script or the subsequent filtering, 
please E-Mail: ariane.hofmann@bsse.ethz.ch


