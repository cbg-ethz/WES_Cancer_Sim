# rank-combination

Usage:

```
/pathTo/bin/Rscript  "rank_combination.R  combined_out_file.txt tool_1.vcf tool_2.vcf ... tool_n.vcf"
```

Command line arguments:

combined_out_file.txt 			= the output file with the combined and ranked variants to be generated
tool_1.vcf tool_2.vcf ... tool_n.vcf	= the vcf files of the tools to be combined. Should be at least two files.

Assumptions:
- The vcf files have a sixth column with the confidence score from the variant callers for ranking.
- In the case of MuTect, the seventh column is considered, which has the "ACCEPT" or "REJECT" label.

