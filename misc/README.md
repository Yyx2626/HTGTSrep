# misc

Some scripts for diversity rarefaction, Venn diagram, and sequence logo

Author: Adam Yongxin Ye @ Boston Children's Hospital / Harvard Medical School


## Scripts

Note: check TODO in scripts, and modify as you like


### CDR3\_diversity\_rarefaction.r

Input: `input_filename = "input.tsv"`
The input tsv file should have header line (colnames) as sample names
Each column can be either the raw sequences (such as CDR3, if character) or the frequencies (if numeric)

Output: `output_prefix`.iNEXT\_DF.tsv , `output_prefix`.iNEXT\_group\_merged.tsv and `output_prefix`.pdf

Demo output: [CDR3_diversity_rarefaction_demo_output.pdf](https://github.com/Yyx2626/HTGTSrep/tree/master/misc/CDR3_diversity_rarefaction_demo_output.pdf)


### clonotype\_diversity\_rarefaction.r

Input: `input_filename = "input.tsv"`
The input tsv file should have header line (colnames) as sample names
Each column can be either the raw sequences (such as CDR3, if character) or the frequencies (if numeric)

Output: `output_prefix`.tsv and `output_prefix`.pdf

Demo output: [clonotype_diversity_rarefaction_demo_output.pdf](https://github.com/Yyx2626/HTGTSrep/tree/master/misc/clonotype_diversity_rarefaction_demo_output.pdf)


### eulerr\_Venn.r

Input: `input_filename = "input.tsv"`
The input tsv file should have header line (colnames) as sample names
Each column can be either the raw sequences

Output: `output_prefix`.pdf
Demo output: [eulerr_Venn_demo_output.pdf](https://github.com/Yyx2626/HTGTSrep/tree/master/misc/eulerr_Venn_demo_output.pdf)


### yyx\_hideRef\_ggseqlogo.r

install the required R packages: tidyverse ggplot2 ggseqlogo readxl

modify the Excel file 'yyx\_hideRef\_ggseqlogo\_inupt.xlsx':  
I only use the second column, and regard the 1st line as reference sequence

change the `output_date` and `working_path` variables as you like

finally run this R script in R or Rscript

Input: yyx\_hideRef\_ggseqlogo\_inupt.xlsx

Output: yyx\_hideRef\_ggseqlogo\_output.`output_date`.pdf

Demo output: [yyx_hideRef_ggseqlogo_output.20230103.pdf](https://github.com/Yyx2626/HTGTSrep/tree/master/misc/yyx_hideRef_ggseqlogo_output.20230103.pdf)

