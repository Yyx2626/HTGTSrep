
options(stringsAsFactors=FALSE)

library(tidyverse)
library(eulerr)


`%.%` = function(x,y) paste0(x,y)


output_prefix = "eulerr_Venn.20220810"   # TODO: change to the output filename prefix you like
## output files would be output_prefix.pdf

input_filename = "input.tsv"   # TODO: change to your input filename
input = read_tsv(input_filename, col_names=TRUE)
## input tsv file should have header line (colnames) as sample names
## each column should be the raw sequences

smpl_colidxes = 1:ncol(input)   # TODO: change to the column indexes for your input file

sequence_list = list()
for(si in 1:length(smpl_colidxes)){
	sequence_list = c(sequence_list, list(na.omit(input[[smpl_colidxes[si]]])))
}
names(sequence_list) = colnames(input)[smpl_colidxes]


pdf(output_prefix %.% ".pdf", width=6, height=6)
par(mar=c(3,3,3,3))

euler(lapply(sequence_list  , unique)) %>% plot(quantities=TRUE) %>% print
euler(lapply(sequence_list  , unique)) %>% plot(labels=FALSE) %>% print

dev.off()
