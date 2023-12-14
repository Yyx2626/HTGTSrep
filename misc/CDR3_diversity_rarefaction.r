
options(stringsAsFactors=FALSE)

library(tidyverse)
library(iNEXT)

`%.%` = function(x,y) paste0(x,y)


output_prefix = "rarefaction.20220808"   # TODO: change to the output filename prefix you like
## output files would be output_prefix.iNEXT_DF.tsv , output_prefix.iNEXT_group_merged.tsv and output_prefix.pdf

input_filename = "input.tsv"   # TODO: change to your input filename
input = read_tsv(input_filename, col_names=TRUE)
## input tsv file should have header line (colnames) as sample names
## each column can be either the raw sequences (such as CDR3, if character) or the frequencies (if numeric)

smpl_colidxes = 1:ncol(input)   # TODO: change to the column indexes for your input file

smpl_groups = c("Ctrl", "Ctrl", "Ctrl", "Test", "Test", "Test")   # TODO: change to the actual group label for each sample in smpl_colidxes

stopifnot(length(smpl_colidxes)==length(smpl_groups))

group_vec = unique(smple_groups)

rarefaction_target_size_vec = c(1000000,500000,200000,100000,50000,20000,10000,5000,2000,1000,500,200,100,50,20,10,5,2,1)   # TODO: change to the proper sizes you like


### call iNEXT to calculate rarefaction

set.seed(1234567)

iNEXT_list = list()
system.time({
for(j in smpl_colidxes){
	now_freq_vec = input[[j]]
	if(is.character(now_freq_vec)){
		now_freq_vec = as.numeric(table(now_freq_vec))
	}
	iNEXT_list[[colnames(input)[j]]] = iNEXT(now_freq_vec, size=rarefaction_target_size_vec)
}
})


### rearrange rarefaction results
tmp_list = iNEXT_list
for(si in 1:length(smpl_colidxes)){
	j = smpl_colidxes[si]
	tmp_list[[colnames(input)[j]]]$smpl = colnames(input)[j]
	tmp_list[[colnames(input)[j]]]$grp = smpl_groups[si]
}
iNEXT_DF = tmp_list %>% bind_rows
dim(iNEXT_DF)

iNEXT_DF = iNEXT_DF %>% rename(N=m) %>% mutate(uniqPct_mean=qD/N, uniqPct_UB=qD.UCL/N, uniqPct_LB=qD.LCL/N)
dim(iNEXT_DF)

iNEXT_DF %>% write_tsv(output_prefix %.% ".iNEXT_DF.tsv")


grp_merged_DF = iNEXT_DF %>% group_by(grp, N) %>% summarize(Nsimu=n(), Nsmpl=(smpl %>% unique %>% length), uniqPct_mean=mean(uniqPct_mean), uniqPct_UB=max(uniqPct_UB), uniqPct_LB=min(uniqPct_LB), any_extrapolated=any(method=="extrapolated"))
dim(grp_merged_DF)   # 84  8

grp_merged_DF %>% write_tsv(output_prefix %.% ".iNEXT_group_merged.tsv")



### plot rarefaction curves (y=unique=qD/m ~ x=m)

grp_merged_DF$grp = factor(grp_merged_DF$grp, levels=group_vec)

selected_grp_vec = group_vec
group_colors = c("blue", "red")   # TODO: change to the colors you like

stopifnot(length(group_colors)==length(selected_grp_vec))

pdf(output_prefix %.% ".pdf", width=6, height=3.5)   # TODO: change width and height if you like
y_at = seq(0,1,by=0.25)

p = grp_merged_DF %>% filter(grp %in% selected_grp_vec) %>% ggplot(aes(x=N, y=uniqPct_mean, ymax=uniqPct_UB, ymin=uniqPct_LB, color=grp)) + geom_line(size=1, linetype="dashed") + geom_errorbar(width=0.1) + geom_line(data=grp_merged_DF %>% filter(grp %in% selected_grp_vec, !any_extrapolated), size=1, linetype="solid") + ggtitle("") + theme_classic() + scale_x_log10(breaks=10^(0:4)) + scale_y_continuous(breaks=y_at, labels=sprintf("%d", y_at*100)) + scale_color_manual(breaks=selected_grp_vec, values=grp_colors) + coord_cartesian(xlim=c(1,1e6), ylim=c(0,1.04), expand=FALSE) + xlab("Total CDR3 (N)") + ylab("Frequency of\nunique CDR3 (%)")
print(p)

p = grp_merged_DF %>% filter(grp %in% selected_grp_vec, !any_extrapolated) %>% ggplot(aes(x=N, y=uniqPct_mean, ymax=uniqPct_UB, ymin=uniqPct_LB, color=grp)) + geom_errorbar(width=0.1) + geom_line(size=1, linetype="solid") + ggtitle("") + theme_classic() + scale_x_log10(breaks=10^(0:4)) + scale_y_continuous(breaks=y_at, labels=sprintf("%d", y_at*100)) + scale_color_manual(breaks=selected_grp_vec, values=grp_colors) + coord_cartesian(xlim=c(1,1e5), ylim=c(0,1.04), expand=FALSE) + xlab("Total CDR3 (N)") + ylab("Frequency of\nunique CDR3 (%)")
print(p)

dev.off()

