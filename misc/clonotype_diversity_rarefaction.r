
options(stringsAsFactors=FALSE)

library(tidyverse)
library(iNEXT)

`%.%` = function(x,y) paste0(x,y)


output_prefix = "rarefaction.20220322"   # TODO: change to the output filename prefix you like
## output files would be output_prefix.tsv and output_prefix.pdf

input_filename = "input.tsv"   # TODO: change to your input filename
input = read_tsv(input_filename, col_names=TRUE)
## input tsv file should have header line (colnames) as sample names
## each column can be either the raw sequences (such as CDR3, if character) or the frequencies (if numeric)

smpl_colidxes = 1:ncol(input)   # TODO: change to the column indexes for your input file

smpl_groups = c("Ctrl", "Ctrl", "Ctrl", "Test", "Test", "Test")   # TODO: change to the actual group label for each sample in smpl_colidxes

stopifnot(length(smpl_colidxes)==length(smpl_groups))

group_vec = unique(smple_groups)


rarefaction_target_size_vec = seq(100,5000,by=100)   # TODO: change to the proper sizes you like

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

qD_DF = data.frame(m=rarefaction_target_size_vec)
for(j in smpl_colidxes){
	qD_DF = qD_DF %>% merge(iNEXT_list[[colnames(input)[j]]]$iNextEst %>% select(m, qD) %>% filter(m %in% qD_DF$m))
	colnames(qD_DF)[ncol(qD_DF)] = colnames(input)[j] %.% ".qD"
}
dim(qD_DF)

for(gi in 1:length(group_vec)){
	now_group = group_vec[gi]
	j = which(smpl_groups==now_group) + 1
	if(j <= 0){
		cat("Warning: something wrong, no samples belong to group " %.% now_group %.% " ?!\n")
		next
	}
	now_dat_DF = qD_DF[,j]
	qD_DF[[now_group %.% ".N"]] = ncol(now_dat_DF)
	qD_DF[[now_group %.% ".qD_mean"]] = apply(now_dat_DF, 1, mean)
	qD_DF[[now_group %.% ".qD_sd"]] = apply(now_dat_DF, 1, sd)
	qD_DF[[now_group %.% ".qD_sem"]] = qD_DF[[now_group %.% ".qD_sd"]] / sqrt(qD_DF[[now_group %.% ".N"]])
}
dim(qD_DF)   # 50 25

for(gi in 1:(length(group_vec)-1)){
	now_group = group_vec[gi]
	j = which(smpl_groups==now_group) + 1
	if(j <= 0){
		cat("Warning: something wrong, no samples belong to group " %.% now_group %.% " ?!\n")
		next
	}
	now_dat_DF = qD_DF[,j]
	for(gj in (gi+1):length(group_vec)){
		now_group_2 = group_vec[gj]
		j = which(smpl_groups==now_group_2) + 1
		if(j <= 0){
			cat("Warning: something wrong, no samples belong to group " %.% now_group_2 %.% " ?!\n")
			next
		}
		now_dat_DF_2 = qD_DF[,j]
		qD_DF[[now_group %.% "_vs_" %.% now_group_2 %.% ".tTest_pVal"]] = NA
		for(i in 1:nrow(qD_DF)){
			qD_DF[[now_group %.% "_vs_" %.% now_group_2 %.% ".tTest_pVal"]][i] = t.test(unlist(now_dat_DF[i,]), unlist(now_dat_DF_2[i,]))$p.value
		}
	}
}
dim(qD_DF)   # 50 37

qD_DF$ANOVA_pVal = NA
for(gi in 1:(length(group_vec)-1)){
	now_group = group_vec[gi]
	for(gj in (gi+1):length(group_vec)){
		now_group_2 = group_vec[gj]
		qD_DF[[now_group %.% "_vs_" %.% now_group_2 %.% ".TukeyHSD_pVal"]] = NA
	}
}
for(i in 1:nrow(qD_DF)){
	tmpDF = data.frame(y=unlist(qD_DF[i, (1:length(smpl_colidxes))+1]))
	tmpDF$x = smpl_groups
	tmpfit = aov(y~x, data=tmpDF)
	qD_DF$ANOVA_pVal[i] = summary(tmpfit)[[1]][["Pr(>F)"]][1]
	tmpTukeyHSD = TukeyHSD(tmpfit)
	for(gi in 1:(length(group_vec)-1)){
		now_group = group_vec[gi]
		for(gj in (gi+1):length(group_vec)){
			now_group_2 = group_vec[gj]
			qD_DF[[now_group %.% "_vs_" %.% now_group_2 %.% ".TukeyHSD_pVal"]][i] = tmpTukeyHSD$x[now_group %.% "-" %.% now_group_2, "p adj"]
		}
	}
}
dim(qD_DF)   # 50 37

qD_DF %>% write_tsv(output_prefix %.% ".tsv")



### plot rarefaction curves (y=unique=qD/m ~ x=m)

group_vec
group_colors = c("blue", "red")   # TODO: change to the colors you like
group_shade_colors = hsv(c(0.67,0), 0.3, 1)   # TODO: change to the colors you like

stopifnot(length(group_colors)==length(group_vec))
stopifnot(length(group_shade_colors)==length(group_vec))

pdf(output_prefix %.% ".pdf", width=8, height=8)   # TODO: change width and height if you like
par(mar=c(3,3,3,1)+.1, mgp=c(2,0.8,0))

plot(0, type="n", xlim=c(0,max(rarefaction_target_size_vec)), ylim=c(0,1), xlab="Number of reads", ylab="Number of clonotypes / number of reads", main="Rarefaction curves")   # TODO: change xlab and ylab according to your input, such as Frequency of unique CDR3 (%) ~ Total CDR3 (N)
for(gi in 1:length(group_vec)){
	now_group = group_vec[gi]
	x = qD_DF$m
	y = qD_DF[[now_group %.% ".qD_mean"]] / x
	y_sem = qD_DF[[now_group %.% ".qD_sem"]] / x
	y_ub = y + y_sem
	y_lb = y - y_sem
	polygon(c(x, rev(x)), c(y_ub, rev(y_lb)), col=group_shade_colors[gi], border=NA)
	lines(x, y, col=group_colors[gi], lwd=2)
}

dev.off()

