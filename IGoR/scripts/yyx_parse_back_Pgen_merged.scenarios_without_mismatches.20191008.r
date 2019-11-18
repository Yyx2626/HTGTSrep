
#### Usage: Rscript this.r <query_mean_log10_Pgen> <title> <output_prefix> <Pgen_merged_filenames> ...

args = c("-6.50983782728612", "NP2", "scenarios_without_mismatches.parsed_back_Pgen_merged_20191008/debug", "run_igor_PP1_5_NP1_2_20190925/scenarios_without_mismatches/Pgen_merged.back_71_4models.20191008.tsv", "run_igor_PP1_5_NP1_2_20190925/scenarios_without_mismatches/Pgen_merged.back_72_4models.20191008.tsv")
args = c("-9.58015253617941", "PP11", "scenarios_without_mismatches.parsed_back_Pgen_merged_20191008/debug", "run_igor_PP6_13_20190925/scenarios_without_mismatches/Pgen_merged.back_6_4models.20191008.tsv")
args = c("-9.88303666161822", "PP12", "parsed_back_Pgen_merged_20191008/debug", "run_igor_PP6_13_20190925/Pgen_merged.20191002/back_7_4models.tsv")
if(!interactive()){
	args = commandArgs(TRUE)
}

#query_part = as.integer(args[1])
#query_seq_idx = as.integer(args[2])
query_log10_Pgen = as.numeric(args[1])
title = args[2]
output_prefix = args[3]


# setwd("/Users/yyx/Documents/work/huan/IGoR")
options(stringsAsFactors=FALSE)


`%.%` = function(x,y) paste0(x,y)

join = function(sep, vec, ...) paste(collapse=sep, c(vec, ...))
cat0 = function(...) cat(sep="", ...)

echo_str <- function(x, sep=" =\t", collapse=", ") deparse(substitute(x)) %.% sep %.% join(collapse, x)
echo <- function(x, sep=" =\t", collapse=", ") cat0(deparse(substitute(x)), sep, join(collapse, x), "\n");


final_DF = character(0)
query_DF = character(0)
#query_log10_Pgen = log10(query_Pgen)

for(Pgen_merged_filename in args[4:length(args)]){
	echo(Pgen_merged_filename)
	input = read.delim(Pgen_merged_filename)
	echo(dim(input))
	
	if(nrow(input) <= 0){
		next
	}

	Pgen_DF = input[,-1]
	log10_Pgen_DF = log10(Pgen_DF)
	names(log10_Pgen_DF) = "log10_" %.% names(log10_Pgen_DF)

	log10_Pgen_DF$mean = rowMeans(log10_Pgen_DF)

	log10_Pgen_DF$mean[log10_Pgen_DF$mean < -20] = -20
#	log10_Pgen_DF$mean[is.nan(log10_Pgen_DF$mean)] = -21

	names(log10_Pgen_DF)[names(log10_Pgen_DF)=="mean"] = "mean_log10_Pgen"
	output_DF = cbind(input, log10_Pgen_DF)
	output_DF$exp10_mean_log10_Pgen = 10^output_DF$mean_log10_Pgen
	echo(dim(output_DF))

#	if(Pgen_merged_filename == "run_igor_20190308/scenarios_without_mismatches/Pgen_merged.back_" %.% query_part %.% "_4models.20181207.tsv"){   # hard-coded, bad
#		query_rowidx = which(input$seq_index == query_seq_idx)
#		query_DF = output_DF[query_rowidx, ]
#		query_log10_Pgen = query_DF$mean_log10_Pgen[1]
#	}

	final_DF = rbind(final_DF, output_DF)
}
echo(dim(final_DF))

print(query_DF)

if(!is.nan(query_log10_Pgen)){
	gt_Pgen_count = sum(final_DF$mean_log10_Pgen > query_log10_Pgen, na.rm=TRUE)
	echo(gt_Pgen_count)
	ge_Pgen_count = sum(final_DF$mean_log10_Pgen >= query_log10_Pgen, na.rm=TRUE)
	echo(ge_Pgen_count)
	output_rowidx = which(final_DF$mean_log10_Pgen >= query_log10_Pgen - 1)
	ge_PgenDiv10_count = length(output_rowidx)
	echo(ge_PgenDiv10_count)
	output_rowidx = which(final_DF$mean_log10_Pgen >= query_log10_Pgen - 2)
	ge_PgenDiv100_count = length(output_rowidx)
	echo(ge_PgenDiv100_count)
}


pdf(output_prefix %.% ".pdf", width=5.5, height=5.5)
par(mar=c(3,3,3,1)+.1, mgp=c(2,0.8,0))

hist_result = hist(final_DF$mean_log10_Pgen, breaks=seq(-20.5, -4.5, by=1), plot=FALSE)
y_max = max(hist_result$density[-1])
plot(hist_result, freq=FALSE, xlab="mean log10(Pgen)", main=title %.% " back-translated CDR3 sequences", ylim=c(0,y_max))
rug(query_log10_Pgen)
#mtext("original CDR3 seq_index = " %.%  query_seq_idx %.% " in part back_" %.% query_part)
#my_text = "P(Pgen==NaN)\n= " %.% nan_Pgen_count %.% "/" %.% nrow(final_DF) %.% " = " %.% format(nan_Pgen_count / nrow(final_DF), digits=3)
if(!is.nan(query_log10_Pgen)){
my_text = "P(mean log10(Pgen) > " %.% format(query_log10_Pgen, digits=3) %.% ")\n"
my_text = my_text %.% "= " %.% gt_Pgen_count %.% "/" %.% nrow(final_DF) %.% " = " %.% format(gt_Pgen_count / nrow(final_DF), digits=3)
my_text = my_text %.% "\nP(mean log10(Pgen) >= " %.% format(query_log10_Pgen, digits=3) %.% ")\n"
my_text = my_text %.% "= " %.% ge_Pgen_count %.% "/" %.% nrow(final_DF) %.% " = " %.% format(ge_Pgen_count / nrow(final_DF), digits=3)
my_text = my_text %.% "\nP(mean log10(Pgen) >= " %.% format(query_log10_Pgen-1, digits=3) %.% ")\n"
my_text = my_text %.% "= " %.% ge_PgenDiv10_count %.% "/" %.% nrow(final_DF) %.% " = " %.% format(ge_PgenDiv10_count / nrow(final_DF), digits=3)
my_text = my_text %.% "\nP(mean log10(Pgen) >= " %.% format(query_log10_Pgen-2, digits=3) %.% ")\n"
my_text = my_text %.% "= " %.% ge_PgenDiv100_count %.% "/" %.% nrow(final_DF) %.% " = " %.% format(ge_PgenDiv100_count / nrow(final_DF), digits=3)
text(-4.5, y_max, adj=c(1,1), my_text, cex=0.8)
}
#text(-4.5, y_max, adj=c(1,1), my_text, cex=0.8)

dev.off()

if(!is.nan(query_log10_Pgen)){
output_DF = final_DF[output_rowidx, ]
output_DF = output_DF[order(-output_DF$mean_log10_Pgen, output_DF$seq_index),]

if(nrow(output_DF) > 600){
	output_DF = output_DF[1:600,]
}

write.table(output_DF, file=output_prefix %.% ".tsv", sep="\t", quote=FALSE, row.names=FALSE)
}

